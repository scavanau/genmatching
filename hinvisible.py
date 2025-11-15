"""Modules for all H->invisible analyses."""

# pylint: disable=R0801,C0302
import os
from typing import Any
from collections import OrderedDict
from utils.pyrat.logger import initialize_logger
from utils.cms.year_run_utils import CampaignInfoMap
from modules.helpers.jet_utils import btagging_syst_variations
from modules.helpers.jet_utils import get_jec_sources, get_jer_sources
from framework.combined_module import build_combined_module_config
from framework.module import get_nano_version_from_datatier

logger = initialize_logger()
Config = dict[str, dict[str, Any]]


def create_module_config(name: str, **kwargs: Any) -> Config:
    """Create a configuration for a single module.

    Args:
        name (str): Name of the module.
        **kwargs (Any): Additional keyword arguments for the module.

    Returns:
        Config: Configuration dictionary for the module.
    """
    module = kwargs.get("module", name)
    return {name: {"name": name, "module": module, **kwargs}}


def get_jec_variations(campaign: str) -> dict[str, dict[str, str]]:
    """Generate a dictionary of JEC variations."""
    year = CampaignInfoMap().get_year(campaign=campaign)

    sources = ["Nominal"] + get_jec_sources(mode="total", year=year) + get_jec_sources(mode="regrouped", year=year) + get_jer_sources(mode="", year=year)
    sources = ["Nominal"]  # by default avoid jec unc to save time
    directions = ["up", "down", ""]
    jec_variations = {}

    for source in sources:
        is_jer = "jer" in source
        for direction in directions:
            if (source == "Nominal" and direction == "") or (source != "Nominal" and direction != ""):
                collection = f"jec_{source}_{direction}" if source != "Nominal" else f"jec_{source}"
                jec_variations[collection] = {
                    "jec_source": "Nominal" if is_jer else source,
                    "jec_direction": "" if is_jer else direction,
                    "jer_source": source if is_jer else "Nominal",
                    "jer_direction": direction if is_jer else "",
                }

    return jec_variations


def get_btagging_infos(campaign: str, datatier: str) -> dict[str, Any]:
    """Add btagging working points."""
    if get_nano_version_from_datatier(datatier) >= 14:
        tagger_type = "UParTAK4"
        branch = "btagUParTAK4B"
    else:
        tagger_type = "robustParticleTransformer"
        branch = "btagRobustParTAK4B"
    measurements = {"l": "light", "c": "mujets", "b": "mujets"}
    btagging_infos = {
        "type": tagger_type,
        "branch": branch,
        "workingpoint": "M",
        "measurements": measurements,
    }
    btagging_infos["variations"] = {
        meas_type: btagging_syst_variations(campaign=campaign, meas_type=meas_type, split="year") for meas_type in set(measurements.values())
    }

    return btagging_infos


def get_lepton_workingpoints() -> dict[str, list[str]]:
    """Return the mapping of working points to lepton object types."""
    working_points: dict[str, list[str]] = {
        "tight": ["Muon", "Electron", "Photon"],
        "loose": ["Tau"],
        "loose_not_tight": ["Muon", "Photon"],
        "veto_not_tight": ["Electron"],
    }

    return working_points


def get_general_cleaning(campaign: str) -> Config:
    """Generate the general cleaning configuration for H->invisible analyses."""
    configs = OrderedDict()
    configs.update(create_module_config(name="GenBosonFinder"))
    configs.update(create_module_config(name="GenLeptonFinder"))
    configs.update(create_module_config(name="NLOCorrections"))
    configs.update(create_module_config(name="MCWeights", list_of_weights=["weight_generator", "weight_pu", "weight_pdf", "weight_scale"]))
    configs.update(create_module_config(name="LumiMask"))
    configs.update(create_module_config(name="NoiseFilterSelection", selection_name="Noise_filter"))

    # Add JEC variations for MET
    for jec_collection, jecs in get_jec_variations(campaign=campaign).items():
        configs.update(create_module_config(name=f"MET_{jec_collection}", module="MET", collection=jec_collection, jecs=jecs))

    # Add cleaning selections
    cleaning_selections = {
        "bad_ecal_crystal_22ee": {"met_flavor": "jec_Nominal_TypeIPuppiMET"},
        "jet_veto_maps": {},
    }
    configs.update(create_module_config(name="DetectorMitigation", cleaning_selections=cleaning_selections))

    return configs


def get_lepton_cleaner(lepton: str, collection: str, working_point: str, **kwargs: Any) -> Config:
    """Generate configuration for lepton cleaning modules based on the lepton type and working point.

    Args:
        lepton (str): Type of lepton, e.g., "Electron", "Muon", "Photon", or "Tau".
        collection (str): Lepton collection name.
        working_point (str): Working point, e.g., "loose", "tight", etc.
        **kwargs: Additional keyword arguments, e.g., `clean_against` for specific lepton types.
    """
    configs = OrderedDict()

    # Define lepton-specific parameters
    lepton_cleaner_map = {
        "Electron": {
            "min_pt": 10 if working_point in ["veto", "loose"] else 40,
            "max_eta": 2.5,
            "lepton_id": "cutBased",
            "lepton_id_cut": working_point,
            "additional_variables": ["charge"],
            "clean_against": kwargs.get("clean_against", []),
        },
        "Muon": {
            "min_pt": 10 if working_point == "loose" else 20,
            "max_eta": 2.4,
            "lepton_id": f"{working_point}Id",
            "iso_variable": "pfRelIso04_all",
            "iso_cut": 0.25 if working_point == "loose" else 0.15,
            "additional_variables": ["charge"] + [f"{wp}Id" for wp in ["tight", "loose"]],
            "clean_against": kwargs.get("clean_against", []),
        },
        "Photon": {
            "min_pt": 20 if working_point == "loose" else 230,
            "max_eta": 2.5 if working_point == "loose" else 1.4442,
            "lepton_id": "cutBased",
            "lepton_id_cut": working_point,
            "additional_ids": {"electronVeto": working_point == "tight"},
            "clean_against": kwargs.get("clean_against", []),
        },
        "Tau": {
            "min_pt": 20,
            "max_eta": 2.5,
            "dz": 0.2,
            "additional_ids": {
                "idDecayModeOldDMs": working_point,  # TODO: Not fully implemented
                # "idDeepTau2018v2p5VSjet": "VVLoose", # TODO: Adjust based on year
                # "idDeepTau2017v2p1VSjet": "VVLoose",
                "idDeepTau2018v2p5VSjet": "Loose",  # TODO: Adjust based on year
                # "idDeepTau2018v2p5VSe": "Loose",
                # "idDeepTau2018v2p5VSmu": "Loose",
                # chargedIso
            },
            "clean_against": kwargs.get("clean_against", []),
        },
    }

    if lepton not in lepton_cleaner_map:
        logger.critical(f"Unsupported lepton type: {lepton}. Supported types are: {list(lepton_cleaner_map.keys())}", exception_cls=ValueError)

    infos: dict[str, Any] = {"collection": collection, "object_name": lepton, **lepton_cleaner_map[lepton]}

    configs.update(create_module_config(name=f"{collection}{lepton}Cleaner", module="LeptonCleaner", **infos))

    return configs


def get_jet_cleaner(campaign: str, jet: str, collection: str, datatier: str, **kwargs: Any) -> Config:
    """Generate configuration for jet cleaning modules based on the jet type and collection.

    Args:
        campaign (str): Data-taking campaign (e.g., 'UL16preVFP').
        datatier (str): NanoAOD Version (e.g., 'NanoAODv12').
        jet (str): Type of jet, e.g., "Jet" or "FatJet".
        collection (str): Jet collection name, e.g., "inclusive", "central".
        **kwargs: Additional keyword arguments, e.g., `clean_against` for specific objects to clean against.
    """
    configs = OrderedDict()

    eta_map = {
        "inclusive": 4.7,
        "central": 2.5,
        "veto": 4.7,
    }

    jet_cleaner_map = {
        "Jet": {
            "min_pt": 40 if collection == "veto" else 20,
            "max_eta": eta_map[collection],
            "clean_against": kwargs.get("clean_against", []),
            "additional_variables": ["chHEF", "neHEF", "chEmEF", "neEmEF", "muEF", get_btagging_infos(campaign=campaign, datatier=datatier)["branch"]],
            "additional_variables_mc_only": ["hadronFlavour"],
        },
        "SubJet": {
            "min_pt": 20,
            "max_eta": eta_map[collection],
            "clean_against": kwargs.get("clean_against", []),
        },
        "FatJet": {
            "min_pt": 250,
            "max_eta": eta_map[collection],
            "additional_variables": ["msoftdrop", "particleNetWithMass_WvsQCD"],
            # "pnet_md_tagger_cuts": {"Xqq": {"working_point": None, "cut": 0.9}},
        },
    }

    if jet not in jet_cleaner_map:
        logger.critical(f"Unsupported jet type: {jet}. Supported types are: {list(jet_cleaner_map.keys())}", exception_cls=ValueError)

    infos = {
        "object_name": jet,
        "flavor": jet,
        "jet_id": "tight",
        **jet_cleaner_map[jet],
        "jecs": {},
        "min_dr": 0.4 if jet != "FatJet" else 0.8,
    }

    if jet == "SubJet":
        infos["jet_id"] = None
    if get_nano_version_from_datatier(datatier) >= 14:
        infos["jet_id"] = None

    jec_variations = get_jec_variations(campaign=campaign)
    for jec_collection, jecs in jec_variations.items():
        infos["collection"] = f"{collection}_{jec_collection}"
        if collection == "inclusive":
            infos["jecs"] = jecs
            infos["jer"] = True
        else:
            infos["object_name"] = f"cleaned_inclusive_{jec_collection}_{jet}"
        if collection == "veto":
            infos["id_selections"] = [
                {"id_var": "chEmEF", "is_flag": False, "min_value": 0.5},
                {"id_var": "muEF", "is_flag": False, "min_value": 0.5},
            ]

        name = f"{collection}{jet}Cleaner_{jec_collection}"
        configs.update(create_module_config(name=name, module="JetCleaner", **infos))

    return configs


def get_scale_factors(correction: str, working_point: str, event_filter: Any, do_veto_eff: bool = False, sample_type: str = "") -> Config:
    """Generate scale factor configurations for a given correction and working point.

    Args:
        correction (str): Type of correction, e.g., "muon_Z", "electron", "photon", etc.
        working_point (str): Working point, e.g., "tight", "loose", etc.
        event_filter (Any): Event filter to apply for the scale factors.
        do_veto_eff (bool): Apply veto efficiency (1-SF).
        sample_type (str): Sample type to apply the scale factors to, e.g. "W", "Z", "ttbar", etc.
    """
    # Scale factor configurations mapped by correction and working point
    # tag names are following CMS combine conventions.
    sfs_map: Config = {
        "muon_Z": {
            "tight": {
                "correction_names": {
                    "id": "NUM_TightID_DEN_TrackerMuons",
                    "iso": "NUM_TightPFIso_DEN_TightID",
                },
            },
            "loose": {
                "force_minpt": 15,
                "correction_names": {
                    "veto_id": "NUM_LooseID_DEN_TrackerMuons",
                    "veto_iso": "NUM_LoosePFIso_DEN_LooseID",
                },
            },
        },
        "electron": {
            "tight": {"tag": "id", "working_point": "Tight", "do_high_pt": True},
            "veto": {"tag": "veto", "working_point": "Veto"},
        },
        "photon": {
            "tight": {"tag": "id", "working_point": "Tight", "do_high_pt": True},
            "loose": {"tag": "loose_id", "working_point": "Loose"},
        },
        "tau": {
            "loose": {"tag": "veto", "working_point": "Loose"},
        },
    }

    if correction not in sfs_map or working_point not in sfs_map[correction]:
        logger.critical(f"Invalid correction '{correction}' or working point '{working_point}'. Please verify the inputs.", exception_cls=ValueError)

    sfs = {
        "correction": correction,
        "event_filter": event_filter,
        "scale_factors": sfs_map[correction][working_point],
        "do_veto_eff": do_veto_eff,
        "sample_type": sample_type,
    }

    return sfs


def get_trigger_module(campaign: str, trigger: str, channels: list[str], apply_sf: bool = True) -> Config:
    """Generate configuration for trigger modules.

    Args:
        campaign (str): Data-taking campaign (e.g., 'UL16preVFP').
        trigger (str): Trigger type, e.g., "MET", "SingleElectron", or "SinglePhoton".
        channels (list[str]): Channels to run on (VBF, Monojet, MonoV).
        apply_sf (bool): Whether to include scale factors in the configuration.
    """
    configs = OrderedDict()
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18", "2022", "2022EE", "2023", "2023BPix", "2024"]
    jec_variations = get_jec_variations(campaign=campaign)

    for jec_collection in jec_variations:
        trigger_infos: dict[str, Any] = {
            "MET": {
                "trigger_selections": {
                    year: [
                        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                    ]
                    for year in years
                },
                "scale_factor": {
                    "name": "trigger_met",
                    "veto_selections": ["SingleElectron_trigger_selection"],
                    "correction_dict": {},
                },
            },
            "SingleElectron": {
                "trigger_selections": {
                    year: [
                        "HLT_Ele32_WPTight_Gsf",
                        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
                        "HLT_Photon200",
                    ]
                    for year in years
                },
                "scale_factor": {
                    "name": "trigger_electron",
                    "correction_dict": {
                        f"diElectron_CR_{jec_collection}": {
                            "json_name": "electron_trigger_sf",
                            "selections": [f"diElectron_CR_{jec_collection}_recoil_selection"],
                            "objects": ["diElectron_lepton_1", "diElectron_lepton_2"],
                            "correlate_scale_factor": True,
                        },
                        f"singleElectron_CR_{jec_collection}": {
                            "json_name": "electron_trigger_sf",
                            "selections": [f"singleElectron_CR_{jec_collection}_recoil_selection"],
                            "objects": [f"WFinder_singleElectron_{jec_collection}_lepton"],
                        },
                    },
                },
            },
            "SinglePhoton": {
                "trigger_selections": {year: ["HLT_Photon200"] for year in years},
                "scale_factor": {
                    "name": "trigger_photon",
                    "correction_dict": {
                        f"singlePhoton_CR_{jec_collection}": {
                            "json_name": "photon_trigger_sf",
                            "selections": [f"singlePhoton_CR_{jec_collection}_recoil_selection"],
                            "objects": ["cleaned_tight_Photon"],
                        },
                    },
                },
            },
        }

        for channel in channels:
            jet_flavor = {
                "VBF": "inclusive",
                "Monojet": "inclusive",
                "MonoV": "inclusive",
            }[channel]

            is_vbf = "VBF" == channel

            # Use 1D sf (vs METNoMu) for Monojet and MonoV, 2D sf (vs METNoMu and mjj) for VBF
            json_name = "METNoMu_2D_trigger_sf" if is_vbf else "METNoMu_trigger_sf"
            for region in ["diMuon_CR", "singleMuon_CR", "SR"]:
                variables = {"pt": (f"{region}_{jec_collection}_recoil", "pt")}

                if is_vbf:
                    variables.update({"mjj": (f"VBF_{jet_flavor}_{jec_collection}_dijet", "mass")})

                trigger_infos["MET"]["scale_factor"]["correction_dict"].update(
                    {
                        f"{channel}_{region}_{jec_collection}": {
                            "json_name": json_name,
                            "selections": [f"{region}_{jec_collection}_recoil_selection", f"{channel}_{jet_flavor}_{jec_collection}_selection"],
                            "variables_map": variables,
                        },
                    }
                )

        # Modify scale factors if not needed
        if not apply_sf:
            for val in trigger_infos.values():
                del val["scale_factor"]

        if trigger not in trigger_infos:
            logger.critical(f"Unsupported trigger type: {trigger}. Supported types are: {list(trigger_infos.keys())}", exception_cls=ValueError)

        infos: dict[str, Any] = {
            "selection_name": f"{trigger}_trigger",
            "enable_single_trigger_sf": True,
            **trigger_infos[trigger],
        }

        if trigger == "MET":
            infos["trigger_selections"].update({year: ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"] for year in ["UL16preVFP", "UL16postVFP", "UL17"]})

        elif trigger == "SingleElectron":
            infos["trigger_selections"].update(
                {
                    "UL16preVFP": ["HLT_Ele27_WPTight_Gsf", "HLT_Photon175"],
                    "UL16postVFP": ["HLT_Ele27_WPTight_Gsf", "HLT_Photon175"],
                    "UL17": ["HLT_Ele35_WPTight_Gsf", "HLT_Photon200"],
                }
            )

        elif trigger == "SinglePhoton":
            infos["trigger_selections"].update({year: ["HLT_Photon175"] for year in ["UL16preVFP", "UL16postVFP"]})

        configs.update(create_module_config(name=f"{trigger}_TriggerSelection", module="TriggerSelection", **infos))
        # configs.update(create_module_config(name=f"{trigger}_TriggerSelection_{jec_collection}", module="TriggerSelection", **infos)) TODO

    return configs


def add_lepton_cleaners(configs: Config, clean_against: list[str]) -> None:
    """Add configurations for lepton cleaning."""
    lepton_wps = get_lepton_workingpoints()
    for collection, leptons in lepton_wps.items():
        working_point = collection.replace("_not_tight", "")
        for lepton in leptons:
            # Skip cleaning for tight muons and electrons (they are reference objects)
            should_clean = not (working_point == "tight" and lepton in {"Muon", "Electron"})
            _clean_against = clean_against if should_clean else []
            if working_point == "tight" and lepton == "Photon":
                _clean_against = [x for x in _clean_against if "Photon" not in x]
            configs.update(get_lepton_cleaner(lepton=lepton, collection=collection, working_point=working_point, clean_against=_clean_against))


def add_jet_cleaners(configs: Config, clean_against: list[str], campaign: str, datatier: str) -> None:
    """Add configurations for jet cleaning."""
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="Jet", collection="inclusive", clean_against=clean_against))
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="Jet", collection="central", clean_against=clean_against))
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="Jet", collection="veto", clean_against=clean_against))
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="SubJet", collection="inclusive", clean_against=clean_against))
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="FatJet", collection="inclusive", clean_against=clean_against))
    configs.update(get_jet_cleaner(campaign=campaign, datatier=datatier, jet="FatJet", collection="central", clean_against=clean_against))


def add_dilepton_selections(configs: Config) -> None:
    """Add configurations for dilepton selections."""
    dilepton_common = {
        "module": "DileptonSelection",
        "z_mass_min": 60,
        "z_mass_max": 120,
        "min_pt_lep2": 15,
    }
    dilepton_configs = {
        "DiMuonSelection": {
            "flavor": "Muon",
            "leptons": "cleaned_tight_Muon",
            "min_pt_lep1": 20,
        },
        "DiElectronSelection": {
            "flavor": "Electron",
            "leptons": "cleaned_tight_Electron",
            "min_pt_lep1": 40,
        },
    }
    for name, params in dilepton_configs.items():
        configs.update(create_module_config(name=name, **dilepton_common, **params))


def add_single_lepton_selections(configs: Config, campaign: str) -> None:
    """Add configurations for single-lepton selections."""
    single_lepton_common = {
        "module": "WFinder",
        "min_w_mt": 0.0,
        "max_w_mt": 160.0,
    }
    single_lepton_infos = {
        "Muon": {
            "leptons": "cleaned_tight_Muon",
            "min_pt_lep": 20,
        },
        "Electron": {
            "leptons": "cleaned_tight_Electron",
            "min_pt_lep": 40,
            "min_met_pt": 100,
        },
    }
    for lepton, infos in single_lepton_infos.items():
        flavor = f"single{lepton}"
        for jec_collection in get_jec_variations(campaign=campaign):
            configs.update(
                create_module_config(
                    name=f"WFinder_{flavor}_{jec_collection}",
                    flavor=flavor,
                    collection=jec_collection,
                    met_flavor=f"{jec_collection}_TypeIPuppiMET",
                    **single_lepton_common,
                    **infos,
                )
            )


def add_lepton_scale_factors(configs: Config) -> None:
    """Add configurations for lepton SFs."""
    leptons_to_correct = {
        "diMuon_lepton_1": get_scale_factors(correction="muon_Z", working_point="tight", event_filter="diMuon_selection"),
        "diMuon_lepton_2": get_scale_factors(correction="muon_Z", working_point="tight", event_filter="diMuon_selection"),
        "WFinder_singleMuon_jec_Nominal_lepton": get_scale_factors(
            correction="muon_Z", working_point="tight", event_filter="WFinder_singleMuon_jec_Nominal_selection"
        ),
        "WFinder_singleElectron_jec_Nominal_lepton": get_scale_factors(
            correction="electron", working_point="tight", event_filter="WFinder_singleElectron_jec_Nominal_selection"
        ),
        "diElectron_lepton_1": get_scale_factors(correction="electron", working_point="tight", event_filter="diElectron_selection"),
        "diElectron_lepton_2": get_scale_factors(correction="electron", working_point="tight", event_filter="diElectron_selection"),
        "cleaned_tight_Photon": get_scale_factors(correction="photon", working_point="tight", event_filter="singlePhoton_selection"),
        "cleaned_veto_not_tight_Electron": get_scale_factors(
            correction="electron", working_point="veto", event_filter="SR_jec_Nominal_selection", do_veto_eff=True
        ),
        "cleaned_loose_not_tight_Muon": get_scale_factors(
            correction="muon_Z", working_point="loose", event_filter="SR_jec_Nominal_selection", do_veto_eff=True
        ),
        "cleaned_loose_Tau": get_scale_factors(correction="tau", working_point="loose", event_filter="SR_jec_Nominal_selection", do_veto_eff=True),
    }

    for lepton, infos in leptons_to_correct.items():
        configs.update(
            create_module_config(
                name=f"{lepton}_SF",
                module="LeptonScaleFactor",
                leptons=lepton,
                **infos,
            )
        )


def add_recoil_selections(configs: Config, campaign: str) -> None:
    """Add configurations for recoil selections."""
    region_leptons_map = {
        "diMuon_CR": {"leptons": ["diMuon_lepton_1", "diMuon_lepton_2"]},
        "diElectron_CR": {"leptons": ["diElectron_lepton_1", "diElectron_lepton_2"]},
        "singlePhoton_CR": {"leptons": ["cleaned_tight_Photon"]},
        "singleMuon_CR": {"leptons": ["WFinder_singleMuon_jec_Nominal_lepton"]},
        "singleElectron_CR": {"leptons": ["WFinder_singleElectron_jec_Nominal_lepton"]},
        "SR": {},
        "QCD_CR": {},
        "QCD_CR_closure_0p2": {},
        "QCD_CR_closure_0p3": {},
    }
    jec_variations = get_jec_variations(campaign=campaign)

    for region, leptons in region_leptons_map.items():
        leptons = leptons if isinstance(leptons, dict) else {}
        if region == "QCD_CR":
            dphi_cuts = {"min_dphi_high": 0.5}
        elif region == "QCD_CR_closure_0p2":
            dphi_cuts = {"min_dphi_low": 0.2, "min_dphi_high": 0.5}
        elif region == "QCD_CR_closure_0p3":
            dphi_cuts = {"min_dphi_low": 0.3, "min_dphi_high": 0.5}
        else:  # For SR and all other CRs
            dphi_cuts = {"min_dphi_low": 0.5}
        for jec_collection in jec_variations:
            configs.update(
                create_module_config(
                    name=f"Recoil_{region}_{jec_collection}",
                    module="RecoilSelection",
                    collection=f"{region}_{jec_collection}",
                    met_flavor=f"{jec_collection}_TypeIPuppiMET",
                    min_pt_recoil=250,
                    do_dphi_selection=True,
                    **dphi_cuts,
                    jets=f"cleaned_inclusive_{jec_collection}_Jet",
                    **leptons,
                )
            )

    configs.update(
        create_module_config(
            name="FilterRecoilSelection",
            module="TriggerSelection",
            selection_name="FilterRecoil",
            trigger_selections={
                year: [f"{region}_{jec_collection}_recoil_selection" for region in region_leptons_map for jec_collection in jec_variations]
                for year in CampaignInfoMap().get_allowed_years()
            },
        )
    )


def add_object_count_selections(configs: Config, campaign: str) -> None:
    """Add configurations for object veto selections."""
    region_objects_count_map = {}

    obj_counts: list[str] = []
    lepton_wps = get_lepton_workingpoints()
    for working_point, leptons in lepton_wps.items():
        for lepton in leptons:
            obj_counts.append(f"{working_point}_{lepton}")
    jec_variations = get_jec_variations(campaign=campaign)
    for jec_collection in jec_variations:
        obj_counts.append(f"central_{jec_collection}_Jet_btagging")
    for obj in obj_counts:
        region_objects_count_map[f"{obj}_exact_count_0"] = {"objects": f"cleaned_{obj}", "exact_count": 0}

    for region, count in region_objects_count_map.items():
        configs.update(
            create_module_config(
                name=f"ObjectCountFilter_{region}",
                module="ObjectCountFilter",
                collection=region,
                **count,
            )
        )


def add_single_photon_selections(configs: Config) -> None:
    """Add configurations for single-photon selections."""
    region = "singlePhoton"
    configs.update(
        create_module_config(
            name=f"ObjectCountFilter_{region}",
            module="ObjectCountFilter",
            collection=region,
            objects="cleaned_tight_Photon",
            exact_count=1,
        )
    )


def add_control_region_selections(configs: Config, campaign: str) -> None:
    """Add configurations for control region selections."""
    jec_variations = get_jec_variations(campaign=campaign)
    years = CampaignInfoMap().get_allowed_years()

    configs.update(
        create_module_config(
            name="tight_objects_veto",
            module="TriggerSelection",
            selection_name="tight_objects_veto",
            trigger_ands={year: [f"tight_{flavor}_exact_count_0_selection" for flavor in get_lepton_workingpoints()["tight"]] for year in years},
        )
    )

    for jec_collection in jec_variations:
        region_selections_map: dict[str, list[str]] = {
            "diMuon_CR": ["diMuon_selection"],
            "diElectron_CR": ["diElectron_selection"],
            "singlePhoton_CR": ["singlePhoton_selection"],
            "singleMuon_CR": [f"WFinder_singleMuon_{jec_collection}_selection"],
            "singleElectron_CR": [f"WFinder_singleElectron_{jec_collection}_selection"],
            "SR": ["tight_objects_veto_selection"],
            "QCD_CR": ["tight_objects_veto_selection"],
            "QCD_CR_closure_0p2": ["tight_objects_veto_selection"],
            "QCD_CR_closure_0p3": ["tight_objects_veto_selection"],
        }
        for region, selections in region_selections_map.items():
            selection_name = f"{region}_{jec_collection}"
            trigger_ands = {year: selections + [f"{selection_name}_recoil_selection"] for year in years}
            configs.update(create_module_config(name=selection_name, module="TriggerSelection", selection_name=selection_name, trigger_ands=trigger_ands))


def add_trigger_modules(configs: Config, campaign: str, channels: list[str]) -> None:
    """Add configurations for trigger modules."""
    trigger_selections = ["SingleElectron", "SinglePhoton", "MET"]
    for trigger in trigger_selections:
        configs.update(get_trigger_module(campaign=campaign, trigger=trigger, channels=channels))

    configs.update(
        create_module_config(
            name="FilterTriggerSelection",
            module="TriggerSelection",
            selection_name="FilterTrigger",
            trigger_selections={year: [f"{trigger}_trigger_selection" for trigger in trigger_selections] for year in CampaignInfoMap().get_allowed_years()},
        )
    )


def add_channel_selections(configs: Config, campaign: str, channels: list[str]) -> None:
    """Add configurations for channel selections."""
    run_vbf = "VBF" in channels
    run_monojet = "Monojet" in channels
    channels_settings: dict[str, dict[str, Any]] = {
        "VBF": {
            "module": "VBFSelection",
            "min_mjj": 200,
            "min_pt_jet1": 80,
            "min_pt_jet2": 40,
            "jets": "cleaned_JETS_Jet",
            "veto_list": [],
        },
        "MonoV": {
            "module": "MonojetSelection",
            "min_pt_jet": 250,
            "max_eta_jet": 2.5,
            "min_msoftdrop": 0.0,
            "jets": "cleaned_JETS_FatJet",
            "veto_list": ["VBF_JETS_selection"] if run_vbf else [],
        },
        "Monojet": {
            "module": "MonojetSelection",
            "min_pt_jet": 100,
            "max_eta_jet": 2.5,
            "jets": "cleaned_JETS_Jet",
            "veto_list": ["VBF_JETS_selection"] if run_vbf else [],
            # "veto_list": ["VBF_JETS_selection", "MonoV_JETS_selection"] #TODO
        },
    }

    jec_variations = get_jec_variations(campaign=campaign)
    channel_selections = {}
    for jec_collection in jec_variations:
        for coll in ["inclusive", "central"]:
            jets = f"{coll}_{jec_collection}"
            for channel in channels:
                settings = channels_settings[channel].copy()
                settings["jets"] = settings["jets"].replace("JETS", jets)
                settings["veto_list"] = [sel.replace("JETS", jets) for sel in settings["veto_list"]]
                channel_selections[f"{channel}_{jets}"] = {**settings}

    for mode, config in channel_selections.items():
        configs.update(create_module_config(name=mode, collection=mode, **config))

    configs.update(
        create_module_config(
            name="Skeleton",
            collection="test",
            fatjets="MonoV_central_jec_Nominal_Fatjets",
        )
    )
    
    if run_monojet:
        configs.update(create_module_config(name="PhotonPurity", jets="Monojet_inclusive_jec_Nominal_jets", photons="Photon"))
        channel_selections["PhotonPurity"] = {"module": "PhotonPurity", "name": "PhotonPurity"}

    regions = [
        "diMuon_CR",
        "diElectron_CR",
        "singlePhoton_CR",
        "singleMuon_CR",
        "singleElectron_CR",
        "SR",
        "QCD_CR",
        "QCD_CR_closure_0p2",
        "QCD_CR_closure_0p3",
    ]
    analysis_channels = list(channel_selections.keys())
    for jec_collection in jec_variations:
        for region in regions:
            analysis_channels.append(f"{region}_{jec_collection}")

    configs.update(
        create_module_config(
            name="FilterChannel",
            module="TriggerSelection",
            selection_name="FilterChannel",
            trigger_selections={year: [f"{mode}_selection" for mode in analysis_channels] for year in CampaignInfoMap().get_allowed_years()},
        )
    )


def add_prefiring_weights(configs: Config, campaign: str, channels: list[str]) -> None:
    """Add configurations for prefiring weights."""
    run_monojet = "Monojet" in channels
    if not run_monojet:
        return
    jec_variations = get_jec_variations(campaign=campaign)
    jet_collections_map = {f"inclusive_{jec_collection}": f"cleaned_inclusive_{jec_collection}_Jet" for jec_collection in jec_variations}
    for jet_collection in jet_collections_map.keys():
        configs.update(
            create_module_config(
                name="PrefiringWeights",
                variable=f"Monojet_{jet_collection}_jets",
                correction="prefiring_sf_Jet",
                event_filter=f"Monojet_{jet_collection}_selection",
            )
        )
    # TODO prefiring weight is getting duplicated for every single jec


def all_dnn_interference(configs: Config, campaign: str) -> None:
    """Add configurations for dnn interference."""
    jec_variations = get_jec_variations(campaign=campaign)
    clean_against = ["cleaned_tight_Photon"]
    for channel in ["Muon", "Electron"]:
        clean_against.append(f"di{channel}_lepton_1")
        clean_against.append(f"di{channel}_lepton_2")
        for jec in jec_variations:
            clean_against.append(f"WFinder_single{channel}_{jec}_lepton")

    configs.update(
        create_module_config(
            name="PFCandCleaner",
            object_name="PFCands",
            min_pt=0.2,
            max_eta=5.2,
            min_dr=0.1,
            additional_variables=["fromPV", "pdgId", "charge", "puppiWeight"],
            pt_pdg_id_cleaning=[
                {"pdgId": 22, "pt": 1},
                {"pdgId": 130, "pt": 1},
                {"pdgId": 1, "pt": 3},
                {"pdgId": 2, "pt": 3},
                {"pdgId": 0, "pt": 3},
            ],
            clean_against=clean_against,
            from_pv=3,
        )
    )

    event_selections = [
        f"{sel}_selection"
        for sel in [
            "lumiMask",
            "Noise_filter",
            "DetectorMitigation",
            "FilterTrigger",
            "FilterRecoil",
            "FilterChannel",
        ]
    ]
    models = ["model_ops12", "best_model_PT_01", "best_model_PN_01"]

    model_path = os.path.join(str(os.getenv("PYRAT_PATH")), "hinvisible", "data", "dnn_models")
    for model in models:
        transformer_name = os.path.join(model_path, f"quantile_transformer_{model}.pkl")
        is_pnet = "PT" not in model
        points_feature = ["eta", "phi"] if is_pnet else ["pt", "eta", "phi", "energy"]
        feature_list = [
            "pt",
            "eta",
            "phi",
            "energy",
            "pdgId",
            "charge",
            "puppiWeight",
        ]
        if "model_ops12" in model:
            feature_list.extend(["energy_log", "pt_log"])
        configs.update(
            create_module_config(
                name=f"VBFTagger_{model}",
                module="VBFTagger",
                collection=model,
                input_variable="cleaned_PFCands",
                model_name=os.path.join(model_path, f"{model}.onnx"),
                transformer_name=transformer_name,
                event_selections=event_selections,
                feature_list=feature_list,
                points_feature=points_feature,
                use_pf_mask=is_pnet,
                pfcand_before_variables="model_ops12" in model,
            )
        )


def add_btagging_weights(configs: Config, campaign: str, datatier: str) -> None:
    """Add configurations for btagging weights."""
    jec_variations = get_jec_variations(campaign=campaign)
    btagging_info = get_btagging_infos(campaign=campaign, datatier=datatier)
    for jec_collection in jec_variations:
        jets = f"cleaned_central_{jec_collection}_Jet"
        configs.update(
            create_module_config(
                name=f"JetTagger_{jets}",
                module="JetTagger",
                jets=jets,
                flavor="Jet",
                jet_btagging=btagging_info,
            )
        )
        configs.update(
            create_module_config(
                name=f"JetTaggingWeights_{jets}",
                module="JetTaggingWeights",
                jets=jets,
                jet_btagging=btagging_info,
                collection=jec_collection,
            )
        )


def add_lhebugfix_weights(configs: Config) -> None:
    """Add temporary bug fixing reweighting modules before NLO PT binned samples are produced."""
    gjets_json_path = {
        "amcatnlo": os.path.join(str(os.getenv("PYRAT_PATH")), "hinvisible", "data", "gjets_bugfix", "amcatnlo_to_amcatnlo.json"),
        "madgraph": os.path.join(str(os.getenv("PYRAT_PATH")), "hinvisible", "data", "gjets_bugfix", "madgraph_to_sherpa.json"),
    }
    configs.update(create_module_config(name="LHEBugFix", bug_fixes=["gjets", "dyto2l", "dyto2nu", "wtolnu"], gjets_json=gjets_json_path))


def add_cleaning_cuts(configs: Config, campaign: str) -> None:
    """Add cleaning cuts between PUPPI and Calo MET."""
    jec_variations = get_jec_variations(campaign=campaign)
    regions = [
        "diMuon_CR",
        "diElectron_CR",
        "singlePhoton_CR",
        "singleMuon_CR",
        "singleElectron_CR",
        "SR",
        "QCD_CR",
        "QCD_CR_closure_0p2",
        "QCD_CR_closure_0p3",
    ]
    for region in regions:
        for jec_collection in jec_variations:
            configs.update(
                create_module_config(
                    name=f"METCleaningFilters_{region}_{jec_collection}",
                    module="METCleaningFilters",
                    collection=f"{region}_{jec_collection}",
                    met_flavor=f"{jec_collection}_TypeIPuppiMET",
                    recoil_flavor=f"{region}_{jec_collection}_recoil",
                    rel_diff_thrs={"tight": 0.4, "medium": 0.5, "loose": 0.6},
                    dphi_thrs={"tight": 1.5, "medium": 2.0, "loose": 2.2},
                )
            )

    configs.update(
        create_module_config(
            name="MCSpikeKiller",
            sanity_checks=[
                {"num": "leading_reco_jet_pt", "den": "gen_ht", "max_value": 1.0},
                {"num": "leading_reco_jet_pt", "den": "qscale", "max_value": 2.3},
            ],
            jets="cleaned_inclusive_jec_Nominal_Jet",
        )
    )


def add_photon_purity_weights(configs: Config) -> None:
    """Return photon pt-dependent photon purity weights."""
    configs.update(
        create_module_config(
            name="PhotonPurityWeights",
            module="PhotonPurityWeights",
            collection="tight",
            photons="cleaned_tight_Photon",
            formula="exp_plus_offset",
            parameters={
                "nominal": {"amp_scale": 172.6476, "exp_decay": 0.0385, "offset": 0.0296},
                "up": {"amp_scale": 215.8095, "exp_decay": 0.0385, "offset": 0.037},
                "down": {"amp_scale": 129.4857, "exp_decay": 0.0385, "offset": 0.0222},
            },
        )
    )


def create_configs(campaign: str, datatier: str) -> Config:
    """Create the configuration dictionary (OrderedDict) for all modules available in the analysis."""
    configs = OrderedDict()
    channels = ["VBF", "Monojet", "MonoV"]
    channels = ["MonoV"]
    clean_against = ["cleaned_tight_Electron", "cleaned_tight_Muon", "cleaned_tight_Photon"]

    # Add configurations
    configs.update(get_general_cleaning(campaign=campaign))
    add_lepton_cleaners(configs=configs, clean_against=clean_against)
    add_jet_cleaners(configs=configs, clean_against=clean_against, campaign=campaign, datatier=datatier)
    add_btagging_weights(configs=configs, campaign=campaign, datatier=datatier)
    add_dilepton_selections(configs=configs)
    add_single_lepton_selections(configs=configs, campaign=campaign)
    add_single_photon_selections(configs=configs)
    add_recoil_selections(configs=configs, campaign=campaign)
    add_object_count_selections(configs=configs, campaign=campaign)
    add_control_region_selections(configs=configs, campaign=campaign)
    add_channel_selections(configs=configs, campaign=campaign, channels=channels)
    add_trigger_modules(configs=configs, campaign=campaign, channels=channels)
    add_lepton_scale_factors(configs=configs)
    add_cleaning_cuts(configs=configs, campaign=campaign)
    add_prefiring_weights(configs=configs, campaign=campaign, channels=channels)
    all_dnn_interference(configs=configs, campaign=campaign)
    add_lhebugfix_weights(configs=configs)  # TODO: temporary bug fixing module
    add_photon_purity_weights(configs=configs)

    objects_to_keep = ["LHEPart"]
    configs.update(create_module_config(name="ObjectsExtractor", module="ObjectsExtractor", objects_to_keep=objects_to_keep))

    # Add the final combined module configuration
    configs["FullAnalysis"] = build_combined_module_config(
        configs=configs,
        apply_filter=True,
        filter_list=["FilterTrigger_selection", "FilterRecoil_selection", "FilterChannel_selection"],
        # variables_to_store=["weight*", "*recoil*", "*selection*"], #TODO simplify output for jec variations
    )

    return configs
