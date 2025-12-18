from typing import Any, Optional
import numpy as np
import awkward as ak  # type: ignore
from framework.module import Module
from framework.events import Events
from modules.helpers.genlevel_utils import ParticleData
from modules.helpers.objects_utils import extract_object, get_all_pairs, add_4vectors, build_awkward_array, add_variable_to_object, sort_objects
from modules.helpers.selection_utils import pt_eta_selection
from modules.helpers.modules_utils import get_sample_info
from modules.helpers.overlap_utils import has_overlap


class GenHadbosonFinder(Module):
    """Module for identifying generated bosons (Z, W, photon) at the generator level.

    Attributes:
        input_filename (str): The name of the input file.
        finder_mode (Optional[str]): Mode for boson identification (e.g., "Z", "W", "photon").
    """

    def __init__(self, cfg: dict[str, Any], events: Optional[Events] = None):
        """Initialize the GenBosonFinder module."""
        super().__init__(cfg, events=events)
        self.sample: str = str(self.cfg.get("sample", None))
        self.fatjets: str = self.cfg.get("fatjets", "")
        self.validate_parameters()

    def call(self, events: Events) -> dict[str, ak.Array]:
        """Identify and extract generated bosons from event data."""
        map_n_bosons = {
            "WW_": 2,
            "WZ_": 2,
            "ZZ_": 2,
            "Wto2Q": 1,
            "Zto2Q": 1,
            "TTtoLNu2Q": 1,
            "TTto4Q": 2,
            "WGto2QG": 1,
            "ZGto2QG": 1,
            "TbarQto2Q": 1,
            "TQbarto2Q": 1,
        }
        n_bosons = 0
        for name, count in map_n_bosons.items():
            if self.sample.startswith(name):
                n_bosons = count
                break
        cleaned_fatjets = extract_object(events=events, obj_name=self.fatjets)
        empty_mask = ak.zeros_like(cleaned_fatjets.pt, dtype=bool)
        out = {
            "matched_fatjets": cleaned_fatjets[empty_mask],
            "unmatched_fatjets": cleaned_fatjets,
            "matched_selection": ak.any(empty_mask, axis=1)
        }
        if n_bosons == 0:
            return out

        gen_particles = extract_object(events=events, obj_name="GenPart", variables=["pdgId", "status"])
        gen_bosons = gen_particles[(gen_particles.status == 62) | (gen_particles.status == 22)]
        gen_bosons = gen_bosons[gen_bosons.pt > 20]
        gen_bosons = gen_bosons[(abs(gen_bosons.pdgId) == 24) | (abs(gen_bosons.pdgId) == 23)]
        overlap = has_overlap(obj_toclean=cleaned_fatjets, clean_against=gen_bosons, max_dr=0.4)

        out["matched_fatjets"] = cleaned_fatjets[overlap]
        out["unmatched_fatjets"] = cleaned_fatjets[~overlap]
        out["matched_selection"] = ak.any(overlap, axis=1)
        return out

    def validate_parameters(self) -> None:
        """Validate the module's parameters."""
        pass
