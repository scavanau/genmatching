from typing import Any, Optional
import awkward as ak # type: ignore
from framework.module import Module
from framework.events import Events
from modules.helpers.objects_utils import extract_object, get_object
from modules.helpers.overlap_utils import calculate_delta_r, has_overlap

class Skeleton(Module):
    """Module for selecting monojet events based on jet kinematics.

    Attributes:
        collection (str): Prefix for output keys to distinguish collections.
        fatjets (str): Name of the cleaned Fatjet collection.
    """
    def __init__(self, cfg: dict[str, Any], events: Optional[Events] = None):
        """Initialize the Skeleton module."""
        super().__init__(cfg, events=events)
        self.fatjets: str = self.cfg.get("fatjets", "")

    def call(self, events: Events) -> dict[str, Any]:
        """Perform skeleton selection and return the selected jets."""
        cleaned_fatjets = extract_object(events=events, obj_name=self.fatjets)
        genpart = extract_object(events=events, obj_name="GenPart", variables=["pdgId", "status"])

        # Example selection: select Z bosons with status 62 from genpart collection
        genpart = genpart[(abs(genpart.pdgId) == 23) & (genpart.status == 62)]

        # Select fatjets that overlap with gen bosons
        overlap = has_overlap(obj_toclean=cleaned_fatjets, clean_against=genpart, max_dr=0.4)
        matched_fatjets = cleaned_fatjets[overlap]
        mismatched_fatjets = cleaned_fatjets[~overlap]

        print(f"Total number of fatjets: {ak.sum(ak.num(cleaned_fatjets))}")
        print(f"Number of tagged fatjets: {ak.sum(ak.num(matched_fatjets))}")
        print(f"Number of mistagged fatjets: {ak.sum(ak.num(mismatched_fatjets))}")

        print(f"cleaned_fatjet_collection: {cleaned_fatjets}")
        print(f"tagged_fatjet_collection: {matched_fatjets}")
        print(f"mistagged_fatjet_collection: {mismatched_fatjets}")
        print(f"genpart_collection: {genpart}")
        
        return {
            f"Tagged_fatjet_collection": matched_fatjets,
            f"Mistagged_fatjet_collection": mismatched_fatjets
        }
