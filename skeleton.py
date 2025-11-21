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
        tagged_fatjets = cleaned_fatjets[overlap]
        mistagged_fatjets = cleaned_fatjets[~overlap]
        
        return {
            f"Tagged_fatjet_collection": tagged_fatjets,
            f"Mistagged_fatjet_collection": mistagged_fatjets,
            f"Genpart_collection": genpart,
        }
