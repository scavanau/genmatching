import uproot, ROOT, vector, os, glob
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np; np.seterr(all="ignore")
from array import array
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

nevents = []
pass_fatjet = []
pass_z = []
pass_dr = []
pteta_cut = []

pT_bins = [200, 250, 300, 350, 400, 500, 600, 700, 800]
pt_edges = array('d', pT_bins)

def make_hists():
    hists = {
        "h_pt":      ROOT.TH1F("h_pt", "FatJet p_{T};p_{T} [GeV];Events", 40, 0, 1000),
        "h_den_pt":  ROOT.TH1F("h_den_pt", "Denominator;p_{T} [GeV];Events", len(pT_bins)-1, pt_edges),
        "h_num_pt":  ROOT.TH1F("h_num_pt", "Numerator;p_{T} [GeV];Events",   len(pT_bins)-1, pt_edges),
        "h_den_msd": ROOT.TH1F("h_den_msoftdrop", "Denominator; msoftdrop [GeV];Events", 40, 0, 200),
        "h_num_msd": ROOT.TH1F("h_num_msoftdrop", "Numerator;   msoftdrop [GeV];Events", 40, 0, 200),
        "h_msd":     ROOT.TH1F("h_msd", "Msoftdrop;M_{softdrop};Events", 40, 0, 200),
        "h_den_dr":  ROOT.TH1F("h_den_dr", "All pairs;#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0),
        "h_num_dr":  ROOT.TH1F("h_num_dr", "Matched (pass #DeltaR);#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0),
        "h_dr":      ROOT.TH1F("h_dr", "Matched #DeltaR(FatJet, Z);#DeltaR;Events", 40, 0.0, 2.0),
    }
    for h in hists.values():
        h.SetDirectory(0)
    return hists

def process_file(input_file):
    h = make_hists()
    with uproot.open(input_file) as f:
        events = f["Events"]
        nevents.append(events.num_entries)
        
        fat_pt  = events["FatJet_pt"].array()
        fat_eta = events["FatJet_eta"].array()
        fat_phi = events["FatJet_phi"].array()
        fat_msd = events["FatJet_msoftdrop"].array()

        gen_pdg = events["GenPart_pdgId"].array()
        gen_eta = events["GenPart_eta"].array()
        gen_phi = events["GenPart_phi"].array()

    nonempty = ak.num(fat_pt) > 0
    fat_pt, fat_eta, fat_phi, fat_msd = fat_pt[nonempty], fat_eta[nonempty], fat_phi[nonempty], fat_msd[nonempty]
    gen_pdg, gen_eta, gen_phi = gen_pdg[nonempty], gen_eta[nonempty], gen_phi[nonempty]

    jet_mask = (fat_pt > 250) & (abs(fat_eta) < 2.4)
    fat_pt  = fat_pt[jet_mask]
    fat_eta = fat_eta[jet_mask]
    fat_phi = fat_phi[jet_mask]
    fat_msd = fat_msd[jet_mask]

    pteta_cut.append(jet_mask)

    hasjet = ak.num(fat_pt) > 0
    fat_pt  = fat_pt[hasjet]
    fat_eta = fat_eta[hasjet]
    fat_phi = fat_phi[hasjet]
    fat_msd = fat_msd[hasjet]
    gen_pdg = gen_pdg[hasjet]
    gen_eta = gen_eta[hasjet]
    gen_phi = gen_phi[hasjet]

    lead_pt  = fat_pt[:, 0]
    lead_eta = fat_eta[:, 0]
    lead_phi = fat_phi[:, 0]
    lead_msd = fat_msd[:, 0]

    pass_fatjet.append(hasjet)
    den_pt = ak.to_numpy(lead_pt)
    den_msd = ak.to_numpy(lead_msd)
        
    is_z = gen_pdg == 23
    z_eta = gen_eta[is_z]
    z_phi = gen_phi[is_z] 

    deta = lead_eta - z_eta
    dphi = np.abs(lead_phi - z_phi)
    dphi = ak.where(dphi > np.pi, 2*np.pi - dphi, dphi)
    dr = np.sqrt(deta**2 + dphi**2)
    
    best_dr = ak.min(dr, axis=1)
    
    matching = best_dr < 0.4

    pt       = (ak.to_numpy(lead_pt[matching]))
    msd      = (ak.to_numpy(lead_msd[matching]))
    num_pt   = (ak.to_numpy(lead_pt[matching]))
    num_msd  = (ak.to_numpy(lead_msd[matching]))
    num_dr   = (ak.to_numpy(best_dr[matching]))
    den_dr   = (ak.to_numpy(best_dr))
    dr_pass  = (ak.to_numpy(best_dr[matching]))

    pass_dr.append(matching)

    for val in pt:       h["h_pt"].Fill(val)
    for val in msd:      h["h_msd"].Fill(val)
    for val in num_dr:   h["h_num_dr"].Fill(val)
    for val in den_dr:   h["h_den_dr"].Fill(val)
    for val in num_pt:   h["h_num_pt"].Fill(val)
    for val in num_msd:  h["h_num_msd"].Fill(val)
    for val in den_msd:  h["h_den_msd"].Fill(val)
    for val in dr_pass:  h["h_dr"].Fill(val)

    output_file = ROOT.TFile(f"output_uproot/hist_{input_file.split('/')[-1]}", "RECREATE")
    for hist in h.values():
        hist.Write()
    output_file.Close()
    
def main():
    eos_path = "/eos/user/s/scavanau/Mono-V/Genmatching/Files"
    files = sorted(glob.glob(os.path.join(eos_path, "20250822_130953_Run3Summer22NanoAODv12_*.root")))[:32]
    for f in files:
        process_file(f)

    total_events    = sum(nevents)
    total_pteta     = sum(np.count_nonzero(x) for x in pteta_cut)
    total_fatjet    = sum(np.count_nonzero(x) for x in pass_fatjet)
    total_genmatch  = sum(np.count_nonzero(x) for x in pass_dr)

    print("============================")
    print(f"  Total events:    {total_events}")
    print(f"  Passing pt eta:  {total_pteta}")
    print(f"  Leading fatjet:  {total_fatjet}")
    print(f"  Genmatched:      {total_genmatch}")
    print("============================")

if __name__ == "__main__":
    main()