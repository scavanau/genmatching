import ROOT, os, glob
import matplotlib.pyplot as plt
import numpy as np
from array import array
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

pT_bins = [200, 250, 300, 350, 400, 500, 600, 700, 800]
def pt_bin_index(x, edges):
    for k in range(len(edges) - 1):
        if edges[k] <= x < edges[k+1]:
            return k
    return None
    
def process_file(input_file):
    tf = ROOT.TFile.Open(input_file)
    itree = tf.Get("Events")
    nEvents = itree.GetEntries()

    h_den_msd = ROOT.TH1F("h_den_msoftdrop", "Denominator; msoftdrop [GeV];Events", 40, 0, 200)
    h_num_msd = ROOT.TH1F("h_num_msoftdrop", "Numerator;   msoftdrop [GeV];Events", 40, 0, 200)

    h_dr     = ROOT.TH1F("h_dr", "Matched #DeltaR(FatJet, Z);#DeltaR;Events", 40, 0.0, 2.0)
    h_den_dr = ROOT.TH1F("h_den_dr", "All pairs;#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0)
    h_num_dr = ROOT.TH1F("h_num_dr", "Matched (pass #DeltaR);#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0)

    h_pt     = ROOT.TH1F("h_pt", "FatJet p_{T};p_{T} [GeV];Events", 40, 0, 1000)
    pt_edges = array('d', pT_bins)
    h_den_pt = ROOT.TH1F("h_den_pt", "Denominator;p_{T} [GeV];Events", len(pT_bins)-1, pt_edges)
    h_num_pt = ROOT.TH1F("h_num_pt", "Numerator;p_{T} [GeV];Events",   len(pT_bins)-1, pt_edges)

    h_den_msd_bypt = []
    h_num_msd_bypt = []
    for i in range(len(pT_bins)-1):
        lo, hi = pT_bins[i], pT_bins[i+1]
        h_den_msd_bypt.append(
            ROOT.TH1F(f"h_den_msoftdrop_pt{i}", f"Denominator msoftdrop; msoftdrop [GeV];Events ({lo}<p_{{T}}<{hi} GeV)", 40, 0, 200)
        )
        h_num_msd_bypt.append(
            ROOT.TH1F(f"h_num_msoftdrop_pt{i}", f"Numerator msoftdrop; msoftdrop [GeV];Events ({lo}<p_{{T}}<{hi} GeV)", 40, 0, 200)
        )
    for h in [h_pt, h_dr, h_den_msd, h_num_msd, h_den_dr, h_num_dr, h_den_pt, h_num_pt, *h_den_msd_bypt, *h_num_msd_bypt]:
        h.SetDirectory(0)
        
    for ev in range(nEvents):
        itree.GetEntry(ev)
        if itree.nFatJet == 0:
            continue
        pt_jet = itree.FatJet_pt[0]
        eta_jet = itree.FatJet_eta[0]
        phi_jet = itree.FatJet_phi[0]
        mass_jet = itree.FatJet_msoftdrop[0]

        if not (pt_jet > 250 and abs(eta_jet) < 2.4):
            continue
            
        ibin = pt_bin_index(pt_jet, pT_bins)
        h_den_msd.Fill(mass_jet)
        h_den_pt.Fill(pt_jet)
        if ibin is not None:
            h_den_msd_bypt[ibin].Fill(mass_jet)
        
        for i in range(itree.nGenPart):
            if not itree.GenPart_pdgId[i] == 23:
                continue
            dphi = abs(phi_jet - itree.GenPart_phi[i])
            if dphi > np.pi:
              dphi = 2*np.pi - dphi
            deta = eta_jet - itree.GenPart_eta[i]
            dr   = (deta*deta + dphi*dphi) ** 0.5

            h_den_dr.Fill(dr)
            
            if dr < 0.4:
                h_num_msd.Fill(mass_jet)
                h_num_dr.Fill(dr)
                h_num_pt.Fill(pt_jet)
                h_dr.Fill(dr)
                h_pt.Fill(pt_jet)
                if ibin is not None:
                    h_num_msd_bypt[ibin].Fill(mass_jet)
                break
                
    tf.Close()
    output_file = ROOT.TFile(f"output/hist_{input_file.split('/')[-1]}", "RECREATE")
    h_pt.Write()
    h_den_pt.Write(); h_num_pt.Write()
    h_den_msd.Write(); h_num_msd.Write()
    h_den_dr.Write();  h_num_dr.Write(); h_dr.Write()
    for h in h_den_msd_bypt + h_num_msd_bypt:
        h.Write()
    output_file.Close()

def main():
    eos_path = "/eos/user/s/scavanau/Mono-V/Genmatching/Files"
    files = sorted(glob.glob(os.path.join(eos_path, "20250822_130953_Run3Summer22NanoAODv12_*.root")))[:31]

    for f in files:
        process_file(f)

if __name__ == "__main__":
    main()
