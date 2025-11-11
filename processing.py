import ROOT, os, glob
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
    
def process_file(input_file):
    tf = ROOT.TFile.Open(input_file)
    itree = tf.Get("Events")
    nEvents = itree.GetEntries()
    nevents.append(nEvents)

    h_msd = ROOT.TH1F("h_msd", "Msoftdrop; msoftdrop [GeV];Events", 40, 0, 200)
    h_den_msd = ROOT.TH1F("h_den_msd", "Denominator; msoftdrop [GeV];Events", 40, 0, 200)
    h_num_msd = ROOT.TH1F("h_num_msd", "Numerator;   msoftdrop [GeV];Events", 40, 0, 200)

    h_dr     = ROOT.TH1F("h_dr", "Matched #DeltaR(FatJet, Z);#DeltaR;Events", 40, 0.0, 2.0)
    h_den_dr = ROOT.TH1F("h_den_dr", "All pairs;#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0)
    h_num_dr = ROOT.TH1F("h_num_dr", "Matched (pass #DeltaR);#DeltaR(FatJet,Z);Events", 40, 0.0, 2.0)

    h_pt     = ROOT.TH1F("h_pt", "FatJet p_{T};p_{T} [GeV];Events", 40, 0, 1000)
    h_den_pt = ROOT.TH1F("h_den_pt", "Denominator;p_{T} [GeV];Events", len(pT_bins)-1, pt_edges)
    h_num_pt = ROOT.TH1F("h_num_pt", "Numerator;p_{T} [GeV];Events",   len(pT_bins)-1, pt_edges)

    for h in [h_msd, h_pt, h_dr, h_den_msd, h_num_msd, h_den_dr, h_num_dr, h_den_pt, h_num_pt]:
        h.SetDirectory(0)
        
    for ev in range(nEvents):
        itree.GetEntry(ev)
        if itree.nFatJet == 0:
            continue

        jet_pass = False
        pt_jet, eta_jet, phi_jet, mass_jet = 0, 0, 0, 0
        for j in range(itree.nFatJet):
            pt_jet = itree.FatJet_pt[j]
            eta_jet = itree.FatJet_eta[j]
            phi_jet = itree.FatJet_phi[j]
            mass_jet = itree.FatJet_msoftdrop[j]
            
            if not (pt_jet > 250 and abs(eta_jet) < 2.4):
                continue
            pt_jet, eta_jet, phi_jet, mass_jet = pt_jet, eta_jet, phi_jet, mass_jet
            jet_pass = True
            
            pteta_cut.append(ev)
            break
        if not jet_pass:
            continue
            
        pass_fatjet.append(ev)
        
        h_den_msd.Fill(mass_jet)
        h_den_pt.Fill(pt_jet)
        
        for i in range(itree.nGenPart):
            if not itree.GenPart_pdgId[i] == 23:
                continue
            pass_z.append(ev)
            
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
                h_msd.Fill(mass_jet)

                pass_dr.append(ev)
                break
                
    tf.Close()
    output_file = ROOT.TFile(f"output/hist_{input_file.split('/')[-1]}", "RECREATE")
    h_pt.Write(); h_den_pt.Write(); h_num_pt.Write()
    h_msd.Write(); h_den_msd.Write(); h_num_msd.Write()
    h_dr.Write(); h_den_dr.Write();  h_num_dr.Write()
    output_file.Close()


def main():
    eos_path = "/eos/user/s/scavanau/Mono-V/Genmatching/Files"
    files = sorted(glob.glob(os.path.join(eos_path, "20250822_130953_Run3Summer22NanoAODv12_*.root")))[:32]

    for f in files:
        process_file(f)

    total_events    = sum(nevents)
    total_pteta     = len(pteta_cut)
    total_fatjet    = len(pass_fatjet)
    total_genmatch  = len(pass_dr)

    print("============================")
    print(f"  Total events:    {total_events}")
    print(f"  Passing pt eta:  {total_pteta}")
    print(f"  Leading fatjet:  {total_fatjet}")
    print(f"  Genmatched:      {total_genmatch}")
    print("============================")


if __name__ == "__main__":
    main()