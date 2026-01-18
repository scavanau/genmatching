#!/usr/bin/env python3
import ROOT
from array import array
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def plot(obj, frame_hist, title, ytitle, pdf_name, ymin, ymax, xmin, xmax, draw_opt):
    c = ROOT.TCanvas("c", "", 600, 600)
    c.SetGrid(1, 1)

    frame = frame_hist.Clone("frame")
    frame.Reset()
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("FatJet p_{T} [GeV]")
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetXaxis().SetRangeUser(xmin, xmax)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.SetTitle(title)
    frame.Draw()

    obj.SetLineColor(ROOT.kBlue + 1)
    obj.SetMarkerColor(ROOT.kBlue + 1)
    obj.SetLineWidth(2)
    obj.SetMarkerStyle(20)
    obj.SetMarkerSize(0.8)

    obj.Draw(draw_opt)

    # save locally
    c.SaveAs(pdf_name)
    # save to EOS
    fname = os.path.basename(pdf_name)
    alt = os.path.join("/eos/user/s/scavanau/Mono-V/efficiency/Plots/23/spho", fname)
    c.SaveAs(alt)


def main():
    # singleMuon_CR_MonoV   singleElectron_CR_MonoV   singlePhoton_CR_MonoV   diElectron_CR_MonoV    diMuon_CR_MonoV
    all_file    = "/eos/user/s/scavanau/Mono-V/efficiency/0/Run3Summer23/NanoAODv12/singlePhoton_CR_MonoV/plots/root/monov_FatJet_pt.root"
    tagged_file = "/eos/user/s/scavanau/Mono-V/efficiency/0p5/Run3Summer23/NanoAODv12/singlePhoton_CR_MonoV/plots/root/monov_FatJet_pt.root"
    hist_name   = "data_obs"

    f_all = ROOT.TFile(all_file)
    f_tag = ROOT.TFile(tagged_file)

    # ---------------- data histos ----------------
    h_all = f_all.Get(hist_name)
    h_tag = f_tag.Get(hist_name)

    # ---------------- MC summed histos ----------------
    h_mc_all = h_all.Clone("h_all_mc")
    h_mc_all.Reset()
    h_mc_all.SetDirectory(0)

    h_mc_tag = h_tag.Clone("h_tag_mc")
    h_mc_tag.Reset()
    h_mc_tag.SetDirectory(0)

    mc_names = []
    for key in f_all.GetListOfKeys():
        name = key.GetName()
        if name == hist_name:
            continue  # skip data
        if name.lower().startswith("stack"):
            continue  # skip any stack

        # MC in 'all' file (required)
        h_all_mc = f_all.Get(name)
        if not isinstance(h_all_mc, ROOT.TH1):
            continue

        mc_names.append(name)
        h_mc_all.Add(h_all_mc)

        # MC in 'tag' file (optional)
        h_tag_mc = f_tag.Get(name)
        if isinstance(h_tag_mc, ROOT.TH1):
            h_mc_tag.Add(h_tag_mc)
        else:
            # if missing in tag file, treat as zero tagged events
            print(f"WARNING: '{name}' not found (or not TH1) in tagged file – counted as 0 tagged")

    print("MC processes in this control region:", mc_names)

    # ---------------- rebin to coarse bins ----------------
    edges = array('d', [200.0, 300.0, 400.0, 800.0])
    nbins = len(edges) - 1

    h_all_rb    = h_all.Rebin(nbins,    "h_all_rb",    edges)
    h_tag_rb    = h_tag.Rebin(nbins,    "h_tag_rb",    edges)
    h_mc_all_rb = h_mc_all.Rebin(nbins, "h_mc_all_rb", edges)
    h_mc_tag_rb = h_mc_tag.Rebin(nbins, "h_mc_tag_rb", edges)

    # ---------------- efficiencies ----------------
    eff_data = ROOT.TEfficiency(h_tag_rb,    h_all_rb)
    eff_mc   = ROOT.TEfficiency(h_mc_tag_rb, h_mc_all_rb)
    eff_data.SetName("h_eff_data")
    eff_mc.SetName("h_eff_mc")

    # ---------------- scale factors ----------------
    h_sf = h_all_rb.Clone("h_sf")
    h_sf.Reset()
    h_sf.SetDirectory(0)

    print("\nScale factors (Data/MC):")
    for i in range(1, h_sf.GetNbinsX() + 1):
        eps_d = eff_data.GetEfficiency(i)
        eps_m = eff_mc.GetEfficiency(i)
        sf = 0.0
        err_sf = 0.0
        if eps_m > 0 and eps_d > 0:
            err_d = max(eff_data.GetEfficiencyErrorUp(i),
                        eff_data.GetEfficiencyErrorLow(i))
            err_m = max(eff_mc.GetEfficiencyErrorUp(i),
                        eff_mc.GetEfficiencyErrorLow(i))
            sf = eps_d / eps_m
            rel2 = 0.0
            if eps_d > 0: rel2 += (err_d / eps_d)**2
            if eps_m > 0: rel2 += (err_m / eps_m)**2
            if rel2 > 0: err_sf = sf * rel2**0.5

        h_sf.SetBinContent(i, sf)
        h_sf.SetBinError(i, err_sf)

        low  = h_sf.GetXaxis().GetBinLowEdge(i)
        high = h_sf.GetXaxis().GetBinUpEdge(i)
        print(f"  {low:.0f}–{high:.0f} GeV : {sf:.3f} ± {err_sf:.3f}")

    # ---------------- plotting ----------------
    base = os.path.basename(tagged_file).replace(".root", "")
    sample_title = tagged_file.split("/")[-4]

    pdf_data = f"eff_output/eff_data_{base}.pdf"
    pdf_mc   = f"eff_output/eff_mc_{base}.pdf"
    pdf_sf   = f"eff_output/sf_{base}.pdf"

    plot(eff_data, h_all_rb, sample_title + " Data", "Efficiency",
         pdf_data, 0.0, 1.1, 200.0, 800.0, "P SAME")
    plot(eff_mc,   h_all_rb, sample_title + " MC", "Efficiency",
         pdf_mc,   0.0, 1.1, 200.0, 800.0, "P SAME")
    plot(h_sf,     h_all_rb, sample_title + " SF", "Scale factor (Data/MC)",
         pdf_sf,   0.0, 2.0, 200.0, 800.0, "E1 SAME")

    # ---------------- save to ROOT file ----------------
    out_name = f"eff_output/eff_{base}.root"
    out = ROOT.TFile(out_name, "RECREATE")
    h_all_rb.Write("h_all")
    h_tag_rb.Write("h_tag")
    h_mc_all_rb.Write("h_all_mc")
    h_mc_tag_rb.Write("h_tag_mc")
    eff_data.Write("h_eff_data")
    eff_mc.Write("h_eff_mc")
    h_sf.Write("h_sf")
    out.Close()


if __name__ == "__main__":
    main()
