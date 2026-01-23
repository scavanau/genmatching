#!/usr/bin/env python3
import ROOT
import os
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)

# -------------------------
# output directory
# -------------------------
outdir = "outputs"
os.makedirs(outdir, exist_ok=True)

# -------------------------
# input file
# -------------------------
f = ROOT.TFile.Open(
    "/eos/user/s/scavanau/Mono-V/efficiency/0/Run3Summer22/"
    "NanoAODv12/singleMuon_CR_MonoV/plots/root/monov_FatJet_pt.root"
)

# -------------------------
# load histograms
# -------------------------
h_data = f.Get("data_obs")

h_dy   = f.Get("QCD_Zll")
h_w    = f.Get("QCD_W")
h_ewkz = f.Get("EWK_Zll")
h_ewkw = f.Get("EWK_W")
h_qcd  = f.Get("QCDJets")

h_top = f.Get("Top")
h_ww  = f.Get("WW")
h_wz  = f.Get("WZ")
h_zz  = f.Get("ZZ")

# -------------------------
# rebinning definition
# -------------------------
edges = array("d", [200, 300, 400, 500, 600, 700, 800, 1000])
nbins = len(edges) - 1

h_data_rb = h_data.Rebin(nbins, "h_data_rb", edges)
h_dy_rb   = h_dy.Rebin(nbins,   "h_dy_rb", edges)
h_w_rb    = h_w.Rebin(nbins,    "h_w_rb", edges)
h_ewkz_rb = h_ewkz.Rebin(nbins, "h_ewkz_rb", edges)
h_ewkw_rb = h_ewkw.Rebin(nbins, "h_ewkw_rb", edges)
h_qcd_rb  = h_qcd.Rebin(nbins,  "h_qcd_rb", edges)

h_top_rb = h_top.Rebin(nbins, "h_top_rb", edges)
h_ww_rb  = h_ww.Rebin(nbins,  "h_ww_rb", edges)
h_wz_rb  = h_wz.Rebin(nbins,  "h_wz_rb", edges)
h_zz_rb  = h_zz.Rebin(nbins,  "h_zz_rb", edges)

# detach from file
for h in [
    h_data_rb, h_dy_rb, h_w_rb, h_ewkz_rb, h_ewkw_rb, h_qcd_rb,
    h_top_rb, h_ww_rb, h_wz_rb, h_zz_rb
]:
    h.SetDirectory(0)

# -------------------------
# data minus non-V MC
# -------------------------
h_data_est = h_data_rb.Clone("h_data_minus_nonV")
h_data_est.SetDirectory(0)

h_data_est.Add(h_dy_rb,   -1)
h_data_est.Add(h_w_rb,    -1)
h_data_est.Add(h_ewkz_rb, -1)
h_data_est.Add(h_ewkw_rb, -1)
h_data_est.Add(h_qcd_rb,  -1)

# -------------------------
# sum of remaining MC
# -------------------------
h_other_mc = h_top_rb.Clone("h_other_mc")
h_other_mc.SetDirectory(0)

h_other_mc.Add(h_ww_rb)
h_other_mc.Add(h_wz_rb)
h_other_mc.Add(h_zz_rb)

# -------------------------
# plot: data − non-V MC
# -------------------------
c1 = ROOT.TCanvas("c1", "", 600, 600)
c1.SetGrid()

h_data_est.SetTitle(
    "singleMuon CR: data - QCD Zll - QCD W - EWK Zll - EWK W - QCD Jets"
)
h_data_est.GetXaxis().SetTitle("FatJet p_{T} [GeV]")
h_data_est.GetYaxis().SetTitle("Events")
h_data_est.SetMarkerStyle(20)
h_data_est.SetLineWidth(2)

h_data_est.Draw("E1")
c1.SaveAs(os.path.join(outdir, "data_minus_nonV.pdf"))

# -------------------------
# plot: other MC sum
# -------------------------
c2 = ROOT.TCanvas("c2", "", 600, 600)
c2.SetGrid()

h_other_mc.SetTitle("singleMuon CR: Top + WW + WZ + ZZ")
h_other_mc.GetXaxis().SetTitle("FatJet p_{T} [GeV]")
h_other_mc.GetYaxis().SetTitle("Events")
h_other_mc.SetMarkerStyle(20)
h_other_mc.SetLineWidth(2)

h_other_mc.Draw("E1")
c2.SaveAs(os.path.join(outdir, "other_MC_sum.pdf"))

# -------------------------
# save ROOT output
# -------------------------
out = ROOT.TFile(os.path.join(outdir, "histograms.root"), "RECREATE")
h_data_est.Write()
h_other_mc.Write()
out.Close()

print("Done. Outputs written to:", outdir)
