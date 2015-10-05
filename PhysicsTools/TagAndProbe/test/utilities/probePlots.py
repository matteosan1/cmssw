import ROOT

ROOT.gROOT.SetBatch(True)

vars = ["mass", "probe_Pho_sieie","probe_Pho_sieip","probe_Pho_s4", "probe_Pho_r9", "probe_Pho_hoe",
        "probe_Pho_phoIso", "probe_Pho_chIso"]
f1 = ROOT.TFile("flashgg/TnP_DY_flashgg_v3.root")
t1 = f1.Get("PhotonToRECO/fitter_tree")
h1 = []
t1.Draw(vars[0]+">>"+vars[0]+"(60, 60, 120)", "(passingID==1 && (probe_sc_abseta>1.5&&tag_sc_abseta>1.5))*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[0]+""))
t1.Draw(vars[1]+">>"+vars[1]+"(100, 0, 0.05)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[1]+""))
t1.Draw(vars[2]+">>"+vars[2]+"(100, -0.0005, 0.0005)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[2]+""))
t1.Draw(vars[3]+">>"+vars[3]+"(100, 0, 1.)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[3]+""))
t1.Draw(vars[4]+">>"+vars[4]+"(100, 0, 1.)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[4]+""))
t1.Draw(vars[5]+">>"+vars[5]+"(100, 0, 0.5)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[5]+""))
t1.Draw(vars[6]+">>"+vars[6]+"(100, 0, 5)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[6]+""))
t1.Draw(vars[7]+">>"+vars[7]+"(100, 0, 1)", "(passingID==1 && probe_sc_abseta>1.5)*totWeight", "GOFF")
h1.append(ROOT.gDirectory.Get(vars[7]+""))

f2 = ROOT.TFile("flashgg/TnP_Run2015C_27Sep.root")
t2 = f2.Get("PhotonToRECO/fitter_tree")
h2 = []
t2.Draw(vars[0]+">>"+vars[0]+"Data(60, 60, 120)", "passingID==1 && (probe_sc_abseta>1.5&&tag_sc_abseta>1.5)", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[0]+"Data"))
t2.Draw(vars[1]+">>"+vars[1]+"Data(100, 0, 0.05)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[1]+"Data"))
t2.Draw(vars[2]+">>"+vars[2]+"Data(100, -0.0005, 0.0005)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[2]+"Data"))
t2.Draw(vars[3]+">>"+vars[3]+"Data(100, 0, 1.)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[3]+"Data"))
t2.Draw(vars[4]+">>"+vars[4]+"Data(100, 0, 1.)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[4]+"Data"))
t2.Draw(vars[5]+">>"+vars[5]+"Data(100, 0, 0.5)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[5]+"Data"))
t2.Draw(vars[6]+">>"+vars[6]+"Data(100, 0, 5)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[6]+"Data"))
t2.Draw(vars[7]+">>"+vars[7]+"Data(100, 0, 1)", "passingID==1 && probe_sc_abseta>1.5", "GOFF")
h2.append(ROOT.gDirectory.Get(vars[7]+"Data"))

c = []
out = ROOT.TFile("out_ee.root", "recreate")
for i in xrange(len(h1)):
    c.append(ROOT.TCanvas("c"+str(i), ""))
    h2[i].Draw("PE")
    h1[i].Draw("HISTSAME")
    h1[i].SetFillColor(ROOT.kRed)
    h2[i].Draw("PESAME")
    h2[i].SetMarkerStyle(20)
    h1[i].Scale(h2[i].Integral()/h1[i].Integral())
    c[-1].Write()

out.Close()
