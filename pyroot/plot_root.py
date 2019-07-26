import ROOT

ROOT.gStyle.SetOptStat(0)

fv = ROOT.TFile("vacuum.root")
tv = fv.Get("artie")
#tv.Print()

fa = ROOT.TFile("argon.root")
ta = fa.Get("artie")
ta.Print()

c1 = ROOT.TCanvas( 'c1', '', 200, 10, 700, 500 )

hev = ROOT.TH1F("hev","",100,40,80)
hea = ROOT.TH1F("hea","",100,40,80)

tv.Draw("1000*gen_energy>>hev", "arrival_time>0")
ta.Draw("1000*gen_energy>>hea", "arrival_time>0")

hev.SetMinimum(0)
hev.Draw("ep")
hea.Draw("epSAME")
hev.SetLineColor(2)
hea.SetLineColor(4)
hev.SetXTitle("Neutron Energy (keV)")

l = ROOT.TLegend(0.8,0.5,0.9,0.6)
l.SetBorderSize(0)
l.AddEntry(hea,"argon","ep")
l.AddEntry(hev,"vacuum","ep")
l.Draw()

c1.Modified()
c1.Update()


input("Press Enter to continue...")
