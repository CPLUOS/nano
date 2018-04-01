import ROOT

#fdir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/"
fdir = "/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180125_131546/0000/"

h = ROOT.TH1D("dr", "dr", 100, 0, 3)
h2 = ROOT.TH1D("mass sig", "mass sig", 100, 1.7, 2)
h3 = ROOT.TH1D("mass bkg", "mass bkg", 100, 1.7, 2)

h_genPt = ROOT.TH1D("gen pt", "gen pt", 100, 0, 100)
h_recoPt = ROOT.TH1D("reco pt", "reco pt", 100, 0, 100)
h_genEta = ROOT.TH1D("gen eta", "gen eta", 100, -3, 3)
h_recoEta = ROOT.TH1D("reco eta", "reco eta", 100, -3, 3)

#hMP = ROOT.TH1D("genmotherIdx", "genmotherIdx", 100, 0, 100)

for ifile in range(1, 2):
    f = ROOT.TFile(fdir+"nanoAOD_%d.root"%ifile)
    print ifile
    for iev, e in enumerate(f.Events):
        if iev%100 == 0 : print "iev: ", iev
        #print "iev: ", iev
        for i in xrange(e.ncmeson):
            if e.cmeson_pdgId[i] != 421: continue
            if e.cmeson_pt[i] < 5: continue
            if abs(e.cmeson_eta[i]) > 2.4: continue
            h3.Fill(e.cmeson_mass[i])

        for i in xrange(e.nGenPart):
            if abs(e.GenPart_pdgId[i]) != 421: continue
            if e.GenPart_pt[i] < 5: continue
            if abs(e.GenPart_eta[i]) > 2.4: continue
            if abs(e.GenPart_status[i]) == 1: continue
            tlv_gen  = ROOT.TLorentzVector()
            tlv_gen.SetPtEtaPhiM(e.GenPart_pt[i], e.GenPart_eta[i], e.GenPart_phi[i], 1.86484)
            dr = 999
            gim = e.GenPart_genPartIdxMother[i]
            
            for m in xrange(e.nGenPart):
                if gim == -1: break
                if abs(e.GenPart_pdgId[gim]) == 6:
                    for j in xrange(e.ncmeson):
                        if e.cmeson_pdgId[j] != 421: continue
                        if e.cmeson_pt[j] < 5: continue
                        if abs(e.cmeson_eta[j]) > 2.4: continue
                        
                        tlv_cme = ROOT.TLorentzVector()
                        tlv_cme.SetPtEtaPhiM(e.cmeson_pt[j], e.cmeson_eta[j], e.cmeson_phi[j], e.cmeson_mass[j])
                        if tlv_gen.DeltaR(tlv_cme) < dr:
                            mindrIdx = j
                            dr = tlv_gen.DeltaR(tlv_cme)
                    if dr == 999: continue
                    h.Fill(dr)
                    if dr < 0.1:
                        h2.Fill(e.cmeson_mass[j])
                        
                        h_recoPt.Fill(e.cmeson_pt[mindrIdx])
                        h_recoEta.Fill(e.cmeson_eta[mindrIdx])
                    h_genPt.Fill(e.GenPart_pt[i])
                    h_genEta.Fill(e.GenPart_eta[i])
                    break
                if abs(e.GenPart_pdgId[gim]) != 6:
                    ppd = e.GenPart_genPartIdxMother[gim]
                    gim = ppd
        
c = ROOT.TCanvas()
h.Draw()
c.Print("dr.png")

#hMP.Draw()
#c.Print("genmotherId.png")

h3.Draw()
h2.Draw("same")
h2.SetFillColor(2)
h3.SetFillColor(4)
h2.SetLineColor(2)
h3.SetLineColor(4)
h3.SetMinimum(0)
c.Print("m.png")

h_genPt.Draw()
h_recoPt.Draw("same")
h_genPt.SetLineColor(2)
c.Print("pt.png")

h_genEta.Draw()
h_recoEta.Draw("same")
h_genEta.SetLineColor(2)
c.Print("eta.png")
