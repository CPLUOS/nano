import ROOT, os

filedir = "./"
#filedir = "/cms/ldap_home/tt8888tt/CMSSW_10_0_0_pre2/src/nano/nanoAOD/prod/"

f = ROOT.TFile(filedir+"nanoAOD.root")
tree = f.Get("Events")

cgen = 0; creco = 0;
for iev, t in enumerate(tree):
    tree.GetEntry(iev)
    for j in range(t.ngenHadron):
        if abs(t.genHadron_pdgId[j]) != 3122: continue
        #if not t.genHadron_isGenHadFromTop[j]: continue
        if abs(t.genHadron_isGenHadFromTsb[j]) != 3: continue
        #if abs(t.genHadron_isGenHadFromTsb[j]) != 3 and abs(t.genHadron_isGenHadFromTsb[j]) != 5: continue
        cgen += 1

    for i in range(t.nhad):
        if abs(t.had_pdgId[i]) != 3122: continue
        if t.hadTruth_nMatched[i] != 2: continue
        if abs(t.hadTruth_isHadFromTsb[i]) != 3: continue
        #if not t.hadTruth_isHadFromTop[i]: continue
        creco+=1
    #print cgen, creco

print float(creco)/cgen*100


cks = 0; cs = 0;
for i, t in enumerate(f.Events):
  t.GetEntry(i)
  for j in range(t.ngenHadron):
    if t.genHadron_pdgId[j]!=310: continue
    if t.genHadron_isGenHadFromTsb[j]==5: cks+=1; break;
  for j in range(t.ngenHadron):
    if t.genHadron_pdgId[j]!=310: continue
    if t.genHadron_isGenHadFromTsb[j]==-5: cks+=1; break;
  for g in range(t.nGenPart):
    if t.GenPart_status[g] != 23: continue
    if abs(t.GenPart_pdgId[g]) != 5: continue
    cs +=1
  if i == 100000: break
print cks, cs, cks/float(cs)
