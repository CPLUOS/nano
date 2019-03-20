import ROOT, math
hlist = { "bfc": {"tt":dict(), "hh":dict()}, "afc" : {"tt":dict(), "hh":dict()} }
limit_list = { 'll_pt': { "bfc":[60,0,400], "afc":[60,0,400] },
               'bb_pt': { "bfc":[60,0,500], "afc":[60,0,500] }, 
               'missing_et': { "bfc":[60,0,400], "afc":[60,0,400]}, 
               'missing_et_phi': { "bfc":[60,-4,4], "afc":[60,-4,4] }, 
               'bbll_mass': { "bfc":[100,0,700], "afc":[100,0,700] }, 
               "bb_mass": { "bfc":[100,0,400], "afc":[100,50,400] }, 
               "ll_mass": { "bfc":[100,0,400], "afc":[100,0,150] },
               "basic_MT2_332_bbll" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_blbl" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_b" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_l" : { "bfc":[80,0,500], "afc":[80,0,200] },
               "bb_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "ll_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bbll_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bl_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bl_min_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "mT" : { "bfc":[75,0,300], "afc":[75,0,300] },
               "tmva_bdtg_output" : { "bfc":[50,-1.2,1.2], "afc":[50,-1.2,1.2] },
               "hig_top_2d" : { "bfc":[100,-10,10], "afc":[100,-10,10]}
              }
def getHistogram(title,binning,lower_limit,upper_limit):
    return ROOT.TH1D(title,title,binning,lower_limit,upper_limit)

for key in ["bfc","afc"]:
    for key2 in ["hh","tt"]:
        for key3 in limit_list.keys():
            if key3 is "hig_top_2d":
                hlist[key][key2][key3] = ROOT.TH2F(key2+" "+key+" "+key3.replace("_"," "),key2+" "+key+" "+key3.replace("_"," "),limit_list[key3][key][0],limit_list[key3][key][1],limit_list[key3][key][2],limit_list[key3][key][0],limit_list[key3][key][1],limit_list[key3][key][2])
            else:
                hlist[key][key2][key3] = getHistogram(key2+" "+key+" "+key3.replace("_"," "),limit_list[key3][key][0],limit_list[key3][key][1],limit_list[key3][key][2])



tt_c = ROOT.TChain("events")
tt_c.Add("TT.root")
hh_c = ROOT.TChain("events")
hh_c.Add("HH_SM.root")
sample_list = {"tt":tt_c,"hh":hh_c}

f = ROOT.TFile("plots_1.root","RECREATE")

#for sample in ["tt", "hh", "hh_B6", "hh_B11"]:
for sample in ["tt", "hh"]:
    c = sample_list[sample]
    for ie, e in enumerate(c):
        cut = ""
        if e.step == 4 and e.tmva_bdtg_output > 0:
            cut = "afc"
        else: cut = "bfc"
        hlist[cut][sample]['hig_top_2d'].Fill(e.higgsness, e.topness)
        hlist[cut][sample]['ll_pt'].Fill(e.ll.Pt())
        hlist[cut][sample]['bb_pt'].Fill(e.bb.Pt())
        hlist[cut][sample]["missing_et_phi"].Fill(e.MET.Phi())
        hlist[cut][sample]["missing_et"].Fill(e.MET.E())
        hlist[cut][sample]["bbll_mass"].Fill(e.bbll.M())
        hlist[cut][sample]["bb_mass"].Fill(e.bb.M())
        hlist[cut][sample]["ll_mass"].Fill(e.ll.M())
        hlist[cut][sample]["basic_MT2_332_bbll"].Fill(e.basic_MT2_332_bbll)
        hlist[cut][sample]["basic_MT2_332_blbl"].Fill(e.basic_MT2_332_blbl)
        hlist[cut][sample]["basic_MT2_332_b"].Fill(e.basic_MT2_332_b)
        hlist[cut][sample]["basic_MT2_332_l"].Fill(e.basic_MT2_332_l)
        hlist[cut][sample]["bb_deltaR"].Fill(e.bb_deltaR)
        hlist[cut][sample]["ll_deltaR"].Fill(e.ll_deltaR)
        hlist[cut][sample]["bbll_deltaR"].Fill(e.bbll_deltaR)
        hlist[cut][sample]["bl_min_deltaR"].Fill(e.bl_min_deltaR)
        hlist[cut][sample]["tmva_bdtg_output"].Fill(e.tmva_bdtg_output)
        hlist[cut][sample]["mT"].Fill(e.mT)
        for bld in e.bl_deltaR:
            hlist[cut][sample]["bl_deltaR"].Fill(bld)

num_bg = hlist["afc"]["tt"]["hig_top_2d"].GetEntries()
num_bg = num_bg * 322627 / 14529280
num_sg = hlist["afc"]["hh"]["hig_top_2d"].GetEntries()
num_sg = num_sg * 4. / 299999
sensitivity = num_sg/math.sqrt(num_bg+num_sg)
print "number of background : " + str(num_bg) + "\n"
print "number of signal : " + str(num_sg) + "\n"
print "sensitivity : " + str(sensitivity) + "\n"

for key in ["bfc","afc"]:
    for key2 in ["tt", "hh"]:
    #for key2 in ["tt", "hh", "hh_B6", "hh_B11"]:
        for key3 in hlist[key][key2].keys():
            #integ = hlist[key][key2][key3].Integral()
            #if integ == 0:
            #    print "%s %s %s has 0 integral" % (key, key2, key3)
            #else:
            #    hlist[key][key2][key3].Scale(1/integ)
            hlist[key][key2][key3].GetXaxis().SetTitle("log H")
            hlist[key][key2][key3].GetYaxis().SetTitle("log T")
            hlist[key][key2][key3].Write()

f.Close()
