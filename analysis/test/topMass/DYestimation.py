import json, os
from nano.analysis.histoHelper import *


def getK(rootfileDir, tname, cut, rdfilelist):
	dycut = "(step1==1)*"
	rfname = rootfileDir + rdfilelist[2-1] +".root"
	rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', [2,0,2], 'tri', dycut+'(%s && channel==2 && step2 ==0)'%(cut))
	rfname = rootfileDir + rdfilelist[3-1] +".root"
	rd_mm_in = makeTH1(rfname, tname,'rd_mm_in', [2,0,2], 'tri', dycut+'(%s && channel==3 && step2 ==0)'%(cut))

	kMM = math.sqrt(rd_mm_in.Integral()/rd_ee_in.Integral())/2.
	kEE = math.sqrt(rd_ee_in.Integral()/rd_mm_in.Integral())/2.
	print 'kMM', kMM
	print 'kEE', kEE

	return kMM, kEE

def	getDYFactor(rfname, front_rfname, scale, front_scale, rootfileDir, tname, cut, weight, step, kMM, kEE, rdfilelist):
	dyest = []
	print 'estimation step... ', step
	dycut = ""
	if step == 1: dycut = "(step1==1)*"
	if step == 2: dycut = "(step1==1)*"
	if step == 3: dycut = "(step1==1)*(step3==1)*"
	if step == 4: dycut = "(step1==1)*(step3==1)*(step4==1)*"
	if step >= 5: dycut = "(step1==1)*(step3==1)*(step4==1)*(step5==1)*"

	mc_ee_in = makeTH1(rfname,tname,"mc_ee_in", [2,0,2], 'tri', dycut+'(%s && channel==2 && step2==0)*(%s)'%(cut,weight), scale)
	mc_mm_in = makeTH1(rfname,tname,"mc_mm_in", [2,0,2], 'tri', dycut+'(%s && channel==3 && step2==0)*(%s)'%(cut,weight), scale)
	mc_ee_out = makeTH1(rfname,tname,"mc_ee_in", [2,0,2], 'tri', dycut+'(%s && channel==2 && step2==1)*(%s)'%(cut,weight), scale)
	mc_mm_out = makeTH1(rfname,tname,"mc_mm_in", [2,0,2], 'tri', dycut+'(%s && channel==3 && step2==1)*(%s)'%(cut,weight), scale)

	mc_ee_in.Add(makeTH1(front_rfname,tname,"mc_ee_in", [2,0,2], 'tri', dycut+'(%s && channel==2 && step2==0)*(%s)'%(cut,weight), front_scale))
	mc_mm_in.Add(makeTH1(front_rfname,tname,"mc_mm_in", [2,0,2], 'tri', dycut+'(%s && channel==3 && step2==0)*(%s)'%(cut,weight), front_scale))
	mc_ee_out.Add(makeTH1(front_rfname,tname,"mc_ee_in", [2,0,2], 'tri', dycut+'(%s && channel==2 && step2==1)*(%s)'%(cut,weight), front_scale))
	mc_mm_out.Add(makeTH1(front_rfname,tname,"mc_mm_in", [2,0,2], 'tri', dycut+'(%s && channel==3 && step2==1)*(%s)'%(cut,weight), front_scale))
	
	rfname = rootfileDir+rdfilelist[1-1]+".root"
	rd_em_in = makeTH1(rfname, tname,'rd_em_in', [2,0,2], 'tri', dycut+'(%s && channel==1 && ((dilep.M() > 76) && (dilep.M() < 106)))'%(cut))
	rfname = rootfileDir + rdfilelist[2-1] +".root"
	rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', [2,0,2], 'tri', dycut+'(%s && channel==2 && step2 ==0)'%(cut))
	rfname = rootfileDir + rdfilelist[3-1] +".root"
	rd_mm_in = makeTH1(rfname, tname,'rd_mm_in', [2,0,2], 'tri', dycut+'(%s && channel==3 && step2 ==0)'%(cut))

	dyest = drellYanEstimation(mc_ee_in.Integral(), mc_ee_out.Integral(), mc_mm_in.Integral(), mc_mm_out.Integral(), rd_ee_in.Integral(), rd_mm_in.Integral(), rd_em_in.Integral(), kMM, kEE)
	return dyest

def printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist):
	dyratio = [[0 for x in range(7)] for x in range(4)]
	kMM, kEE = getK(rootfileDir, tname, cut, rdfilelist)

	rfname = rootfileDir + 'DYJets' +".root"
	data = findDataSet('DYJets', datasets)
	scale = datalumi*data["xsec"]
	tfile = ROOT.TFile(rfname)
	wentries = tfile.Get("nevents").Integral()
	scale = scale/wentries

	front_rfname = rootfileDir + 'DYJets_10to50' +".root"
	front_data = findDataSet('DYJets_10to50', datasets)
	front_scale = datalumi*front_data["xsec"]
	front_wentries = getWeightedEntries(front_rfname, tname, "tri", weight)
	front_scale = scale/wentries

	for step in range(1,7):
		dyratio[1][step] = 1.
		dyest = getDYFactor(rfname, front_rfname, scale, front_scale, rootfileDir, tname, cut, weight, step, kMM, kEE, rdfilelist)

		print "DY estimation for s", step, "ee =",dyest[0], "mm =",dyest[1]   
		dyratio[2][step] = dyest[0]
		dyratio[3][step] = dyest[1]
	f=open('DYFactor.json','w')
	json.dump(dyratio, f)
	f.close()

