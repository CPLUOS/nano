#./DYdraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x-name> -y <y_name> -f <f_name> -d <dolog> 

std_cut     = 'Dilep.M()>60&&Step>6'

BDT_cut     = ''

weight      = 'weight*mueffweight*btagweight'

json_used   = 'Golden'

x_name_l    = ["Invariant Mass [GeV]", "Transverse Momentum [GeV]"]

plotvar_l   = ["dilep.M\(\)", "dilep.Pt\(\)"]

binset_l    = ["[300,0,300]", "[250,0,300]", "[200,0,300]", "[150,0,300]"]

MCat_l      = ["VBFT_M", "GGFT_M", "VBFL_M", "JetT_M", "JetL_M"]

PTCat_l     = ["VBFT_Pt", "GGFT_Pt", "VBFL_Pt", "JetT_Pt", "JetL_Pt"] 

lst         = []

User_input = raw_input(" How would you like it listed? valid, FL, FH, SL : ") 
if User_input == "valid":

    std_cut = 'Dilep.M()>60&&charge==1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [75,50,200] -p 'Dilep.M()' -x 'm(\mu^+\mu^-)(GeV)' -y 'Event' -f 'Dilept_M' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    std_cut = 'Dilep.M()>12&&charge==1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [60,0,300] -p 'Dilep.Pt()' -x 'p_T(\mu^+\mu^-)(GeV)' -y 'Event' -f 'Dilept_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu1.Pt()' -x '\[p_T(\mu_1)(GeV)\]' -y 'Event' -f 'lept1_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu1.Eta()' -x '\[\eta(\mu_1)\]' -y 'Event' -f 'lept1_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu2.Pt()' -x '\[p_t(\mu_2)\](GeV)' -y 'Event' -f 'lept2_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu2.Eta()' -x '\[\eta(\mu_2)\]' -y 'Event' -f 'lept2_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'njet' -x 'nJets' -y 'Event' -f 'NJets' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nbjet' -x 'nBJets' -y 'Event' -f 'NBMed' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nelec' -x 'nElec' -y 'Event' -f 'nElec' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nmuon' -x 'nMuon' -y 'Event' -f 'nMu' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nlep' -x 'nLep' -y 'Event' -f 'nLep' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,300] -p 'Met' -x '\[E_T^miss\](GeV)' -y 'Event' -f 'MET' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0.8484,1.05] -p 'CSVv2' -x 'CSVv2-Medium' -y 'Event' -f 'CSVv2' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-3.2,3.2] -p 'Met_phi' -x 'MET_phi' -y 'Event' -f 'MET_phi' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   

if User_input == "FL":

    std_cut = 'Dilep.M()>60&&FL==1&&njet>=1&&nbjet>=1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [75,50,200] -p 'Dilep.M()' -x 'm(\mu^+\mu^-)(GeV)' -y 'Event' -f 'FL_Dilept_M' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    std_cut = 'Dilep.M()>12&&FL==1&&njet>=1&&nbjet>=1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [60,0,300] -p 'Dilep.Pt()' -x 'p_T(\mu^+\mu^-)(GeV)' -y 'Event' -f 'FL_Dilept_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu1.Pt()' -x '\[p_T(\mu_1)(GeV)\]' -y 'Event' -f 'FL_lept1_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu1.Eta()' -x '\[\eta(\mu_1)\]' -y 'Event' -f 'FL_lept1_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu2.Pt()' -x '\[p_t(\mu_2)\](GeV)' -y 'Event' -f 'FL_lept2_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu2.Eta()' -x '\[\eta(\mu_2)\]' -y 'Event' -f 'FL_lept2_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'njet' -x 'nJets' -y 'Event' -f 'FL_NJets' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nbjet' -x 'nBJets' -y 'Event' -f 'FL_NBMed' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nelec' -x 'nElec' -y 'Event' -f 'FL_nElec' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nmuon' -x 'nMuon' -y 'Event' -f 'FL_nExtraMu' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nlep' -x 'nLep' -y 'Event' -f 'FL_nLep' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,300] -p 'Met' -x '\[E_T^miss\](GeV)' -y 'Event' -f 'FL_MET' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0.8484,1.05] -p 'CSVv2' -x 'CSVv2-Medium' -y 'Event' -f 'FL_CSVv2' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-3.2,3.2] -p 'Met_phi' -x 'MET_phi' -y 'Event' -f 'FL_MET_phi' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   

if User_input == "FH":
   for i in [2,3,4]: 
    std_cut = 'Dilep.M()>60&&FH%s==1'%i
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [75,50,200] -p 'Dilep.M()' -x 'm(\mu^+\mu^-)(GeV)' -y 'Event' -f 'FH%s_Dilept_M' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    std_cut = 'Dilep.M()>12&&FH%s==1'%i
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [75,0,300] -p 'Dilep.Pt()' -x 'p_T(\mu^+\mu^-)(GeV)' -y 'Event' -f 'FH%s_Dilept_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu1.Pt()' -x '\[p_T(\mu_1)(GeV)\]' -y 'Event' -f 'FH%s_lept1_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu1.Eta()' -x '\[\eta(\mu_1)\]' -y 'Event' -f 'FH%s_lept1_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu2.Pt()' -x '\[p_t(\mu_2)\](GeV)' -y 'Event' -f 'FH%s_lept2_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu2.Eta()' -x '\[\eta(\mu_2)\]' -y 'Event' -f 'FH%s_lept2_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'njet' -x 'nJets' -y 'Event' -f 'FH%s_NJets' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nbjet' -x 'nBJets' -y 'Event' -f 'FH%s_NBMJet' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nelec' -x 'nElec' -y 'Event' -f 'FH%s_nElec' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nmuon' -x 'nMuon' -y 'Event' -f 'FH%s_nExtraMu' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nlep' -x 'nLep' -y 'Event' -f 'FH%s_nLep' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,300] -p 'Met' -x '\[E_T^miss\](GeV)' -y 'Event' -f 'FH%s_MET' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0.8484,1.05] -p 'CSVv2' -x 'CSVv2-Medium' -y 'Event' -f 'FH%s_CSVv2' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-3.2,3.2] -p 'Met_phi' -x 'MET_phi' -y 'Event' -f 'FH%s_MET_phi' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight, i)
    lst.append(cmd)   

if User_input == "SL":

    std_cut = 'Dilep.M()>60&&SL==1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [75,50,200] -p 'Dilep.M()' -x 'm(\mu^+\mu^-)(GeV)' -y 'Event' -f 'SL_Dilept_M' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    std_cut = 'Dilep.M()>12&&SL==1'
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [60,0,300] -p 'Dilep.Pt()' -x 'p_T(\mu^+\mu^-)(GeV)' -y 'Event' -f 'SL_Dilept_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu1.Pt()' -x '\[p_T(\mu_1)(GeV)\]' -y 'Event' -f 'SL_lept1_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu1.Eta()' -x '\[\eta(\mu_1)\]' -y 'Event' -f 'SL_lept1_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'Mu2.Pt()' -x '\[p_t(\mu_2)\](GeV)' -y 'Event' -f 'SL_lept2_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Mu2.Eta()' -x '\[\eta(\mu_2)\]' -y 'Event' -f 'SL_lept2_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'njet' -x 'nJets' -y 'Event' -f 'SL_NJets' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nbjet' -x 'nBJets' -y 'Event' -f 'SL_NBMed' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nelec' -x 'nElec' -y 'Event' -f 'SL_nElec' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nmuon' -x 'nMuon' -y 'Event' -f 'SL_nExtraMu' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'nlep' -x 'nLep' -y 'Event' -f 'SL_nLep' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,300] -p 'Met' -x '\[E_T^miss\](GeV)' -y 'Event' -f 'SL_MET' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0.8484,1.05] -p 'CSVv2' -x 'CSVv2-Medium' -y 'Event' -f 'SL_CSVv2' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
    cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-3.2,3.2] -p 'Met_phi' -x 'MET_phi' -y 'Event' -f 'SL_MET_phi' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
    lst.append(cmd)   
for l in lst:
    print l


  #  cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'jet1.Pt()' -x '\[p_t(j\_1)\](GeV)' -y 'Event' -f 'jet1_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
  #  lst.append(cmd)   
  #  cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-5,5] -p 'jet1.Eta()' -x '\[\eta(j\_1)\]' -y 'Event' -f 'jet1a_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
  #  lst.append(cmd)   
  #  cmd ="./tthDraw.py -c '%s' -w '%s' -b [45,0,200] -p 'jet2.Pt()' -x '\[p_t(j\_2)\](GeV)' -y 'Event' -f 'jet2_Pt' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
  #  lst.append(cmd)   
  #  cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,-5,5] -p 'jet2.Eta()' -x '\[\eta(j\_2)\]' -y 'Event' -f 'jet2_eta' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
  #  lst.append(cmd)   
  # if User_input =="BDT":
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [150,50,200] -p 'Dilep.M()' -x 'm(\mu^+\mu^-)(GeV)' -y 'Event' -f 'Dilept_M' -j 'Golden' -d 'True' > /dev/null &"%(std_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [60,0,300] -p 'Dilep.Pt()' -x 'p_t(\mu^+\mu^-)(GeV)' -y 'Event' -f 'Dilept_Pt' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #   #  cmd ="./tthDraw.py -c '%s' -w '%s' -b [52,0,150] -p 'met' -x '\[E_T^miss\](GeV)' -y 'Event' -f 'MET' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #   #  lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,-2.5,2.5] -p 'Dilep.Eta()' -x '\eta(\mu^+\mu^-)' -y 'Event' -f 'Dilept_Eta' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [12,0,5] -p 'Muetadiff' -x '|\Delta\eta(\mu^+\mu^-)|' -y 'Event' -f 'Mu_DeltaEta' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [25,0,3.2] -p 'Muphidiff' -x '|\Delta\phi(\mu^+\mu^-)|' -y 'Event' -f 'Mu_DeltaPhi' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'central_jets' -x 'nJetsCent' -y 'Event' -f 'nJets_Cent' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,0,2000] -p 'DijetM1' -x '\[m(jj\_1)\]' -y 'Event' -f 'Dijet_M1' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,10] -p 'DijetEta1' -x '\[|\Delta\eta(jj\_1)|\]' -y 'Event' -f 'Dijet_Eta1' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [11,0,11] -p 'forward_jets' -x 'nJetsFwd' -y 'Event' -f 'nJets_Fwd' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [50,0,2000] -p 'DijetM2' -x '\[m(jj\_2)\]' -y 'Event' -f 'Dijet_M2' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #     lst.append(cmd)   
  #     cmd ="./tthDraw.py -c '%s' -w '%s' -b [30,0,10] -p 'DijetEta2' -x '\[|\Delta\eta(jj\_2)|\]' -y 'Event' -f 'Dijet_Eta2' -j 'Golden' -d 'True' > /dev/null &"%(BDT_cut, weight)
  #  
  #  lst.append(cmd)   
