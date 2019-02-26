Dilep_m=("[60,20,320]" "dilep.M()" "Invariant_Mass_(ll)_[GeV]")
Dilep_pt=("[80,20,420]" "dilep.Pt()" "Transverse_Momentum_(ll)_[GeV]")
Dilep_eta=("[80,-4,4]" "dilep.Eta()" "Eta_(ll)")
Dilep_phi=("[100,-5,5]" "dilep.Phi()" "Phi_(ll)")
MET_pt=("[60,20,320]" "MET_pt" "Missing_Energy_[GeV]")
MET_sumEt=("[80,200,5000]" "MET_sumEt" "Summed_Transverse_Energy_[GeV]")
nJet=("[30,0,30]" "nJet" "Jet_Multiplicity")

had_d=("[30,0,30]" "d" "had_d")
had_pt=("[60,0,100]" "pt" "Transverse_Momentum_(had)_[GeV]")
had_eta=("[80,-4,4]" "eta" "Eta_(had)")
had_phi=("[80,-4,4]" "phi" "Phi_(had)")
had_m=("[60,0.25,0.85]" "mass" "Invariant_Mass_(had)_[GeV]")
had_angleXY=("[100,-1,1]" "angleXY" "AngleXY_(had)")
had_chi2=("[60,0,10]" "chi2" "Chi2_(had)")
had_dca=("[60,0,2]" "dca" "DCA_(had)")

jet_cmult=("[60,0,60]" "cmult" "Charged_Multiplicity")
jet_nmult=("[40,0,40]" "nmult" "Neutral_Multiplicity")
jet_pt=("[100,20,1020]" "pt" "Transverse_Momentum_(jet)_[GeV]")
jet_eta=("[50,-2.5,2.5]" "eta" "Eta_(jet)")
jet_phi=("[70,-3.5,3.5]" "phi" "Phi_(jet)")
jet_m=("[200,0,200]" "mass" "Invariant_Mass_(jet)_[GeV]")
jet_cx1=("[100,0,1]" "c_x1" "Charged_top_1_pT_Ratio_(jet)")
jet_cx2=("[100,0,1]" "c_x2" "Charged_top_2_pT_Ratio_(jet)")
jet_cx3=("[100,0,1]" "c_x3" "Charged_top_3_pT_Ratio_(jet)")
jet_nx1=("[100,0,1]" "n_x1" "Neutral_top_1_pT_Ratio_(jet)")
jet_nx2=("[100,0,1]" "n_x2" "Neutral_top_2_pT_Ratio_(jet)")
jet_nx3=("[100,0,1]" "n_x3" "Neutral_top_3_pT_Ratio_(jet)")
jet_axis1=("[35,0,0.35]" "axis1" "Major_Axis_(jet)")
jet_axis2=("[20,0,0.2]" "axis2" "Minor_Axis_(jet)")

jet_KS_d=("[30,0,30]" "KS_d" "had_d(KS_in_jet)")
jet_KS_pt=("[60,0,100]" "KS_pt" "Transverse_Momentum_(KS_in_jet)_[GeV]")
jet_KS_eta=("[80,-4,4]" "KS_eta" "Eta_(KS_in_jet)")
jet_KS_phi=("[80,-4,4]" "KS_phi" "Phi_(KS_in_jet)")
jet_KS_m=("[60,0.25,0.85]" "KS_mass" "Invariant_Mass_(KS_in_jet)_[GeV]")
jet_KS_angleXY=("[100,-1,1]" "KS_angleXY" "AngleXY_(KS_in_jet)")
jet_KS_chi2=("[60,0,10]" "KS_chi2" "Chi2_(KS_in_jet)")
jet_KS_dca=("[60,0,2]" "KS_dca" "DCA_(KS_in_jet)")

./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${Dilep_pt[0]} -p ${Dilep_pt[1]} -x ${Dilep_pt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${Dilep_m[0]} -p ${Dilep_m[1]} -x ${Dilep_m[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${Dilep_eta[0]} -p ${Dilep_eta[1]} -x ${Dilep_eta[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${Dilep_phi[0]} -p ${Dilep_phi[1]} -x ${Dilep_phi[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${MET_pt[0]} -p ${MET_pt[1]} -x ${MET_pt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${MET_sumEt[0]} -p ${MET_sumEt[1]} -x ${MET_sumEt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "event" -c "passedEvent" -w "genweight*puweight*mueffweight*tri" -b ${nJet[0]} -p ${nJet[1]} -x ${nJet[2]} -y "Events" -j "Golden" -d "True"

./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_pt[0]} -p ${had_pt[1]} -x ${had_pt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_m[0]} -p ${had_m[1]} -x ${had_m[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_eta[0]} -p ${had_eta[1]} -x ${had_eta[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_phi[0]} -p ${had_phi[1]} -x ${had_phi[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_d[0]} -p ${had_d[1]} -x ${had_d[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_angleXY[0]} -p ${had_angleXY[1]} -x ${had_angleXY[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_chi2[0]} -p ${had_chi2[1]} -x ${had_chi2[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_had" -c "" -w "genweight*puweight*mueffweight*tri" -b ${had_dca[0]} -p ${had_dca[1]} -x ${had_dca[2]} -y "Events" -j "Golden" -d "True"

#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_cmult[0]} -p ${jet_cmult[1]} -x ${jet_cmult[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_nmult[0]} -p ${jet_nmult[1]} -x ${jet_nmult[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_pt[0]} -p ${jet_pt[1]} -x ${jet_pt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_eta[0]} -p ${jet_eta[1]} -x ${jet_eta[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_phi[0]} -p ${jet_phi[1]} -x ${jet_phi[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_m[0]} -p ${jet_m[1]} -x ${jet_m[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_cx1[0]} -p ${jet_cx1[1]} -x ${jet_cx1[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_cx2[0]} -p ${jet_cx2[1]} -x ${jet_cx2[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_cx3[0]} -p ${jet_cx3[1]} -x ${jet_cx3[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_nx1[0]} -p ${jet_nx1[1]} -x ${jet_nx1[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_nx2[0]} -p ${jet_nx2[1]} -x ${jet_nx2[2]} -y "Events" -j "Golden" -d "True"
#./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_nx3[0]} -p ${jet_nx3[1]} -x ${jet_nx3[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_axis1[0]} -p ${jet_axis1[1]} -x ${jet_axis1[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_axis2[0]} -p ${jet_axis2[1]} -x ${jet_axis2[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_pt[0]} -p ${jet_KS_pt[1]} -x ${jet_KS_pt[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_m[0]} -p ${jet_KS_m[1]} -x ${jet_KS_m[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_eta[0]} -p ${jet_KS_eta[1]} -x ${jet_KS_eta[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_phi[0]} -p ${jet_KS_phi[1]} -x ${jet_KS_phi[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_d[0]} -p ${jet_KS_d[1]} -x ${jet_KS_d[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_angleXY[0]} -p ${jet_KS_angleXY[1]} -x ${jet_KS_angleXY[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_chi2[0]} -p ${jet_KS_chi2[1]} -x ${jet_KS_chi2[2]} -y "Events" -j "Golden" -d "True"
./MCDraw.py -t "MVA_jet" -c "" -w "genweight*puweight*mueffweight*tri" -b ${jet_KS_dca[0]} -p ${jet_KS_dca[1]} -x ${jet_KS_dca[2]} -y "Events" -j "Golden" -d "True"




