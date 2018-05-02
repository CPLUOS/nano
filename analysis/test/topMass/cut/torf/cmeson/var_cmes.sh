



#for meson in d0, jpsi
#do


sed -i 's/whatcm/d0/g' var_cm.py
sed -i 's/number/421/g' var_cm.py

#sed 's/var/lxySig/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/lxy/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 1))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/lxyE/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/l3D/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 2))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/l3DE/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dca/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, .2))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/mass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 2., 4.))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/mass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 1.7, 2.))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/angleXY/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.05, 1.05))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/angleXYZ/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.05, 1.05))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/jetDR/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 40, 0, .4))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/legDR/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 40, 0, .8))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/diffMass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 35, 0.135, 0.17))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/eta/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, -3, 3))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/pt/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/phi/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 40, -4, 4))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/chi2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
#
sed 's/var/jet_btagCMVA/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.05, 1.05))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/jet_btagCSVV2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -15, 5))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/jet_btagDeepB/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -2.5, 2.5))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/jet_btagDeepC/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.5, 1.))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau1_chi2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau1_ipsigXY/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau1_ipsigZ/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau1_nHits/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau1_pt/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 20))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau2_chi2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau2_ipsigXY/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau2_ipsigZ/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau2_nHits/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dau2_pt/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 20))/g' all_cm.py && cmsRun all_cm.py

sed -i 's/d0/whatcm/g' var_cm.py
sed -i 's/421/number/g' var_cm.py
#done
