



#for meson in d0, jpsi
#do


sed -i 's/whatcm/jpsi/g' var_cm.py
sed -i 's/number/443/g' var_cm.py

sed 's/var/lxySig/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, .5))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/lxy/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 1.))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/l3D/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 1.))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/dca/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, .2))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/mass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 2., 4.))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/mass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 1.6, 2.2))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/angleXY/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.05, 1.05))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/angleXYZ/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -1.05, 1.05))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/jetDR/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 40, 0, .4))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/legDR/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 40, 0, .8))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/eta/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, -3, 3))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/trk_normalizedChi2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 100, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/trk_nHits/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/trk_pt/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 100, 0, 20))/g' all_cm.py && cmsRun all_cm.py

sed 's/var/trk_ipsigXY/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
sed 's/var/trk_ipsigZ/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 10))/g' all_cm.py && cmsRun all_cm.py
#
#sed 's/var/pt/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 100))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/mass/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, 0, 2.5))/g' all_cm.py && cmsRun all_cm.py
#sed 's/var/chi2/g' var_cm.py > all_cm.py && sed -i 's/, 50, 0, .5))/, 50, -150, 50))/g' all_cm.py && cmsRun all_cm.py

sed -i 's/jpsi/whatcm/g' var_cm.py
sed -i 's/443/number/g' var_cm.py
#done
