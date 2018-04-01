sed 's/var/lxy/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 50, 0, 1.))/g' all.py && cmsRun all.py
sed 's/var/l3D/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 50, 0, 1.))/g' all.py && cmsRun all.py
sed 's/var/dca/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 50, 0, .2))/g' all.py && cmsRun all.py
sed 's/var/mass/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 50, 1.6, 3.4))/g' all.py && cmsRun all.py
#sed 's/var/diffMass/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 60, 0.12, 0.18))/g' all.py && cmsRun all.py
sed 's/var/eta/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 60, -3, 3))/g' all.py && cmsRun all.py
sed 's/var/angleXY/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 42, -1.05, 1.05))/g' all.py && cmsRun all.py
sed 's/var/angleXYZ/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 42, -1.05, 1.05))/g' all.py && cmsRun all.py
sed 's/var/jetDR/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 40, 0, .4))/g' all.py && cmsRun all.py
sed 's/var/legDR/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 40, 0, .8))/g' all.py && cmsRun all.py
sed 's/var/trk_normalizedChi2/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 100, 0, 10))/g' all.py && cmsRun all.py
sed 's/var/trk_nHits/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 60, 0, 30))/g' all.py && cmsRun all.py
sed 's/var/trk_pt/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 100, 0, 20))/g' all.py && cmsRun all.py
sed 's/var/trk_ipsigXY/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 100, 0, 5))/g' all.py && cmsRun all.py
sed 's/var/trk_ipsigZ/g' var.py > all.py && sed -i 's/, 50, 0, .5))/, 100, 0, 5))/g' all.py && cmsRun all.py


