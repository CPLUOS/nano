json_used='Golden'
cut='passedEvent'
weight='genweight*puweight*eleffweight*mueffweight*tri'
plotvar_dilep_m='dilep.M()'
plotvar_dilep_pt='dilep.Pt()'
binning="[60, 20, 320]"
x_name_dilep_m='Invariant mass(ll) [GeV]'
x_name_dilep_pt='Transverse momentum(ll) [GeV]'
y_name='Events'
dolog='True'


./vtsDraw.py -c $cut -w $weight -b $binning -p $plotvar_dilep_m -x $x_name_dilep_m -y $y_name -j $json_used -d $dolog
./vtsDraw.py -c $cut -w $weight -b $binning -p $plotvar_dilep_pt -x $x_name_dilep_pt -y $y_name -j $json_used -d $dolog
