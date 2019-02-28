
for num in {0..1}
do
  python limit_calc.py pp_combined_JKS_BDT_highest $num
  python limit_calc.py pp_combined_JKS_BDT_highBDT $num
  python limit_calc.py pp_combined_J_BDT_highest $num
  python limit_calc.py pp_combined_J_BDT_highBDT $num

  python limit_calc.py pp_s_vs_b_JKS_BDT_highest $num
  python limit_calc.py pp_s_vs_b_JKS_BDT_highBDT $num
  python limit_calc.py pp_s_vs_b_J_BDT_highest $num
  python limit_calc.py pp_s_vs_b_J_BDT_highBDT $num

  python limit_calc.py JKS_wo_KS_highest $num
  python limit_calc.py JKS_wo_KS_highBDT $num
  python limit_calc.py JKS_wo_KS_s_vs_b_highest $num
  python limit_calc.py JKS_wo_KS_s_vs_b_highBDT $num

  python limit_calc.py JKS_wo_Jet_highest $num
  python limit_calc.py JKS_wo_Jet_highBDT $num
  python limit_calc.py JKS_wo_Jet_s_vs_b_highest $num
  python limit_calc.py JKS_wo_Jet_s_vs_b_highBDT $num

  python limit_calc.py all_highest $num
  python limit_calc.py all_highBDT $num
  python limit_calc.py all_s_vs_b_highest $num
  python limit_calc.py all_s_vs_b_highBDT $num

  python limit_calc.py all_wo_KS_highest $num
  python limit_calc.py all_wo_KS_highBDT $num
  python limit_calc.py all_wo_KS_s_vs_b_highest $num
  python limit_calc.py all_wo_KS_s_vs_b_highBDT $num
done
