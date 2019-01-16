
for num in [0.1] do;
  python limit_calc_sim.py pp_combined_JKS_BDT_highest pp_combined_J_BDT_highest $num
  python limit_calc_sim.py pp_combined_JKS_BDT_highBDT pp_combined_J_BDT_highBDT $num
  python limit_calc_sim.py pp_s_vs_b_JKS_BDT_highest pp_s_vs_b_J_BDT_highest $num
  python limit_calc_sim.py pp_s_vs_b_JKS_BDT_highBDT pp_s_vs_b_J_BDT_highBDT $num
done;
