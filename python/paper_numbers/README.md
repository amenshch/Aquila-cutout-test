# Scripts behind the numbers of hires2_correction.tex

Each script takes the directory of hicorrect2.py (with the constants of the adopted law) as its first argument and prints the numbers that the paper quotes. The comments starting with "% NUMBERS" in the LaTeX file say which script and which input each number comes from.

| Script | Inputs | Paper |
|---|---|---|
| law_checks.py | attenuation_law.txt, attenuation_law-leave_one_cylinder_out.txt | Sect. 3.2: departures of spheres and cylinders, leave-one-cylinder change of the factor |
| cylinder_table.py | the seven hires/*_isolated directories | Table tab:cyl, Sect. 5.2 |
| simsky_numbers.py | simsky6_RT_enh_cirrdepth=2_nocrater, simsky6_RT_enh_cirronly, simsky6.w40.cores.cores.txt | Sect. 5.3 |
| aquila_numbers.py | hires/Aquila_cutout_g33 | Sect. 6 |
| inclination_errors.py | (hicorrect2.py only) | Sect. 7, inclination |
| superposition_numbers.py | superposition_test.py (next to hicorrect2.py) | Table tab:concentr, Appendix B |
| field_factor_errors.py | hires/harvest_calib3_G0=*.txt.hist.npz | Table tab:fieldR, Sect. 7 radiation field |

Not here, because they are the user's own scripts: fit_attenuation_law.py (the law, Sect. 3.2), validate_calib3.py (Table tab:valgrid, Sect. 5.1), diagnose_shielding_branches.py (the outermost-shell excess of 2.4 K).

After a change of the law: refit, update the constants in hicorrect2.py, hicorrect2.for and make_hires2_figures12.py, rerun hicorrect2 on the cylinders, the simulated sky and Aquila, rerun validate_calib3.py, then run the scripts above and update the numbers next to each "% NUMBERS" comment.
