echo "no constraints"
python scripts/count_stem_basepair_cml.py test_no_constraints/ --output_csv no_constraints_bulge_len.csv
echo "restrict_loop"
python scripts/count_stem_basepair_cml.py test_restrict_loop/ --output_csv restrict_loop_bulge_len.csv
echo "restrict_lower_stem"
python scripts/count_stem_basepair_cml.py test_restrict_lower_stem/ --output_csv restrict_lower_stem_bulge_len.csv
echo "restrict_lower_stem_loop"
python scripts/count_stem_basepair_cml.py test_restrict_lower_stem_loop/ --output_csv restrict_lower_stem_loop_bulge_len.csv
echo "Done"