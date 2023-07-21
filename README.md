# Proteomic changes induced by longevity promoting interventions in mice
Code for processing and analyzing data described in the manuscript "Proteomic changes induced by longevity promoting interventions in mice" (2023).


## File order:

1.) /00_Format_Data/lc_proteomics_m005_formatting.r

R code for formatting available peptide data for downstream analyses.

2.) /01_Protein_Rollup/lc_proteomics_m005_rollup.do

STATA code for rolling up measured peptide into protein level measurements.

2.) /02_Normalization/lc_proteomics_m005_normalize.do

STATA code for normalizing protein measurements.

3.) /03_Hurdle_Model/lc_proteomics_m005_hurdle_model.do

STATA code for hurdle and linear models to estimate effects of longevity on protein measurements.
