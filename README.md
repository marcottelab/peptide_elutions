# peptide_elutions
Discovering proteolysis from peptide elutions

for f in data_files/BCM]*elut_long.csv; do echo "Rscript scripts/fit_mixedmodels.R
 --elut_long $f --outfile ${f%.csv}_mixed_models.csv --numgauss 2" ; done > mixed_models_COMMANDS.sh
