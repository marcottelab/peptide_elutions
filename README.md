# peptide_elutions

Analysis accompanying "Alternative proteoforms and proteoform-dependent assemblies in humans and plants"

Claire D. McWhite, Wisath Sae-Lee, Anna L. Mallam, Nicolas A. Gort Freitas, and Edward M. Marcotte

### Input data

Four column csv, describing the which fraction each peptide is found in, and its amount. 

ExperimentID     | FractionID                 | Peptide                           | PeptideCount
---------------- | -------------------------- | --------------------------------- | ------------
Hemolysate_IEX_1 | Hemolysate_IEX_06_10032017 | ACANPAAGSVILLENLR                 | 3.0
Hemolysate_IEX_1 | Hemolysate_IEX_06_10032017 | ADGLAVIGVLMK                      | 33.0
Hemolysate_IEX_1 | Hemolysate_IEX_07_10032017 | ADGLAVIGVLMKVGEANPK               | 2.0
Hemolysate_IEX_1 | Hemolysate_IEX_07_10032017 | ADGLAVIGVLMKVGEANPKLQKVLDALQAIK   | 4.0
Hemolysate_IEX_1 | Hemolysate_IEX_07_10032017 | ADGLAVIGVLMKVGEANPKLQKVLDALQAIKTK | 3.0

Alternately, begin with a wide csv of peptide identifications x fractionID

Hemolysate_IEX_1.wide.csv
Peptide | Hemolysate_IEX_07 | Hemolysate_IEX_08 | Hemolysate_IEX_09_10032017
--------|-------------------|-------------------|---------------------------


### Processing

### Image creation

### Scoring

### Visual caterogizing







for f in data_files/BCM]*elut_long.csv; do echo "Rscript scripts/fit_mixedmodels.R
 --elut_long $f --outfile ${f%.csv}_mixed_models.csv --numgauss 2" ; done > mixed_models_COMMANDS.sh
