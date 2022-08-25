# peptide_elutions

Analysis accompanying "Alternative proteoforms and proteoform-dependent assemblies in humans and plants"

Claire D. McWhite, Wisath Sae-Lee, Anna L. Mallam, Nicolas A. Gort Freitas, and Edward M. Marcotte

***TO-DO***
A short paragraph to explain the pipeline and how it works.

### Input data

***TO-DO***
Simplify input data format since peptide_identification.R will reformat input file to short_tidy form anyway.

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
***TO-DO***
modify peptide_identification.R to be able to take result file from both MSFragger and MSBlender
write up and give example of how to run the script

Process peptide files for Gaussian fitting 

### Scoring
1. Idenitfy peaks from peptide elution profile of each protein in a fractionation experiment using Gaussian Mixture Model. Multiple peaks suggest the existence of proteoforms or intact proteins eluting with different binding partners. 
Script: Gaussian_fitting.R
Input file:
![alt text](https://user-images.githubusercontent.com/32718019/186776934-acc71510-69cb-474d-91ae-63d62fa8c032.png)

Example:

Expected result:


2. Calculate terminal bias score in order to prioritize proteins to inspect manually.
-make .R file for Line 1097 - Line 1149

### Visual caterogizing
***TO-DO***
-clean up elution_viewer.Rmd

Visualize proteins with terminal bias score>2 through a Shiny app
A section to explain how to look for proteoform visually.



for f in data_files/BCM]*elut_long.csv; do echo "Rscript scripts/fit_mixedmodels.R
 --elut_long $f --outfile ${f%.csv}_mixed_models.csv --numgauss 2" ; done > mixed_models_COMMANDS.sh


# Peptide identification and formatting
## 1) Format result files from MSFragger using format_MSFragger_files.py
	python3 /stor/home/mwsaelee/peptide_elutions-master/scripts/format_MSFragger_files.py --root_folder '/stor/work/Marcotte/MS/processed/Lumos_Marcotte/Momo/menstrual_blood_sup_SEC_1/fragger' --fractionation_name 'MB_sup_SEC1' --output_file '/stor/home/mwsaelee/peptide_elutions-master/scripts/outputfile.csv'


## 2) Format result files from MSBlender by following the steps below
### These steps use scripts from two other github repositories

- github.com/marcottelab/MS_grouped_lookup
- github.com/marcottelab/protein_complex_maps


## Key steps

### Get fractionation peptides into one file - consolidate identified peptides from multiple experiments into a single file by combining all pepcount files together:
### Create combined table of all experimentally identified peptides
#### In this case, input files are the output from github.com/marcottelab/run_msblender and github.com/marcottelab/msblender

python2 /path/to/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files example_files/*pep*count*1 --output_filename consolidated_pepcounts/example.pepcount --fraction_name_from_filename --msblender_format --spectral_count_type TotalCount --pepcount

**Input tsvs of outputs of the marcottelab run_msblender/msblender pipeline:**

filename1: example_experiment_fraction10.pep_count_mFDRpsm001 

#Peptide FDR: 0.1
#PepSeq TotalCount      example_experiment_fraction10
MEATAK       5       5
ELVISR   3       3

filename2: example_experiment_fraction11.pep_count_mFDRpsm001 
#Peptide FDR: 0.1
#PepSeq TotalCount      example_experiment_fraction11
EAPEPTIDE 6	6	
ELVISR	1	1


**Output tsv combined from all fractions:**
filename: consolidated_pepcounts/example.pepcount
	example_experiment_fraction10	example_experiment_fraction11
MEATAK	5	0
EAPEPTIDE 0	6
ELVISR	3	1


### Convert .pepcount table into into tidy formatted csv: 
python2 /path/to/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution consolidated_pepcounts/example.pepcount --outfile consolidated_pepcounts/example.pepcount.tidy --firstcol Peptide --valuecol PeptideCount --experiment_id example_experiment

**Input is output of previous step:**

**Output csv:**
filename:consolidate_pepcounts/example.pepcount.tidy
ExperimentID,FractionID,Peptide,PeptideCount
OP_Corn_20181009,Corn_SEC_08a_102018,AAAAAGGGLFPMPDPK,1.0
example_experiment,example_experiment_fraction10,MEATAK,5
example_experiment,example_experiment_fraction11,EAPEPTIDE,6
example_experiment,example_experiment_fraction10,ELVISR,3
example_experiment,example_experiment_fraction11,ELVISR,1

### From a proteome fasta file, generate possible trypsinized peptide, allowing up to two missed cleavages

python /path/to/MS_grouped_lookup/scripts/proteome_utils/trypsin.py -i example_files/example.fasta -o protein_identification/example_peptides.csv -m 2 -p TRUE

**Input fasta file:**

>prot1
MEATAKEAPEPTIDE
>prot2
ELVISRLIVES
>prot3
MAKELVISR

**Output csv of possible peptides:**

ProteinID,Peptide,Start,End
prot1,MEATAK,1,6
prot1,MEATAKEAPEPTIDE,1,15
prot1,EAPEPTIDE,7,15
prot2,ELVISR,1,6
prot2,ELVISRLIVES,1,11
prot2,LIVES,7,11
prot3,ELVISR,4,9
prot3,MAKELVISR,1,9


### Reduce these possible peptides to only ones that are unique to an individual protein
python /path/to/MS_grouped_lookup/scripts/lookup_utils/define_grouping.py --peptides example_files/example_peptides.csv --output_basename example_unique_peptides.csv

**Input is output of previous step**

**Output csv of only unique peptides:**
filename:example_unique_peptides.csv
ProteinID,Peptide
prot1,MEATAK,1,6
prot1,MEATAKEAPEPTIDE,1,15
prot1,EAPEPTIDE,7,15
prot2,ELVISRLIVES,1,11
prot2,LIVES,7,11
prot3,MAKELVISR,1,9

### With files of experimentally observed peptides, and protein-unique peptides, can now identify proteins from each experiment


Rscript scripts/peptide_identification.R --elut_wide example_file/example.pepcount --peps example_files/example_unique_peptides.csv


**Input: Output of previous steps**

**Outputs:**


1. pepcount.annot.long.tidy 
Peptide,ExperimentID,experiment_order,experiment_name,ProteinID,Start,End,FractionID,pepcount,FractionOrder,totfracts
MFFESR,soybn_OP_Soy_sprout_SEC_20172110,19,SOYBN_SEC1,tr|I1KSR1|I1KSR1_SOYBN,1,6,Soy_sprout_SEC_26_1a_10212017,1,24,67

2. pepcount.annot.short.tidy
Same as above file, but with peptides with zero observations removed

3. fraction_order.csv: Preserves order of each fraction in the experiment
FractionOrder,FractionID
1,fractionid_Soy_sprout_SEC_03_1a_10212017
2,fractionid_Soy_sprout_SEC_04_1a_10212017
3,fractionid_Soy_sprout_SEC_05_1a_10222017

