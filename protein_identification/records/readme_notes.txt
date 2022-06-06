### template commands can be found at: 
/project/cmcwhite/data/peptide_elutions/protein_identification/records

### generate in silico trypsinized peptides: 
python /project/cmcwhite/pipelines/MS_grouped_lookup/scripts/proteome_utils/trypsin.py -i /project/vyqtdang/proteomes_eggNOG/human/uniprot-proteome_Homosapiens.fasta -o peptides/human_uniprot-proteome_human_reviewed_peptides.csv -m 2 TRUE

### generate in silico trypsinized peptides WITH positions: 
python /project/cmcwhite/pipelines/MS_grouped_lookup/scripts/proteome_utils/trypsin.py -i /project/vyqtdang/proteomes_eggNOG/pig/uniprot-proteome_Susscrofa.fasta -o peptides_DB/pig_uniprot-proteome_peptides_positions.csv -m 2 -p TRUE


### get fractionation peptides into one file - consolidate identified peptides from multiple experiments into a single file by combining all pepcount files together:
python2 /project/cmcwhite/github/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /MS.processed/Lumos_Marcotte/Vy/Human_23555_SEC_20210828/results/*/output/*pep*count*1 --output_filename consolidated_pepcounts/Human_23555_SEC_20210828.pepcount --fraction_name_from_filename --msblender_format --spectral_count_type TotalCount --pepcount

### get .pepcount into tidy format: 
python2 /stor/home/vyqtdang/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution consolidated_pepcounts/Human_23555_SEC_20210828.pepcount --outfile consolidated_pepcounts/Human_23555_SEC_20210828.pepcount.tidy --firstcol Peptide --valuecol PeptideCount --experiment_id Human_23555_SEC_20210828

### get unique trypsinized peptides from whole proteome:
python /project/cmcwhite/pipelines/MS_grouped_lookup/scripts/lookup_utils/define_grouping.py --peptides peptides/human_uniprot-proteome_human_reviewed_peptides.csv --output_basename human_uniprot-proteome_human_reviewed_unique_peptides.csv
