#All can be done with parallel


#Digest proteomes
record_trypsin.sh

#Sort to eggnog
record_emapper.sh

#Format eggnog
record_reformat.sh

#Assign to groups
record_grouping2.sh

#Get fractionation weighted proteins counts (for comparison)
record_msblender2elution_proteins.sh
#Tidy fractionation weighted counts
record_elutionTidy_proteins.sh

#Get fractionation peptides into one file
record_msblender2elution_peptides.sh
#Tidy the fractionation peptides to prepare for lookup
record_elutionTidy_peptides.sh

#Look up according to groups
record_lookup_plants.sh

#Made wide version of elution matrix
record_widen_elution.sh
