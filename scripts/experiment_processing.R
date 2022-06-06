library(tidyverse)
library(argparse)

source('scripts/fasta_utils.R')

parser <- ArgumentParser(description='preprocess peptides for a species to prepare to use to interpret experiment peptides')
parser$add_argument('--proteome', dest='proteome', action='store', required=TRUE,
    help='A .fasta file containing the proteome of the species')

parser$add_argument('--all_peps', dest='all_peps', action='store', required=TRUE,
    help='File of all peptides for the species')

parser$add_argument('--protuniq_peps', dest='protuniq_peps', action='store', required=TRUE,
    help='File of all peptides unique to a protein')

parser$add_argument('--orthouniq_peps', dest='orthouniq_peps', action='store', required=TRUE,
    help='File of all peptides unique to an orthogroup')

parser$add_argument('--close_orthology', dest='close_orthology', action='store', required=TRUE,
    help='Mapping of each protein to the veNOG or virNOG orthology')

parser$add_argument('--inclusive_orthology', dest='inclusive_orthology', action='store', required=TRUE,
    help='Mapping of each protein to the euNOG orthology. spec.euNOG.mapping')

parser$add_argument('--spec', dest='spec', action='store', required=TRUE,
    help='Which species is being processed (for outfile naming)')


args = parser$parse_args()

proteome <- args$proteome
elution_pattern <- args$elution_pattern
euNOG_orthology <- args$inclusive_orthology
close_orthology <- args$close_orthology
all_peps <- args$all_peps
protuniq_peps <- args$protuniq_peps
orthouniq_peps <- args$orthouniq_peps
spec <- args$spec


protein_unique_peptides <- read_csv(args$protuniq_peps) %>%
    mutate(Peptide = str_replace_all(Peptide, "[I|L]", "J")) %>%
               mutate(status = "protein_unique")

protein_unique_peptides %>%
    write_csv(paste0("data/",spec, "_protein_unique_peptides.csv"))

orthogroup_unique_peptides <- read_csv(args$orthouniq_peps) %>%
    mutate(Peptide = str_replace_all(Peptide, "[I|L]", "J")) %>%
               mutate(status = "orthogroup_unique")

orthogroup_unique_peptides %>%
    write_csv(paste0("data/",spec, "_orthogroup_unique_peptides.csv"))


protein_all_peptides <- read_csv(all_peps) %>%
               mutate(Peptide = str_replace_all(Peptide, "[I|L]", "J"))

protein_all_peptides %>%
    write_csv(paste0("data/",spec, "_protein_all_peptides.csv"))

close_orthology <- read_delim(close_orthology, delim = "\t")
euNOG_orthology <- read_delim(euNOG_orthology, delim = "\t") 


# Derive proteins where each orthogroup unique peptide could come from
orthogroup_unique_peptides_identified <- orthogroup_unique_peptides %>%
  full_join(close_orthology, by = c("ID")) 

#The orthology is missing proteins that aren't sorted into orthogroups.
#Add them back in
orthogroup_unique_peptides_identified$ProteinID[is.na(orthogroup_unique_peptides_identified$ProteinID)] <- 
  orthogroup_unique_peptides_identified$ID[is.na(orthogroup_unique_peptides_identified$ProteinID)]


#Combine list of all peptides with uniquenss information
#Keep peptides that are at least orthogroup-unique 

peps <- protein_all_peptides %>%
      left_join(protein_unique_peptides, by = c("ProteinID", "Peptide")) %>%
      inner_join(orthogroup_unique_peptides_identified, by = c("ProteinID", "Peptide"),
                                 suffix = c("", "_secondary_uniqueness"))


# If a peptide isn't protein unique, but is orthogroup unique, give it orthogroup-unique status. 
# If protein unique and orthogroup unique, keep as protein unique
peps$status[is.na(peps$status)] <- peps$status_secondary_uniqueness[is.na(peps$status)]

# Add on euNOG information for finding conserved events
peps <- peps %>% select(-status_secondary_uniqueness, -ID) %>%
  left_join(euNOG_orthology, by = c("ProteinID"))


peps %>%  %>%
    write_csv(paste0("data/",spec, "_annotated_peptides.csv"))


