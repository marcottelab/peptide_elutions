library(tidyverse)
library(argparse)


#This script is used to format peptide file from fractionation experiment for peak fitting.

parser <- ArgumentParser(description='format peptide file from fractionation experiment for peak fitting.')

parser$add_argument('--elut_wide_file', dest='elut_wide_file', action='store', required=TRUE,
    help='input file in wide format')

parser$add_argument('--fraction_order', dest='fraction_order', action='store', required=TRUE,
    help='fraction order file')

parser$add_argument('--peps', dest='peps', action='store', required=TRUE,
    help='file containing info on in-silico digested peptides from trypsin.py')

parser$add_argument('--seqlen', dest='seqlen', action='store', required=TRUE,
    help='proteome file for the organism of interest (tsv format)')

parser$add_argument('--spec', dest='spec', action='store', required=TRUE,
    help='species of the sample. text format')

parser$add_argument('--output_file', dest='output_file', action='store', required=TRUE,
    help='ouput filename')

args = parser$parse_args()


elut_wide <- read_delim(args$elut_wide_file,  delim = ",") 
fraction_order <- read_delim(args$fraction_order,  delim = ",")
peps <- read_delim(args$peps,  delim = ",")
col_name <- unlist(strsplit(fraction_order$FractionID[1], split = '_' ))[1]

print("Identify peptides")
elut_wide_identified <- elut_wide %>%
  inner_join(peps, by = "Peptide")


print("Identify proteins with no peptides")
proteins_w_unique_peps <- elut_wide_identified %>%
  select(ProteinID, ExperimentID) %>%
  unique


print("Convert wide to tidy")
elut_long <- elut_wide_identified %>%
  pivot_longer(cols = starts_with(col_name),
               names_to = "FractionID",
               values_to = "pepcount") %>%
  mutate(pepcount = as.numeric(pepcount)) %>%
  left_join(fraction_order, by = "FractionID") %>%
  mutate(FractionID = str_replace(FractionID, "^fractionid_", "")) %>%
  filter(!is.na(pepcount)) # Unneeded
  

seq_len <- read_delim(args$seqlen,  delim = "\t")
seq_len_formatted <- seq_len %>% 
  select(c('Entry', 'Length'))
colnames(seq_len_formatted)[colnames(seq_len_formatted) == 'Entry'] <- 'ID'
colnames(seq_len_formatted)[colnames(seq_len_formatted) == 'Length'] <- 'seqlen'



elut_short <- elut_long %>% 
  filter(pepcount > 0) %>%
  mutate(totfracts = c(nrow(fraction_order))) %>%
  mutate(experiment_name = c(unique(elut_long$ExperimentID))) %>%
  mutate(Protein_ID=paste(ProteinID)) %>% 
  mutate(spec = c(args$spec)) %>% 
  mutate(experiment_order = c(1)) %>% 
  separate(col=Protein_ID, c('sp', 'ID', 'gene'), sep="\\|") %>% 
  select(-c('sp', 'gene')) %>% 
  inner_join(seq_len_formatted, by = "ID")


print("writing result file")
elut_short %>% write_csv(args$output_file)

