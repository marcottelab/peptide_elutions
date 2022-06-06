library(tidyverse)
library(argparse)

source("scripts/fasta_utils.R")

parser <- ArgumentParser(description='Get lengths of all proteins in a proteome')

parser$add_argument('--proteome', dest='proteome', action='store', required=TRUE,
    help='A proteome is fasta format')
parser$add_argument('--spec', dest='spec', action='store', required=TRUE,
    help='Name outfile by species code')



args = parser$parse_args()

proteome_file <- args$proteome 
spec <- args$spec

protein_lengths <- read_fasta(proteome_file) %>%
  mutate(seqlen = str_length(Sequence)) %>%
  select(ID, seqlen) %>%
  rename(ProteinID = ID)

protein_lengths %>% 
  write_csv(paste0("data/", spec, "_protein_lengths.csv"))

