library(tidyverse)
library(argparse)
library(furrr)
library(purrr)
library(future)
plan("multisession")
#This script is used to calculate terminal bias score for each protein in a fractionation experiment. Only works with R>4.0.

parser <- ArgumentParser(description='Calculate terminal bias score for each protein')

parser$add_argument('--input_file', dest='input_file', action='store', required=TRUE,
    help='input file from peptide_identification.R')

parser$add_argument('--peaks', dest='peaks', action='store', required=TRUE,
    help='file containing peaks from Gaussian_fitting.R')

parser$add_argument('--output_file', dest='output_file', action='store', required=TRUE,
    help='output file ****add more info here*****')


args = parser$parse_args()

print("load functions for terminal bias calculation")

library(tidyverse)
#library(argparse)

elut_wide_file <- '/Users/mwsaelee/Desktop/peptide_elutions-master/test/pivot_test.csv'
peps_file <- '/Users/mwsaelee/Desktop/peptide_elutions-master/test/uniprot_human_digested.csv'
fraction_order_file  <- '/Users/mwsaelee/Desktop/peptide_elutions-master/test/fraction_order_test.csv'
  

elut_wide <- read_delim(elut_wide_file,  delim = ",") 
fraction_order <- read_delim(fraction_order_file,  delim = ",")
peps <- read_delim(peps_file,  delim = ",")

print("Identify peptides")
elut_wide_identified <- elut_wide %>%
  inner_join(peps, by = "Peptide")

print("Identify proteins with no peptides")
proteins_w_unique_peps <- elut_wide_identified %>%
  select(ProteinID, ExperimentID) %>%
  unique

print("Convert wide to tidy")
elut_long <- elut_wide_identified %>%
  pivot_longer(cols = starts_with("MB_Sup_HSA"),
               names_to = "FractionID",
               values_to = "pepcount") %>%
  mutate(pepcount = as.numeric(pepcount)) %>%
  left_join(fraction_order, by = "FractionID") %>%
  mutate(FractionID = str_replace(FractionID, "^fractionid_", "")) %>%
  filter(!is.na(pepcount)) # Unneeded

#fraction_order %>% write_csv(paste0(elut_wide_file, "fraction_order.csv"))

elut_long %>% write_csv(paste0(elut_wide_file, ".annot.long.tidy"))

total_fractions <- c(nrow(fraction_order))


#10x smaller
elut_short <- elut_long %>% 
  filter(pepcount > 0) %>%
  mutate(totfracts = total_fractions)




elut_short %>% write_csv(paste0(elut_wide_file, ".annot.short.tidy"))



print("load input file")
short_tidy_unique <- read_csv(args$input_file)
peaks  <- read_csv(args$peaks)

#Calculate for terminal bias 

#join input table (info of peptides) to resulting peaks from Gaussian fitting
short_tidy_unique_peaks <-  short_tidy_unique %>%
     left_join(peaks,  by = c("ProteinID","experiment_name", "FractionOrder"))

#Get proteins that have at least 10 unique peptides
numpeps <- short_tidy_unique_peaks %>%
    select(ProteinID, Peptide) %>% unique %>%
    group_by(ProteinID) %>%
       summarize(numpeps = n()) %>% ungroup 

numpeps_sel <- numpeps %>% filter(numpeps >= 9)


#Remove proteins with low peptide coverage (<9 unique peptides)
short_tidy_unique_peaks_expanded <- short_tidy_unique_peaks %>%
  semi_join(numpeps_sel) %>% # Get rid of proteins with <9 unique peptides per an experiment
     
  select(-ExperimentID, -spec, -ID, -experiment_order) %>% # Unused columns
  select(-FractionID, -pepcount) %>% unique %>%
   mutate(pepid = paste0(Start, "-", End)) %>%
   mutate(presence = "observed") %>%
  
    group_by(ProteinID) %>%
     complete(nesting(pepid, Peptide), nesting(experiment_name, peak)) %>%
   ungroup %>%
     separate(pepid, into = c("Start", "End"), sep = "-", convert = TRUE, remove = FALSE) %>% # get start and end for filled in peptides
    mutate(presence = replace_na(presence, "missing")) # Filled in peptides get zero presence

# Get peaks that have at least 7 peptides spread across fractions
enough_cov_peaks <- short_tidy_unique_peaks %>%
   filter(!is.na(peak)) %>%
    group_by(ProteinID, peak, experiment_name) %>%
    mutate(n = n()) %>%
  ungroup %>%
  filter(n >= 7) %>%
  select(-n)

short_tidy_unique_peaks_expanded_fragments  <- short_tidy_unique_peaks_expanded %>%
  filter(!is.na(peak)) %>%
  semi_join(enough_cov_peaks, by = c("ProteinID", "peak", "experiment_name"))

bps <- short_tidy_unique_peaks_expanded_fragments %>%
  arrange(Start, End)%>%
  split(list(.$ProteinID, .$experiment_name, .$peak)) %>%
  future_map_dfr(~possibly_get_bps(.))


print("writing result file")
bps %>% write_csv(args$output_file)

