library(tidyverse)
library(argparse)
library(furrr)
library(purrr)
library(future)
library(argparse)
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

get_terminal_ranges <- function(oneexp){
  firstobs <- oneexp %>%
    filter(presence == "observed") %>%
    filter(Start == min(Start)) %>% pull(Start) %>% head(1)
  lastobs <- oneexp %>%
    filter(presence == "observed") %>%
    filter(Start == max(Start)) %>% pull(Start) %>% head(1)
  
  # observed = in this experiment, possible = seen at all across experiments / detectable
  # Look from first possible peptide to last observed peptide
  missing_n <- oneexp %>%
    filter(Start <= lastobs)
  # Look from first observed peptide to last possible peptide
  missing_c <- oneexp %>%
    filter(Start >= firstobs)
  return(list(missing_n$Start, missing_c$Start))
}

get_missing_obs_counts <- function(start, oneexp){
  
  df <-oneexp %>%
    group_by(ProteinID, peak, experiment_name) %>%
    mutate(fragment_end = case_when(Start < {{ start }}  ~ "nterm",
                                    Start >= {{ start }} ~ "cterm")
    ) %>% ungroup %>%
    group_by(ProteinID, peak, experiment_name, fragment_end) %>%
    dplyr::count(presence) %>% # Overwritted by seqinr
    ungroup %>%
    mutate(n= n + 1) %>%
    pivot_wider(names_from = c(presence, fragment_end), values_from = n, values_fill = 1) %>%
    mutate(bp_start = {{ start }})
  return(df)
}


get_bps <- function(oneexp){
  # Terminal bias calculation
  missing <- get_terminal_ranges(oneexp)
  missing_n <- missing[[1]]
  missing_c <- missing[[2]]
  
  print(missing_n)
  # C terminal fragments are missing n termini
  c_fragments <- map_dfr(missing_n, ~get_missing_obs_counts(., oneexp)) 
  
  c_fragments$terminal <- "cterm_extension"
  # N terminal fragments are missing c termini
  n_fragments <- map_dfr(missing_c, ~get_missing_obs_counts(., oneexp)) 
    
   n_fragments$terminal <-"nterm_extension"
  
  
  fragments <- bind_rows(n_fragments, c_fragments) %>%
    mutate_at(vars(matches("_[nc]term")), replace_na, 1) 
  
  fragments <- fragments %>%
    
    mutate(       ratio_cterm = observed_cterm/(observed_cterm + missing_cterm),
                  ratio_nterm = observed_nterm/(observed_nterm + missing_nterm)
    ) %>% 
    #select(-breakpointleft_pepid, -breakpointright_pepid) %>%
    mutate(fc = ratio_nterm/ratio_cterm) %>%
    mutate(log2fc = log2(fc)) %>%
    mutate(abs_log2fc = abs(log2fc))
  
  return(fragments)
  
}


possibly_get_bps <- possibly(get_bps, otherwise = data.frame())





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

