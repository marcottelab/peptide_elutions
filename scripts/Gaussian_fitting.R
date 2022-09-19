library(tidyverse)
library(argparse)
library(furrr)
library(purrr)
library(future)
library(AdaptGauss)
plan("multisession")
#This script is used to find number of elution peaks for each protein in a fractionation experiment. Only works with R>4.0.

parser <- ArgumentParser(description='Find number of elution peaks for each protein in a fractionation experiment. Only works with R>4.0.')

parser$add_argument('--input_file', dest='input_file', action='store', required=TRUE,
    help='input file from peptide_identification.R')

parser$add_argument('--simple_AdapGauss', dest='simple_AdapGauss', action='store', required=TRUE,
    help='file location for simple_AdaptGauss.R')

parser$add_argument('--output_file', dest='output_file', action='store', required=TRUE,
    help='output file has 4 columns: ProteinID, experiment_name, FractionOrder, peak')


args = parser$parse_args()

source(args$simple_AdapGauss)

print("load input file")
short_tidy_unique <- read_csv(args$input_file)


print("find peaks")
peaks <-  
  short_tidy_unique %>%
  filter(pepcount > 0) %>%
  split(list(.$ProteinID, .$experiment_name)) %>%
  future_map_dfr(~possibly_get_boundaries(.x), .id = "id", .progress = TRUE) %>% 
  separate(id, into = c("ProteinID", "experiment_name"), sep = "[.]")

print("writing result file")
peaks %>% write_csv(args$output_file)

