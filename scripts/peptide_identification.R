library(tidyverse)
library(argparse)

Rscript scripts/peptide_identification.R --elut_wide example_file/example_experiment.pepcount --peps example_files/example_unique_peptides.csv

parser <- ArgumentParser(description='Identify peptide -> protein in an experiments')

parser$add_argument('--elut_wide', dest='elut_wide', action='store', required=TRUE,
    help='Wide elution table')

parser$add_argument('--peps', dest='peps', action='store', required=TRUE,
    help='Output of define_grouping.py')


args = parser$parse_args()

elut_wide_file <- args$elut_wide
peps <- read_csv(args$peps)

# Load elution
elut_wide <- read_delim(elut_wide_file,  delim = "\t") %>%
        #read_delim(elut_wide_file, col_names = FALSE, delim = "\t") %>%
        mutate(ExperimentID = str_replace({{ elut_wide_file }}, ".pepcount", "")) %>%
        mutate(ExperimentID = str_replace(ExperimentID, "data/", "")) %>%
         rename(Peptide = X1) %>%
        filter(!is.na(Peptide)) %>%
        mutate(Peptide = str_replace_all(Peptide, "[I|L]", "J")) %>%
  rename_with(., function(x){paste0("fractionid_", x)},  c(-ExperimentID,-Peptide))

fraction_order <- elut_wide %>% head(0) %>%
  select(starts_with("fractionid_")) %>%
  colnames() %>%
  as_tibble(rownames = "FractionOrder") %>%
  mutate(FractionOrder = as.integer(FractionOrder)) %>%
  rename(FractionID = value)



print("Identify peptides")
elut_wide_identified <- elut_wide %>%
  inner_join(peps, by = "Peptide")

print("Identify proteins with no peptides")
proteins_w_unique_peps <- elut_wide_identified %>%
  select(ProteinID, experiment_name) %>%
  unique

print("Remove proteins with no unique peptides")
# Each protein has to have at least one unique peptide
  elut_wide_identified <- elut_wide_identified %>%
    inner_join(proteins_w_unique_peps, by = c("ProteinID", "experiment_name"))

print("Convert wide to tidy")
elut_long <- elut_wide_identified %>%
  pivot_longer(cols = starts_with("fractionid_"),
               names_to = "FractionID",
               values_to = "pepcount") %>%
  mutate(pepcount = as.numeric(pepcount)) %>%
  left_join(fraction_order, by = "FractionID") %>%
  mutate(FractionID = str_replace(FractionID, "^fractionid_", "")) %>%
  filter(!is.na(pepcount)) # Unneeded

fraction_order %>% write_csv(paste0(elut_wide_file, "fraction_order.csv"))

elut_long %>% write_csv(paste0(elut_wide_file, ".annot.long.tidy"))

total_fractions <- nrow(fraction_order)


#10x smaller
elut_short <- elut_long %>% 
  filter(pepcount > 0) %>%
  mutate(totfracts = total_fractions)

elut_short %>% write_csv(paste0(elut_wide_file, ".annot.short.tidy"))




