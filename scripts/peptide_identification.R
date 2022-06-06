library(tidyverse)
library(argparse)


parser <- ArgumentParser(description='Identify peptides in an experiment, after running peptide_preprocessing.R for a species')

parser$add_argument('--elut_wide', dest='elut_wide', action='store', required=TRUE,
    help='Wide elution table')

parser$add_argument('--peps', dest='peps', action='store', required=TRUE,
    help='Output of peptide_preprocessing.R')

parser$add_argument('--protein_lengths', dest='protein_lengths', action='store', required=TRUE,
    help='How long each protein is, two column csv, header ProteinID,seqlen')
parser$add_argument('--exp_meta', dest='exp_meta', action='store', required=TRUE,
    help='csv containing ExperimentID, experiment_name')


args = parser$parse_args()

elut_wide_file <- args$elut_wide
peps <- read_csv(args$peps)
protein_lengths <- read_csv(args$protein_length)
exp_meta <- args$exp_meta

#elut_wide_file <-  "data/Anna_HEK293T_IEX_mixed_bed_RNASE_050117.pepcount"
#peps <-  read_csv("data/human_annotated_peptides.csv")
#exp_meta<- "annotation_files/experiment_meta.csv"
#protein_lengths <- read_csv("data/human_protein_lengths.csv")

exp_meta <- read_csv(exp_meta) %>%
            rowid_to_column(var = "experiment_order")
print("Read elution")

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
  left_join(exp_meta, by = "ExperimentID") %>%
  inner_join(peps, by = "Peptide") %>%
  arrange(Start, End) %>%
  left_join(protein_lengths, by = "ProteinID")

print("Identify proteins with no peptides")
proteins_w_unique_peps <- elut_wide_identified %>%
  filter(status == "protein_unique") %>%
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

#total_fractions <- elut_long %>%
#  select(experiment_name, FractionID) %>%
#  unique %>%
#  group_by(experiment_name) %>%
#     summarize(totfracs = max(FractionID))

#fraction_order %>% write_csv(paste0(elut_wide_file, "fraction_order.csv"))
total_fractions <- nrow(fraction_order)


#10x smaller
elut_short <- elut_long %>% 
  filter(pepcount > 0) %>%
  mutate(totfracts = total_fractions)
  #left_join(total_fractions, by = "experiment_name")

elut_short %>% write_csv(paste0(elut_wide_file, ".annot.short.tidy"))




