library(dplyr)
library(tidyr)
library(readr)
library(broom)
library(purrr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the data
data <- read_csv(input_file)

# Your analysis here (adjust as needed)
result <- data %>%
  group_by(ProteinID, experiment_name, peak, midpoint) %>%
  nest() %>%
  mutate(test_result = map(data, ~broom::tidy(fisher.test(xtabs(n ~ pep_half + presence, data = .x))))) %>%
  unnest(test_result)

# Save the result
write_csv(result, output_file)

