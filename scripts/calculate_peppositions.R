library(tidyverse)

read_fasta <- function(fasta_filename, annot = FALSE){
    fasta <- seqinr::read.fasta(fasta_filename, as.string = TRUE)

    # Convert seqinr SeqFastadna object to data.frame
    fasta_df <- fasta %>%
                   sapply(function(x){x[1:length(x)]}) %>%
                   as.data.frame %>%
                   broom::fix_data_frame(newcol = "ID", newnames = "Sequence")

    if(annot == TRUE){
        annot_df <- getAnnot(fasta) %>%
                         sapply(function(x){x[1:length(x)]}) %>%
                         as.data.frame() %>%
                         broom::fix_data_frame(newnames = "Annot")

        fasta_df <- cbind(fasta_df, annot_df)
    }
    return(fasta_df)
}


# This one is ~0.2 ms faster
cov_columns <- function(df, peptide_column = Peptide, sequence_column = Sequence){
  
  # Replace with regex for redundant I/L in peptide
  # Get start and end position of peptide in protein sequence
  # Separate out start and end to two columns
  # Add on a column for length of the protein
  
  df_col <- df %>% 
    mutate(pepregex = toupper(str_replace_all({{ peptide_column }}, "[I|J|L]", "(I|J|L)"))) %>% 
    rowwise() %>%
        mutate(positions = paste0(str_locate(toupper({{ sequence_column }}), pepregex ), collapse=",")) %>% 
    
    separate(positions, into = c("start", "end"), remove = TRUE) %>% 
    mutate(start = as.integer(start), end = as.integer(end)) %>% #not necessary if could do unlist better
    select(-pepregex) %>%
    mutate(length = str_length({{ sequence_column }})) 
  
  return(df_col)
}

pep_data_path <- "protein_identification/peptides"
# This for getting order info
protein_all_pep_filenames <- dir(pep_data_path, pattern = "*peptides.csv") # Get filenames
protein_all_pep_filenames <-protein_all_pep_filenames[!grepl("unique", protein_all_pep_filenames)]
protein_all_pep_filenames <-protein_all_pep_filenames[!grepl("virNOG", protein_all_pep_filenames)]


protein_all_peptides <- tibble(filename = protein_all_pep_filenames) %>%
               mutate(file_contents = map(filename, ~ read_csv(file.path(pep_data_path, .)))) %>%
               unnest() %>%
               mutate(Peptide = str_replace_all(Peptide, "[I|L]", "J")) %>%
               mutate(spec = str_extract(filename,"^.....")) %>%
               select(-filename)



# Please work

data_path <- "protein_identification/proteomes"
proteome_files <- dir(data_path, pattern = "*.fasta$") # Get filenames
proteomes <- tibble(filename = proteome_files) %>%
               mutate(file_contents = map(filename, ~read_fasta(file.path(data_path, .)))) %>%
  
               unnest()

protein_all_peptides_sequences <- protein_all_peptides %>%
  left_join(proteomes, by = c("ProteinID"= "ID"))

print("Finding peptide positions")
peptide_locations <- cov_columns(protein_all_peptides_sequences)
peptide_locations %>% write_csv("peptide_locations.csv")
