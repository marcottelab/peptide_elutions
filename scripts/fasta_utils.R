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



complete_fxn <- function(df_sel, df_meta){
# Complete missing observations of peptides in fractions


  df_complete <- df_sel %>%
           group_by(ExperimentName) %>%
              bind_rows(df_meta) %>% #adding on so we have all the Fraction IDs

              complete(Peptide, nesting(FractionID, ExperimentName, ExperimentName_order, ExperimentID), fill = list(PeptideCount = 0.0)) %>%
              ungroup %>%
           filter(!is.na(Peptide))

   return(df_complete)

}


# This can be done better?
cov_columns2 <- function(df, peptide_column="Peptide", sequence_column="Sequence"){

    df$pepregex <- mapply(gsub, pattern = "[I|J|L]",
                       replacement = "(I|J|L)", df[,peptide_column]) %>% as.vector()
    #create Length vector, to be appended to he original data frame afterwards

    df <- df %>% mutate(Length = str_length(Sequence))

    getStartorEnd <- function(peptide_column, sequence_column, s=1){

        #s=1 for start position, s=2 for end position.
        comparison <- stringr::str_locate_all(sequence_column[1],coll(peptide_column[1], ignore_case = TRUE))[[1]][,s][1]
        return(comparison)
    }
    df$Start <- mapply(getStartorEnd, df$pepregex, df$Sequence, 1)
    df$End <- mapply(getStartorEnd, df$pepregex, df$Sequence, 2)

    df$pepregex <- NULL

    return(df)
}

# This one is ~0.2 ms faster
# Deprecated
# Now calculating during trypsin.py step (-p TRUE)
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

















