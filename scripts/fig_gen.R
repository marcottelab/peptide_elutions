library(tidyverse)
library(cowplot)
library(broom)
library(patchwork)
library(argparse)

source("scripts/theme_utils.R")
source("scripts/saved_ggplots.R")


parser <- ArgumentParser(description='Make many sparklines')
parser$add_argument('--tidy_psms', dest='tidy_psms', action='store', required=TRUE,
    help='File of complete PSMs in tidy form, output of Rscript peptide_identification.R')
parser$add_argument('--exp', dest='exps', action='store',
    help='Filename of experiment_names, one per line')


parser$add_argument('--IDs', dest='id_file', action='store',
    help='Filename of IDs, one per line')
parser$add_argument('--ID', dest='ID', action='store',
    help='One ID')
parser$add_argument('--rel_widths', action='store', default=10,
    help='How much space to give panels vs peptide ladder 1:x')
parser$add_argument('--height', action = 'store', default=2,
    help='Figure height')
parser$add_argument('--width', action = 'store', default=15,
    help='Figure width')


args = parser$parse_args()


palette <- rep(c("#0072B2","#E69F00","#009E24","#FF0000", "#979797","#5530AA"), 8)


theme_set(theme_cowplot_consistent_text(font_size = 8))

helper2 <- function(elut_short, proteome, protein = NULL, exact = FALSE, exps = NULL, rel_widths = c(1,4), figwidth = 12, figheight = 4, units = "in", outfile_base = "images2020/", axis_tick_spacing = 400, image_type = "png"){
  
  
  if(!is.null(protein)){
    if(exact==TRUE){
      
       elut_short <- elut_short %>% 
         filter(ProteinID == protein)
    }
    else {
       elut_short <- elut_short %>% 
         filter(grepl(protein, ProteinID))
    }
   }
  
   heatmap_plot <- elut_short %>%
   pep_heatmap_fxn(.) + theme(axis.text.x = element_text(angle = 0)) 


   coverage_plot <- elut_short %>%
   pepcov_fxn(.)

    
     


    
    final <-
      heatmap_plot +   
      coverage_plot + 
      plot_layout(widths = rel_widths) 

    orthogroup <- elut_short %>%
      


    filename <- paste0(outfile_base, image_type)
   
    print(filename)
    
    final %>% ggsave(filename, .,device = image_type, width = figwidth, height = figheight,  units = units)
 
    #final %>% ggsave(pdf_filename, ., device = "pdf", width = figwidth, height = figheight, units = units)
    return(final)
}


elut_long <- read_csv(args$tidy_psms)
proteome <- read_csv(args$proteome)
protein <- args$ID
outfile <- paste0("figures/genfig_", protein)
width <- as.integer(args$width)
height <- as.integer(args$height)
rel_widths <- c(1, as.integer(args$rel_widths))
print(width)
print(height)
print(rel_widths)


helper2(elut_long, proteome,  protein = protein, exact = FALSE,
                 #exps = c( "HEMO_SEC1", "HEK_SEC1", "HEKR_IEX1", "HEMO_IEX1", "HEMO_IEX2"), 
                 rel_widths = rel_widths, figwidth = width, figheight = height, units = "in", outfile_base = outfile, axis_tick_spacing = 400 )
