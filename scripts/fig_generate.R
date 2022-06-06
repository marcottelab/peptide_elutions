library(tidyverse)
library(ggrepel)
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
#parser$add_argument('--proteome', dest='proteome', action='store',
#                    help='Two column proteome, with')

parser$add_argument('--domaindis', dest='domaindis', action='store',
                    help='Annotation file of domain and disordered locations')

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

helper2 <- function(elut_long, domaindis, protein = NULL, exact = FALSE, exps = NULL, rel_widths = c(1,4), figwidth = 12, figheight = 4, units = "in", outfile_base = "images2020/", axis_tick_spacing = 400, image_type = "png", rescale_pepcounts = TRUE){
  
  width <- 5
  rel_widths <- c(1,4,1)
  rel_widths <- c(1, as.integer(args$rel_widths))
  domaindis <- read_tsv("data_files/domaindis_setup.txt")
  elut_long <- read_csv("data_files/short_tidy_unique.csv")
  #protein  <- "A0A2K3CN31_CHLRE"  
  #protein <- "TRIPB_HUMAN"
  #protein <- "Y4554_ARATH"
  #protein <- "CDC73_HUMAN"
  exact <- FALSE
  image_type <- "png"
  
  

    
  
  if(!is.null(protein)){
    if(exact==TRUE){
      
       elut_short <- elut_long %>% 
         filter(ProteinID == protein)
    }
    else {
       elut_short <- elut_long %>% 
         dplyr::filter(grepl({{ protein }}, ProteinID))
    }
  }
   orthogroup <- elut_short$ID[[1]]
   entry_name <- elut_short %>% separate(ProteinID, into = c(NA, NA, "Entry_name"), sep= "[|]") %>%
    pull(Entry_name) %>% unique()
  
   
   if(rescale_pepcounts == TRUE){
     elut_short <-elut_short %>%
       group_by(experiment_name, Peptide) %>%
          mutate(pepcount = pepcount/max(pepcount, na.rm = TRUE)) %>% ungroup()
   }
   #msfragger_cov %>%  filter(grepl(protein, ProteinID)) %>% 
   #  select(filename, Peptide, start, end )
   
   
   
   
   
   domaindis_short <- domaindis %>% filter(Entry_name == {{entry_name}})
  
   
   
   print("make_plots")
   heatmap_plot <- elut_short %>%
      pep_heatmap_fxn(., pepmod = FALSE) + theme(axis.text.x = element_text(angle = 0)) 
  
   domaindis_plot <- domaindis_short %>%
     plot_domain_dis(.)

   coverage_plot <- elut_short %>%
            pepcov_fxn(., axis_tick_spacing = 400)

    print("plotsmade")
    
    layout <- "
    ##C
    ABC
    ##C
    "
    
    final <-    
      coverage_plot + 
      heatmap_plot+ 
      domaindis_plot +
      plot_layout(design = layout, widths =rel_widths, heights = c(0.0001,1,0.0001)) 

    filename <- paste0(outfile_base, "-", orthogroup, "-", entry_name, ".", image_type)
   
    print(filename)
    
    final %>% ggsave(filename, .,device = image_type, width = figwidth, height = figheight,  units = units)
    #final %>% ggsave(filename, .,device = image_type, width = 12, height = 4,  units = units)
    
    #final %>% ggsave(pdf_filename, ., device = "pdf", width = figwidth, height = figheight, units = units)
    return(final)
}

#Rscript scripts/fig_generate.R 
#--ID A0A2K3CN31_CHLRE 
#--proteome data_files/chlre_proteome.csv 
#--tidy_psms data_files/elut_long.csv 
#--rel_widths 4 
#--width 5



elut_long <- read_csv(args$tidy_psms)
#proteome <- read_csv(args$proteome)
protein <- args$ID
domaindis <- read_tsv(args$domaindis)
outfile <- "genfig"
width <- as.integer(args$width)
height <- 4
rel_widths <- c(1, as.integer(args$rel_widths),1)
print(width)
print(height)
print(rel_widths)


#Figures_to_make

helper2(elut_long, protein = protein, domaindis = domaindis, exact = FALSE,
                 #exps = c( "HEMO_SEC1", "HEK_SEC1", "HEKR_IEX1", "HEMO_IEX1", "HEMO_IEX2"), 
                 rel_widths = rel_widths, figwidth = width, figheight = height, units = "in", outfile_base = outfile, axis_tick_spacing = 400 )
