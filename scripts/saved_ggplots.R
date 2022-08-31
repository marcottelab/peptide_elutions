
# This function is for plotting one protein at a time
b_w_plot <- function(short_tidy_unique_sel, pubtheme = TRUE){
  
  numprots <-  short_tidy_unique_sel %>% pull(ProteinID) %>% unique %>% length()
  print(numprots)
  plt <- short_tidy_unique_sel %>%
    group_by(Peptide, experiment_name) %>%
    mutate(pepcount = pepcount/ max(pepcount)) %>%
    ungroup %>%
    ggplot(aes( x = FractionOrder , y = fct_rev(fct_reorder(Peptide, Start)), fill = pepcount)) +
    
    geom_tile() +
    geom_blank(aes(x = totfracs)) +  
    facet_grid(ProteinID  ~ experiment_name , scales = "free_y", drop = FALSE) +
    scale_fill_gradient(low = "#FFFFC6", high = "#3F1CC6", na.value = '#FFFFC6' ) +
    scale_x_continuous(limits = c(0, NA),breaks = seq(0, 1000, 30) ) +
    scale_y_discrete(expand = c(0,0)) + #This gets rid of the bottom grey space to help with alignment
    
    theme(panel.background = element_rect(fill = "#FFFFC6"),
          panel.spacing = unit(0.2, "lines"), 
          axis.line.x = element_line(color = "black"),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),  
          axis.text.x = element_text(color = "black", angle = 45,vjust =1 , hjust = 1 ),
          legend.position = "none",
          strip.background = element_rect(fill = "white"))
  
  if(pubtheme == TRUE){
    plt <- plt +
      theme_nothing() + 
      theme( panel.background = element_rect(fill = "#FFFFC6"),
             panel.spacing = unit(0.25, "lines"),
             strip.text = element_text(), 
             plot.margin = unit(margin(0.5,0.5,0.5,0.5), "lines"))
    
  }
  
  return(plt)
  
}

pep_heatmap_fxn <- function(elut_short, pepmods = FALSE){

  hm_plot <- elut_short %>%
    arrange(Start, End) %>%
  ggplot(aes( x = FractionOrder , y =fct_rev(fct_inorder(Peptide)))) +
    geom_tile(aes(fill = pepcount)) +
    #geom_tile(fill = "#6F236B") +
    geom_blank(aes(x = totfracs))+
    facet_wrap(~experiment_order + experiment_name, nrow = 1, scales = "free_x") +

    scale_fill_gradient(low = "#bea9c2", high = "#073bf2", na.value= "#FFFFBF") +
    scale_x_continuous(limits = c(0, NA),breaks = seq(0, 1000, 30), expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) + #This gets rid of the bottom grey space to help with alignment

    theme(panel.background = element_rect(fill = "#FFFFBF", colour = "black"),
          panel.spacing = unit(0.2, "lines"),
          #axis.line.x = element_line(color = "black"),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
         # panel.border = element_line(colour = "black"),
          axis.text.x = element_text(color = "black", angle = 45,vjust =1 , hjust = 1 ),
          legend.position = "none",
          strip.background = element_rect(fill = "white"))

  if(pepmods == TRUE){
    
    hm_plot <- hm_plot + 
      geom_point(aes(shape = source)) +
      scale_shape_manual(values = c(1,4)) +
      theme(legend.position = "none")
    
    
  }
  
  return(hm_plot)

}

pep_point_fxn <- function(elut_short, pepmods = TRUE){
  
  hm_plot <- elut_short %>%
    arrange(Start, End) %>%
    ggplot(aes( x = FractionOrder , y =fct_rev(fct_inorder(Peptide)))) +
    geom_point(aes(color = pepcount)) +
    #geom_tile(fill = "#5F236B") +
    geom_blank(aes(x = totfracs))+
    facet_wrap(~experiment_order + experiment_name, nrow = 1, scales = "free_x") +
    
    scale_color_gradient(low = "#bea9c2", high = "#073bf2", na.value= "#FFFFBF") +
    scale_x_continuous(limits = c(0, NA),breaks = seq(0, 1000, 30), expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) + #This gets rid of the bottom grey space to help with alignment
    
    theme(panel.background = element_rect(fill = "#FFFFBF", colour = "black"),
          panel.spacing = unit(0.2, "lines"),
          #axis.line.x = element_line(color = "black"),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          # panel.border = element_line(colour = "black"),
          axis.text.x = element_text(color = "black", angle = 45,vjust =1 , hjust = 1 ),
          legend.position = "none",
          strip.background = element_rect(fill = "white"))
  
  if(pepmods == TRUE){
    
    hm_plot <- hm_plot + 
      geom_point(aes(shape = source)) +
      scale_shape_manual(values = c(1,4, 6)) +
      theme(legend.position = "none")
    
    
  }
  
  return(hm_plot)
  
}


plot_domain_dis <- function(domaindis_setup ){

  domaindis_setup$set <- factor(domaindis_setup$set, levels = c("Disordered", "Domains"))
  
  
  domain_dis_plot <- domaindis_setup %>%
  ggplot() +
  geom_rect(fill = "grey90", aes(xmin = 0, xmax = seqlen, ymin = 0.1, ymax = 0.8)) +
  geom_rect(aes(xmin = ProtBegin, xmax = ProtEnd, ymin = ymin, ymax = ymax, fill = set)) +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
 # scale_fill_manual(values = c("orange", "cadetblue4")) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  scale_x_reverse(expand = c(0,0)) +
  #scale_y_continuous()
  geom_text_repel(y = 1, ylim = c(1,15), size = 2.5, aes(x = ProtBegin, label = Annotation)) +
  theme(legend.position = "none") +
    scale_fill_manual(values = c("Disordered" = "orange", "Domains" = "cadetblue4", "Transmembrane" = "purple"))
  

  return(domain_dis_plot)

}




pepcov_fxn <- function(elut_short, axis_tick_spacing = 100, min_graphic_details = FALSE, pubtheme = TRUE){
  
  pepcov <- elut_short %>%
    
    select(Peptide, Start, End, seqlen) %>%
    
    unique %>%
    arrange(Start, End) %>%
    mutate(label = paste0(Start, "-", End, "-", Peptide)) %>%
    ggplot(aes(fill = "grey90", x = fct_rev(fct_inorder(label)))) +
    geom_boxplot(stat = "identity", aes(ymin = Start, ymax = End, lower = Start, upper = End, middle = Start), size = 0) +
    scale_fill_manual(values = c("black", "grey50"))  +
    theme(strip.background = element_rect(fill = "white")) +
    scale_x_discrete(position = "top") + #Really y right
    coord_flip()  +
    background_grid(major = "xy", minor = "none", size.major = 0.8) + # Really x
    scale_y_continuous(breaks = seq(0, 100000, axis_tick_spacing), labels = seq(0, 100000, axis_tick_spacing), minor_breaks = seq(100, 100000, 100), expand = c(0,0)) + # Really x
    expand_limits(y = 0) + #Really x
    facet_wrap(~"") +
    
    theme(axis.text.y = element_text(size = 6), # Really y
          axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),
          axis.title.x =  element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),

          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          
          
          legend.position = 'none')
  
  if(min_graphic_details == TRUE){
    pepcov <- pepcov + theme_nothing()
    
  }
  
  if(pubtheme == TRUE){
    pepcov <- pepcov + theme(axis.text.y = element_blank(), # Really y
                             axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),
                             axis.title.y = element_blank())
    
  }
  
  return(pepcov)
  
}



pepcov_fxn_lightly_deprecate <- function(elut_short, axis_tick_spacing = 100, min_graphic_details = FALSE, pubtheme = TRUE){

  pepcov <- elut_short %>%

   select(Peptide, Start, End, seqlen, status) %>%

   unique %>%
    mutate(label = paste0(Start, "-", End, "-", Peptide)) %>%
   ggplot(aes(fill = status, x = fct_rev(fct_inorder(label)))) +
      geom_boxplot(stat = "identity", aes(ymin = Start, ymax = End, lower = Start, upper = End, middle = Start), size = 0) +
  scale_fill_manual(values = c("black", "grey50"))  +
  theme(strip.background = element_rect(fill = "white")) +
          scale_x_discrete(position = "top") + #Really y right
          coord_flip()  +
          background_grid(major = "xy", minor = "none", size.major = 0.8) + # Really x
          scale_y_continuous(breaks = seq(0, 100000, axis_tick_spacing), labels = seq(0, 100000, axis_tick_spacing), minor_breaks = seq(100, 100000, 100), expand = c(0,0)) + # Really x
          expand_limits(y = 0) + #Really x
          facet_wrap(~"") +

    theme(axis.text.y = element_text(size = 6), # Really y
        axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),
                axis.title.x =  element_blank(),
        axis.title.y = element_blank(),
                        axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),


                legend.position = 'none')

  if(min_graphic_details == TRUE){
      pepcov <- pepcov + theme_nothing()

  }

  if(pubtheme == TRUE){
      pepcov <- pepcov + theme(axis.text.y = element_blank(), # Really y
                axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),
                axis.title.y = element_blank())

  }

  return(pepcov)

}



# Unused
pepcov_fxn_deprecated <- function(df_cov, axis_tick_spacing = 400, min_graphic_details = FALSE, pubtheme = TRUE){
   pepcov <- ggplot(data = df_cov, aes(x = fct_reorder(Peptide, desc(start)))) +
          geom_boxplot(stat = "identity", fill="black", aes(ymin = start, ymax = end, lower = start, upper = end, middle = start), size = 0) +
          theme(strip.background = element_rect(fill = "white")) +
          scale_x_discrete(position = "top") + #Really y right
          coord_flip()  +
          background_grid(major = "xy", minor = "none", size.major = 0.8) + # Really x
          scale_y_continuous(breaks = seq(0, 100000, axis_tick_spacing), labels = seq(0, 100000, axis_tick_spacing), minor_breaks = seq(100, 100000, 100), expand = c(0,0)) + # Really x
          expand_limits(y = 0) + #Really x
          facet_wrap(~"") + #Force top spacing +
          background_grid(major = "x", minor = "x", colour.major = "grey80", colour.minor = "grey80", size.major = 0.4, size.minor = 0.4) +
          theme(panel.background = element_blank()) +
          ylab("Residue") + #

  theme(axis.text.y = element_text(size = 6), # Really y
        axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),

        plot.margin = unit(c(0, 0, 0, 0), "cm"))

  if(min_graphic_details == TRUE){
      pepcov <- pepcov + theme_nothing()

  }

  if(pubtheme == TRUE){
      pepcov <- pepcov + theme(axis.text.y = element_blank(), # Really y
                axis.text.x = element_text(angle = 45,vjust =1 , hjust = 1 ),
                axis.title.y = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = 'none')
  }

  return(pepcov)

}


# Unused
pep_sparkline_fxn <- function(df_full){

  sparkplot <- ggplot(df_full, aes( x = FractionID_order , y =fct_reorder(Peptide, desc(start)), height = PeptideCount, group = Peptide)) +
   # theme(strip.background = element_rect(fill = "white") +
   geom_ridgeline(fill= "#E69F00", scale = 1, color = "grey50", size = 0.1) +#, min_height = 0.001) +
   #geom_ridgeline(fill= "grey40", scale = 1.5, color = "grey80") +
    theme_sparkline_peps() +
    #theme(axis.ticks.x = element_line(color = "black"))+
    theme(axis.ticks.y = element_blank()) +
    scale_y_discrete(expand = c(0,0)) + #This gets rid of the bottom grey space to help with alignment
    facet_wrap(~fct_reorder(ExperimentName, ExperimentName_order), scales = "free_x", nrow = 1) +
          theme(plot.background =  element_blank()) +
         # panel_border() +
         theme(panel.spacing = unit(0.5, "lines")) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.text.x.top = element_blank()) +
    theme(axis.ticks.x.top = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.x = element_text(colour = "black")) +
    scale_x_continuous(breaks = seq(0, 1000, 20), expand = c(0,0)) +
    theme(axis.text.x = element_text(color = "black", angle = 45,vjust =1 , hjust = 1 )) +
    xlab("Elution position")
  return(sparkplot)

}
