theme_cowplot_consistent_text <- function (font_size = 8) {

  theme_cowplot() %+replace%
       theme(strip.text = element_text(size = font_size),
             axis.text = element_text(colour = "black", size = font_size),
             plot.title = element_text(size = font_size),
             axis.title = element_text(size = font_size),
             legend.text = element_text(size = font_size),
             legend.title = element_text(size = font_size),
             legend.key.size = unit(0.5, "lines"),
             axis.title.x = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), vjust = 1),
             axis.text.x = element_text(margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0), vjust = 1),
             axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), angle = 90, vjust = 1),
             axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0), hjust = 1)
            )

}


theme_sparkline_peps <- function () {
    theme_cowplot_consistent_text() %+replace%
        theme(
           #axis.text.x = element_blank(),
           axis.text.y = element_text(vjust = -1),
          # axis.ticks.x = element_blank(),
           axis.line.x = element_blank(),
           #axis.title.x = element_blank(),
           axis.title.y = element_blank(),

           panel.spacing = unit(0.5, "lines"),
           panel.background = element_rect(fill = "white"),
           plot.background=element_rect(fill="grey90"),
           ##strip.text = element_text(angle = 45, hjust = 0),
           #strip.text = element_blank(),
           strip.background = element_rect(fill = NA, color = "black"),
           axis.line.y = element_blank(),
           plot.margin = unit(c(0, 0, 0, 0), "cm")#,
           #panel.border = element_rect(color = "grey85", fill = NA, linetype = 1, size =0.5)#,
            # axis.line.x.top = element_line(color = "red", size = 0.5, linetype = 1)
          )
}
