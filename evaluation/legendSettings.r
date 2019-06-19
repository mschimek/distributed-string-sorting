
colours = c("kurpicz" = "coral3", "hQuick" = "burlywood4", "MS_LCP_S" = "cornflowerblue", "MS_NoLCP_S" = "black", "PD_NoGolomb" = "darkolivegreen4",
                                                   "MS_LCP_C" = "darkorange2", "MS_NoLCP_C" = "darkgoldenrod3", "PD_SeqGolomb" = "cyan4", "MS_Simple_S" = "aquamarine4")
shapes = c("kurpicz" = 1, "hQuick" = 2, "MS_LCP_S" = 3, "MS_NoLCP_S" = 4, "PD_NoGolomb" = 5, 
                                                   "MS_LCP_C" = 6, "MS_NoLCP_C" = 7, "PD_SeqGolomb" = 2, "MS_Simple_S" = 1)
linetypes = c("kurpicz" = "solid", "hQuick" = "twodash", "MS_LCP_S" = "dotted", "MS_NoLCP_S" = "dashed", "PD_NoGolomb" = "longdash",
              "MS_LCP_C" = "dotdash", "MS_NoLCP_C" = "dotted", "PD_SeqGolomb" = "twodash", "MS_Simple_S" = "solid")


addSettings <- function(plot) {
  plot <- plot + scale_colour_manual(values = colours) 
  plot <- plot + scale_shape_manual(values = shapes)
  plot <- plot + scale_linetype_manual(values = linetypes)
  return(plot)
}

addBox <- function(plot) {
  plot <- plot + theme(legend.box.background = element_rect(colour = "black"))
  return(plot)
}

getLegend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
