
colours = c("fkss" = "coral3", "hQuick" = "burlywood4", "msLCPS" = "cornflowerblue", "msNoLCPS" = "black", "pdNoGolombC" = "darkolivegreen4","pdGolombC" = "deeppink", "pdNoGolombS"="deepskyblue3", "msLCPC" = "darkorange2", "msNoLCPC" = "darkgoldenrod3", "pdGolombS" = "cyan4", "msSimpleS" = "aquamarine4")
shapes = c("fkss" = 1, "hQuick" = 2, "msLCPS" = 3, "msNoLCPS" = 4, "pdNoGolombS" = 5, 
                                                   "msLCPC" = 6, "msNoLCPC" = 7, "pdGolombS" = 2, "msSimpleS" = 1, "pdGolombC"=2, "pdNoGolombC" = 3)
linetypes = c("fkss" = "solid", "hQuick" = "twodash", "msLCPS" = "dotted", "msNoLCPS" = "dashed", "pdNoGolombS" = "longdash",
              "msLCPC" = "dotdash", "msNoLCPC" = "dotted", "pdGolombS" = "twodash", "msSimpleS" = "solid", "pdNoGolombC" = "solid", "pdGolombC" = "dotted")


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
