
colours = c("fkss" = "red", "hQuick" = "burlywood4", "ms2LCPS" = "blue", "msLCPS" = "black", "pdNoGolombC" = "darkolivegreen4","pdC" = "deeppink", "pdNoGolombS"="deepskyblue3", "ms2LCPC" = "darkorange1", "msLCPC" = "darkgoldenrod3", "pdS" = "purple", "msBase" = "darkgreen")
#colours = c("fkss" = "coral3", "hQuick" = "burlywood4", "ms2LCPS" = "blue", "msLCPS" = "black", "pdNoGolombC" = "darkolivegreen4","pdC" = "deeppink", "pdNoGolombS"="deepskyblue3", "ms2LCPC" = "darkorange2", "msLCPC" = "darkgoldenrod3", "pdS" = "cyan4", "msBase" = "darkgreen")
shapes = c("fkss" = 1, "hQuick" = 2, "ms2LCPS" = 3, "msLCPS" = 4, "pdNoGolombS" = 5, 
                                                   "ms2LCPC" = 6, "msLCPC" = 7, "pdS" = 2, "msBase" = 1, "pdC"=2, "pdNoGolombC" = 3)
linetypes = c("fkss" = "solid", "hQuick" = "twodash", "ms2LCPS" = "solid", "msLCPS" = "dashed", "pdNoGolombS" = "longdash",
              "ms2LCPC" = "dotdash", "msLCPC" = "solid", "pdS" = "twodash", "msBase" = "solid", "pdNoGolombC" = "solid", "pdC" = "longdash")


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
