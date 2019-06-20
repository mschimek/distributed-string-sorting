library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)
library(latex2exp)
library(cowplot)
library(gridExtra)
library(xtable)
library(dplyr)
source("legendSettings.r")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

colTypeSpec = cols(numberProcessors = col_integer(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())



pathToJSON = args[[1]]
pdfname = args[[2]]
jsonMetaObject <- read_json(pathToJSON, simplifyDataFrame=TRUE)
jsonObject <- jsonMetaObject["data"][[1]]
title <- jsonMetaObject["title"]
isD2N <- TRUE == (jsonMetaObject["isDToN"])
print(str(jsonObject))
length <- length(jsonObject)
data <- vector("list", length)
filters <- vector("list", length)
print(paste("number datasets: ", length))
for (i in c(1:length)){
  print(i)
  curPath = jsonObject[[i]]$path
  filename = paste(curPath, "/data.txt", sep="")
  data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
  data[[i]] <- filter(data[[i]], iteration != 0)
  filters[[i]] <- jsonObject[[i]]$filter
  print(filters[[i]])
}
numberExpandedDatasets <- function() {
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  return(counter)
}

lineplot <- function(datasets, operation_, type_ = "maxTime", title = " ", work = FALSE) {
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  print(paste("counter: ", counter))
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    print(paste("numberFilters: ", numberFilters))
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      localSets[[counter]] <- filter(data[[i]], operation == operation_, type == type_)
      #print(localSets[[counter]])
      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      #print(localSets[[counter]])
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], numberProcessors, dToNRatio,  operation, type, name, iteration, value)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
    }
    set$numberProcessors <- as.factor(set$numberProcessors)
    set$value <- set$value / 10^9
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  operation, type, name)
    valueMean <- summarise(groupByIterations, sd = sd(value), value = mean(value, rm.na = TRUE))
    
  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, group = name, colour = name, shape = name, linetype = name))
  plot <- plot + ylab("time [sec]")
  plot <- plot + xlab("PEs")
  plot <- plot + theme_light()
  plot <- plot + geom_point(position=position_dodge(width=0.1))
  plot <- plot + geom_line(position=position_dodge(width=0.1))

  #plot <- plot + theme(legend.direction = "horizontal")
plot <- plot + theme(legend.box.background = element_rect(colour = "black"))
plot <- plot + theme(legend.position = "bottom")
  #plot <- plot + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #                 position=position_dodge(.1))

  #plot <- plot + scale_y_continuous(trans='log2')
  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #plot <- plot + theme(legend.position = "none")
  if (isD2N) {
  plot <- plot + facet_wrap(~dToNRatio, labeller = label_both, nrow=1)
  } 
  plot <- plot + ggtitle(title)
  return (plot)

}

lineplotMemory <- function(datasets, type_ = "number", title = " ", work = FALSE) {
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  print(paste("counter: ", counter))
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    print(paste("numberFilters: ", numberFilters))
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      localSets[[counter]] <- filter(data[[i]], type == type_, rawCommunication == 1, phase != "none")
      #print(localSets[[counter]])
      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      #print(localSets[[counter]])
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], numberProcessors, dToNRatio,  phase, name, iteration, value, size)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
    }
    set$value <- set$value / (set$size * set$numberProcessors) 
    set$numberProcessors <- as.factor(set$numberProcessors)
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  name, phase)
    valueMean <- summarise(groupByIterations,  value = mean(value, rm.na = TRUE))
    value <- ungroup(valueMean)
    valueMean <- group_by(valueMean, numberProcessors, dToNRatio,  name)
    valueMean <- summarise(valueMean, value = sum(value)) 
    print(valueMean,n = 50)
    
  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, group = name, colour = name, shape = name, linetype = name))
  plot <- plot + ylab("sent bytes per string")
  plot <- plot + xlab("PEs")
  plot <- plot + theme_light()
  plot <- plot + geom_point(position=position_dodge(width=0.1))
  plot <- plot + geom_line(position=position_dodge(width=0.1))

  #plot <- plot + theme(legend.direction = "horizontal")
plot <- plot + theme(legend.box.background = element_rect(colour = "black"))
plot <- plot + theme(legend.position = "bottom")
  #plot <- plot + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #                 position=position_dodge(.1))

  #plot <- plot + scale_y_continuous(trans='log2')
  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #plot <- plot + theme(legend.position = "none")
  if (isD2N) {
  plot <- plot + facet_wrap(~dToNRatio, labeller = label_both, nrow=1)
  } 
  plot <- plot + ggtitle(title)
  return (plot)

}

stackedBarPlot <- function(datasets, dToNRatio_, operations_, type_ = "maxTime", title = " ", work = FALSE) {
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  print(paste("counter: ", counter))
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    print(paste("numberFilters: ", numberFilters))
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      if (isD2N) {
        localSets[[counter]] <- filter(data[[i]], operation %in% operations_, type == type_, dToNRatio == dToNRatio_)
      } else {
        localSets[[counter]] <- filter(data[[i]], operation %in% operations_, type == type_)
      }
      #print(localSets[[counter]])
      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      #print(localSets[[counter]])
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], numberProcessors, dToNRatio,  operation, type, name, iteration, value)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
    }
    set$numberProcessors <- as.factor(set$numberProcessors)
    set$value <- set$value / 10^9
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  operation, type, name)
    valueMean <- summarise(groupByIterations , value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean, mapping = aes(x = name, y = value, fill=operation))
  plot <- plot + ylab("time [sec]")
  plot <- plot + xlab("")
  plot <- plot + theme_light()
  plot <- plot + geom_bar(stat="identity")

  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(plot.subtitle=element_text(hjust=0.5, size=12))
  plot <- plot + facet_wrap(~numberProcessors, labeller = label_value, nrow=1) 
  plot <- plot + ggtitle(title)
  plot <- plot + labs(subtitle = "number of PEs")
  return (plot)

}
stackedBarPlotMemory <- function(datasets, dToNRatio_, type_ = "maxTime", title = " ", work = FALSE) {
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  print(paste("counter: ", counter))
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    print(paste("numberFilters: ", numberFilters))
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      if (isD2N) {
        localSets[[counter]] <- filter(data[[i]], type == type_, dToNRatio == dToNRatio_, rawCommunication == 1, phase != "none")
      } else {
        localSets[[counter]] <- filter(data[[i]],  type == type_, rawCommunication == 1, phase != "none")
      }
      #print(localSets[[counter]])
      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      #print(localSets[[counter]])
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], numberProcessors, phase, dToNRatio, name, iteration, value)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
    }
    set$numberProcessors <- as.factor(set$numberProcessors)
    set$value <- set$value / 10^9
    print(set)
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  name, phase)
    valueMean <- summarise(groupByIterations , value = mean(value))
    valueMean <- ungroup(valueMean)
    valueMean <- select(valueMean, numberProcessors, dToNRatio, name, value, phase)
    valueMean <- ungroup(valueMean)
    print(valueMean, n = 100 )
    valueMean <- group_by(valueMean, dToNRatio, numberProcessors,  name)
    valueMean <- summarise(valueMean , value = sum(value))
    print(valueMean, n = 100 )

  plot <- ggplot(data = valueMean, mapping = aes(x = name, y = value, fill=name))
  plot <- plot + ylab("sent bytes per string [GB]")
  plot <- plot + xlab("")
  plot <- plot + theme_light()
  plot <- plot + geom_bar(stat="identity")

  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(plot.subtitle=element_text(hjust=0.5, size=12))
  plot <- plot + facet_wrap(~numberProcessors, labeller = label_value, nrow=1) 
  plot <- plot + ggtitle(title)
  plot <- plot + labs(subtitle = "number of PEs")
  return (plot)

}
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

pdf(paste("./plots/", pdfname, ".pdf",sep=""), width=10, height=9)
operations = c("sort_splitter","prefix_decompression", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "sort_locally", "bloomfilter_overall")
l <- lineplot(c(1:length(data)), "sorting_overall", "maxTime", title)
s <- stackedBarPlot(c(1:length(data)), dToNRatio_ = 0.5, operations_ = operations, "maxTime", title)
m <- stackedBarPlotMemory(c(1:length(data)), dToNRatio_ = 1.0,  "number", title)
print("memline")
memLine <- lineplotMemory(c(1:length(data)),  "number", title)
memLine <- addSettings(memLine)
l <- addSettings(l)
l <- l + theme(legend.direction = "horizontal")
l <- l + theme(legend.box.background = element_rect(colour = "black"))
l <- l + theme(legend.title = element_blank()) 
l
legend <- get_legend(l)
l <- l + theme(legend.position = "none")
memLine <- memLine + theme(legend.position = "none")
memLine <- memLine + ggtitle("")

top <- plot_grid(l,ncol = 1)
middle <- plot_grid(memLine, ncol=1)
legend <- plot_grid(NULL, legend, NULL, ncol=3, rel_widths=c(0.25, 0.5, 0.25))
plot_grid(top, middle,   legend, nrow = 3, rel_heights = c(1, 1, 0.2))
#plot_grid(s, l, labels = c("s", "l"), ncol = 2, nrow = 1)
#grid.arrange(s,l, legend, ncol=2, nrow = 2,layout_matrix = rbind(c(1,2), c(3,3)),
#widths = c(2.7, 2.7), heights = c(2.2, 0.4))
#if (isD2N) {
#  for (d in unique(data[[i]]$dToNRatio)) {
#    print(stackedBarPlot(c(1:length(data)), dToNRatio_ = d, operations_ = operations, "maxTime", paste(title, " with D/N-ratio: ", d)))
#    print(stackedBarPlotMemory(c(1:length(data)), dToNRatio_ = d, "number", paste("Communication Volume: ", title, " with D/N-ratio: ", d)))
#  }
#} 


