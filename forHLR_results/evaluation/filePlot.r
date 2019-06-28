library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)
library(latex2exp)
library(cowplot)
library(gridExtra)
library(xtable)
library(ggpubr)


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


limit1 <- as.numeric(args[[5]])
limit2 <- as.numeric(args[[6]])
limit3 <- as.numeric(args[[7]])

pathToJSON = args[[1]]
pdfTitle = args[[2]]
plots <- read_json(pathToJSON, simplifyDataFrame=TRUE)

numberPlots <- length(plots)
print(str(plots))
plots.data <- vector("list", numberPlots)
plots.filters <- vector("list", numberPlots)
plots.title <- vector("list", numberPlots)

for (k in c(1:numberPlots)) {
  jsonMetaObject <- plots[[k]]
  print(str(jsonMetaObject))
  jsonObject <- jsonMetaObject["data"][[1]]
  plots.title[[k]] <- jsonMetaObject["title"]
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
  #jsonMetaObject <- plots["CommonCrawlUnique"][[1]]
  print(str(jsonMetaObject))
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
    print(filename)
    data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
    data[[i]] <- filter(data[[i]], iteration != 0)
    filters[[i]] <- jsonObject[[i]]$filter
    print(filters[[i]])
  }
  plots.data[[k]] <- data
  plots.filters[[k]] <- filters
}


divideByPE320 <- function(vec) {
    valueBase <- filter(vec, numberProcessors == 320)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
    vec$value <- (valueBase / vec$value) 
  print(vec)
    return(vec)
}

#jsonMetaObject <- plots["CommonCrawlUnique"][[1]]
#print(str(jsonMetaObject))
#jsonObject <- jsonMetaObject["data"][[1]]
#title <- jsonMetaObject["title"]
#isD2N <- TRUE == (jsonMetaObject["isDToN"])
#print(str(jsonObject))
#length <- length(jsonObject)
#data <- vector("list", length)
#filters <- vector("list", length)
#print(paste("number datasets: ", length))
#for (i in c(1:length)){
#  print(i)
#  curPath = jsonObject[[i]]$path
#  filename = paste(curPath, "/data.txt", sep="")
#  data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
#  data[[i]] <- filter(data[[i]], iteration != 0)
#  filters[[i]] <- jsonObject[[i]]$filter
#  print(filters[[i]])
#}
lineplotSpeedup <- function(data, filters, datasets, operation_, type_ = "maxTime", title = " ", work = FALSE) {
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
    valueMean <- summarise(groupByIterations, value = mean(value, rm.na = TRUE))
    valueMean <- ungroup(valueMean)
    valueMean <- group_by(valueMean, name)
    speedup <- do(valueMean, divideByPE320(.))
    speedup$numberProcessors <- as.factor(speedup$numberProcessors)

  plot <- ggplot(data = speedup, mapping = aes(x = numberProcessors, y = value, group = name, colour = name, shape = name, linetype = name))
  plot <- plot + ylab("speed up")
  plot <- plot + xlab("PEs")
  plot <- plot + theme_light()
  plot <- plot + geom_point(position=position_dodge(width=0.1))
  plot <- plot + geom_line(position=position_dodge(width=0.1))

  #plot <- plot + theme(legend.direction = "horizontal")
plot <- plot + theme(legend.box.background = element_rect(colour = "black"))
  #plot <- plot + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #                 position=position_dodge(.1))

  #plot <- plot + scale_y_continuous(trans='log2')
  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(title)
  return (plot)

}
lineplot <- function(data, filters, datasets, operation_, type_ = "maxTime", title = " ", blacklist = c()) {
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
    set <- filter(set, !(name %in% blacklist))
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
  #plot <- plot + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #                 position=position_dodge(.1))

  #plot <- plot + scale_y_continuous(trans='log2')
  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(strip.background =element_rect(fill="white"))
  plot <- plot + theme(strip.text = element_text(colour = 'black'))
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(title)
  return (plot)

}
stackedBarPlot <- function(data, filters, datasets, dToNRatio_, operations_, type_ = "maxTime", title = " ", work = FALSE) {
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
      isD2N <- FALSE
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
    print(set)
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

lineplotMemory <- function(data, filters, datasets, title = " ", size, blacklist) {
  set <- "dummy"
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      localSets[[counter]] <- filter(data[[i]], type == "number", rawCommunication == 1, phase != "none")
      print(localSets[[counter]])
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
    #set$value <- set$value / (set$size * set$numberProcessors) 
    set$numberProcessors <- as.factor(set$numberProcessors)
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  name, phase)
    valueMean <- summarise(groupByIterations,  value = mean(value, rm.na = TRUE))
    value <- ungroup(valueMean)
    valueMean <- group_by(valueMean, numberProcessors, dToNRatio,  name)
    valueMean <- summarise(valueMean, value = sum(value)) 
    print(valueMean,n = 50)
    divisor <- as.double(size)
  valueMean$value <- valueMean$value / divisor
  valueMean <- filter(valueMean, !(name %in% blacklist))
    
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

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

pdf(paste("./plots/", pdfTitle, ".pdf",sep=""), width=8, height=11)
operations = c("sort_splitter","prefix_decompression", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "sort_locally", "bloomfilter_overall")
#bar <- stackedBarPlot(plots.data[[4]], plots.filters[[4]],  c(1:length(plots.data[[4]])), dToNRatio_ = 0.5, operations_ = operations,  "maxTime", plots.title[[4]])

plot.1 <- lineplot(plots.data[[1]], plots.filters[[1]],  c(1:length(plots.data[[1]])), "sorting_overall", "maxTime", plots.title[[1]])
plot.2 <- lineplot(plots.data[[2]], plots.filters[[2]],  c(1:length(plots.data[[2]])), "sorting_overall", "maxTime", plots.title[[2]])
plot.7 <- lineplot(plots.data[[1]], plots.filters[[1]],  c(1:length(plots.data[[1]])), "sorting_overall", "maxTime", plots.title[[1]], c("kurpicz", "hQuick"))
plot.8 <- lineplot(plots.data[[2]], plots.filters[[2]],  c(1:length(plots.data[[2]])), "sorting_overall", "maxTime", plots.title[[2]], c("kurpicz", "hQuick"))
plot.3 <- lineplotMemory(plots.data[[1]], plots.filters[[1]], c(1:length(plots.data[[1]])), plots.title[[1]], as.double(args[[3]]), c("kurpicz", "hQuick"))
plot.4 <- lineplotMemory(plots.data[[2]], plots.filters[[2]], c(1:length(plots.data[[2]])), plots.title[[2]], as.double(args[[4]]), c("kurpicz", "hQuick"))
plot.5 <- lineplotSpeedup(plots.data[[1]], plots.filters[[1]],  c(1:length(plots.data[[1]])), "sorting_overall", "maxTime", plots.title[[1]])
plot.6 <- lineplotSpeedup(plots.data[[2]], plots.filters[[2]],  c(1:length(plots.data[[2]])), "sorting_overall", "maxTime", plots.title[[2]])
#plot.5 <- lineplot(plots.data[[5]], plots.filters[[5]],  c(1:length(plots.data[[5]])), "sorting_overall", "maxTime", plots.title[[5]])
#plot.6 <- lineplot(plots.data[[6]], plots.filters[[6]],  c(1:length(plots.data[[6]])), "sorting_overall", "maxTime", plots.title[[6]])


plot.1 <- addSettings(plot.1)
plot.1 <- plot.1 + theme(legend.direction = "horizontal")
plot.1 <- plot.1 + theme(legend.box.background = element_rect(colour = "black"))
plot.1 <- plot.1 + theme(legend.title = element_blank()) 
#plot.1
plot.2 <- addSettings(plot.2)
#plot.2
plot.3 <- addSettings(plot.3)
#plot.3
plot.4 <- addSettings(plot.4)
plot.5 <- addSettings(plot.5)
plot.6 <- addSettings(plot.6)
plot.7 <- addSettings(plot.7)
plot.8 <- addSettings(plot.8)
#plot.4
#plot.5 <- addSettings(plot.5)
#plot.5
#plot.6 <- addSettings(plot.6)
#plot.6

#plot.3 <- plot.3 + guides(shape = guide_legend(title = "algorithm", title.position = "top", ncol = 4, nrow = 2)) 
#plot.3 <- plot.3 + guides(colour = guide_legend(title = "algorithm", title.position = "top", ncol = 4, nrow = 2)) 
#plot.3 <- plot.3 + guides(linetype = guide_legend(title = "algorithm", title.position = "top", ncol = 4, nrow = 2)) 
legend <- getLegend(plot.1)
plot.1 <- plot.1 + theme(legend.position = "none") + expand_limits(y=c(0, limit1))
plot.2 <- plot.2 + theme(legend.position = "none") + expand_limits(y=c(0, limit1))
plot.3 <- plot.3 + theme(legend.position = "none") + expand_limits(y=c(0,limit2)) + ggtitle("")
plot.4 <- plot.4 + theme(legend.position = "none")+ expand_limits(y=c(0, limit2)) + ggtitle("")
plot.7 <- plot.7 + theme(legend.position = "none") +expand_limits(y=c(0, limit3)) + ggtitle("")
plot.8 <- plot.8 + theme(legend.position = "none") +expand_limits(y=c(0, limit3))+ ggtitle("")
#plot.1 <- plotk1 + theme(legend.position = "none") + expand_limits(y=c(0, 40))
#plot.2 <- plot.2 + theme(legend.position = "none", axis.title.x = element_blank()) + expand_limits(y=0)
#plot.3 <- plot.3 + theme(legend.position = "none") + expand_limits(y=c(0, 40))
#plot.4 <- plot.4 + theme(legend.position = "none") + expand_limits(y=0)
#plot.5 <- plot.5 + theme(legend.position = "none") + expand_limits(y=c(0, 40))
#plot.6 <- plot.6 + theme(legend.position = "none") + expand_limits(y=0)
top <- plot_grid(plot.1, plot.2, ncol = 2)
middle <- plot_grid(plot.7, plot.8, ncol = 2)
bottom <- plot_grid(plot.3, plot.4, ncol=2)

legend <- plot_grid(NULL, legend, NULL, ncol=3, rel_widths=c(0.25, 0.5, 0.25))
plot_grid(top, middle, bottom,  legend, nrow = 4, rel_heights = c(1, 1, 1, 0.2))
plot.5 
plot.6 
