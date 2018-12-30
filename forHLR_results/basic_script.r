library(ggplot2)
library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)
#
## test if there is at least one argument: if not, return an error
#if (length(args) != 1) {
#  stop("At least one argument must be supplied (input file)", call.=FALSE)
#}

#print(args[[1]])

columns <- c("numberProcessors",
             "samplePolicy",
             "iteration",
             "size",
             "operation",
             "type",
             "time")

colTypeSpec = cols(numberProcessors = col_integer(),
       samplePolicy = col_character(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       time = col_double())

# "Constants"
POINTSIZE = 0.05

#allData <- read.table("./output.txt", comment.char = "#", col.names = columns)
allData <- read_delim(file = "2018_12_28_H23_48/data.txt", "|", col_types = colTypeSpec, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1 <- filter(allData, iteration != 1, Timer != "Timer")

availableProcessorCounts <- unique(allData$numberProcessors)

#functions
scatterPlotOnOperation <- function(data_, operation_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_point(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

boxPlotOnOperation <- function(data_, operation_, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_boxplot(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type)) 
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title) 
  return(plot)
}

computeYLimitsBoxplot <- function(data_) {
  return(boxplot.stats(data)$stats[c(1,5)])
} 

computeYLimitsBoxplot_v <- function(data_) {
  min <- 999999999
  max <- 0
  for (name in colnames(data_)) {
    yLimits <- computeYLimitsBoxplot(data_$name)
    if (min > yLimits[0]) {
      min <- yLimits[0]
    }
    if (max < yLimits[1]) {
      max <- yLimits[1]
    }
  }
  c(min, max)
}



boxPlotOnOperationProcsAreFacets <- function(data_, operation_, types_, title) {
  filteredData <- filter(data_, operation == operation_, type %in% types_)
  plot <- ggplot(data = filteredData) 
  plot <- plot + geom_boxplot(mapping = aes(x = samplePolicy, y = time, colour = type)) 
  plot <- plot + facet_wrap(~ numberProcessors, nrow = 1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title) 
  return(plot)
}

boxPlotOnOperationProcsAreFacets_ <- function(data_, operation_, types_, title, ylimits) {
  plot <- boxPlotOnOperationProcsAreFacets(data_, operation_, types_, title)
  plot <- plot + coord_cartesian(ylim = ylimits)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxis <- function(data_, operation_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_))
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = time, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxisType <- function(data_, operation_, type_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_, type %in% type_))
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = time, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}




pdf(paste("plots.pdf",sep=""), width=10, height=5)

smallStringSet <- filter(allDataWithoutIt1, size == 1000000)
bigStringSet <- filter(allDataWithoutIt1, size == 5000000)

#All-to-All Strings
operation = "all_to_all_strings"
filterTime = c("avgTime", "maxTime", "minTime")
filterLoss = c("avgLoss", "maxLoss", "minLoss")
scatterPlotOnOperation(smallStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 1000000")
scatterPlotOnOperation(bigStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 5000000")
boxPlotOnOperation(smallStringSet, operation, "All-to-All Strings, StringSetSize = 1000000")
boxPlotOnOperation(bigStringSet, operation, "All-to-All Strings, StringSetSize = 5000000")
scatterPlotOnOperationProcsOnXAxis(smallStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 1000000")
scatterPlotOnOperationProcsOnXAxisType(smallStringSet, operation,
                                       filterTime, POINTSIZE, "All-to-All Strings, Time only, StringSetSize = 1000000")
scatterPlotOnOperationProcsOnXAxisType(smallStringSet, operation,
                                       filterLoss, POINTSIZE, "All-to-All Strings, Loss only, StringSetSize = 1000000")
boxPlotOnOperationProcsAreFacets(smallStringSet, operation, filterTime, "All-to-All Strings, StringSetSize = 1000000")
boxPlotOnOperationProcsAreFacets_(filter(smallStringSet, numberProcessors != 32), operation, filterTime, "All-to-All Strings, ZoomIn, StringSetSize = 1000000", c(250000000, 420000000))
#Sort String locally
scatterPlotOnOperation(smallStringSet, "sort_locally", POINTSIZE, "Sort locally, StringSetSize = 1000000")
scatterPlotOnOperation(bigStringSet, "sort_locally", POINTSIZE, "Sort locally, StringSetSize = 5000000")
boxPlotOnOperation(smallStringSet, "sort_locally", "Sort Strings locally, StringSetSize = 1000000")
boxPlotOnOperation(bigStringSet, "sort_locally", "Sort Strings locally, StringSetSize = 5000000")

#Merge Sequences
scatterPlotOnOperation(smallStringSet, "merge_ranges", POINTSIZE, "Merge Sequences, StringSetSize = 1000000")
scatterPlotOnOperation(bigStringSet, "merge_ranges", POINTSIZE, "Merge Sequences, StringSetSize = 5000000")
boxPlotOnOperation(smallStringSet, "merge_ranges", "Merge Sequences, StringSetSize = 1000000")
boxPlotOnOperation(bigStringSet, "merge_ranges", "Merge Sequences, StringSetSize = 5000000")

#Overall
barPlot <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation != operation_, type == type_, size == size_)
  group <- group_by(filteredData, numberProcessors, samplePolicy, size, operation, type)
  timeMean <- summarise(group, time = mean(time, rm.na = TRUE))
  timeMean$size <- as.factor(timeMean$size)
  plot <- ggplot(data = timeMean)
  plot <- plot + geom_bar(mapping = aes(x = samplePolicy, y = time, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ numberProcessors, labeller = label_both, nrow=1)
  plot <- plot + ggtitle(title)

  return(plot)
}

barPlot(allDataWithoutIt1, operation_ = "sorting_overall", type_ = "avgTime", size_ = 1000000, "AvgTime and StringSetSize = 1000000")
barPlot(allDataWithoutIt1, operation_ = "sorting_overall", type = "avgTime", size = 5000000, "AvgTime and StringSetSize = 5000000")

#plot1
