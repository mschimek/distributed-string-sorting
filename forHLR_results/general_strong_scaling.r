library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}
dirName <- args[[1]]
columns <- c("numberProcessors",
             "samplePolicy",
             "iteration",
             "size",
             "operation",
             "type",
             "value")

colTypeSpec = cols(numberProcessors = col_integer(),
       samplePolicy = col_character(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())

# "Constants"
POINTSIZE = 0.05
filename <- paste(args[[1]], "/data.txt", sep="")
#allData <- read.table("./output.txt", comment.char = "#", col.names = columns)
allData <- read_delim(file = filename, "|", col_types = colTypeSpec, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1 <- filter(allData, iteration != 1, Timer != "Timer")

availableProcessorCounts <- unique(allData$numberProcessors)

#functions
scatterPlotOnOperation <- function(data_, operation_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_point(mapping = aes(x = interaction(numberProcessors, type), y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

boxPlotOnOperation <- function(data_, operation_, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_boxplot(mapping = aes(x = interaction(numberProcessors, type), y = value, colour = type)) 
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



boxPlotOnOperationProcsAreFacets_ <- function(data_, operation_, types_, title, ylimits) {
  plot <- boxPlotOnOperationProcsAreFacets(data_, operation_, types_, title)
  plot <- plot + coord_cartesian(ylim = ylimits)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxis <- function(data_, operation_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_))
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxisType <- function(data_, operation_, type_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_, type %in% type_))
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}

overallBarPlot <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation != operation_, type == type_, size == size_, Timer != "Timer")
  group <- group_by(filteredData, numberProcessors, samplePolicy, size, operation, type)
  timeMean <- summarise(group, time = mean(value, rm.na = TRUE))
  timeMean$size <- as.factor(timeMean$size)
  plot <- ggplot(data = timeMean)
  plot <- plot + geom_bar(mapping = aes(x = samplePolicy, y = time, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ numberProcessors, labeller = label_both, nrow=1)
  plot <- plot + ggtitle(title)

  return(plot)
}

boxPlotOnOperationProcsAreFacets <- function(data_, operation_, type_, size_, title) {
  filteredData <- filter(data_, operation == operation_, type == type_, size == size_, Timer != "Timer")
  print(filteredData)
  plot <- ggplot(data = filteredData) 
  plot <- plot + geom_boxplot(mapping = aes(x = samplePolicy, y = value, colour = operation)) 
  plot <- plot + facet_wrap(~ numberProcessors, nrow = 1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
#scatterPlotOnOperation(smallStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 1000000")
#scatterPlotOnOperation(bigStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 5000000")
#boxPlotOnOperation(smallStringSet, operation, "All-to-All Strings, StringSetSize = 1000000")
#boxPlotOnOperation(bigStringSet, operation, "All-to-All Strings, StringSetSize = 5000000")
#scatterPlotOnOperationProcsOnXAxis(smallStringSet, operation, POINTSIZE, "All-to-All Strings, StringSetSize = 1000000")
#scatterPlotOnOperationProcsOnXAxisType(smallStringSet, operation,
#                                       filterTime, POINTSIZE, "All-to-All Strings, Time only, StringSetSize = 1000000")
#scatterPlotOnOperationProcsOnXAxisType(smallStringSet, operation,
#                                       filterLoss, POINTSIZE, "All-to-All Strings, Loss only, StringSetSize = 1000000")
#boxPlotOnOperationProcsAreFacets(smallStringSet, operation, filterTime, "All-to-All Strings, StringSetSize = 1000000")
#boxPlotOnOperationProcsAreFacets_(filter(smallStringSet, numberProcessors != 32), operation, filterTime, "All-to-All Strings, ZoomIn, StringSetSize = 1000000", c(250000000, 420000000))
##Sort String locally
#scatterPlotOnOperation(smallStringSet, "sort_locally", POINTSIZE, "Sort locally, StringSetSize = 1000000")
#scatterPlotOnOperation(bigStringSet, "sort_locally", POINTSIZE, "Sort locally, StringSetSize = 5000000")
#boxPlotOnOperation(smallStringSet, "sort_locally", "Sort Strings locally, StringSetSize = 1000000")
#boxPlotOnOperation(bigStringSet, "sort_locally", "Sort Strings locally, StringSetSize = 5000000")
#
##Merge Sequences
#scatterPlotOnOperation(smallStringSet, "merge_ranges", POINTSIZE, "Merge Sequences, StringSetSize = 1000000")
#scatterPlotOnOperation(bigStringSet, "merge_ranges", POINTSIZE, "Merge Sequences, StringSetSize = 5000000")
#boxPlotOnOperation(smallStringSet, "merge_ranges", "Merge Sequences, StringSetSize = 1000000")
#boxPlotOnOperation(bigStringSet, "merge_ranges", "Merge Sequences, StringSetSize = 5000000")
#
#Overall


dirNameWithoutSlash <- str_sub(dirName, start = 1, end = -2)
pdf(paste(dirNameWithoutSlash, "_plots_strong_scaling", ".pdf",sep=""), width=10, height=5)
overallBarPlot(allDataWithoutIt1, operation_ = "sorting_overall", type_ = "avgTime", size_ = 1000000, "AvgTime and StringSetSize = 1000000")
overallBarPlot(allDataWithoutIt1, operation_ = "sorting_overall", type = "avgTime", size = 5000000, "AvgTime and StringSetSize = 5000000")

boxPlotOnOperationProcsAreFacets(allDataWithoutIt1, operation_ = "merge_ranges", type = "avgTime", size = 5000000, "AvgTime and StringSetSize = 5000000")
boxPlotOnOperationProcsAreFacets(allDataWithoutIt1, operation_ = "all_to_all_strings", type = "avgTime", size = 5000000, "AvgTime and StringSetSize = 5000000")
boxPlotOnOperationProcsAreFacets(allDataWithoutIt1, operation_ = "sort_locally", type = "avgTime", size = 5000000, "AvgTime and StringSetSize = 5000000")
