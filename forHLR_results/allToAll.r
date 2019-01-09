library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(args[[1]])

colTypeSpec = cols(numberProcessors = col_integer(),
       samplePolicy = col_character(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())

# "Constants"
POINTSIZE = 0.1

filename = paste(args[[1]], "/data.txt", sep="")
allData <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1_Timer <- filter(allData, iteration != 1, Timer == "Timer")
allDataWithoutIt1_EmptyTimer <- filter(allData, iteration != 1, Timer == "EmptyTimer")

availableProcessorCounts <- unique(allData$numberProcessors)
availableByteEncoders <- unique(allData$ByteEncoder)

scatterPlot <- function(data_, operations_, type_, pointSize, title) {
plot <- ggplot(data = filter(data_, operation %in% operations_, type ==  type_))
  plot <- plot + geom_point(mapping = aes(x = operation, y = value, colour = ByteEncoder),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ MPIAllToAllRoutine)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

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



boxPlotOnOperationProcsAreFacets <- function(data_, operation_, types_, title) {
  filteredData <- filter(data_, operation == operation_, type %in% types_)
  plot <- ggplot(data = filteredData) 
  plot <- plot + geom_boxplot(mapping = aes(x = samplePolicy, y = value, colour = type)) 
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
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxisType <- function(data_, operation_, type_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_, type %in% type_))
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = type),
                            size = pointSize, position = "jitter"

  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}


scatterPlotAllProcessors <- function(data_, operations_, type_) {
  procs <- unique(allData$numberProcessors)
  procs <- sort(procs)
  for (i in procs) {
    print(i)
    print(scatterPlot(filter(data_, numberProcessors == i), operations_, type_, 0.5, paste(i, " procs type of measurement ", type_, sep = "")))
  }
}

pureDirName <- str_sub(args, start = 1, end = -2)
pdf(paste(pureDirName, "_plots_alltoall.pdf",sep=""), width=10, height=5)

operations <- c("all_to_all_strings")
data <- allDataWithoutIt1_EmptyTimer
scatterPlotAllProcessors(data, operations, "avgTime")
scatterPlotAllProcessors(data, operations, "maxTime")
scatterPlotAllProcessors(data, operations, "avgLoss")
scatterPlotAllProcessors(data, operations, "maxLoss")

operations <- c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi")
data <- allDataWithoutIt1_Timer
scatterPlotAllProcessors(data, operations, "avgTime")
scatterPlotAllProcessors(data, operations, "maxTime")
scatterPlotAllProcessors(data, operations, "avgLoss")
scatterPlotAllProcessors(data, operations, "maxLoss")
#
#plot1 <- scatter(filter(allDataWithoutIt1_Timer, numberProcessors == 2), c("all_to_all_strings"), 0.5, " 2 procs")
#plot2 <- scatter(filter(allDataWithoutIt1_Timer, numberProcessors == 4), c("all_to_all_strings"), 0.5, " 4 procs")
#plot3 <- scatter(filter(allDataWithoutIt1_Timer, numberProcessors == 8), c("all_to_all_strings"), 0.5, " 8 procs")
#plot4 <- scatter(filter(allDataWithoutIt1_Timer, numberProcessors == 32), c("all_to_all_strings"), 0.5, " 32 procs")
#plot1
#plot2
#plot3
#plot4
##plot5 <- scatterPlot(filter(allDataWithoutIt1_EmptyTimer, numberProcessors == 2), c("all_to_all_strings"), 0.5, " 2 procs")
##plot6 <- scatterPlot(filter(allDataWithoutIt1_EmptyTimer, numberProcessors == 4), c("all_to_all_strings"), 0.5, " 4 procs")
##plot7 <- scatterPlot(filter(allDataWithoutIt1_EmptyTimer, numberProcessors == 8), c("all_to_all_strings"), 0.5, " 8 procs")
##plot8 <- scatterPlot(filter(allDataWithoutIt1_EmptyTimer, numberProcessors == 32), c("all_to_all_strings"), 0.5, " 32 procs")
##plot5
##plot6
##plot7
##plot8
##scatterPlot(filter(allDataWithoutIt1_Timer, numberProcessors == 2), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "2")
##scatterPlot(filter(allDataWithoutIt1_Timer, numberProcessors == 4), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "4")
##scatterPlot(filter(allDataWithoutIt1_Timer, numberProcessors == 8), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "8")
##scatterPlot(filter(allDataWithoutIt1_Timer, numberProcessors == 32), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "32")
