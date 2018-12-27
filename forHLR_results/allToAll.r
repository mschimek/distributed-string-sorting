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
POINTSIZE = 0.1

#allData <- read.table("./output.txt", comment.char = "#", col.names = columns)
allData <- read_delim(file = "2018_12_27_H13_38/data.txt", delim = "|", col_types = colTypeSpec, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1_intern <- filter(allData, iteration != 1, Timer == "Timer")
allDataWithoutIt1_noIntern <- filter(allData, iteration != 1, Timer == "EmptyTimer")

availableProcessorCounts <- unique(allData$numberProcessors)
availableByteEncoders <- unique(allData$ByteEncoder)

scatter <- function(data_, operations_, pointSize, title) {
plot <- ggplot(data = filter(data_, operation %in% operations_, type %in% c("avgTime")))
  plot <- plot + geom_point(mapping = aes(x = operation, y = time, colour = ByteEncoder),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ MPIAllToAllRoutine)
  plot <- plot + ggtitle(title)
  return(plot)
}
plot1 <- scatter(allDataWithoutIt1_intern, c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "blabla")

plot2 <- scatter(filter(allDataWithoutIt1_intern, numberProcessors == 2), c("all_to_all_strings"), 0.5, "blabla")
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




pdf(paste("plots_alltoall.pdf",sep=""), width=10, height=5)
plot1 <- scatter(filter(allDataWithoutIt1_intern, numberProcessors == 2), c("all_to_all_strings"), 0.5, " 2 procs")
plot2 <- scatter(filter(allDataWithoutIt1_intern, numberProcessors == 4), c("all_to_all_strings"), 0.5, " 4 procs")
plot3 <- scatter(filter(allDataWithoutIt1_intern, numberProcessors == 8), c("all_to_all_strings"), 0.5, " 8 procs")
plot1
plot2
plot3
plot4 <- scatter(filter(allDataWithoutIt1_noIntern, numberProcessors == 2), c("all_to_all_strings"), 0.5, " 2 procs")
plot5 <- scatter(filter(allDataWithoutIt1_noIntern, numberProcessors == 4), c("all_to_all_strings"), 0.5, " 4 procs")
plot6 <- scatter(filter(allDataWithoutIt1_noIntern, numberProcessors == 8), c("all_to_all_strings"), 0.5, " 8 procs")
plot4
plot5
plot6
scatter(filter(allDataWithoutIt1_intern, numberProcessors == 2), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "blabla")
scatter(filter(allDataWithoutIt1_intern, numberProcessors == 4), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "blabla")
scatter(filter(allDataWithoutIt1_intern, numberProcessors == 8), c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi"), 0.5, "blabla")
