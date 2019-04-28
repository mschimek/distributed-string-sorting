library(ggplot2)
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(args[[1]])

colTypeSpec = cols(numberOfProcessor = col_integer(),
       iteration = col_integer(),
       operation = col_character(),
       type = col_character(),
       value = col_double())

# "Constants"
POINTSIZE = 0.1
globalSize <- 5000000

filename = paste(args[[1]], "/data.txt", sep="")
allData <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
allData <- filter(allData, iteration != 0)
allData$numberOfProcessor <- as.factor(allData$numberOfProcessor)

isD2N <- TRUE


scatterPlot <- function(data_, operations_, type_, pointSize, title) {
  print(operations_)
  data_ <- filter(data_, operation %in% operations_, type ==  type_)
  print(data_)
  print("halloJK")
plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = operation, y = value, colour = ByteEncoder),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio )
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotOnOperation <- function(data_, operation_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_point(mapping = aes(x = interaction(numberOfProcessor, type), y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

boxPlotOnOperation <- function(data_, operation_, title) {
  plot <- ggplot(data = filter(data_, operation == operation_)) 
  plot <- plot + geom_boxplot(mapping = aes(x = interaction(numberOfProcessor, type), y = value, colour = type)) 
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
  plot <- plot + facet_wrap(~ numberOfProcessor, nrow = 1)
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
  plot <- plot + geom_point(mapping = aes(x = numberOfProcessor, y = value, colour = type),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotOnOperationProcsOnXAxisType <- function(data_, operation_, type_, pointSize, title) {
  plot <- ggplot(data = filter(data_, operation == operation_, type %in% type_))
  plot <- plot + geom_point(mapping = aes(x = numberOfProcessor, y = value, colour = type),
                            size = pointSize, position = "jitter")

  plot <- plot + facet_wrap(~ samplePolicy)
  plot <- plot + ggtitle(title)
  return(plot)
}


scatterPlotAllProcessors <- function(data_, operations_, type_) {
  procs <- unique(allData$numberOfProcessor)
  procs <- sort(procs)
  for (i in procs) {
    print(i)
    print(scatterPlot(filter(data_, numberOfProcessor == i), operations_, type_, 0.5, paste(i, " procs type of measurement ", type_, sep = "")))
  }
}

barPlot <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

# Groups data such that a group consists of all operations in operations_
# values of each group are added up
sumOperations <- function(data_, operations_) {
  filteredData <- filter(data_, operation %in% operations_)
  valueSum <- group_by(filteredData, iteration, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, type) %>% summarise(value = sum(value, rm.na = TRUE))
  valueSum$operation <- "unified"
  return(valueSum)
}

# calculate speedup
divideByPE0 <- function(vec) {
  valueBase <- filter(vec, numberOfProcessor == 1)$value
  vec$numberOfProcessor <- as.numeric(as.character(vec$numberOfProcessor))
  vec$value <- (valueBase / vec$value) 
  return(vec)
}

# calculate efficiency
divideByPE0AndPENumber <- function(vec) {
  valueBase <- filter(vec, numberOfProcessor == 1)$value
  vec$numberOfProcessor <- as.numeric(as.character(vec$numberOfProcessor))
  vec$value <- (valueBase / vec$value) / vec$numberOfProcessor
  return(vec)
}

#Group: Mean of all iterations of operation_ -> all means are grouped such that a group consists of all numberOfProcessors with which operation_ is run
applyProcedurePerGroup <- function(data_, operation_, type_, size_, title = " ", procedure = divideByPE0) {
  filteredData <- filter(data_, operation == operation_, type == type_, size == size_)
  
  #Not sure whether the following two lines are really needed
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  filteredData$size <- as.factor(filteredData$size)

  # Calculate mean of each group
  group <- group_by(filteredData, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))

  # new Grouping: each group must contain multiple processor numbers -> want to apply procedure per PE
  speedup <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  speedup <- do(speedup, procedure(.))

  plot <- ggplot(data = speedup, mapping = aes(x = numberOfProcessor, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

speedup <- function(data_, operation_, type_, size_, title = " ") {
  applyProcedurePerGroup(data_, operation_, type_, size_, title, procedure = divideByPE0)
 }

efficiency <- function(data_, operation_, type_, size_, title = " ") {
  applyProcedurePerGroup(data_, operation_, type_, size_, title, procedure = divideByPE0AndPENumber)
 }

speedup_sum <- function(data_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  filteredData$size <- as.factor(filteredData$size)

  operations <- unique(filteredData$operation)
  valueSum <- sumOperations(filteredData, operations)
  return (speedup(valueSum, "unified", type_, size_, title))
}

efficiency_sum <- function(data_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  filteredData$size <- as.factor(filteredData$size)

  operations <- unique(filteredData$operation)
  valueSum <- sumOperations(filteredData, operations)
  return (efficiency(valueSum, "unified", type_, size_, title) )
}

sumData <- function(data_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)

  #Sum over all operations of one iteration
  operations <- unique(filteredData$operation)
  valueSum <- sumOperations(filteredData, operations)

  #Mean over the iterations
  valueMean <- group_by(valueSum, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueMean <- summarise(valueMean, value = mean(value, rm.na = TRUE))

  plot <- ggplot(data = valueMean, mapping = aes(x = numberOfProcessor, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}


barPlot <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(numberOfProcessor ~ samplePolicy, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

numBytesSent <- function(data_) {
  data_ = filter(data_, operation == "bytes_sent", size == 1000000)
  data_$dToNRatio <- as.factor(data_$dToNRatio)

  group <- group_by(data_, numberOfProcessor, samplePolicy, size, dToNRatio, ByteEncoder)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueSum <- summarise(group, value = sum(value, rm.na = TRUE))

  plot <- ggplot(data = valueSum)
  plot <- plot + geom_point(mapping = aes(x = ByteEncoder, y = value), size = 0.1, position = "jitter")
  plot <- plot + facet_wrap(dToNRatio ~ numberOfProcessor)
  plot <- plot + ggtitle("sum bytes sent")
  return(plot)
}

numAllgatherBytesSent <- function(data_) {
  data_ = filter(data_, operation == "allgather_splitters_bytes_sent", size == globalSize)
  data_$dToNRatio <- as.factor(data_$dToNRatio)
  data_$value <- as.factor(data_$value)

  #group <- group_by(data_, numberOfProcessor, samplePolicy, size, dToNRatio, ByteEncoder)
  #valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  #valueSum <- summarise(group, value = sum(value, rm.na = TRUE))
  #print(valueSum)

  print(data_)
  plot <- ggplot(data = data_)
  plot <- plot + geom_point(mapping = aes(x = ByteEncoder, y = value, colour = ByteEncoder), size = 0.1, position = "jitter")
  plot <- plot + facet_wrap(dToNRatio ~ numberOfProcessor)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle("Bytes Sent in allGather splitters")
  return(plot)
}

communicationVolume <- function(data, title = " ") {
  allPhases <- unique(data$phase)
  onlyRawCommunication <- filter(data, rawCommunication == 1, phase != "none")
  onlyRawCommunication$numberOfProcessor <- as.factor(onlyRawCommunication$numberOfProcessor)
  processorGroup <- group_by(onlyRawCommunication, numberOfProcessor, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, StringSet, iteration, size, strongScaling, phase, operation) 

  summedValues <- summarise(processorGroup, value = sum(value, rm.na = TRUE))

  meanOfIteration <- group_by(summedValues, numberOfProcessor, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, StringSet, size, phase, operation)
  
  valueMean <- summarise(meanOfIteration, value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = phase, y = value, fill = numberOfProcessor, colour = operation), position="dodge", stat="identity")
  if (isD2N) {
  plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(ByteEncoder ~ samplePolicy)
  }
  plot <- plot + ylab("# sent bytes")
  plot <- plot + ggtitle(title)
  plot
}

numberPlot <- function(data, operation_, title = " ") {
  filteredData <- filter(data, operation == operation_)
  processorGroup <- group_by(filteredData, numberOfProcessor, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, StringSet, iteration, size, strongScaling, phase, operation) 

  summedValues <- summarise(processorGroup, value = sum(value, rm.na = TRUE))

  meanOfIteration <- group_by(summedValues, numberOfProcessor, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, StringSet, size, phase)
  
  valueMean <- summarise(meanOfIteration, value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberOfProcessor, y = value), position="dodge", stat="identity")
  if (isD2N) {
  plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(ByteEncoder ~ samplePolicy)
  }
  plot <- plot + ggtitle(title)
  plot
}
barPlotWhitelist <- function(data_, operations_, type_, title = " ", work = FALSE) {
  filteredData <- filter(data_, operation %in% operations_, type == type_)
#  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
#  sumInternIteration <- group_by(filteredData, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
#  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
#  group <- group_by(sumInternIteration, numberOfProcessor, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
#  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
#  valueMean$size <- as.factor(valueMean$size)
#  if (work) {
#  valueMean$numberOfProcessor <- as.character(valueMean$numberOfProcessor)
#  valueMean$numberOfProcessor <- as.numeric(valueMean$numberOfProcessor)
#  valueMean$value <- valueMean$value * valueMean$numberOfProcessor
#  valueMean$numberOfProcessor <- as.factor(valueMean$numberOfProcessor)
  sumInternIteration <- group_by(filteredData, numberOfProcessor, operation, type, iteration)
  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, numberOfProcessor,  operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  if (work) {
  valueMean$numberOfProcessor <- as.character(valueMean$numberOfProcessor)
  valueMean$numberOfProcessor <- as.numeric(valueMean$numberOfProcessor)
  valueMean$value <- valueMean$value * valueMean$numberOfProcessor
  valueMean$numberOfProcessor <- as.factor(valueMean$numberOfProcessor)
  }
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberOfProcessor, y = value, fill = operation), stat="identity")
  if (work) {
    plot <- plot + ylab("work in nano sec * #procs")
  } else {
  plot <- plot + ylab("time in  nano sec")
  }
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlot <- function(data_, operations_, type_, pointSize, title) {
  print(operations_)
  data_ <- filter(data_, operation %in% operations_, type ==  type_)
  print(data_)
plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = numberOfProcessor, y = value, colour = operation),
                            size = pointSize, position = "jitter")
  if (isD2N) {
  plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(ByteEncoder ~ samplePolicy)
  }
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}

pureDirName <- str_sub(args, start = 1, end = -2)
pdf(paste(pureDirName, "_plots_prefixCompression.pdf",sep=""), width=10, height=5)

#communicationVolume(allData, "communication volume ")
#scatterPlot(allData, c("num_received_chars"), "number", 0.1, "Number of recv. characters")
#scatterPlot(allData, c("num_recv_strings"), "number", 0.1, "Number of recv. strings")
#numberPlot(allData, "charactersInSet", "number characters to sort")

 print("1")
 operations <- c("distributed_sort", "allgatherv", "local_sort")
 barPlotWhitelist(data_ = allData, operations_ = operations, type_ = "maxTime", title = "overall max Time")
 print("1.1")
 barPlotWhitelist(data_ = allData, operations_ = operations, type_ = "avgTime", title = "overall avg Time")
 barPlotWhitelist(data_ = allData, operations_ = operations, type_ = "maxTime", title = "overall Work", TRUE)
 barPlotWhitelist(data_ = allData, operations_ = operations, type_ = "maxLoss", title = "overall avg Loss")
 barPlotWhitelist(data_ = allData, operations_ = operations, type_ = "maxLoss", title = "overall avg Loss Work", TRUE)
 print("2")
