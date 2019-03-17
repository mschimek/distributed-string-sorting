library(ggplot2)
library(tidyverse)
library(magrittr)
library(stringr)

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
globalSize <- 5000000

filename = paste(args[[1]], "/data.txt", sep="")
allData <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1_Timer <- filter(allData,  iteration != 0)
allDataWithoutIt1_EmptyTimer <- filter(allData,  iteration != 0)

availableProcessorCounts <- unique(allData$numberProcessors)
availableByteEncoders <- unique(allData$ByteEncoder)

scatterPlot <- function(data_, operations_, type_, pointSize, title) {
  print(operations_)
  data_ <- filter(data_, operation %in% operations_, type ==  type_)
  print(data_)
plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = operation, y = value, colour = ByteEncoder),
                            size = pointSize, position = "jitter")
  plot <- plot + facet_wrap(~ dToNRatio )
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
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

communicationVolume <- function(data, title = " ") {
  allPhases <- unique(data$phase)
  onlyRawCommunication <- filter(data, rawCommunication == 1, phase != "none")
  onlyRawCommunication$numberProcessors <- as.factor(onlyRawCommunication$numberProcessors)
  processorGroup <- group_by(onlyRawCommunication, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, iteration, size, strongScaling, phase, operation) 

  summedValues <- summarise(processorGroup, value = sum(value, rm.na = TRUE))

  meanOfIteration <- group_by(summedValues, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, phase, operation)
  
  valueMean <- summarise(meanOfIteration, value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = phase, y = value, fill = numberProcessors, colour = operation), position="dodge", stat="identity")
  plot <- plot + facet_wrap(~ GolombEncoding)
  plot <- plot + ylab("# sent bytes")
  plot <- plot + ggtitle(title)
  plot
}

numberPlot <- function(data, operation_, title = " ") {
  filteredData <- filter(data, operation == operation_)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, iteration, size, strongScaling, phase, operation) 

  summedValues <- summarise(processorGroup, value = sum(value, rm.na = TRUE))

  meanOfIteration <- group_by(summedValues, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, phase)
  
  valueMean <- summarise(meanOfIteration, value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value), position="dodge", stat="identity")
  plot <- plot + facet_wrap(dToNRatio ~ GolombEncoding)
  plot <- plot + ggtitle(title)
  plot
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


scatterPlotAllProcessors <- function(data_, operations_, type_) {
  procs <- unique(allData$numberProcessors)
  procs <- sort(procs)
  for (i in procs) {
    print(i)
    print(scatterPlot(filter(data_, numberProcessors == i), operations_, type_, 0.5, paste(i, " procs type of measurement ", type_, sep = "")))
  }
}

#barPlot <- function(data_, operations_, type_, size_, title = " ") {
#  print("hallo")
#  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
#  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
#  sumInternIteration <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
#  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
#  group <- group_by(sumInternIteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
#  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
#  valueMean$size <- as.factor(valueMean$size)
#  plot <- ggplot(data = valueMean)
#  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
#  plot <- plot + facet_wrap(~numberProcessors, labeller = label_both, nrow=1)
#  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#  plot <- plot + ggtitle("blabla")
#  return(plot)
#}

# Groups data such that a group consists of all operations in operations_
# values of each group are added up
sumOperations <- function(data_, operations_) {
  filteredData <- filter(data_, operation %in% operations_)
  valueSum <- group_by(filteredData, iteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type) %>% summarise(value = sum(value, rm.na = TRUE))
  valueSum$operation <- "unified"
  return(valueSum)
}

# calculate speedup
divideByPE0 <- function(vec) {
  valueBase <- filter(vec, numberProcessors == 1)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
  vec$value <- (valueBase / vec$value) 
  return(vec)
}

# calculate efficiency
divideByPE0AndPENumber <- function(vec) {
  valueBase <- filter(vec, numberProcessors == 1)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
  vec$value <- (valueBase / vec$value) / vec$numberProcessors
  return(vec)
}

#Group: Mean of all iterations of operation_ -> all means are grouped such that a group consists of all numberOfProcessors with which operation_ is run
applyProcedurePerGroup <- function(data_, operation_, type_, size_, title = " ", procedure = divideByPE0) {
  filteredData <- filter(data_, operation == operation_, type == type_, size == size_)
  
  #Not sure whether the following two lines are really needed
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  filteredData$size <- as.factor(filteredData$size)

  # Calculate mean of each group
  group <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))

  # new Grouping: each group must contain multiple processor numbers -> want to apply procedure per PE
  speedup <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  speedup <- do(speedup, procedure(.))

  plot <- ggplot(data = speedup, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
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
  valueMean <- group_by(valueSum, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueMean <- summarise(valueMean, value = mean(value, rm.na = TRUE))

  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}


barPlot <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  sumInternIteration <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotWhitelist <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation %in% operations_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomfilterWhitelist <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation %in% operations_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, iteration, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  # sum rounds
  sumMean <- summarise(processorGroup, value = sum(value, rm.na = TRUE))
  processorGroup <- group_by(sumMean, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  # mean over iteration 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ GolombEncoding, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomfilterRoundWhitelist <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation %in% operations_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  filteredData$round <- as.factor(filteredData$round)
  
  processorGroup <- group_by(filteredData, round, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  # mean over iteration 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation, colour = round), position="dodge", stat="identity")
  plot <- plot + facet_wrap(~ GolombEncoding, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomFilterPerInternIteration <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)

  remainingOperations <- unique(filteredData$operation)
  bloomfilterOperations <- remainingOperations[startsWith(remainingOperations, "bloomfilter")]
  filteredData <- filter(filteredData, operation %in% bloomfilterOperations)
  sumInternIteration <- filteredData
  print(sumInternIteration)
  #sumInternIteration <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
  #sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, round, GolombEncoding, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(GolombEncoding ~ round, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomFilter <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  
  #let pass only bloomfilter related operations
  remainingOperations <- unique(filteredData$operation)
  bloomfilterOperations <- remainingOperations[startsWith(remainingOperations, "bloomfilter_sendEncoded") | startsWith(remainingOperations, "bloomfilter_golomb")]
  filteredData <- filter(filteredData, operation %in% bloomfilterOperations)

  sumInternIteration <- group_by(filteredData, numberProcessors, GolombEncoding, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, numberProcessors, GolombEncoding, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~GolombEncoding, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomFilter_Number <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  
  #let pass only bloomfilter related operations
  remainingOperations <- unique(filteredData$operation)
  bloomfilterOperations <- remainingOperations[startsWith(remainingOperations, "bloomfilter_")]
  filteredData <- filter(filteredData, operation %in% bloomfilterOperations)

  sumInternIteration <- group_by(filteredData, numberProcessors, GolombEncoding, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, iteration)
  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, numberProcessors, GolombEncoding, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~GolombEncoding, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

numBytesSent <- function(data_) {
  data_ = filter(data_, operation == "bytes_sent", size == 1000000)
  data_$dToNRatio <- as.factor(data_$dToNRatio)

  group <- group_by(data_, numberProcessors, samplePolicy, size, dToNRatio, ByteEncoder)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueSum <- summarise(group, value = sum(value, rm.na = TRUE))

  plot <- ggplot(data = valueSum)
  plot <- plot + geom_point(mapping = aes(x = ByteEncoder, y = value), size = 0.1, position = "jitter")
  plot <- plot + facet_wrap(dToNRatio ~ numberProcessors)
  plot <- plot + ggtitle("sum bytes sent")
  return(plot)
}

numAllgatherBytesSent <- function(data_) {
  data_ = filter(data_, operation == "allgather_splitters_bytes_sent", size == globalSize)
  data_$dToNRatio <- as.factor(data_$dToNRatio)
  data_$value <- as.factor(data_$value)

  #group <- group_by(data_, numberProcessors, samplePolicy, size, dToNRatio, ByteEncoder)
  #valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  #valueSum <- summarise(group, value = sum(value, rm.na = TRUE))
  #print(valueSum)

  print(data_)
  plot <- ggplot(data = data_)
  plot <- plot + geom_point(mapping = aes(x = ByteEncoder, y = value, colour = ByteEncoder), size = 0.1, position = "jitter")
  plot <- plot + facet_wrap(dToNRatio ~ numberProcessors)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle("Bytes Sent in allGather splitters")
  return(plot)
}

numRecvHashValues <- function(data_) {
  data <- filter(data_, operation == "bloomfilter_numberCandidates", type == "number", size == globalSize, iteration == 1, round == 2)
  print(data$numberProcessors)
  print(data$value)
  data_ = filter(data_, operation == "bloomfilter_recvHashValues", type == "number", size == globalSize, iteration == 1)
  data_$iteration <- as.factor(data_$iteration)
  data_$round <- as.factor(data_$round)
  print(unique(data_$round))

  plot <- ggplot(data = data_)
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = round, shape = iteration), size = 0.5, position = "jitter")
  #plot <- plot + facet_wrap(dToNRatio ~ numberProcessors)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle("number recv hash values per PE")
  return(plot)
}

pureDirName <- str_sub(args, start = 1, end = -2)
pdf(paste(pureDirName, "_plots_prefixDoubling.pdf",sep=""), width=10, height=5)
data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)

dToNRatios <- unique(data$dToNRatio)
sort(dToNRatios, TRUE)
for (dToNRatio_ in sort(dToNRatios)) {
  localData  <- filter(data, dToNRatio == dToNRatio_)
  title <- paste("Communication Volume for dToN = ", dToNRatio_)
  print(communicationVolume(localData, title))
  print(numberPlot(localData, "localL", "Sum of local L"))
}

numberPlot(data, "string_exchange_bytes_sent", "bytes sent in string exchange")
numberPlot(data, "localL", "Sum of local L")
numberPlot(data, "localD", "Sum of Local D")
numberPlot(data, "charactersInSet", "Characters to sort")
numberPlot(data, "bloomfilter_sentEncodedValues", " encoded values sent")
numberPlot(data, "bloomfilter_recvHashValues", " recvHash values ")
#operations <- c("all_to_all_strings")
#data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
#operations <- c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi")
#data <- filter(allDataWithoutIt1_Timer, size == globalSize)
# operations <- c("prefix_decompression")
#  
  data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
 barPlot(data_ = filter(data, operation == "sorting_overall"), operations_ = c(), type_ = "avgTime", size_ = globalSize, title = "overall avg Time")
 barPlot(data_ = filter(data, operation == "sorting_overall"), operations_ = c(), type_ = "avgLoss", size_ = globalSize, title = "overall avg Loss")

 operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall", "sort_locally")
 barPlotWhitelist(data_ = data, operations_ = operations, type_ = "avgTime", size_ = globalSize, title = "overall avg Time")
 barPlotWhitelist(data_ = data, operations_ = operations, type_ = "maxTime", size_ = globalSize, title = "overall max Time")
 barPlotWhitelist(data_ = data, operations_ = operations, type_ = "avgLoss", size_ = globalSize, title = "overall avg Loss")
 barPlotWhitelist(data_ = data, operations_ = operations, type_ = "maxLoss", size_ = globalSize, title = "overall max Loss")

 operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall")
 GolombEncodings <- unique(data$GolombEncoding)
 for (i in GolombEncodings) {
   localData <- filter(data, GolombEncoding == i)
   operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall", "sort_locally")
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgTime", size_ = globalSize, title = paste("overall avg Time with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", size_ = globalSize, title = paste("overall max Time with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgLoss", size_ = globalSize, title = paste("overall avg Loss with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", size_ = globalSize, title = paste("overall max Loss with GolombEncoding: ", i)))
   operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall")
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgTime", size_ = globalSize, title = paste("overall avg Time with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", size_ = globalSize, title = paste("overall max Time with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgLoss", size_ = globalSize, title = paste("overall avg Loss with GolombEncoding: ", i)))
   print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", size_ = globalSize, title = paste("overall max Loss with GolombEncoding: ", i)))
 }


 #Bloomfilter internal

   bloomfilterOuterOps <- c("bloomfilter_generateHashStringIndices", "bloomfilter_sortHashStringIndices", "bloomfilter_indicesOfLocalDuplicates", "bloomfilter_ReducedHashStringIndices", "bloomfilter_sendHashStringIndices", "bloomfilter_addPEIndex","bloomfilter_getIndices", "bloomfilter_setDepth")    


 dToNRatios <- unique(data$dToNRatio)
 for (i in dToNRatios) {
   localData <- filter(data, dToNRatio == i)
   print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "avgTime", size = globalSize, title = paste("Bloomfilter outer ops overall avg Time with dToN: ", i)))
   print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", size = globalSize, title = paste("Bloomfilter outer ops overall max Time with dToN: ", i)))
   print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "avgLoss", size = globalSize, title = paste("Bloomfilter outer ops overall avg Loss with dToN: ", i)))
   print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", size = globalSize, title = paste("Bloomfilter outer ops overall max Loss with dToN: ", i)))
}
 
procs <- unique(data$numberProcessors)
procs <- procs[procs != 1] 
print(procs)
 for (i in dToNRatios) {
   for ( j in procs) {

   localData <- filter(data, dToNRatio == i, numberProcessors == j, GolombEncoding == "sequentialGolombEncoding")
   print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "avgTime", size = globalSize, title = paste("Bloomfilter outer ops overall avg Time with dToN: ", i)))
   print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", size = globalSize, title = paste("Bloomfilter outer ops overall max Time with dToN: ", i)))
   print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "avgLoss", size = globalSize, title = paste("Bloomfilter outer ops overall avg Loss with dToN: ", i)))
   print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", size = globalSize, title = paste("Bloomfilter outer ops overall max Loss with dToN: ", i)))
   }
}


#  barPlot(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "sorting_overall"), operations_ = c(), type_ = "maxTime", size_ = globalSize, title = "overall max Time")
#  barPlot(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "bloomfilter_overall"), operations_ = c(), type_ = "avgTime", size_ = globalSize, title = "overall Time")
#  barPlot(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "sort_locally"), operations_ = c(), type_ = "avgTime", size_ = globalSize, title = "overall Time")
# allDataWithoutIt1_EmptyTimer$numberProcessors <- as.character(allDataWithoutIt1_EmptyTimer$numberProcessors)
#  allDataWithoutIt1_EmptyTimer$numberProcessors <- as.numeric(allDataWithoutIt1_EmptyTimer$numberProcessors)
#  numRecvHashValues(data_ = filter(allDataWithoutIt1_EmptyTimer, numberProcessors <= 16))
#  numRecvHashValues(data_ = filter(allDataWithoutIt1_EmptyTimer, numberProcessors > 16))
#  atLeast2Processors <- filter(allDataWithoutIt1_EmptyTimer, !(numberProcessors %in% c(1)))
#  atLeast2Processors$numberProcessors <- as.factor(atLeast2Processors$numberProcessors)
#  filterOperation <- c("sorting_overall", "bloomfilter_overall", "bloomfilter_totalLoop", "bloomfilter_filterTotal", "bloomfilter_findDuplicatesOverallIntern", "bloomfilter_findDuplicatesOverallIntern_Barrier")
#  barPlotBloomFilter(data_ = atLeast2Processors, 
#                     operations_ = filterOperation, 
#                     type_ = "avgTime", 
#                     size_ = globalSize, 
#                     title = "Bloomfilter internal")
#  barPlotBloomFilter_Number(data_ = atLeast2Processors, 
#                     operations_ = filterOperation, 
#                     type_ = "number", 
#                     size_ = globalSize, 
#                     title = "Bloomfilter internal numbers")
#
#  barPlotBloomFilter(data_ = atLeast2Processors, 
#                     operations_ = filterOperation, 
#                     type_ = "minTime", 
#                     size_ = globalSize, 
#                     title = "Bloomfilter internal min time")
#  barPlotBloomFilter(data_ = atLeast2Processors, 
#                     operations_ = filterOperation, 
#                     type_ = "maxTime", 
#                     size_ = globalSize, 
#                     title = "Bloomfilter internal max")
#
#
#  print("hallo")
#  barPlot(data_ = atLeast2Processors, operations_ = c("sorting_overall", "bloomfilter_overall"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
#  print("hallo1.5")
#  for (np in sort(unique(atLeast2Processors$numberProcessors))) {
#    print(np)
#    plot <- barPlotBloomFilterPerInternIteration(data_ = filter(atLeast2Processors, numberProcessors == np), 
#                                           operations_ = filterOperation, 
#                                           type_ = "avgTime", 
#                                           size_ = globalSize, 
#                                           title = paste("Bloomfilter internal with avgTime", np, " PEs"))
#    print(plot)
#    plot <- barPlotBloomFilterPerInternIteration(data_ = filter(atLeast2Processors, numberProcessors == np), 
#                                                 operations_ = filterOperation, 
#                                                 type_ = "maxTime", 
#                                                 size_ = globalSize, 
#                                                 title = paste("Bloomfilter internal max Time with ", np, " PEs"))
#
#    print(plot)
#    plot <- barPlotBloomFilterPerInternIteration(data_ = filter(atLeast2Processors, numberProcessors == np), 
#                                                 operations_ = filterOperation, 
#                                                 type_ = "minTime", 
#                                                 size_ = globalSize, 
#                                                 title = paste("Bloomfilter internal min Time with ", np, " PEs"))
#
#    print(plot)
#    plot <- barPlotBloomFilterPerInternIteration(data_ = filter(atLeast2Processors, numberProcessors == np), 
#                                                 operations_ = c("sorting_overall"), 
#                                                 type_ = "number", 
#                                                 size_ = globalSize, 
#                                                 title = paste("Bloomfilter internal with ", np, " PEs"))
#    print(plot)
#
#  }
#  print("hallo2")
#  operations <- unique(allDataWithoutIt1_EmptyTimer$operation)
#  barrierOperations <- operations[endsWith(operations, "Barrier")]
#  bloomfilterOperation <- operations[startsWith(operations, "bloomfilter")]
#  dataWithoutBarriers <- filter(atLeast2Processors, !(operation %in% barrierOperations))
#  bloomfilterDataWithoutBarriers <- filter(dataWithoutBarriers, (operation %in% bloomfilterOperation))
#      barPlot(data_ = dataWithoutBarriers,  operations_ = c("sorting_overall", "bloomfilter_overall", "prefix_decompression", "allgatherv_test_before", "bloomfilter_findDuplicatesOverallIntern"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
#  barPlot(data_ = bloomfilterDataWithoutBarriers, operations_ = c("sorting_overall", "bloomfilter_overall", "prefix_decompression", "sort_locally", "allgatherv_test_before"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #
#  barPlot(data_ = dataWithoutBarriers, operations_ = c("sorting_overall", "bloomfilter_overall", "prefix_decompression", "sort_locally", "all_to_all_strings", "merge_ranges", "compute_interval_sizes"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #
#  barPlot(data_ = dataWithoutBarriers, operations_ = c("choose_splitters","allgatherv_test_before", "bloomfilter_overall", "sorting_overall", "prefix_decompression", "sort_locally", "all_to_all_strings", "merge_ranges", "compute_interval_sizes"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #
#
