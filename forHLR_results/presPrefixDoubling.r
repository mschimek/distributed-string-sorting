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


isD2N <- length(unique(allDataWithoutIt1_Timer$dToNRatio)) > 1

availableProcessorCounts <- unique(allData$numberProcessors)
availableByteEncoders <- unique(allData$ByteEncoder)

scatterPlot <- function(data_, operations_, type_, pointSize, title) {
  print(operations_)
  data_ <- filter(data_, operation %in% operations_, type ==  type_)
  print(data_)
plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = operation),
                            size = pointSize, position = "jitter")
  if (isD2N) {
  plot <- plot + facet_wrap(~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(~ samplePolicy)
  }
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotPerRound <- function(data_, operations_, type_, pointSize, title) {
  data_ <- filter(data_, operation %in% operations_, type ==  type_)
  data_$round <- as.factor(data_$round)
plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = round),
                            size = pointSize, position = "jitter")
  if (isD2N) {
  plot <- plot + facet_wrap(GolombEncoding ~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(GolombEncoding ~ samplePolicy)
  }
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
  print(isD2N)
  if (isD2N) {
  plot <- plot + facet_wrap(GolombEncoding ~ dToNRatio)
  } else {
  plot <- plot + facet_wrap(GolombEncoding ~ samplePolicy)
  }
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
  plot <- plot + facet_wrap(GolombEncoding ~ samplePolicy)
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


barPlot <- function(data_, operations_, type_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_)
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

barPlotWhitelist <- function(data_, operations_, type_, title = " ", cpuTime = FALSE) {
  filteredData <- filter(data_, operation %in% operations_, type == type_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  if (cpuTime) {
  valueMean$numberProcessors <- as.character(valueMean$numberProcessors)
  valueMean$numberProcessors <- as.numeric(valueMean$numberProcessors)
  valueMean$value <- valueMean$value * valueMean$numberProcessors
  valueMean$numberProcessors <- as.factor(valueMean$numberProcessors)
  }
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  if (isD2N) {
  plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  } else {
  plot <- plot + facet_wrap(~ samplePolicy, labeller = label_both, nrow=1)
  }
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (cpuTime) {
    plot <- plot + ylab("work in nano sec * #procs")
  } else {
  plot <- plot + ylab("time in  nano sec")
  }
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotWhitelistMinMax <- function(data_, operations_, types_, title = " ", cpuTime = FALSE) {
  filteredData <- filter(data_, operation %in% operations_, type %in% types_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, type, strongScaling, phase, operation) 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  if (cpuTime) {
  valueMean$numberProcessors <- as.character(valueMean$numberProcessors)
  valueMean$numberProcessors <- as.numeric(valueMean$numberProcessors)
  valueMean$value <- valueMean$value * valueMean$numberProcessors
  valueMean$numberProcessors <- as.factor(valueMean$numberProcessors)
  }
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = type, y = value, fill = operation), stat="identity")
  if (isD2N) {
  plot <- plot + facet_wrap(numberProcessors ~ dToNRatio, labeller = label_both, nrow=1)
  } else {
  plot <- plot + facet_wrap(numberProcessors ~ samplePolicy, labeller = label_both, nrow=1)
  }
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ylab("time in  nano sec")
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomfilterWhitelist <- function(data_, operations_, type_, title = " ", work = FALSE) {
  filteredData <- filter(data_, operation %in% operations_, type == type_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  processorGroup <- group_by(filteredData, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, iteration, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  # sum rounds
  sumMean <- summarise(processorGroup, value = sum(value, rm.na = TRUE))
  processorGroup <- group_by(sumMean, numberProcessors, samplePolicy, StringGenerator, dToNRatio, stringLength, MPIAllToAllRoutine,  ByteEncoder, GolombEncoding, StringSet, size, strongScaling, phase, operation) 
  # mean over iteration 
  valueMean <- summarise(processorGroup, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  if (work) {
  valueMean$numberProcessors <- as.character(valueMean$numberProcessors)
  valueMean$numberProcessors <- as.numeric(valueMean$numberProcessors)
  valueMean$value <- valueMean$value * valueMean$numberProcessors
  valueMean$numberProcessors <- as.factor(valueMean$numberProcessors)
  }
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ GolombEncoding, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (work) {
    plot <- plot + ylab("work in nano sec * #procs")
  } else {
  plot <- plot + ylab("time in  nano sec")
  }
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlotBloomfilterRoundWhitelist <- function(data_, operations_, type_, title = " ") {
  filteredData <- filter(data_, operation %in% operations_, type == type_)
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

barPlotBloomFilterPerInternIteration <- function(data_, operations_, type_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_)
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

barPlotBloomFilter <- function(data_, operations_, type_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_)
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

barPlotBloomFilter_Number <- function(data_, operations_, type_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_)
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
  data_ = filter(data_, operation == "allgather_splitters_bytes_sent")#, size == globalSize)
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
  data <- filter(data_, operation == "bloomfilter_numberCandidates", type == "number", iteration == 1, round == 2)
  print(data$numberProcessors)
  print(data$value)
  data_ = filter(data_, operation == "bloomfilter_recvHashValues", type == "number", iteration == 1)
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
data <- filter(allDataWithoutIt1_EmptyTimer)#, size == globalSize)

dToNRatios <- unique(data$dToNRatio)
sort(dToNRatios, TRUE)
for (dToNRatio_ in sort(dToNRatios)) {
  localData  <- filter(data, dToNRatio == dToNRatio_)
  title <- "Communication Volume" 
  if (length(dToNRatios) > 1) {
  title <- paste("Communication Volume for dToN = ", dToNRatio_)
  }
  print(communicationVolume(localData, title))
}

#scatterPlot(data, c("num_received_chars"), "number", 0.1, "Number of recv. characters")
#scatterPlot(data, c("num_recv_strings"), "number", 0.1, "Number of recv. strings")
#scatterPlot(data, c("bloomfilter_sentEncodedValues"), "number", 0.1, "sent encoded values")

#scatterPlotPerRound(filter(data, GolombEncoding == "noGolombEncoding"), "bloomfilter_numberCandidates", "number", 0.1, " number candidates per Round ")
#scatterPlotPerRound(filter(data, GolombEncoding == "sequentialGolombEncoding"), "bloomfilter_numberCandidates", "number", 0.1, " number candidates per Round ")
#scatterPlotPerRound(filter(data, GolombEncoding == "sequentialGolombEncoding"), "bloomfilter_sentEncodedValues", "number", 0.1, " sent encoded values per Round ")
#scatterPlotPerRound(filter(data, GolombEncoding == "GolombPipelined"), "bloomfilter_numberCandidates", "number", 0.1, " number candidates per Round ")
#scatterPlotPerRound(filter(data, GolombEncoding == "noGolombEncoding"), "bloomfilter_recvHashValues", "number", 0.1, " number recv hash values per Round noGolombEncoding")
#scatterPlotPerRound(filter(data, GolombEncoding == "sequentialGolombEncoding"), "bloomfilter_recvHashValues", "number", 0.1, " number recv hash values per Round sequentialGolombEncoding")
#numberPlot(data, "string_exchange_bytes_sent", "bytes sent in string exchange")
#numberPlot(data, "localL", "Sum of local L")
#numberPlot(data, "localD", "Sum of Local D")
#numberPlot(data, "charactersInSet", "Characters to sort")
#numberPlot(data, "bloomfilter_sentEncodedValues", " encoded values sent")
#numberPlot(data, "bloomfilter_recvHashValues", " recvHash values ")
  
  data <- filter(allDataWithoutIt1_EmptyTimer)
 operations <- c("sorting_overall")
 localData <- filter(data, GolombEncoding == "noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgTime", title = "overall avg Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall work noGolombEncoding", TRUE)
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss Work noGolombEncoding", TRUE)
 operations <- c("sort_splitter")
 localData <- filter(data, GolombEncoding == "noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgTime", title = "overall avg Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall work noGolombEncoding", TRUE)
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss Work noGolombEncoding", TRUE)


 operations <- c("sort_splitter","writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall", "sort_locally")
 localData <- filter(data, GolombEncoding == "noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "avgTime", title = "overall avg Time noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall work noGolombEncoding", TRUE)
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss noGolombEncoding")
 barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = "overall Loss Work noGolombEncoding", TRUE)

 operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "choose_splitters", "allgather_splitters", "sample_splitters", "bloomfilter_overall")
 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = paste("overall Time with GolombEncoding: ", "noGolombEncoding")))
 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = paste("overall Work with GolombEncoding: ", "noGolombEncoding"), TRUE))

 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = paste("overall Loss with GolombEncoding: ", "noGolombEncoding")))

 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = paste("overall Loss Work with GolombEncoding: ", "noGolombEncoding"), TRUE))

 operations <- c("writeback_permutation", "merge_ranges", "compute_ranges", "all_to_all_strings", "compute_interval_sizes", "allgather_splitters", "sample_splitters", "bloomfilter_overall")
 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = paste("overall Time with GolombEncoding: ", "noGolombEncoding")))
 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = paste("overall Work with GolombEncoding: ", "noGolombEncoding"), TRUE))

 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = paste("overall Loss with GolombEncoding: ", "noGolombEncoding")))

 print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxLoss", title = paste("overall Loss Work with GolombEncoding: ", "noGolombEncoding"), TRUE))
 #Bloomfilter internal

   bloomfilterFindDupOps <- c("bloomfilter_findDuplicatesSendDups", "bloomfilter_findDuplicatesFind", "bloomfilter_findDuplicatesMerge", "bloomfilter_findDuplicatesSetup")
   bloomfilterOuterOps <- c("bloomfilter_generateHashStringIndices", "bloomfilter_sortHashStringIndices", "bloomfilter_indicesOfLocalDuplicates", "bloomfilter_ReducedHashStringIndices", "bloomfilter_sendHashStringIndices", "bloomfilter_addPEIndex","bloomfilter_getIndices", "bloomfilter_setDepth", "bloomfilter_findDuplicatesOverallIntern")

   bloomfilterOuterOps #<- bloomfilterFindDupOps

 dToNRatios <- unique(data$dToNRatio)
 for (i in dToNRatios) {

   localData <- filter(data, dToNRatio == i, samplePolicy == "IndexedNumStrings")
   if (length(dToNRatios)== 1) {

     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals Time")))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals work"), TRUE))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss")))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss Work"), TRUE))
   } else {
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals Time with dToN: ", i)))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals work with dToN: ", i), TRUE))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss with dToN: ", i)))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss Work with dToN: ", i), TRUE))
   }
 }
   bloomfilterOuterOps <- c("bloomfilter_sendToFilterSetup", "bloomfilter_sendEncodedValuesOverall", "bloomfilter_golombEncoding", "bloomfilter_sendEncodedValues", "bloomfilter_golombDecoding")
 for (i in dToNRatios) {

   localData <- filter(data, dToNRatio == i, samplePolicy == "IndexedNumStrings")
   if (length(dToNRatios)== 1) {

     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals Time")))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals work"), TRUE))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss")))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss Work"), TRUE))
   } else {
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals Time with dToN: ", i)))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals work with dToN: ", i), TRUE))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss with dToN: ", i)))
     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals Loss Work with dToN: ", i), TRUE))
   }
  }
# print("im here 0")
#   bloomfilterOuterOps <- c("bloomfilter_findDuplicatesSetup", "bloomfilter_findDuplicatesMerge", "bloomfilter_findDuplicatesFind", "bloomfilter_findDuplicatesSendDups")
# for (i in dToNRatios) {
#
#   localData <- filter(data, dToNRatio == i, samplePolicy == "IndexedNumStrings")
#   if (length(dToNRatios)== 1) {
#
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time")))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max work"), TRUE))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals max Loss")))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals max Loss Work"), TRUE))
#   } else {
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time with dToN: ", i)))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max work with dToN: ", i), TRUE))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals max Loss with dToN: ", i)))
#     print(barPlotBloomfilterWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxLoss", title = paste("Bloomfilter Internals max Loss Work with dToN: ", i), TRUE))
#   }
# }
#   bloomfilterOuterOps <- c("bloomfilter_sendToFilterSetup", "bloomfilter_sendEncodedValuesOverall", "bloomfilter_golombEncoding", "bloomfilter_sendEncodedValues", "bloomfilter_golombDecoding","bloomfilter_generateHashStringIndices", "bloomfilter_sortHashStringIndices", "bloomfilter_indicesOfLocalDuplicates", "bloomfilter_ReducedHashStringIndices", "bloomfilter_sendHashStringIndices", "bloomfilter_addPEIndex","bloomfilter_getIndices", "bloomfilter_setDepth", "bloomfilter_findDuplicatesOverallIntern")
# procs <- unique(data$numberProcessors)
# procs <- procs[procs != 1] 
# print(procs) 
# for (i in dToNRatios) {
#   for ( j in procs) {
#     text <- i
#     localData <- filter(data, dToNRatio == i, numberProcessors == j, GolombEncoding == "noGolombEncoding", samplePolicy == "IndexedNumStrings")
#     cnumbers <- c("bloomfilter_findDuplicatesSendDups")
#     if (length (dToNRatios) == 1) {
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time" )))
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = cnumbers, type = "number", title = paste("Bloomfilter Internals" )))
#     } else {
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time with dToN: ", text)))
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = cnumbers, type = "number", title = paste("Bloomfilter Internals with dToN: ", text)))
#     }
#     localData <- filter(data, dToNRatio == i, numberProcessors == j, GolombEncoding == "sequentialGolombEncoding", samplePolicy == "IndexedNumStrings")
#     cnumbers <- c("bloomfilter_findDuplicatesSendDups")
#     if (length(dToNRatios) == 1) {
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time" )))
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = cnumbers, type = "number", title = paste("Bloomfilter Internals Time" )))
#     } else {
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = bloomfilterOuterOps, type = "maxTime", title = paste("Bloomfilter Internals max Time with dToN: ", text)))
#       print(barPlotBloomfilterRoundWhitelist(data_ = localData, operations_ = cnumbers, type = "number", title = paste("Bloomfilter Internals  with dToN: ", text)))
#     }
#   }
# }
#
#
