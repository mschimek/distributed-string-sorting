library(ggplot2)
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(args[[1]])


colTypeSpec = cols(numberProcessors = col_integer(),
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
allDataFirstIter <- filter(allData, iteration == 0)
allDataAllButFirstIter <- filter(allData, iteration != 0)
allData$numberProcessors <- as.factor(allData$numberProcessors)

isD2N <- TRUE




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
  print(data_, n = 100)
  sumInternIteration <- group_by(filteredData, numberProcessors, dToNRatio, stringLength, operation, type, iteration)
  sumInternIteration <- summarise(sumInternIteration, value = sum(value, rm.na = TRUE))
  group <- group_by(sumInternIteration, numberProcessors,  operation, dToNRatio, type, stringLength)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$numberProcessors <- as.factor(valueMean$numberProcessors)
  if (work) {
  valueMean$numberProcessors <- as.character(valueMean$numberProcessors)
  valueMean$numberProcessors <- as.numeric(valueMean$numberProcessors)
  valueMean$value <- valueMean$value * valueMean$numberProcessors
  valueMean$numberProcessors <- as.factor(valueMean$numberProcessors)
  }
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(stringLength ~ dToNRatio)
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
  plot <- plot + geom_point(mapping = aes(x = numberProcessors, y = value, colour = operation),
                            size = pointSize)
  if (isD2N) {
  plot <- plot + facet_wrap(~dToNRatio)
  } else {
  plot <- plot + facet_wrap(~dToNRatio)
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
stringLengths <- unique(allData$stringLength)
print(stringLengths)
for (i in stringLengths) {
  localData <- filter(allData, stringLength == i)
  localData <- allDataFirstIter
#  operations <- c("firstBarrier", "secondBarrier") 
#  print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time Initial Barriers"))
  #operations <- c( "Splitter_move_to_pow_of_two_t_Barrier", "Splitter_shuffle_Barrier","Splitter_sortLocally_Barrier", "Splitter_median_select_Barrier", "Splitter_median_select_Barrier", "Splitter_partition_Barrier", "Splitter_exchange_Barrier", "Splitter_merge_Barrier", "Splitter_split_Barrier", "Splitter_baseCase_Barrier")
  #print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time Barriers first Iteration"))
  #operations <- c("Splitter_move_to_pow_of_two_t", "Splitter_shuffle","Splitter_sortLocally", "Splitter_median_select", "Splitter_median_select", "Splitter_partition", "Splitter_exchange", "Splitter_merge", "Splitter_split", "Splitter_baseCase")
  #print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time first Iteration"))
  #operations <- c("sorting_overall")
  #print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time first Iteration"))
  #print("2")
#  print(scatterPlot(localData, c("inbalance"), "number", 0.1, "inbalance"))
}

for (i in stringLengths) {
  #localData <- filter(allData, stringLength == i)
  localData <- allDataAllButFirstIter
#  operations <- c("firstBarrier", "secondBarrier") 
#  print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time Initial Barriers"))
  #operations <- c( "Splitter_move_to_pow_of_two_t_Barrier", "Splitter_shuffle_Barrier","Splitter_sortLocally_Barrier", "Splitter_median_select_Barrier", "Splitter_median_select_Barrier", "Splitter_partition_Barrier", "Splitter_exchange_Barrier", "Splitter_merge_Barrier", "Splitter_split_Barrier", "Splitter_baseCase_Barrier")
  #print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time Barriers, iterations > 0"))
  #operations <- c("Splitter_move_to_pow_of_two_t", "Splitter_shuffle","Splitter_sortLocally", "Splitter_median_select", "Splitter_median_select", "Splitter_partition", "Splitter_exchange", "Splitter_merge", "Splitter_split", "Splitter_baseCase")
  #print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time, iterations > 0"))
  operations <- c("sorting_overall")
  print(barPlotWhitelist(data_ = localData, operations_ = operations, type_ = "maxTime", title = "overall max Time, iterations > 0"))
  print("2")
#  print(scatterPlot(localData, c("inbalance"), "number", 0.1, "inbalance"))
} 
