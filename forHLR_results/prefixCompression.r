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
globalSize <- 5000000

filename = paste(args[[1]], "/data.txt", sep="")
allData <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
allData <- filter(allData)#, iteration == 0)
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1_Timer <- filter(allData,  Timer == "Timer", iteration != 0)
allDataWithoutIt1_EmptyTimer <- filter(allData,  Timer == "EmptyTimer", iteration != 0)

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

barPlot <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

efficiency_sum <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, iteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueSum <- summarise(group, value=sum(value, rm.na = TRUE))
  ungroup(valueSum)
  group <- group_by(valueSum, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  ungroup(valueMean)
  valueMean$size <- as.factor(valueMean$size)
  valueMean <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, type)
  
  myFunc <- function(vec) {
   tmp <- filter(vec, numberProcessors == 1)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
   #vec$numberProcessors <- as.numeric(vec$numberProcessors)
   print(vec)
   print(vec$numberProcessors)
   print(vec$value)
   print(tmp)
   vec$value <- (tmp / vec$value) / (as.numeric(vec$numberProcessors))
   print(vec$value)
   print("-----")
   return(vec)
  } 
  valueMean <- do(valueMean, myFunc(.))
  print(valueMean)

  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  #plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

speedup_sum <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, iteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueSum <- summarise(group, value=sum(value, rm.na = TRUE))
  ungroup(valueSum)
  group <- group_by(valueSum, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  ungroup(valueMean)
  valueMean$size <- as.factor(valueMean$size)
  valueMean <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, type)
  
  myFunc <- function(vec) {
   tmp <- filter(vec, numberProcessors == 1)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
   #vec$numberProcessors <- as.numeric(vec$numberProcessors)
   print(vec)
   print(vec$numberProcessors)
   print(vec$value)
   print(tmp)
   vec$value <- (tmp / vec$value) 
   print(vec$value)
   print("-----")
   return(vec)
  } 
  valueMean <- do(valueMean, myFunc(.))
  print(valueMean)

  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  #plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

sumData <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, iteration, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueSum <- summarise(group, value=sum(value, rm.na = TRUE))
  ungroup(valueSum)
  group <- group_by(valueSum, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  ungroup(valueMean)
  valueMean$size <- as.factor(valueMean$size)
  valueMean <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, type)
  
  myFunc <- function(vec) {
   tmp <- filter(vec, numberProcessors == 1)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
   #vec$numberProcessors <- as.numeric(vec$numberProcessors)
   print(vec)
   print(vec$numberProcessors)
   print(vec$value)
   print(tmp)
   vec$value <- (tmp / vec$value) 
   print(vec$value)
   print("-----")
   return(vec)
  } 
#  valueMean <- do(valueMean, myFunc(.))
  print(valueMean)

  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  #plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}
speedup <- function(data_, operation_, type_, size_, title = " ") {
  filteredData <- filter(data_, operation == operation_, type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  valueMean <- group_by(valueMean, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  
  myFunc <- function(vec) {
   tmp <- filter(vec, numberProcessors == 2)$value
  vec$numberProcessors <- as.numeric(as.character(vec$numberProcessors))
   #vec$numberProcessors <- as.numeric(vec$numberProcessors)
   print(vec)
   print(vec$numberProcessors)
   print(vec$value)
   print(tmp)
   vec$value <- (tmp / vec$value) #/ (as.numeric(vec$numberProcessors))
   print(vec$value)
   print("-----")
   return(vec)
  } 
  valueMean <- do(valueMean, myFunc(.))
  print(valueMean)

  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, colour = dToNRatio, shape=ByteEncoder))
  plot <- plot + geom_line()
  plot <- plot + geom_point()
  #plot <- plot + facet_wrap(~ dToNRatio, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

barPlot <- function(data_, operations_, type_, size_, title = " ") {
  filteredData <- filter(data_, !(operation %in% operations_), type == type_, size == size_)
  filteredData$dToNRatio <- as.factor(filteredData$dToNRatio)
  group <- group_by(filteredData, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  valueMean$size <- as.factor(valueMean$size)
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = ByteEncoder, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(numberProcessors ~ dToNRatio, labeller = label_both, nrow=1)
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

pureDirName <- str_sub(args, start = 1, end = -2)
pdf(paste(pureDirName, "_plots_prefixCompression.pdf",sep=""), width=10, height=5)

efficiency_sum(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "sort_locally"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "efficiency sorting")
efficiency_sum(data_ = filter(allDataWithoutIt1_EmptyTimer, operation != "sorting_overall"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "efficiency overall")
speedup_sum(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "sort_locally"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "speedup sorting")
speedup_sum(data_ = filter(allDataWithoutIt1_EmptyTimer, operation != "sorting_overall"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "speedup overall")
sumData(data_ = filter(allDataWithoutIt1_EmptyTimer, operation == "sort_locally"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "sorting")
sumData(data_ = filter(allDataWithoutIt1_EmptyTimer, operation != "sorting_overall"), operation_ = "blabla", type_ = "avgTime", size_ = globalSize, title = "sum overall")
#numAllgatherBytesSent(filter(allDataWithoutIt1_EmptyTimer, numberProcessors %in% c(2,4) ))
  #numAllgatherBytesSent(filter(allDataWithoutIt1_EmptyTimer, numberProcessors %in% c(8, 16) ))
  #
  #operations <- c("all_to_all_strings")
  #data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
  #scatterPlotAllProcessors(data, operations, "avgTime")
  #
  #
  #operations <- c("allgather_splitters")
  #data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
  #scatterPlotAllProcessors(data, operations, "avgTime")
  #
  #operations <- c("allgatherv_test_before")
  #scatterPlotAllProcessors(data, operations, "avgTime")
#operations <- c("allgather_test_before")
#scatterPlotAllProcessors(data, operations, "avgTime")
#operations <- c("allgatherv_test_after")
#scatterPlotAllProcessors(data, operations, "avgTime")
#operations <- c("allgather_test_after")
#scatterPlotAllProcessors(data, operations, "avgTime")



operations <- c("all_to_all_strings")
data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
#scatterPlotAllProcessors(data, operations, "avgTime")
#scatterPlotAllProcessors(data, operations, "maxTime")
#scatterPlotAllProcessors(data, operations, "avgLoss")
#scatterPlotAllProcessors(data, operations, "maxLoss")

operations <- c("all_to_all_strings_intern_copy", "all_to_all_strings_read", "all_to_all_strings_mpi")
data <- filter(allDataWithoutIt1_Timer, size == globalSize)
#scatterPlotAllProcessors(data, operations, "avgTime")
#scatterPlotAllProcessors(data, operations, "maxTime")
#scatterPlotAllProcessors(data, operations, "avgLoss")
#scatterPlotAllProcessors(data, operations, "maxLoss")

  operations <- c("prefix_decompression")
  
  data <- filter(allDataWithoutIt1_EmptyTimer, size == globalSize)
  #scatterPlotAllProcessors(data, operations, "avgTime")
  
  barPlot(data_ = filter(allDataWithoutIt1_EmptyTimer, dToNRatio == 0.0), operations_ = c("sorting_overall", "prefix_decompression", "all_to_all_strings", "merge_ranges", "compute_interval_sizes"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
  barPlot(data_ = filter(allDataWithoutIt1_EmptyTimer, dToNRatio == 0.0), operations_ = c("sorting_overall"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
  barPlot(data_ = allDataWithoutIt1_EmptyTimer, operations_ = c("sorting_overall"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
  barPlot(data_ = allDataWithoutIt1_EmptyTimer,  operations_ = c("sorting_overall", "prefix_decompression"), type_ = "avgTime", size_ = globalSize, title = "overall Time")
  barPlot(data_ = allDataWithoutIt1_EmptyTimer, operations_ = c("sorting_overall", "prefix_decompression", "sort_locally"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #
  barPlot(data_ = allDataWithoutIt1_EmptyTimer, operations_ = c("sorting_overall", "prefix_decompression", "sort_locally", "all_to_all_strings", "merge_ranges", "compute_interval_sizes"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #
  barPlot(data_ = allDataWithoutIt1_EmptyTimer, operations_ = c("choose_splitters","allgatherv_test_before", "sorting_overall", "prefix_decompression", "sort_locally", "all_to_all_strings", "merge_ranges", "compute_interval_sizes"), type_ = "avgTime", size_ = globalSize, title = "overall Time") #

#numBytesSent(allDataWithoutIt1_Timer)
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
