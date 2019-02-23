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
       sizeInBytesAllgather = col_double(),
       sizeInBytesAllToAll = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())

# "Constants"
POINTSIZE = 0.1
globalSize <- 5000000

filename = paste(args[[1]], "/data.txt", sep="")
allData <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
allDataWithoutIt0 <- filter(allData, iteration !=  0)

barPlot <- function(data_, operations_, type_, title = " ") {
  data_$numberProcessors <- as.factor(data_$numberProcessors)
  data_$sizeInBytesAllToAll <- as.factor(data_$sizeInBytesAllToAll)

  filteredData <- filter(data_, !(operation %in% operations_), type == type_)
  group <- group_by(filteredData, numberProcessors, sizeInBytesAllgather, sizeInBytesAllToAll, operation, type)
  valueMean <- summarise(group, value = mean(value, rm.na = TRUE))
  plot <- ggplot(data = valueMean)
  plot <- plot + geom_bar(mapping = aes(x = sizeInBytesAllToAll, y = value, fill = operation), stat="identity")
  plot <- plot + facet_wrap(~ numberProcessors, labeller = label_both, nrow=1)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + ggtitle(title)
  return(plot)
}

pureDirName <- str_sub(args, start = 1, end = -2)
pdf(paste(pureDirName, "plots_strangeMpiBehaviour",sep=""), width=10, height=5)

operations <- c("")
barPlot(allDataWithoutIt0, operations, "avgTime", "mpi test")
operations <- c("alltoall")
barPlot(allDataWithoutIt0, operations, "avgTime", "mpi test")
operations <- c("alltoall", "allgatherv_1")
barPlot(allDataWithoutIt0, operations, "avgTime", "mpi test")
