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

#allData <- read.table("./output.txt", comment.char = "#", col.names = columns)
allData <- read_table(file = "2018_12_02_H22_35/data.txt", col_types = 
  cols(numberProcessors = col_integer(),
       samplePolicy = col_character(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       time = col_double())
, comment="-")
allData$numberProcessors <- as.factor(allData$numberProcessors)
allDataWithoutIt1 <- filter(allData, iteration != 1)

availableProcessorCounts <- unique(allData$numberProcessors)


#examine all-to-all-strings

allToAllStrings <- filter(allDataWithoutIt1, operation == "all_to_all_strings")

allToAllScatterPlot1000000 <- ggplot(data = filter(allToAllStrings, size == 1000000)) + 
  geom_point(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type), position = "jitter") +
  facet_wrap(~ samplePolicy) +
  ggtitle("All-To-All Strings, StringSetSize = 1000000") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

allToAllBoxPlot1000000 <- ggplot(data = filter(allToAllStrings, size == 1000000)) + 
  geom_boxplot(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type)) +
  facet_wrap(~ samplePolicy) +
  ggtitle("All-To-All Strings, StringSetSize = 1000000") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

allToAllScatterPlot5000000 <- ggplot(data = filter(allToAllStrings, size == 5000000)) + 
  geom_point(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type), position = "jitter") +
  facet_wrap(~ samplePolicy) +
  ggtitle("All-To-All Strings, StringSetSize = 5000000") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

allToAllBoxPlot5000000 <- ggplot(data = filter(allToAllStrings, size == 5000000)) + 
  geom_boxplot(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type)) +
  facet_wrap(~ samplePolicy) +
  ggtitle("All-To-All Strings, StringSetSize = 5000000") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))






plot1 <- ggplot(data=filter(allData, numberProcessors==4)) + 
  geom_boxplot(mapping = aes(x=operation,y=time, color=type)) + 
  facet_wrap(~ samplePolicy)

plot2 <- ggplot(data=filter(allData, numberProcessors==128)) + 
  geom_boxplot(mapping = aes(x=operation,y=time, color=type)) + 
  facet_wrap(~ samplePolicy)

pdf(paste("plots.pdf",sep=""), width=10, height=5)

allToAllScatterPlot1000000
allToAllScatterPlot5000000
allToAllBoxPlot1000000
allToAllBoxPlot5000000

dev.off()
