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

#rawData <- read.table("./output.txt", comment.char = "#", col.names = columns)
rawData <- read_table(file = "2018_12_02/data.txt", col_names = TRUE, comment="-")

availableProcessorCounts <- unique(rawData$numberProcessors)

#examine all-to-all-strings

allToAllStrings <- filter(rawData, operation == "all_to_all_strings", size == 1000000)
allToAllStrings$numberProcessors <- as.factor(allToAllStrings$numberProcessors)

allToAllScatterPlot <- ggplot(data = allToAllStrings) + 
                geom_point(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type), position = "jitter") +
                facet_wrap(~ samplePolicy) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

allToAllBoxPlot <- ggplot(data = allToAllStrings) + 
                geom_boxplot(mapping = aes(x = interaction(numberProcessors, type), y = time, colour = type)) +
                facet_wrap(~ samplePolicy) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))





plot1 <- ggplot(data=filter(rawData, numberProcessors==4)) + 
          geom_boxplot(mapping = aes(x=operation,y=time, color=type)) + 
          facet_wrap(~ samplePolicy)
        
plot2 <- ggplot(data=filter(rawData, numberProcessors==128)) + 
          geom_boxplot(mapping = aes(x=operation,y=time, color=type)) + 
          facet_wrap(~ samplePolicy)
