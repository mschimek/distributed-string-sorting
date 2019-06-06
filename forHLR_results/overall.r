library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

colTypeSpec = cols(numberProcessors = col_integer(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())



pathToJSON = args[[1]]
pdfTitle = args[[2]]
jsonMetaObject <- read_json(pathToJSON, simplifyDataFrame=TRUE)
jsonObject <- jsonMetaObject["data"][[1]]
title <- jsonMetaObject["title"]
isD2N <- TRUE == (jsonMetaObject["isDToN"])
print(str(jsonObject))
length <- length(jsonObject)
data <- vector("list", length)
filters <- vector("list", length)
print(paste("number datasets: ", length))
for (i in c(1:length)){
  print(i)
  curPath = jsonObject[[i]]$path
  filename = paste(curPath, "/data.txt", sep="")
  data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
  data[[i]] <- filter(data[[i]], iteration != 0)
  filters[[i]] <- jsonObject[[i]]$filter
  print(filters[[i]])
}
numberExpandedDatasets <- function() {
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  return(counter)
}

barPlotWhitelist <- function(datasets, operation_, type_ = "maxTime", title = " ", work = FALSE) {
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  print(paste("counter: ", counter))
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    print(paste("numberFilters: ", numberFilters))
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      #print(data[[i]])
      localSets[[counter]] <- filter(data[[i]], operation == operation_, type == type_)
      #print(localSets[[counter]])
      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      #print(localSets[[counter]])
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], numberProcessors, dToNRatio,  operation, type, name, iteration, value)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
    }
    set$numberProcessors <- as.factor(set$numberProcessors)
    groupByIterations <- group_by(set, numberProcessors, dToNRatio,  operation, type, name)
    valueMean <- summarise(groupByIterations, value = mean(value, rm.na = TRUE))
    valueMean$value <- valueMean$value / 10^9
  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, group = name, colour = name, shape = name, linetype = name))
  #plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = name), stat="identity", position="dodge", colour="black", width=0.5)
  plot <- plot + ylab("time [sec]")
  plot <- plot + xlab("PEs")
  plot <- plot + theme_light()
  plot <- plot + geom_point()
  plot <- plot + geom_line()
  plot <- plot + scale_y_continuous(trans='log2')
  plot <- plot + theme(panel.spacing = unit(1.25, "lines"))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5))
  if (isD2N) {
  plot <- plot + facet_wrap(~dToNRatio, labeller = label_both, nrow=1)
  } 
  plot <- plot + ggtitle(title)
  return (plot)

}

pdf(paste("evaluation/", pdfTitle, ".pdf",sep=""), width=10, height=5)
barPlotWhitelist(c(1:length(data)), "sorting_overall", "maxTime", title)
