library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)
library(latex2exp)
library(cowplot)
library(gridExtra)
library(scales)
library(xtable)
source("legendSettings.r")

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
tablename = args[[2]]
plots <- read_json(pathToJSON, simplifyDataFrame=TRUE)

numberPlots <- length(plots)
plots.data <- vector("list", numberPlots)
plots.filters <- vector("list", numberPlots)
plots.title <- vector("list", numberPlots)
plots.size <- vector("list", numberPlots)

for (k in c(1:numberPlots)) {
  jsonMetaObject <- plots[[k]]
  jsonObject <- jsonMetaObject["data"][[1]]
  plots.title[[k]] <- jsonMetaObject["title"]
  plots.size[[k]] <- as.numeric(jsonMetaObject["numStrings"])
  length <- length(jsonObject)
  data <- vector("list", length)
  filters <- vector("list", length)
  for (i in c(1:length)){
    curPath = jsonObject[[i]]$path
    filename = paste(curPath, "/data.txt", sep="")
    data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
    data[[i]] <- filter(data[[i]], iteration != 0)
    filters[[i]] <- jsonObject[[i]]$filter
  }
  #jsonMetaObject <- plots["CommonCrawlUnique"][[1]]
  jsonObject <- jsonMetaObject["data"][[1]]
  title <- jsonMetaObject["title"]
  isD2N <- TRUE == (jsonMetaObject["isDToN"])
  length <- length(jsonObject)
  data <- vector("list", length)
  filters <- vector("list", length)
  for (i in c(1:length)){
    curPath = jsonObject[[i]]$path
    filename = paste(curPath, "/data.txt", sep="")
    data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
    data[[i]] <- filter(data[[i]], iteration != 0)
    filters[[i]] <- jsonObject[[i]]$filter
  }
  plots.data[[k]] <- data
  plots.filters[[k]] <- filters
}
numberExpandedDatasets <- function() {
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  return(counter)
}

createData <- function(datasets, filters,  operations_, type_ = "maxTime", title = " ") {
  set <- "dummy"
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  localSets <- vector("list", counter)
  counter <- 1
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      isD2N <-FALSE

      localSets[[counter]] <- filter(datasets[[i]], operation %in% operations_, type == type_)

      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
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
  set$value <- set$value / 10^9
  groupByIterations <- group_by(set, numberProcessors, dToNRatio,  operation, type, name)
  valueMean <- summarise(groupByIterations , value = mean(value, rm.na = TRUE))

  valueMean <- ungroup(valueMean)
  reduced <- select(valueMean, dToNRatio, value, name, numberProcessors) 
  return(reduced)
}
createDataMemory <- function(datasets, filters, numStrings) {
  set <- "dummy"
  counter <- 0
  for (i in c(1:length(datasets))) {
    numberFilters <- nrow(filters[[i]][[1]])
    counter <- counter + numberFilters
  }
  localSets <- vector("list", counter)
  counter <- 1
  isStrongScaling <- FALSE
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    numberFilters <- nrow(filters[[i]][[1]])
    for (k in c(1:numberFilters)){

      df = filters[[i]][[1]]
      names <- colnames(df)
      isD2N <-FALSE

      localSets[[counter]] <- filter(datasets[[i]], type == "number", rawCommunication == 1, phase != "none")
      if (nrow(localSets[[counter]]) > 0) {
        if (localSets[[counter]]$strongScaling[1] == 1) {
          isStrongScaling <- TRUE
        }
      }


      if (length(names) > 1) {
        for (j in c(2:length(names))){
          localSets[[counter]] <- filter(localSets[[counter]], UQ(as.name(names[j])) == df[k,j])
        }
      }
      localSets[[counter]] <- mutate(localSets[[counter]], name = df[k, 1])
      localSets[[counter]] <- select(localSets[[counter]], size, numberProcessors, dToNRatio,  phase, name, iteration, value)
      if (k == 1 && i == 1){
        set <- localSets[[counter]]
      }
      else {
        set <- rbind(set, localSets[[counter]])
      }
    }
    counter <- counter + 1
  }
  divisor <- 1
  if (!isStrongScaling) {
    divisor <- set$size * set$numberProcessors
  }else {
  divisor <- numStrings 
  }
  set$value <- set$value / (divisor)
  set$numberProcessors <- as.factor(set$numberProcessors)
  groupByIterations <- group_by(set, numberProcessors, dToNRatio,  phase, name)
  valueMean <- summarise(groupByIterations , value = mean(value, rm.na = TRUE))
  valueMean <- group_by(valueMean, numberProcessors, dToNRatio,  name)
  valueMean <- summarise(valueMean , value = sum(value))

  valueMean <- ungroup(valueMean)
  reduced <- select(valueMean, dToNRatio, value, numberProcessors, name) 
  return(reduced)
}

getMin <- function(data, competitors) {
  PEs <- unique(data$numberProcessors)
  numberPEs <- length(PEs)
  minValues <- rep(1000000, numberPEs) 
    for (name_ in competitors) {
      for (k in c(1:numberPEs)) {
        f <- filter(data, name == name_, numberProcessors == PEs[k])
        if (nrow(f) > 0)
          minValues[k] <- min(minValues[k], f$value)
      }
    }
  return(minValues)
}

createSingleTable <-function(tablename, inputData, communicationData, titles) {
  PEs <- unique(inputData[[1]]$numberProcessors)
  numberPEs <- length(PEs)

  sink(paste("./tables/", tablename, ".tex", sep=""))
  header <- paste("\\begin{tabular}{@{}l")
  for (i in c(1:numberPEs)) {
    header <- paste(header, "ll", sep="") 
  }
  header <- paste(header, "@{}}\n",sep="")
  toprule <- paste("\\toprule\n")
  columNames <- paste("PEs & $\\varnothing$ ")
  for (p in PEs) {
    columNames <- paste(columNames," & \\multicolumn{1}{c}{", p, "}", sep="" )
  }
  columNames <- paste(columNames, "\\\\\n")
  bottomRule <- paste("\\bottomrule\\\\\n")
  footer <- paste("\\end{tabular}\n")
  cat(header)
  cat(toprule)
  cat(columNames)
  for (i in c(1:length(inputData))) {
    curInputData <- inputData[[i]]
    curCommunicationData <- communicationData[[i]]
    multicolumn <- paste("& \\multicolumn{", numberPEs, "}{l}{\\textbf{", titles[[i]], "}}\\\\\n", sep="") 
    names <- unique(curInputData$name) 
    cat("\\hline\n")
    cat(multicolumn)
    cat("\\hline\n")
    minValues <- rep(100000, numberPEs) 
    minValuesC <- rep(100000, numberPEs) 
    for (name_ in names) {
      for (k in c(1:numberPEs)) {
        f <- filter(curInputData, name == name_, numberProcessors == PEs[k])
        c <- filter(curCommunicationData, name == name_, numberProcessors == PEs[k])
        if (nrow(f) > 0)
          minValues[k] <- min(minValues[k], f$value)
        if (nrow(c) > 0)
          minValuesC[k] <- min(minValuesC[k], c$value)
      }
    }
    names <- sort(names)
    for (name_ in names) {
        line <- ""
        line <- paste(line, name_)
        competitors <- c("kurpicz", "hQuick")
      minValues <- getMin(curInputData, competitors) 
      avg <- 0
      counter <- 0
      #avg
      for (k in c(1:numberPEs)) {
        f <- filter(curInputData, name == name_, numberProcessors == PEs[k])

        value <- format(round(f$value, 1), nsmall = 2)
        if (nrow(f) == 0) {

        }else {

          if( minValues[k] < 10000){
          value <- minValues[k] / f$value
          avg <- avg + value
          counter <- counter + 1
          }


        }
      }
      avg <- avg / counter
      avg <- format(round(avg, 1), nsmall=1)
        line <- paste(line, " & ",  avg)
      for (k in c(1:numberPEs)) {
        f <- filter(curInputData, name == name_, numberProcessors == PEs[k])
        c <- filter(curCommunicationData, name == name_, numberProcessors == PEs[k])

        value <- format(round(f$value, 1), nsmall = 2)
        commValue <- format(round(c$value, 0), nsmall = 0)
        if (nrow(f) == 0) {

          line <- paste(line, " & ", value, sep="")
        }else {

          runtime <- paste(" & ")
          if( minValues[k] >= 10000)
            runtime <- paste(runtime, "$-$")
          else {
          value <- minValues[k] / f$value
        value <- format(round(value, 1), nsmall = 1)
            runtime <- paste(runtime, value, " ", sep="")
          }

          line <- paste(line, runtime, sep="")   

        }
      }
      line <- paste(line, "\\\\\n")
      cat(line)
    }
  }
  cat(bottomRule)
  cat(footer)
}

numberPlots <- length(plots.data)
data.runtime <- vector("list", numberPlots)
data.communication <- vector("list", numberPlots)
for (i in c(1:numberPlots)) {
       data.runtime[[i]] <- createData(plots.data[[i]], plots.filters[[i]], c("sorting_overall"), "maxTime", title)
       data.communication[[i]] <- createDataMemory(plots.data[[i]], plots.filters[[i]], plots.size[[i]])
}
createSingleTable(tablename, data.runtime, data.communication, plots.title)


