library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)
library(latex2exp)
library(cowplot)
library(gridExtra)
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

createData <- function(datasets, operations_, type_ = "maxTime", title = " ", work = FALSE) {
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
      isD2N <-FALSE

      localSets[[counter]] <- filter(data[[i]], operation %in% operations_, type == type_)

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
  set$value <- set$value / 10^9
  groupByIterations <- group_by(set, numberProcessors, dToNRatio,  operation, type, name)
  valueMean <- summarise(groupByIterations , value = mean(value, rm.na = TRUE))

  valueMean <- ungroup(valueMean)
  reduced <- select(valueMean, dToNRatio, value, name, numberProcessors) 
  return(reduced)
}

createSingleTable <-function(tablename, inputData) {
  dToN <- unique(inputData$dToNRatio)
  print(dToN)
  data <- vector("list", length(dToN))
  for (i in c(1:length(dToN))) {
    data[[i]] <- filter(inputData, dToNRatio == dToN[i])
  }
  PEs <- unique(inputData$numberProcessors)
  numberPEs <- length(PEs)

  sink(paste("./tables/", tablename, ".tex", sep=""))
  header <- paste("\\begin{tabular}{@{}lcccccc@{}}\n")
  toprule <- paste("\\toprule\n")
  columNames <- paste("PEs")
  for (p in PEs) {
    columNames <- paste(columNames," & ", p )
  }
  columNames <- paste(columNames, "\\\\\n")
  bottomRule <- paste("\\bottomrule\\\\\n")
  footer <- paste("\\end{tabular}\n")
  cat(header)
  cat(toprule)
  cat(columNames)
  for (i in c(1:length(data))) {
    multicolumn <- paste("& \\multicolumn{", numberPEs, "}{l}{\\textbf{r = ", dToN[i], "}}\\\\\n", sep="") 
    names <- unique(inputData$name) 
    cat("\\hline\n")
    cat(multicolumn)
    cat("\\hline\n")
    for (name_ in names) {
      line <- ""
      line <- paste(line, name_)
      minValue <- 100000
      for (p in PEs) {
        f <- filter(inputData, name == name_, numberProcessors == p, dToNRatio == dToN[i])
        if (nrow(f) > 0)
          minValue <- min(minValue, f$value)
      }
      for (p in PEs) {
        f <- filter(inputData, name == name_, numberProcessors == p, dToNRatio == dToN[i])
        if (nrow(f) == 0) {

          value <- format(round(f$value, 2), nsmall = 2)
          line <- paste(line, " & ", value, sep="")
        }else {
          if (minValue == f$value) {
            value <- format(round(f$value, 2), nsmall = 2)
            line <- paste(line, " & \\textbf{", value,"}", sep="")
          }else {
            value <- format(round(f$value, 2), nsmall = 2)
            line <- paste(line, " & ", value, sep="")
          }
        }
      }
      line <- paste(line, "\\\\\n")
      cat(line)
    }
  }
  cat(bottomRule)
  cat(footer)
}

l <- createData(c(1:length(data)), c("sorting_overall"), "maxTime", title)
createSingleTable(tablename, l)

