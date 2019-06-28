library(ggplot2)
library(tidyverse)
library(magrittr)
library(jsonlite)
library(latex2exp)
library(cowplot)
library(gridExtra)
library(xtable)
library(ggpubr)

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
print(str(plots))
plots.data <- vector("list", numberPlots)
plots.filters <- vector("list", numberPlots)
plots.title <- vector("list", numberPlots)

for (k in c(1:numberPlots)) {
  jsonMetaObject <- plots[[k]]
  print(str(jsonMetaObject))
  jsonObject <- jsonMetaObject["data"][[1]]
  plots.title[[k]] <- jsonMetaObject["title"]
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
  #jsonMetaObject <- plots["CommonCrawlUnique"][[1]]
  print(str(jsonMetaObject))
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
    print(filename)
    data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
    data[[i]] <- filter(data[[i]], iteration != 0)
    filters[[i]] <- jsonObject[[i]]$filter
    print(filters[[i]])
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

createData <- function(data, filters, datasets, operations_, type_ = "maxTime", title = " ", work = FALSE) {
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
getMin <- function(data, competitors, dToN) {
  PEs <- unique(data$numberProcessors)
  numberPEs <- length(PEs)
  minValues <- rep(1000000, numberPEs) 
    for (name_ in competitors) {
      for (k in c(1:numberPEs)) {
        f <- filter(data, name == name_, numberProcessors == PEs[k], dToNRatio == dToN)
        if (nrow(f) > 0)
          minValues[k] <- min(minValues[k], f$value)
      }
    }
  return(minValues)
}

createSingleTable <-function(tablename, inputData1, inputData2, title1, title2) {
  dToN <- unique(inputData1$dToNRatio)
  print(dToN)
  data1 <- vector("list", length(dToN))
  data2 <- vector("list", length(dToN))
  for (i in c(1:length(dToN))) {
    data1[[i]] <- filter(inputData1, dToNRatio == dToN[i])
    data2[[i]] <- filter(inputData1, dToNRatio == dToN[i])
  }
  print("data1:")
  print(inputData1)
  print("data2:")
  print(inputData2)
  PEs <- unique(inputData1$numberProcessors)
  numberPEs <- length(PEs)

  sink(paste("./tables/", tablename, ".tex", sep=""))
  header <- paste("\\begin{tabular}{@{}l")
  for (i in c(1:numberPEs)) {
    header <- paste(header, "c", sep="")
  }
  header <- paste(header, "@{\\hskip 1cm}", sep="")
  for (i in c(1:numberPEs)) {
    header <- paste(header, "c", sep="")
  }
  header <- paste(header, "@{}}\n")
  subheader <- paste("& \\multicolumn{", numberPEs, "}{c}{\\textbf{", title1,"}} & \\multicolumn{", numberPEs, "}{c}{\\textbf{", title2, "}}\\\\\n",sep="")

  toprule <- paste("\\toprule\n")

  columNames <- paste("PEs")
  for (p in PEs) {
    columNames <- paste(columNames," & ", p )
  }
  for (p in PEs) {
    columNames <- paste(columNames," & ", p )
  }

  columNames <- paste(columNames, "\\\\\n")
  bottomRule <- paste("\\bottomrule\\\\\n")
  footer <- paste("\\end{tabular}\n")
  cat(header)
  cat(toprule)
  cat(subheader)
  cat(columNames)
  for (i in c(1:length(data1))) {
    multicolumn <- paste("\\multicolumn{", 2 * numberPEs + 1, "}{l}{\\textbf{r = ", dToN[i], "}}\\\\\n", sep="") 
    names <- unique(inputData1$name) 
    cat("\\hline\n")
    cat(multicolumn)
    cat("\\hline\n")
    for (name_ in names) {
      line <- ""
      line <- paste(line, name_)
      minValue <- 100000
      competitors <- c("kurpicz", "hQuick")
      minValues <- getMin(inputData1, competitors, dToN[i])
      for (k in c(1:numberPEs)) {
        f <- filter(inputData1, name == name_, numberProcessors == PEs[k], dToNRatio == dToN[i])
        if (nrow(f) == 0) {

          value <- format(round(f$value, 2), nsmall = 2)
          line <- paste(line, " & ", value, sep="")
        }else {
            if (minValues[k] > 100000){
          line <- paste(line, " & $\\infty$", sep="")}
          else {
            value <- minValues[k] / f$value 
            value <- format(round(value, 2), nsmall = 2)
            line <- paste(line, " & ", value, sep="")
          }
          
        }
      }
      minValues <- getMin(inputData2, competitors, dToN[i])
      for (k in c(1:numberPEs)) {
        f <- filter(inputData2, name == name_, numberProcessors == PEs[k], dToNRatio == dToN[i])
        if (nrow(f) == 0) {

          value <- format(round(f$value, 2), nsmall = 2)
          line <- paste(line, " & ", value, sep="")
        }else {
            if (minValues[k] > 100000){
          line <- paste(line, " & $\\infty$", sep="")}
          else {
            value <- minValues[k] / f$value 
            value <- format(round(value, 2), nsmall = 2)
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

data1 <- createData(plots.data[[1]], plots.filters[[1]],  c(1:length(plots.data[[1]])), c("sorting_overall"), "maxTime", plots.title[[1]])
data2 <- createData(plots.data[[2]], plots.filters[[2]],  c(1:length(plots.data[[2]])), c("sorting_overall"), "maxTime", plots.title[[2]])
createSingleTable(tablename, data1, data2, plots.title[[1]], plots.title[[2]])

