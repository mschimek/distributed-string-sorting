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
       samplePolicy = col_character(),
       iteration = col_integer(),
       size = col_double(),
       operation = col_character(),
       type = col_character(),
       value = col_double())



pathToJSON = args[[1]]
jsonObject <- read_json(pathToJSON, simplifyDataFrame=TRUE)
print(str(jsonObject))
length <- length(jsonObject)
data <- vector("list", length)
filters <- vector("list", length)
print(length)
for (i in c(1:length)){
  print(i)
  curPath = jsonObject[[i]]$path
  filename = paste(curPath, "/data.txt", sep="")
  data[[i]] <- read_delim(file = filename, delim = "|", col_types = colTypeSpec, comment="-")
  data[[i]] <- filter(data[[i]], iteration != 0)
  filters[[i]] <- jsonObject[[i]]$filter
}
print(data)
print(filters)
# print(str(jsonObject[[i]]))
#  #filters[i] <- x["filter"][[i]]
#
#print("hier")
