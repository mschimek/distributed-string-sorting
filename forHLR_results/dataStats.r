library(ggplot2)
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(args[[1]])

colTypeSpec = cols(i = col_integer(),
       dist = col_double(),
       length_ = col_double(),
       dups = col_double())

# "Constants"
POINTSIZE = 0.1
globalSize <- 5000000

filename = paste(args[[1]],  sep="")
allData <- read_delim(file = filename, delim = " ", col_types = colTypeSpec, comment="-")

scatterPlot <- function(data_, pointSize, title) {
  plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = i, y = dist),
                            size = pointSize)
  #plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio )
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}
scatterPlotDups <- function(data_, pointSize, title) {
  plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = i, y = dups),
                            size = pointSize)
  #plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio )
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}

scatterPlotLength <- function(data_, pointSize, title) {
  plot <- ggplot(data = data_) 
  plot <- plot + geom_point(mapping = aes(x = i, y = length_),
                            size = pointSize)
  #plot <- plot + facet_wrap(ByteEncoder ~ dToNRatio )
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                       legend.text = element_text(size = 5),
                       legend.title = element_text(size =  6))
  plot <- plot + ggtitle(title)
  return(plot)
}
sum_ <- function(data_) {
  overallSum <- sum(data_$dist)
  d <- data.frame(data_)
  curSum <- 0
  for(index in 1:(nrow(data_) - 1)) {
    d$dist[index] <- overallSum - curSum
    curSum <- curSum + data_$dist[index]
  }
  return(as.tibble(d))
}
intervals <- function(data_) {
  
  maxValue <- max(data_$i)
  print(maxValue)
  d <- data.frame("i" = c(), "dist" = c())
  prevI <- 0
  for(index in c(0, 2,4,8,16,32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144)) {
  cur <- filter(data_, i >= index)

  curSum <- sum(cur$dist)
  curSumDups <- sum(cur$dups)

  d <- rbind(d, data.frame("i" = c(index), "dist" = c(curSum - curSumDups)))
  prevI <- index
    if (index > maxValue)
      return(as.tibble(d))
  }
  return(as.tibble(d))
}


pureDirName <- str_sub(args, start = 1, end = -2)
print(allData)
pdf(paste(pureDirName, "_plots_prefixCompression.pdf",sep=""), width=10, height=5)

scatterPlot(filter(allData, i < 50), POINTSIZE, "i < 50")
scatterPlotDups(filter(allData, i < 50), POINTSIZE, "i < 50")
scatterPlot(filter(allData, i < 200), POINTSIZE, "i < 200")
scatterPlot(filter(allData, i < 2000), POINTSIZE, "i < 2000")
scatterPlot(filter(allData, i < 20000), POINTSIZE, "i < 20000")
small <- filter(allData, i < 2000)
scatterPlot(small, POINTSIZE, "small dist")
scatterPlotDups(small, POINTSIZE, "small dups")
scatterPlot(intervals(small), POINTSIZE, "intervals on sums")
scatterPlot(sum_(small), POINTSIZE, "sum all dists for i > x")


