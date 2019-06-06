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
pdfTitle = args[[2]]
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
  print(filters[[i]])
}

barPlotWhitelist <- function(datasets, operation_, type_ = "maxTime", title = " ", work = FALSE) {
  isD2N <- TRUE
  localSets <- vector("list", length(datasets))
  set <- "dummy"
  print(paste("length data set: ", length(datasets)))
  for (i in c(1:length(datasets))) {
    j <- ncol(filters[[i]][[1]]) 
    df = filters[[i]][[1]]
    print(df)
    names <- colnames(df)
    localSets[[i]] <- filter(data[[i]], operation == operation_, type == type_)
    for (j in c(2:length(names))){
      localSets[[i]] <- filter(localSets[[i]], UQ(as.name(names[j])) == df[1,j], operation == operation_)
    }
    localSets[[i]] <- mutate(localSets[[i]], name = df[1, 1])
    if (i == 1)
      set <- localSets[[i]]
    else {
      set <- rbind(set, localSets[[i]])
    }
   

    }
    set$numberProcessors <- as.factor(set$numberProcessors)
    groupByIterations <- group_by(set, numberProcessors, dToNRatio, samplePolicy, ByteEncoder, size, operation, type, name)
    valueMean <- summarise(groupByIterations, value = mean(value, rm.na = TRUE))
    print(as.data.frame(valueMean))
    valueMean$value <- valueMean$value / 10^9
  plot <- ggplot(data = valueMean, mapping = aes(x = numberProcessors, y = value, group = name, colour = name, shape = name))
  #plot <- plot + geom_bar(mapping = aes(x = numberProcessors, y = value, fill = name), stat="identity", position="dodge", colour="black", width=0.5)
  plot <- plot + ylab("time [sec]")
  plot <- plot + xlab("PEs")
  plot <- plot + theme_light()
  plot <- plot + geom_point()
  plot <- plot + geom_line()
  plot <- plot + theme(panel.margin = unit(-1.25, "lines"))
  if (isD2N) {
  plot <- plot + facet_wrap(~dToNRatio, labeller = label_both, nrow=1)
  } else {
  plot <- plot + facet_wrap(samplePolicy ~ ByteEncoder, labeller = label_both, nrow=1)
  }
  return (plot)

}

pdf(paste("evaluation/", pdfTitle, ".pdf",sep=""), width=10, height=5)
barPlotWhitelist(c(1,2), "sorting_overall")
