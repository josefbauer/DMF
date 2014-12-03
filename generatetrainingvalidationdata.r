generatetrainingvalidationdata <- function (sampledata,datasource) {
#generates training and validation part of the data and calculates user and item averages
#size of the training set and number of batches can be set below for each data set


if (datasource == "tasteprofile") {

trainingsamplesize <- 1200000
numberofbatches <- 4

}


if (datasource == "E-Commerce") {

trainingsamplesize <- 180000
numberofbatches <- 4

}


set.seed(2)
validationindices <- sample(1:length(sampledata[,1]), size = (length(sampledata[,1])-trainingsamplesize), replace = FALSE)


validationmatrix <- sampledata[validationindices,]
validationsamplesize <- length(validationmatrix[,1])

trainingmatrix <- sampledata[setdiff(1:length(sampledata[,1]),validationindices),]


trainingmatrixdataframe <- as.data.frame(trainingmatrix)
colnames(trainingmatrixdataframe) <- c("user","item","rating")

usermeanvalues <- aggregate(trainingmatrixdataframe$"rating", by = list(trainingmatrixdataframe$"user"), FUN = mean)
colnames(usermeanvalues) <- c("user","usermeanvalue")
itemmeanvalues <- aggregate(trainingmatrixdataframe$"rating", by = list(trainingmatrixdataframe$"item"), FUN = mean)
colnames(itemmeanvalues) <- c("item","itemmeanvalue")

usermedianvalues <- aggregate(trainingmatrixdataframe$"rating", by = list(trainingmatrixdataframe$"user"), FUN = median)
colnames(usermedianvalues) <- c("user","usermedianvalue")
itemmedianvalues <- aggregate(trainingmatrixdataframe$"rating", by = list(trainingmatrixdataframe$"item"), FUN = median)
colnames(itemmedianvalues) <- c("item","itemmedianvalue")


trainingdatamean <- mean(trainingmatrix[,3])
trainingdatamedian <- median(trainingmatrix[,3])

if (distribution == "normal") {biasaveragetype <- "mean"} else {biasaveragetype <- "median"}

if (biasaveragetype == "mean") {
predictionadjustment <- trainingdatamean
}
if (biasaveragetype == "median") {
predictionadjustment <- trainingdatamedian
}

if (biasaveragetype == "mean") {
userbiases <- usermeanvalues
itembiases <- itemmeanvalues
}

if (biasaveragetype == "median") {
userbiases <- usermedianvalues
itembiases <- itemmedianvalues
}

userbiases[,2] <- userbiases[,2] - predictionadjustment
itembiases[,2] <- itembiases[,2] - predictionadjustment
colnames(userbiases) <- c("user","userratingbias")
colnames(itembiases) <- c("item","itemratingbias")


numberofusers <- length(unique(sampledata[,1]))
numberofitems <- length(unique(sampledata[,2]))


batchsize <- length(trainingmatrix[,1])/numberofbatches


rm(trainingmatrixdataframe)
rm(usermeanvalues)
rm(itemmeanvalues)
rm(usermedianvalues)
rm(itemmedianvalues)
rm(trainingdatamean)
rm(trainingdatamedian)


localvariablelist <- ls(,envir=environment())

for (variable in localvariablelist) {
  
  eval(parse(text=paste0(variable,"<<-",variable)), envir = environment())
  
}


}

