createuseritemdata <- function (datasource) {

#Taste Profile
if (datasource == "tasteprofile") {
  
sampledata <- read.table(paste0(getwd(),"/data/tasteprofilesample.csv"), header = FALSE, sep = ",", colClasses = "integer")

}

#E-Commerce
if (datasource == "E-Commerce") {

sampledata <- read.table(paste0(getwd(),"/data/E-Commerce.csv"), header = TRUE, colClasses = "character", sep = ",")
colnames(sampledata) <- c("user","item","rating")
sampledata[,1] <- as.integer(sampledata[,1])
sampledata[,2] <- as.integer(sampledata[,2])
sampledata[,3] <- as.integer(sampledata[,3])

}

sampledata <<- sampledata

}

