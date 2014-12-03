processuseritemdata <- function (sampledata,truncatethresholdmin,truncatethresholdmax,excludesparseuseritemsthreshold) {
#generates filtered data in the form (user_id,item_id,rating), all being integer

sampledata <- sampledata[(sampledata[,3] <= truncatethresholdmax)&(sampledata[,3] >= truncatethresholdmin),]

usercounttable <- table(sampledata[,1])
userstoexcludechar <- names(usercounttable[usercounttable <= excludesparseuseritemsthreshold])
itemcounttable <- table(sampledata[,2])
itemstoexcludechar <- names(itemcounttable[itemcounttable <= excludesparseuseritemsthreshold])
sampledata <- sampledata[!is.element(as.character(sampledata[,1]),userstoexcludechar)&!is.element(as.character(sampledata[,2]),itemstoexcludechar),]

sampledata[,1] <- as.integer(as.factor(sampledata[,1]))
sampledata[,2] <- as.integer(as.factor(sampledata[,2]))

sampledata <<- sampledata

}
