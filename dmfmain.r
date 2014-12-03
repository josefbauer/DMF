####SETTINGS####
#rm(list = ls())

#set working directory (where files and data are stored)
setwd("C:/DMF")


#data set to use
datasource <- "tasteprofile"

#distribution templates
#distribution <- "normal"
distribution <- "lognormal"
#distribution <- "poisson"
#distribution <- "gamma"
#distribution <- "pareto"


#number of epochs
maxepoch <- 300

#exclude count values larger than
truncatethresholdmax <- 300

#exclude count values smaller than
truncatethresholdmin <- 3

#exclude users and items with observations less or equal than
excludesparseuseritemsthreshold <- 9


#distribution settings (grid search as well as parameter ranges and initial values)
#grid search: parametervaluematrix is the matrix of parameter values to be iterated through (learning rate, regularization parameter, momentum, number of features, 
#initial value parameter, additional distribution parameter)
#(note that the learning rate specified here is divided by the batch size internally)
#values for each parameter are separated by a comma
#
#note that in the definition of the densities the variable parameter is denoted by theta, vartheta is the hyperparameter, 
#and xval denotes the variable argument at which the density is evaluated


if (distribution == "normal"){
densityexpression <- "1/sqrt(2*pi)*exp(-(xval-theta)^2/2)"
gtosolveexpression <- "theta - expval"
expvaltosolveexpression <- "theta - expval"
useritemfeatmodulationexpression <- "xval"
useritembiasmodulationexpression <- "xval"
#distribution settings (parameter ranges, initial values and grid search)
usebiasmodminmax <- TRUE
linnergmin <- 1
linnergmax <- Inf
useadj0 <- FALSE
initcoeff <- 1
useposinit <- FALSE
parametervaluematrix <- expand.grid(c(50), c(0.25), c(0.4), c(20), c(0), c(1))
}


if (distribution == "lognormal"){
densityexpression <- "1/(xval*sqrt(2*pi)*vartheta)*exp(-(log(xval)-theta)^2/(2*vartheta^2))"
gtosolveexpression <- "theta - log(expval)"
expvaltosolveexpression <- "exp(theta) - expval"
useritemfeatmodulationexpression <- "xval"
useritembiasmodulationexpression <- "log(xval)"
#distribution settings (parameter ranges, initial values and grid search)
usebiasmodminmax <- TRUE
linnergmin <- 1
linnergmax <- Inf
useadj0 <- FALSE
initcoeff <- 0.2
useposinit <- FALSE
parametervaluematrix <- expand.grid(c(500), c(0.01), c(0.4), c(20), c(0), c(0.5))
}


if (distribution == "poisson"){
#unshifted - note that if you want to shift the support of the distribution (e.g., by replacing xval with (xval-3)) as described in the paper, you also need to add the corresponding value in the predictions in the code below
densityexpression <- "theta^xval/gamma(xval+1)*exp(-theta)"
gtosolveexpression <- "theta - expval"
expvaltosolveexpression <- "theta - expval"
useritemfeatmodulationexpression <- "xval"
useritembiasmodulationexpression <- "xval"
#distribution settings (parameter ranges and initial values)
usebiasmodminmax <- FALSE
linnergmin <- 1
linnergmax <- Inf
useadj0 <- TRUE
initcoeff <- 4
useposinit <- TRUE
parametervaluematrix <- expand.grid(c(200,500), c(0.2), c(0.4), c(20), c(0), c(1))
}


if (distribution == "gamma"){
densityexpression <- "1/(gamma(vartheta)*theta^vartheta)*xval^(vartheta-1)*exp(-xval/theta)"
gtosolveexpression <- "theta - expval"
expvaltosolveexpression <- "theta - expval"
useritemfeatmodulationexpression <- "xval"
useritembiasmodulationexpression <- "xval"
#distribution settings (parameter ranges and initial values)
usebiasmodminmax <- FALSE
linnergmin <- 1
linnergmax <- Inf
useadj0 <- TRUE
initcoeff <- 0.01
useposinit <- TRUE
parametervaluematrix <- expand.grid(c(2000), c(0.01), c(0.4), c(20), c(0), c(2))
}


if (distribution == "pareto"){
densityexpression <- "theta*vartheta^theta/xval^(theta+1)"
gtosolveexpression <- "theta - expval"
expvaltosolveexpression <- "theta - expval"
useritemfeatmodulationexpression <- "xval"
useritembiasmodulationexpression <- "xval"
#distribution settings (parameter ranges and initial values)
usebiasmodminmax <- FALSE
linnergmin <- 0.1
linnergmax <- Inf
useadj0 <- TRUE
initcoeff <- 1.2
useposinit <- TRUE
parametervaluematrix <- expand.grid(c(20,100), c(0.15), c(0.4), c(20), c(0), c(3))
}



####COMPUTATION####

library(compiler)

enableJIT(3)

library(inline)

library(SuppDists)

library(VGAM)

#C++ helper function for gradient descent update calculation
cpp_grad_f_b <- "
  Rcpp::NumericVector v(a);
  Rcpp::NumericMatrix M(B);
  Rcpp::NumericMatrix N(C);
  int size_v = v.size();
  for (int i=0; i < size_v; i++){
		M(v[i]-1,_) = M(v[i]-1,_) + N(i,_);
  }
  return M;
"
cpp_grad_f <- cxxfunction(signature(a="numeric", B="numeric", C="numeric"), body = cpp_grad_f_b, plugin = "Rcpp")



source(paste(getwd(),"/autodiffsymfunctions.r",sep=""))

.Jython$exec("from __future__ import division")


generatevectorizedfunction(densityexpression,c("xval","theta"),"fdensityvec")

sympydensityexpressionstring <- convertexpressionrtosympystring(densityexpression)

generatepartialderivative(sympydensityexpressionstring,"theta","fdensitypartialderivthetaexpression",simplify=TRUE)

generatevectorizedfunction(fdensitypartialderivthetaexpression,c("xval","theta"),"fdensitypartialderivthetavec")

sympygtosolveexpressionstring <- convertexpressionrtosympystring(gtosolveexpression)

generateexpressionsolution(sympygtosolveexpressionstring,"theta","gtrafosolutionexpression")

generatevectorizedfunction(gtrafosolutionexpression,c("expval"),"gtrafovec")

sympygexpressionstring <- convertexpressionrtosympystring(gtrafosolutionexpression)

generatepartialderivative(sympygexpressionstring,"expval","gtrafopartialderivxvalexpression",simplify=TRUE)

generatevectorizedfunction(gtrafopartialderivxvalexpression,c("expval"),"gtrafopartialderivxvalvec")

sympyexpvaltosolveexpressionstring <- convertexpressionrtosympystring(expvaltosolveexpression)

generateexpressionsolution(sympyexpvaltosolveexpressionstring,"expval","expvalexpression")

generatevectorizedfunction(expvalexpression,c("theta"),"expvalfunctionvec")

generatelinnergexpression(useritemfeatmodulationexpression,useritembiasmodulationexpression,"xval","linnergexpression")

generatevectorizedfunction(linnergexpression,c("featval","biasval"),"linnergvec")

sympylinnergexpressionstring <- convertexpressionrtosympystring(linnergexpression)

generatepartialderivative(sympylinnergexpressionstring,"featval","linnergpartialderivfeatvalexpression",simplify=TRUE)
generatevectorizedfunction(linnergpartialderivfeatvalexpression,c("featval","biasval"),"linnergpartialderivfeatvalvec")
generatepartialderivative(sympylinnergexpressionstring,"biasval","linnergpartialderivbiasvalexpression",simplify=TRUE)
generatevectorizedfunction(linnergpartialderivbiasvalexpression,c("featval","biasval"),"linnergpartialderivbiasvalvec")

generatealgebraiccomposition("1",convertexpressionrtosympystring(densityexpression),"/","oneodensityexpression")
generatealgebraiccomposition(convertexpressionrtosympystring(oneodensityexpression),convertexpressionrtosympystring(fdensitypartialderivthetaexpression),"*","fdensitypartderivthetadivivedbyfdensityexpression")
generatevectorizedfunction(fdensitypartderivthetadivivedbyfdensityexpression,c("xval","theta"),"fdensitypartderivthetadivivedbyfdensityvec")


source(paste(getwd(),"/createuseritemdata.r",sep=""))

createuseritemdata(datasource)

source(paste(getwd(),"/processuseritemdata.r",sep=""))

processuseritemdata(sampledata,truncatethresholdmin,truncatethresholdmax,excludesparseuseritemsthreshold)

source(paste(getwd(),"/generatetrainingvalidationdata.r",sep=""))

generatetrainingvalidationdata(sampledata,datasource)



minvalidationerrorvalue <- 1000
minvalidationerrorvalueparamlog <- "not yet defined"



for (parindex in 1:length(parametervaluematrix[,1])){

#learning rate
epsilon <- parametervaluematrix[parindex,1]

#regularization parameter
lambda <- parametervaluematrix[parindex,2]

#momentum parameter
momentum <- parametervaluematrix[parindex,3]

numberoffeatures <- parametervaluematrix[parindex,4]

initialvalueparameter <- parametervaluematrix[parindex,5]

vartheta <- parametervaluematrix[parindex,6]


cat("\n\n")
print(sprintf("epsilon: %f  lambda: %f  momentum: %f  number of features: %i", epsilon, lambda, momentum, numberoffeatures))


if (useposinit){
itemfeaturematrix <- sqrt(1/numberoffeatures)*(sqrt(initialvalueparameter) + initcoeff*matrix(abs(rnorm(numberofitems*numberoffeatures, mean = 0, sd = 1)), numberofitems, numberoffeatures))
userfeaturematrix <- sqrt(1/numberoffeatures)*(sqrt(initialvalueparameter) + initcoeff*matrix(abs(rnorm(numberofusers*numberoffeatures, mean = 0, sd = 1)), numberofusers, numberoffeatures))
} else{
itemfeaturematrix <- sqrt(1/numberoffeatures)*(sqrt(initialvalueparameter) + initcoeff*matrix(rnorm(numberofitems*numberoffeatures, mean = 0, sd = 1), numberofitems, numberoffeatures))
userfeaturematrix <- sqrt(1/numberoffeatures)*(sqrt(initialvalueparameter) + initcoeff*matrix(rnorm(numberofusers*numberoffeatures, mean = 0, sd = 1), numberofusers, numberoffeatures))
}

itembiasmatrix <- matrix(0, numberofitems, 1)

userbiasmatrix <- matrix(0, numberofusers, 1)


if (useadj0){
predictionadjustment <- 0
}

itemfeaturematrix_inc <- matrix(0, numberofitems, numberoffeatures)

userfeaturematrix_inc <- matrix(0, numberofusers, numberoffeatures)

itembiasmatrix_inc <- matrix(0, numberofitems, 1)

userbiasmatrix_inc <- matrix(0, numberofusers, 1)

validationerror <- rep(0, maxepoch)



for (epoch in 1:maxepoch){

randomtrainingdatapermutation <- sample(1:trainingsamplesize, size = trainingsamplesize, replace = FALSE)

trainingmatrix <- trainingmatrix[randomtrainingdatapermutation,]


for (batch in 1:numberofbatches){

#computation of the predictions on the training set

	userbatch <- trainingmatrix[((batch-1)*batchsize+1):(batch*batchsize),1]

	itembatch <- trainingmatrix[((batch-1)*batchsize+1):(batch*batchsize),2]
	
	ratingbatch <- trainingmatrix[((batch-1)*batchsize+1):(batch*batchsize),3]

	featdotprodbatch <- rowSums(itemfeaturematrix[itembatch,]*userfeaturematrix[userbatch,])

	#computation of gradients
	
	biasbatch <- predictionadjustment + userbiasmatrix[userbatch,] + itembiasmatrix[itembatch,]

	if (usebiasmodminmax) {
	biasbatch <- pmin(pmax(biasbatch,truncatethresholdmin),truncatethresholdmax)
	}

	linnergvaluevec <- linnergvec(featdotprodbatch,biasbatch,vartheta)

	linnergvaluevec <- pmin(pmax(linnergvaluevec,linnergmin),linnergmax)

	thetavalvec <- gtrafovec(linnergvaluevec,vartheta)
  
	D_useritemcoefficient <- fdensitypartderivthetadivivedbyfdensityvec(ratingbatch,thetavalvec,vartheta)*
									gtrafopartialderivxvalvec(linnergvaluevec,vartheta)
  
	M_base_feat <- kronecker(-D_useritemcoefficient*linnergpartialderivfeatvalvec(featdotprodbatch,biasbatch,vartheta), matrix(1,1,numberoffeatures))
	
	M_item <- M_base_feat*userfeaturematrix[userbatch,] + 2*lambda*itemfeaturematrix[itembatch,]
	
	M_user <- M_base_feat*itemfeaturematrix[itembatch,] + 2*lambda*userfeaturematrix[userbatch,]
  
	M_base_bias <- as.matrix(-D_useritemcoefficient*linnergpartialderivbiasvalvec(featdotprodbatch,biasbatch,vartheta))
	
	M_itembias <- M_base_bias + 2*lambda*itembiasmatrix[itembatch,]
	
	M_userbias <- M_base_bias + 2*lambda*userbiasmatrix[userbatch,]

	M_item[is.na(M_item)|!is.finite(M_item)] <- 0
	M_user[is.na(M_user)|!is.finite(M_user)] <- 0
	M_itembias[is.na(M_itembias)|!is.finite(M_itembias)] <- 0
	M_userbias[is.na(M_userbias)|!is.finite(M_userbias)] <- 0
  
	
	d_itemfeaturematrix <- matrix(0,numberofitems,numberoffeatures)
	
	d_userfeaturematrix <- matrix(0,numberofusers,numberoffeatures)
	
	d_itemfeaturematrix <- cpp_grad_f(itembatch,d_itemfeaturematrix,M_item)
	
	d_userfeaturematrix <- cpp_grad_f(userbatch,d_userfeaturematrix,M_user)
	
	d_itembiasmatrix <- matrix(0,numberofitems,1)
	
	d_userbiasmatrix <- matrix(0,numberofusers,1)
	
	d_itembiasmatrix <- cpp_grad_f(itembatch,d_itembiasmatrix,M_itembias)
	
	d_userbiasmatrix <- cpp_grad_f(userbatch,d_userbiasmatrix,M_userbias)

	
	#update of item features
	itemfeaturematrix_inc <- momentum*itemfeaturematrix_inc + epsilon*(d_itemfeaturematrix/batchsize)
	itemfeaturematrix <- itemfeaturematrix - itemfeaturematrix_inc
	
	#update of user features
	userfeaturematrix_inc <- momentum*userfeaturematrix_inc + epsilon*(d_userfeaturematrix/batchsize)
	userfeaturematrix <- userfeaturematrix - userfeaturematrix_inc
	
	#update of item biases
	itembiasmatrix_inc <- momentum*itembiasmatrix_inc + epsilon*(d_itembiasmatrix/batchsize)
	itembiasmatrix <- itembiasmatrix - itembiasmatrix_inc
	
	#update of user biases
	userbiasmatrix_inc <- momentum*userbiasmatrix_inc + epsilon*(d_userbiasmatrix/batchsize)
	userbiasmatrix <- userbiasmatrix - userbiasmatrix_inc
	
}

	
	#computation of the predictions after parameter updates
    
	biasbatch <- predictionadjustment + userbiasmatrix[userbatch,] + itembiasmatrix[itembatch,]

	if (usebiasmodminmax) {
	biasbatch <- pmin(pmax(biasbatch,truncatethresholdmin),truncatethresholdmax)
	}
  
	featdotprodbatch <- rowSums(itemfeaturematrix[itembatch,]*userfeaturematrix[userbatch,])
  
	linnergvaluevec <- linnergvec(featdotprodbatch,biasbatch,vartheta)

	linnergvaluevec <- pmin(pmax(linnergvaluevec,linnergmin),linnergmax)

	thetavalvec <- gtrafovec(linnergvaluevec,vartheta)

	predictionbatch <- expvalfunctionvec(thetavalvec,vartheta)
	
	if (distribution == "poisson"){
	predictionbatch <- qpois(0.5, thetavalvec)
	}
	
	if (distribution == "gamma"){
	predictionbatch <- qgamma(0.5, scale=thetavalvec, shape=vartheta) 
	}
	
	if (distribution == "pareto"){
	predictionbatch <- qpareto(0.5, vartheta, thetavalvec)
	}

	predictionbatch <- pmin(pmax(predictionbatch,truncatethresholdmin),truncatethresholdmax)


#computation of the predictions on the validation set

	uservalidationbatch <- validationmatrix[,1]
	
	itemvalidationbatch <-  validationmatrix[,2]
	
	ratingvalidationbatch <- validationmatrix[,3]
	
	biasvalidationbatch <- predictionadjustment + userbiasmatrix[uservalidationbatch,] + itembiasmatrix[itemvalidationbatch,]

	featdotprodvalidationbatch <- rowSums(itemfeaturematrix[itemvalidationbatch,]*userfeaturematrix[uservalidationbatch,])

	linnergvaluevec <- linnergvec(featdotprodvalidationbatch,biasvalidationbatch,vartheta) 

	linnergvaluevec <- pmin(pmax(linnergvaluevec,linnergmin),linnergmax)

	thetavalvec <- gtrafovec(linnergvaluevec,vartheta)

	predictionvalidationbatch <- expvalfunctionvec(thetavalvec,vartheta)

	if (distribution == "poisson"){
	predictionvalidationbatch <- qpois(0.5, thetavalvec)
	}
	
	if (distribution == "gamma"){
	predictionvalidationbatch <- qgamma(0.5, scale=thetavalvec, shape=vartheta) 
	}
	
	if (distribution == "pareto"){
	predictionvalidationbatch <- qpareto(0.5, vartheta, thetavalvec)
	}
	
	predictionvalidationbatch <- pmin(pmax(predictionvalidationbatch,truncatethresholdmin),truncatethresholdmax)

	
	#mean absolute error on the validation set
	validationerror[epoch] <- sum(abs(predictionvalidationbatch - ratingvalidationbatch))/validationsamplesize

	if (is.na(validationerror[epoch]) | is.nan(validationerror[epoch]) | !(is.finite(validationerror[epoch]))){

	validationerror[epoch:maxepoch] <- NaN

	break

	}

	if (validationerror[epoch] < minvalidationerrorvalue){
	minvalidationerrorvalue <- validationerror[epoch]
	minvalidationerrorvalueparamlog <- paste("At epoch ",toString(epoch)," for parameter settings "," epsilon = ",toString(epsilon),", lambda = ",toString(lambda),", momentum = ",toString(momentum),", nfeatures = ",toString(numberoffeatures),
		" distparam = ",toString(vartheta)," maxepoch = ",toString(maxepoch)," ",toString(datasource)," truncmax = ",toString(truncatethresholdmax)," truncmin = ",toString(truncatethresholdmin)," excl = ",toString(excludesparseuseritemsthreshold))
	}

	
#error output
print(sprintf("Epoch %i:  Validation Mean Absolute Error: %6.4f", epoch, validationerror[epoch]))


}


print("Minimum Validation Mean Absolute Error achieved:")
print(round(minvalidationerrorvalue,4))
print(minvalidationerrorvalueparamlog)


}
