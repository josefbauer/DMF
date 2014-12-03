library(rSymPy)

sympy("var('xval')")
sympy("var('theta')")
sympy("var('vartheta')")
sympy("var('expval')")
sympy("var('featval')")
sympy("var('biasval')")


convertexpressionrtosympystring <- function (rexpression) {
  
sympyexpression <- gsub("\\^","\\*\\*",rexpression)

}

convertexpressionsympytorstring <- function (sympyexpression) {
  
rexpression <- gsub("\\*\\*","\\^",sympyexpression)
rexpression <- gsub("None","0",sympyexpression)

}



generatepartialderivative <- function (sympyexpressionstring,diffvariable,partialderivativeexpressionname,simplify=TRUE) {

sympyfsymbolicpartialderivative <- sympy(paste0("simplify(diff(",sympyexpressionstring,",",diffvariable,"))"))

if (simplify == FALSE) {sympyfsymbolicpartialderivative <- sympy(paste0("diff(",sympyexpressionstring,",",diffvariable,")"))}

symbolicpartialderivative <- convertexpressionsympytorstring(sympyfsymbolicpartialderivative)

eval(parse(text=paste0(partialderivativeexpressionname," <<- symbolicpartialderivative")), envir = environment())

}



generatevectorizedfunction <- function (functionexpressionrstring,variablestovectorizevec,vectorizedfunctionname) {

functionvectorizedbody <- functionexpressionrstring

for (variablenumber in 1:length(variablestovectorizevec)){
  
varvariable <- paste0("var",variablestovectorizevec[variablenumber])
renameprotectvarvariable <- "protvarnm"
functionvectorizedbody <- gsub(varvariable,renameprotectvarvariable,functionvectorizedbody)
functionvectorizedbody <- gsub(variablestovectorizevec[variablenumber],paste0("argvec[,",variablenumber,"]"),functionvectorizedbody)
functionvectorizedbody <- gsub(renameprotectvarvariable,varvariable,functionvectorizedbody)

}

functionvectorizedbody <- convertexpressionsympytorstring(functionvectorizedbody)

functionvectorized <- 
eval(parse(text=paste0("function (",paste0(variablestovectorizevec,"vec",collapse=","),",","vartheta) {}")), envir = environment())

parsedassignment <- parse(text=paste0("argvec <- cbind(",paste0(variablestovectorizevec,"vec",collapse=","),")"))
body(functionvectorized)[[2]] <- substitute(eval(parsedassignment))

parsedexpression <- parse(text=functionvectorizedbody)

body(functionvectorized)[[3]] <- substitute(eval(parsedexpression))

eval(parse(text=paste0(vectorizedfunctionname," <<- functionvectorized")), envir = environment())

}



generateexpressionsolution <- function (sympyexpressionstring,solutionvariable,solutionexpressionname) {

sympysolve <- sympy(paste0("solve([",sympyexpressionstring,"], [",solutionvariable,"])"))

sympysolve <- gsub("^.*?: ","",sympysolve)

sympysolve <- gsub("}.*?$","",sympysolve)

sympysolve <- convertexpressionsympytorstring(sympysolve)

eval(parse(text=paste0(solutionexpressionname," <<- sympysolve")), envir = environment())

}



generatelinnergexpression <- function (useritemfeatmodulationexpressionrstring,useritembiasmodulationexpressionrstring,targetvariable,outputexpressionname) {

useritemfeatmodulationexpressionrstring <- gsub(targetvariable,"featval",useritemfeatmodulationexpressionrstring)

useritembiasmodulationexpressionrstring <- gsub(targetvariable,"biasval",useritembiasmodulationexpressionrstring)

outputexpressionrstring <- paste0(useritemfeatmodulationexpressionrstring," + ",useritembiasmodulationexpressionrstring)

eval(parse(text=paste0(outputexpressionname," <<- outputexpressionrstring")), envir = environment())

}



generatealgebraiccomposition <- function (sympyexpressionstring1,sympyexpressionstring2,operation,outputexpressionname,simplify=TRUE) {

sympyfsymboliccomposition <- sympy(paste0("simplify(","(",sympyexpressionstring1,")",operation,"(",sympyexpressionstring2,")",")"))

if (simplify == FALSE) {sympyfsymboliccomposition <- sympy(paste0("(",sympyexpressionstring1,")",operation,"(",sympyexpressionstring2,")"))}

sympyfsymboliccomposition <- convertexpressionsympytorstring(sympyfsymboliccomposition)

eval(parse(text=paste0(outputexpressionname," <<- sympyfsymboliccomposition")), envir = environment())

}

