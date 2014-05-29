# <OWNER> = Saip CISS
# <ORGANIZATION> = QueensBridge Quantitative
# <YEAR> = 2014

# LICENSE 
# BSD 3-CLAUSE LICENSE

# Copyright (c) 2014, Saip CISS
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# END OF LICENSE 

randomUniformForest <- function(...) UseMethod("randomUniformForest")
	
importance <- function(object, ...) UseMethod("importance")

randomUniformForest.formula <- function(formula, data = NULL, ...)
{
	if (is.null(data))
	stop ("please provide data")
	data <- fillVariablesNames(data)	
	mf <- model.frame(formula = formula, data = as.data.frame(data))
	x <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
	y <- model.response(mf)
	RUFObject <- randomUniformForest.default(x, Y = y, ...)
	RUFObject$call <- match.call()
	RUFObject$formula <- formula
	class(RUFObject) <- c("randomUniformForest.formula", "randomUniformForest")
	RUFObject
}

print.randomUniformForest <- function(x,...)
{
	object <- x
	cat("Call:\n")
	print(object$call)
	cat("\n")
	cat("Type of random uniform forest: ") 
	if (object$forest$regression) 	{  cat("Regression\n")	}
	else { 	cat("Classification\n") }
	cat("\n")
	print(object$forestParams)
	cat("\n")
	if (!is.null(object$forest$OOB))
	{ 
		cat("Out-of-bag (OOB) evaluation")
		if (!object$forest$regression)
		{
			OOBErrorRate = mean(100*object$forest$pred.error)
			cat("\nOOB estimate of error rate: ", round(OOBErrorRate, 2),"%\n", sep ="")
			cat("OOB error rate bound (with 1% deviation): ",round(OOBErrorRate + OOBErrorRate*(1 - estimatePredictionAccuracy(floor(length(object$forest$OOB.predicts)*0.368))), 2),"%\n", sep="")
			cat("\nOOB confusion matrix:\n")
			colnames(object$forest$OOB)[1:length(object$classes)] = object$classes
			rownames(object$forest$OOB)[1:length(object$classes)] = object$classes
			print(round(object$forest$OOB,4))			
			if ((length(object$classes) == 2) & (rownames(object$forestParams)[1] != "reduceDimension"))
			{	
				cat("\nOOB estimate of AUC: ", round(pROC::auc(as.numeric(object$y), as.numeric(object$forest$OOB.predicts))[[1]], 4), 
				sep = "") 
			}
			cat("\nGeometric mean:", round(gMean(object$forest$OOB),4),"\n")
			if (nrow(object$forest$OOB) > 2)
			{ cat("Geometric mean of the precision:", round(gMean(object$forest$OOB, precision = TRUE),4),"\n") }
			if (!is.null(object$forest$OOB.strengthCorr))
			{
				cat("\nTheorical (Breiman) bounds")
				cat("\nPrediction error (expected to be lower than): ", round(mean(100*rmNA(object$forest$OOB.strengthCorr$PE)),2),"%\n", sep ="")
				cat("Upper bound of prediction error: ", round(mean(100*rmNA(object$forest$OOB.strengthCorr$std.strength^2/object$forest$OOB.strengthCorr$strength^2)),2),"%\n", sep ="")
				cat("Average correlation between trees:", round(mean(round(rmNA(object$forest$OOB.strengthCorr$avg.corr),4)),4),"\n")
				cat("Strength (margin):", round(mean(round(rmNA(object$forest$OOB.strengthCorr$strength),4)),4),"\n")
				cat("Standard deviation of strength:", round(mean(round(rmNA(object$forest$OOB.strengthCorr$std.strength),4)),4),"\n")
			}
		}
		else
		{
			cat("\nMean of squared residuals:", round(mean(object$forest$pred.error),8), "\n")
			cat("Mean squared error bound (experimental):", round(object$forest$pred.error + object$forest$pred.error*(1- estimatePredictionAccuracy(floor(length(object$forest$OOB.predicts)*(1 - as.numeric(as.vector(object$forestParams["subsamplerate",1])))))), 6),"\n")
			if (length(object$forest$percent.varExplained) > 1)
			{	
				varExplained = object$forest$percent.varExplained
				for (i in 1:length(varExplained))
				{	varExplained[i] = paste(object$forest$percent.varExplained[i], "%", sep="")	}
				cat("Variance explained:", varExplained, "\n")	
			}
			else { cat("Variance explained: ", object$forest$percent.varExplained, "%\n", sep = "")	} #percent.varExplained
			cat("\nOOB residuals:\n")
			Residuals = summary(rmNA(object$forest$OOB.predicts - object$y))
			names(Residuals) = c("Min", "1Q", "Median", "Mean", "3Q", "Max")
			print(Residuals)
			cat("Mean of absolute residuals:", sum(abs(object$forest$OOB.predicts - object$y))/length(object$y),"\n")						
			if (!is.null(object$forest$OOB.strengthCorr))
			{
				cat("\nTheorical (Breiman) bounds")
				cat("\nTheorical prediction error:", round(rmNA(object$forest$OOB.strengthCorr$PE.forest),6) ,"\n")
				cat("Upper bound of prediction error:", round(rmNA(object$forest$OOB.strengthCorr$PE.max),6),"\n")
				cat("Mean prediction error of a tree:", round(rmNA(object$forest$OOB.strengthCorr$PE.tree),6),"\n")
				cat("Average correlation between trees residuals:", round(rmNA(object$forest$OOB.strengthCorr$mean.corr),4),"\n")
				residuals_hat = vector(length = ncol(object$forest$OOB.votes))
				residuals_hat <- apply(object$forest$OOB.votes, 2, function(Z) mean(rmInf(Z) - mean(rmNA(object$forest$OOB.predicts))))		
				cat("Expected squared bias (experimental):", round(mean(rmNA(residuals_hat))^2,6),"\n")
			}
		}		
	}  
    if (!is.null(object$errorObject))
	{	
	    cat("\nTest set")
		if (!object$forest$regression)
		{
			cat("\nError rate: ")
			cat(round(100*object$errorObject$error,2), "%\n", sep="")
			cat("\nConfusion matrix:\n")
			print(round(object$errorObject$confusion,4))			
			if (!is.null(object$errorObject$AUC)) 
			{ 
				cat("\nArea Under ROC Curve:", round(object$errorObject$AUC,4)) 
				cat("\nF1 score:", round(fScore(object$errorObject$confusion),4)) 	
			}
			cat("\nGeometric mean:", round(gMean(object$errorObject$confusion),4),"\n")
			if (nrow(object$errorObject$confusion) > 2)
			{ cat("Geometric mean of the precision:", round(gMean(object$errorObject$confusion, precision = TRUE),4),"\n") }
		}
		else
		{	
			cat("\nMean of squared residuals: ", round(object$errorObject$error, 6), "\n", sep="")			
			cat("Variance explained: ", object$errorObject$percent.varExplained, "%\n\n", sep = "")		
			Residuals <- object$errorObject$Residuals
			names(Residuals) = c("Min", "1Q", "Median", "Mean", "3Q", "Max")
			cat("Residuals:\n")
			print(Residuals)
			cat("Mean of absolute residuals: ", round(object$errorObject$meanAbsResiduals, 6), "\n", sep="")							
		} 		
	}
}

summary.randomUniformForest <- function(object, maxVar = 30, border = NA,...)
{	
	object <- filter.object(object)	
	if (!is.null(object$forest$variableImportance))
	{	
		par(las=1) 
		maxChar = floor(2 + max(nchar(object$variablesNames))/2)		
		par(mar=c(5, maxChar + 1,4,2))		
		varImportance1 = varImportance = object$forest$variableImportance		
		if (!object$forest$regression)
		{
			varImportance[,"class"] = object$classes[as.numeric(varImportance[,"class"])]
			varImportance1[,"class"] = varImportance[,"class"]
		}		
		nVar = nrow(varImportance)
		if (nVar > maxVar) { varImportance = varImportance[1:maxVar,] }	
		barplot( varImportance[nrow(varImportance):1,"percent.importance"], horiz = TRUE, col = sort(heat.colors(nrow(varImportance)), decreasing = TRUE), 
		names.arg = varImportance[nrow(varImportance):1,"variables"], border = border,
		xlab = "Relative Influence (%)", main = "Variable Importance based on information gain")
		abline(v = 100/nVar, col ='grey')
		cat("Variables summary: ", "\n")
		print(varImportance1)
		cat("\n")
		cat("Average tree size (number of nodes) summary: ", "\n")
		nodesView = unlist(lapply(object$forest$object,nrow))
		print(floor(summary(nodesView)))
		cat("\n")
		cat("Average Leaf nodes (number of terminal nodes) summary: ", "\n")
		terminalNodesView = unlist(lapply(object$forest$object, function(Z) length(which(Z[,"status"] == -1))))
		print(floor(summary(terminalNodesView)))
		cat("\n")
		cat("Leaf nodes size (number of observations per leaf node) summary: ", "\n")
		print(summary(unlist(lapply(object$forest$object, function(Z) Z[which(Z[,"prediction"] != 0), "nodes"]) )))
		cat("\n")
		cat("Average tree depth :", round(log(mean(nodesView))/log(2),0), "\n")
		cat("\n")
		cat("Theorical tree depth :", round(log(length(object$y), base = 2),0), "\n")
		cat("\n")		
	}
	else
	{ 
		cat("Average tree size (number of nodes) summary: ", "\n")
		nodesView = unlist(lapply(object$forest$object,nrow))
		print(floor(summary(nodesView)))
		cat("\n")
		cat("Average Leaf nodes (number of terminal nodes) summary: ", "\n")
		terminalNodesView = unlist(lapply(object$forest$object, function(Z) length(which(Z[,"status"] == -1))))
		print(floor(summary(terminalNodesView)))
		cat("\n")
		cat("Leaf nodes size (number of observations per leaf node) summary: ", "\n")
		print(summary(unlist(lapply(object$forest$object, function(Z) Z[which(Z[,"prediction"] != 0), "nodes"]) )))
		cat("\n")
		cat("Average tree depth :", round(log(mean(nodesView))/log(2),0), "\n")
		cat("\n")
		cat("Theorical tree depth :", round(log(length(object$y), base = 2),0), "\n")
		cat("\n")		
	}
}

plot.randomUniformForest <- function(x, ...) 
{
	object <- x
	if (!is.null(object$forest$OOB.votes))
	{
		Ytrain = object$y
		OOBMonitoring = object$forest$OOB.votes
		
		ZZ <- monitorOOBError(OOBMonitoring, Ytrain, regression = object$forest$regression, ...)
		
		if (object$forest$regression)
		{
			plot(ZZ, type = 'l', lty=2, xlab = "Trees", ylab ="OOB mean squared error") 
		}
		else
		{	    
			plot(apply(ZZ[,1:3],1, min), type='l', lty=2, col = "green", xlab = "Trees", ylab ="OOB error")
			points(apply(ZZ[,1:3],1, mean), type='l', lty=3)
			points(apply(ZZ[,1:3],1, max), type='l', lty=3, col='red')
		}
		grid()
	}
	else
	{ print("no OOB data to plot")	}
}

getTree.randomUniformForest <- function(object, whichTree, labelVar = TRUE)
{
	if (labelVar)
	{
		Tree = data.frame(object$forest$object[[whichTree]])
		idx = which(Tree[, "split.var"] != 0)
		Tree[idx, "split.var"] = object$variablesNames[Tree[idx, "split.var"]]
		
		return(Tree)
	}
	else
	{  return(object$forest$object[[whichTree]])  }
}
		
predict.randomUniformForest <- function(object, X, 
	type = c("response", "prob", "votes", "confInt", "ranking", "quantile", "truemajority", "all"),
	classcutoff = c(0,0), 
	conf = 0.95,
	whichQuantile = NULL,
	rankingIDs = NULL,
	threads = "auto", 
	parallelpackage = "doParallel", ...) rUniformForestPredict(object, X, type = type, classcutoff = classcutoff, 
									conf = conf,
									whichQuantile = whichQuantile,
									rankingIDs = rankingIDs,
									threads = threads, 
									parallelpackage = parallelpackage, ...)

residualsRandomUniformForest <- function(object, Y = NULL) 
{
	object = filter.object(object)
	if (is.null(Y))
	{
		if (is.null(object$forest$OOB.predicts))
		{	stop("please enable OOB option when computing a random uniform forest") }
		else
		{	
			print("OOB residuals:")
			cat("\n")
			return(object$y - object$forest$OOB.predicts) 
		}
	}
	else
	{	
		if (is.numeric(object))
		{ return(object - Y) }
		else
		{ 
			if (!is.null(object$predictionObject$majority.vote))
			{ return(object$predictionObject$majority.vote - Y) }
			else
			{	stop("Please provide model responses to compute residuals")	}
		}
	}	
}

genericOutput <- function(xtest, ytest, paramsObject, RUF.model, ytrain = NULL, classcutoff = c(0,0))
{
	classes = NULL
	if (!is.null(ytest))
	{
		if (is.factor(ytest)) { YNames = classes = levels(ytest)   }
	}
	else
	{
		if (!is.null(ytrain))
		{
			if (!RUF.model$regression) 
			{ 
				ytrain = as.factor(ytrain) 
				YNames = classes = levels(ytrain)  
			}
		}
		else
		{
			if (!RUF.model$regression)  { YNames = classes = sort(unique(RUF.model$object[[2]][,6])[-1]) }
		}
	}		
	if (!is.null(RUF.model$OOB))
	{
		if (!RUF.model$regression & is.factor(ytest)) { row.names(RUF.model$OOB) = colnames(RUF.model$OOB)[-(length(YNames)+1)] = YNames }
	}
	if (!is.null(xtest)) 
	{ 
		classwtString = as.vector(paramsObject[which(row.names(paramsObject) == "classwt"),1])	
		classwt = FALSE
		if (is.na(as.logical(classwtString))) {  classwt = TRUE }
		
		if (!RUF.model$regression & (as.numeric(classcutoff[2]) != 0)) 
		{  
			classcutoff = c(which(classes == as.character(classcutoff[1])), as.numeric( classcutoff[2]))
		}				
		if (classcutoff[2] != 0 ) {  classcutoff[2] = 0.5/classcutoff[2] }		
		RUF.predicts <- randomUniformForestCore.predict(RUF.model, xtest, pr.classwt = classwt, pr.imbalance = classcutoff)		 
		if (!is.null(ytest)) 
		{ 
			errorObject <- someErrorType(RUF.predicts, factor2vector(ytest)$vector, generic = FALSE, regression = RUF.model$regression)
			if (!RUF.model$regression & is.factor(ytest)) 
			{ row.names(errorObject$confusion) = colnames(errorObject$confusion)[-(length(YNames)+1)] = YNames }			
			RUFObject = list(forest = RUF.model, predictionObject = RUF.predicts, errorObject = errorObject, forestParams = paramsObject, classes = classes) 
		}
		else
		{	RUFObject = list(forest = RUF.model, predictionObject = RUF.predicts, forestParams = paramsObject,  classes = classes) }		
	}
	else
	{   RUFObject = list(forest = RUF.model, forestParams = paramsObject, classes = classes)	}	
	RUFObject
}

randomUniformForest.default <- function(X, Y = NULL, xtest = NULL, ytest = NULL, ntree = 100, 
	mtry = ifelse(bagging,ncol(X),floor(4/3*ncol(X))), 
	nodesize = 1,
	maxnodes = Inf,
    depth = Inf,
    depthcontrol = NULL,	    
	regression = ifelse(is.factor(Y), FALSE, TRUE),
	replace = ifelse(regression,FALSE,TRUE),
	OOB = TRUE,
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
	subsamplerate = ifelse(regression,0.7,1),
    importance = TRUE,	
	bagging = FALSE,
	unsupervised = FALSE, 
	proximities = FALSE,
	classwt = NULL,
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
    rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L2", "L1"),	
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	threads = "auto",
	parallelpackage = "doParallel", 
	...)
{
	{
		if (threads != 1)
		{
			if (sample(3,1) == 3) { rm.tempdir() }
		}
		if (is.null(Y)) { stop("Unsupervised learning is not currently supported. Please provide a response vector.\n") }
		if (ntree < 2) { stop("Please use at least 2 trees for computing forest.\n") }			
		if ( (subsamplerate == 1) &  (replace == FALSE))  { OOB = FALSE }
		if (depth < 3) { stop("Stumps are not allowed. Minimal depth is 3, leading, at most, to 8 leaf nodes.\n") }
		if ( (BreimanBounds & (length(unique(Y)) > 2)) | (BreimanBounds & ntree > 500) ) 
		{ cat("Warning : Breiman Bounds, especially for multiclass problems, are computationnaly intensive.\n") }
		if (maxnodes < 6) { stop("Maximal number of nodes must be above 5.\n") }
		X <- fillVariablesNames(X)        
		getFactors <- which.is.factor(X, count = TRUE)	 
		if (is.data.frame(X)) 
		{ 
			cat("X is a dataframe. String or factors have been converted to numeric values.\n") 
			X <- NAfactor2matrix(X)	
		}				
		if (!is.null(categoricalvariablesidx))
		{
			if (categoricalvariablesidx[1] == "all")
			{	
				factorVariables <- which(getFactors > 0)
				if (length(factorVariables) > 0) { 	categoricalvariablesidx = factorVariables	}
			}
			else
			{
				maxLengthOfCat <- max(apply(X[,categoricalvariablesidx, drop =FALSE], 2, function(Z) length(unique(Z))))
				if (maxLengthOfCat > 10) { getFactors <- which.is.factor(X, maxClasses = maxLengthOfCat, count = TRUE) 	}			
			}
		}
		if ( length(which( (X == Inf) | (X == -Inf) ) ) > 0) 
		{ stop("Inf or -Inf values found in data. Learning can not be done.\nRemove or replace them with NA in order to learn.\n") }
		NAInputs = which(is.na(X), arr.ind = TRUE)[,1]
		matchNA = (length(NAInputs) > 0)		
		if ( (length(NAInputs) > (nrow(X) - 30)) & (na.action[1] == "omit") ) { stop("Too much missing values in data.\n") }		
		if (!regression & !is.factor(Y) ) { Y = as.factor(Y) }
		NALabels = which(is.na(Y))
		matchNALabels = (length(NALabels) > 0)		
		if ( (length(NALabels) > (length(Y) - 30)) & (na.action[1] == "omit") )	{ stop("Too much missing values in responses.\n") }			
		if (matchNA | matchNALabels)
		{
			if (!is.null(Y))
			{
				newY <- as.vector(NAfactor2matrix(matrix(as.numeric(Y))))
				if (na.action[1] != "omit") 
				{ 
					cat("NA found in data. Imputation (fast or accurate, depending on option) is used for missing values.\nIt is strongly recommended to impute values outside modelling, using one of many available models.\n")					
					XY <- na.impute(cbind(X, newY), type = na.action[1])
					newY <- XY[,ncol(X)+1]
					X <- XY[,-(ncol(X)+1)]
					rm(XY)					
					if (is.factor(Y) & matchNALabels)
					{	
						levelsY = levels(Y)
						Y[NALabels] = levelsY[newY[NALabels]]
					}
				}
				else
				{
					if (matchNALabels & matchNA) {	rmRows = unique(c(NAInputs, NALabels))	}
					else
					{
						if (matchNA) {	rmRows = unique(NAInputs)	}						
						if (matchNALabels) 	{	rmRows = unique(NALabels)	}
					}					 
					X = X[-rmRows,]
					Y = Y[-rmRows]
				}
			}
			else
			{
				cat("If accuracy is needed, it is strongly recommended to impute values outside modelling, using one of many available models.\n")
				if (na.action[1] != "omit") { 	X <- na.impute(X, type = na.action[1])	}
				else
				{		
					if (matchNA) { X <- X[-NAInputs,] 	}		
				}
			}
		}		
		if (randomcombination[1] > 0)
		{ 
			randomcombinationString = randomcombination[1]
			L.combination = length(randomcombination)
			if (L.combination%%3 != 0)
			{	weights = round(replicate(length(randomcombination)/2, sample(c(-1,1),1)*runif(1)), 2)	}
			else
			{ 	
				weights = round(randomcombination[(L.combination + 1 - L.combination/3):L.combination],2)
				randomcombination = randomcombination[1:(L.combination - (L.combination/3))]
			}			
			for (i in 2:length(randomcombination)) {  randomcombinationString = paste(randomcombinationString, randomcombination[i], sep=",") }
			for (i in 1:length(weights)) {  randomcombinationString = paste(randomcombinationString, weights[i], sep=",") }
			if (!is.null(xtest)) { xtest <- randomCombination(NAfactor2matrix(xtest), combination = randomcombination, weights = weights) }			
			randomcombinationObject = list(randomcombination, weights)					
			X <- randomCombination(X, combination = randomcombinationObject[[1]], weights = randomcombinationObject[[2]])			
		}
		else
		{   randomcombinationObject = randomcombination }
	}	
	RUF.model <- randomUniformForestCore(X, trainLabels = Y, ntree = ntree, nodeMinSize = nodesize, maxNodes = maxnodes, 
		features = mtry, rf.bootstrap = replace, depth = depth, depthControl = depthcontrol, rf.subagging = subsamplerate, classwt = classwt, classCutOff = classcutoff, rf.overSampling = oversampling, rf.targetClass = targetclass, rf.rebalancedSampling = rebalancedsampling, 
		rf.outputPerturbationSampling = outputperturbationsampling,	rf.randomCombination = randomcombinationObject, 
		rf.randomFeature = randomfeature, rf.x.bagging = bagging, rf.featureSelectionRule = featureselectionrule[1], 
		rf.regression = regression, use.OOB = OOB, BreimanBounds = BreimanBounds, variableImportance = importance, 
		whichCatVariables = categoricalvariablesidx, logX = logX, threads = threads, factors = getFactors, 
		parallelPackage = parallelpackage[1])
	if (!is.null(classwt))
	{ 
		classwtString = classwt[1]
		for (i in 2:length(classwt)) { 	classwtString = paste(classwtString, classwt[i], sep=",") }	
	}	
	if (length(targetclass) > 1)
	{ 
		targetclassString = targetclass[1]
		for (i in 2:length(targetclass)) { 	targetclassString = paste(targetclassString, targetclass[i], sep=",") }	
	}	
	if (length(rebalancedsampling) > 1)
	{ 
		rebalancedsamplingString = rebalancedsampling[1]
		for (i in 2:length(rebalancedsampling)) 
		{ rebalancedsamplingString = paste(rebalancedsamplingString, rebalancedsampling[i], sep= ",") }	
	}	
	if (RUF.model$regression) {	classcutoff = c(0,0)  }	
	if (as.numeric(classcutoff[2]) == 0) {	classcutoffString = FALSE }
	else
	{
		classcutoffString = levels(Y)[which(levels(Y) == as.character(classcutoff[1]))]
		classcutoffString = paste("Class ", classcutoffString, "," , as.numeric(classcutoff[2])*100, "%", sep ="")
	}		
	paramsObject = c(ntree, mtry, nodesize, maxnodes, as.character(replace), bagging, depth, 
		ifelse(is.null(depthcontrol), length(depthcontrol)> 0, depthcontrol), OOB, importance, subsamplerate, 
		ifelse(is.null(classwt), length(classwt) > 0, classwtString), classcutoffString, 
		ifelse((oversampling == 0), (oversampling != 0), oversampling), outputperturbationsampling, 
		ifelse(length(targetclass) > 1, targetclassString, targetclass), 
		ifelse(length(rebalancedsampling) > 1, rebalancedsamplingString, rebalancedsampling), 
		ifelse((randomcombination[1] == 0),(randomcombination != 0), randomcombinationString), randomfeature, 
		ifelse(is.null(categoricalvariablesidx), length(categoricalvariablesidx) > 0, length(categoricalvariablesidx)), featureselectionrule[1])				
	if (RUF.model$regression) 
	{ 
		paramsObject[21] = if (featureselectionrule[1] == "L1") { "Sum of absolute residuals" }  else { "Sum of squared residuals" }
	    if ((paramsObject[5] == "FALSE") & (subsamplerate == 1)) { paramsObject[8] = FALSE }
	}	
	names(paramsObject) = c("ntree", "mtry", "nodesize", "maxnodes", "replace", "bagging", "depth", "depthcontrol", "OOB", "importance", "subsamplerate", "classwt", "classcutoff", "oversampling", "outputperturbationsampling", "targetclass", "rebalancedsampling", "randomcombination", "randomfeature", "categorical variables", "featureselectionrule")						
	paramsObject = as.data.frame(paramsObject)	
	RUF.model$logX = logX	
	if (!regression & !is.null(ytest) ) 
	{ 
		if (!is.factor(ytest))
		{ ytest = as.factor(ytest) }
	}
	RUFObject <- genericOutput(xtest, ytest, paramsObject, RUF.model, ytrain = Y, classcutoff = classcutoff)
	RUFObject$logX = logX
	if (!is.null(Y)) { 	RUFObject$y = Y	 }	
	if (is.null(colnames(X))) 
	{  
		varNames = NULL
		for (i in 1:ncol(X)) { varNames = c(varNames, paste("V", i, sep="")) }
		RUFObject$variablesNames = varNames
	}
	else
	{	RUFObject$variablesNames = colnames(X)	}	
	if (randomcombination[1] > 0)
	{ 	
		tempVarNames = vector(length = length(randomcombination)/2)
		idx = 1
		for (i in 1:(length(randomcombination)/2)) 
		{ 
			tempVarNames[i] = paste("V", randomcombination[idx], "x", randomcombination[idx+1], sep="") 
			idx = idx + 2
		}
		RUFObject$variablesNames = c(RUFObject$variablesNames, tempVarNames)
	}
	if (!is.null(categoricalvariablesidx)) {  RUFObject$categoricalvariables = categoricalvariablesidx }
	RUFObject$call <- match.call()
	class(RUFObject) <- "randomUniformForest"
	RUFObject
}

randomUniformForestCore <- function(trainData, trainLabels = 0, 
	features = floor(4/3*(ncol(trainData))), 
	ntree = 100, 
	nodeMinSize = 1,
	maxNodes = Inf,
	depth = Inf,
    depthControl = NULL,
	splitAt = "random", 
	rf.regression = TRUE, 
	rf.bootstrap = ifelse(rf.regression, FALSE, TRUE), 
	use.OOB = TRUE,
	BreimanBounds = ifelse(use.OOB, TRUE, FALSE),
	rf.subagging = ifelse(rf.regression, 0.7, 1),
	rf.x.bagging = FALSE, 
	classwt = NULL,
	classCutOff = c(0,0),
	rf.overSampling = 0,
	rf.targetClass = -1,
	rf.outputPerturbationSampling = FALSE,
    rf.rebalancedSampling = FALSE,
	rf.featureSelectionRule = c("entropy", "gini", "random", "L2", "L1"),
	variableImportance = TRUE,	
	unsupervised = FALSE, 
	proximities = FALSE,
	rf.randomCombination = 0,
	rf.randomFeature = FALSE,
	whichCatVariables = NULL,
	logX = FALSE,
	factors = NULL,
	threads = "auto",	
	parallelPackage = "doParallel",
	export = c("uniformDecisionTree", "CheckSameValuesInAllAttributes", "CheckSameValuesInLabels", "fullNode", "genericNode", "leafNode", "randomUniformForestCore.predict", "onlineClassify", "overSampling", "predictDecisionTree", "options.filter", "majorityClass", "randomCombination", "randomWhichMax", "which.is.na", "factor2vector", "outputPerturbationSampling", "rmNA", "count.factor", 
	"find.idx", "classifyMatrixCPP", "L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP",
	"entropyInformationGainCPP", "runifMatrixCPP"))
{
	set.seed(sample(ntree,1))	
	if ((rf.randomFeature) & (features == "random") ) { stop("random feature is a special case of mtry = 'random'") }
	if (is.numeric(threads))
	{
		if (threads < 1) { stop("Number of threads must be positive") }
	}
	if (!is.matrix(trainData)) { trainData <- NAfactor2matrix(trainData) }
	if (!is.null(whichCatVariables)) 
	{ 
	  cat("Please also consider Dummies for categorical variables. Seem yet less reliable.\nFormula automatically builds dummies.\n")	
	}
	{
		if (!is.null(classwt)) 
		{ 
			classwt = classwt/sum(classwt) 
			if (is.na(sum(classwt))) { stop ("NA found in class weights.\n") }			
			if (sum(classwt) > 1.1)  { stop ("(inverse) weights do not sum to 1.\n") }
		}		
		if (length(trainLabels) != 1)
		{	
			if (nrow(trainData) != length(trainLabels)) 
			{ stop ("X and Y don't have the same size.\n")	}	
		}		
		if (!is.matrix(trainData))
		{ stop ("X cannot be converted to a matrix. Please provide true matrix not dataframe.") }		
		if (rf.subagging == 0) { rf.subagging = 1 }		
		if (rf.x.bagging) { features = min(features, ncol(trainData))  }
		if ( (rf.subagging < 0.5) & (nrow(trainData) < 300) ) { rf.subagging = 0.5;  cat("Too small output matrix. Subsample rate has been set to 0.5.\n") }		
		if (!is.character(nodeMinSize))
		{
			if ( (nodeMinSize != 1) & (nodeMinSize >  floor(nrow(trainData)/4)) ) { stop("nodeMinSize is too high. Not suitable for a random forest.\n") }	
			if ( (nodeMinSize < 1) ) {  nodeMinSize  = 1; cat("Minimal node size has been set to 1.\n") }
		}		
		if ( features == 1 ) { rf.randomFeature = TRUE }
		if ( features < 1 )  
		{ 
			features = floor(4/3*(ncol(trainData)))  
			bagging = FALSE
			cat("Error setting mtry. Resetting to default values.\n") 
		}		
		if (rf.randomFeature) { variableImportance = FALSE }		
		if (is.factor(trainLabels))
		{ 	
			if ((as.numeric(classCutOff[2]) != 0)) 
			{  	
				classCutOff = c(which(levels(trainLabels) == as.character(classCutOff[1])), 0.5/as.numeric(classCutOff[2])) 
				if (length(classCutOff) == 1)
				{ stop("Label not found. Please provide name of the label instead of its position.\n") }
			}					
			rf.regression = FALSE 
			labelsObject = factor2vector(trainLabels)
			trainLabels = as.numeric(as.vector(labelsObject$vector))
			if (length(unique(trainLabels)) == 1 ) { stop ("Y is a constant value. Thus, learning is not needed.\n")	}			
		}
		else
		{	
			if (rf.regression & (length(unique(trainLabels)) > 32)) 
			{ 
				cat("Regression has been performed.\n")				
				if ((rf.subagging < 1) & (rf.bootstrap == TRUE))
				{ cat("For only accuracy, use option 'subsamplerate = 1' and 'replace = FALSE' \n") }				
				classCutOff =  c(0,0) 				
			}
			else
			{
				if (!rf.regression)
				{
					if ((as.numeric(classCutOff[2]) != 0)) 
					{  	
						classCutOff = c(which(levels(trainLabels) == as.character(classCutOff[1])), 0.5/as.numeric(classCutOff[2])) 
						if (length(classCutOff) == 1)
						{ stop("Label not found. Please provide name of the label instead of its position") }
					}					
					labelsObject = factor2vector(trainLabels)
					trainLabels = as.numeric(as.vector(labelsObject$vector))
					labelsObject$factors = levels(as.factor(trainLabels))
					if (length(unique(trainLabels)) == 1 ) { stop ("Y is a constant value. Thus, learning is not needed\n")	}
				}
				else
				{		
					if (length(unique(trainLabels)) <= 32)  
					{ 
						cat("Regression has been performed but there is less than 32 distinct values.\nRegression option could be set to FALSE, if classification is required.\nIf outputs are factors, conversion is automatically done.\n")
						classCutOff =  c(0,0)
					}
				}
			}
		}							
		if (!rf.regression)
		{
			rf.classes <- as.numeric(levels(as.factor(trainLabels)) )			
			if (as.character(rf.classes[1]) != labelsObject$factors[1]) 
			{	
				cat("Labels", labelsObject$factors, "have been converted to", rf.classes, 
				"for ease of computation and will be used internally as a replacement.\n")	
			}			
			if (!is.null(classwt))
			{ 	
				if (length(classwt) != length(rf.classes)) { stop("Length of class weigths is not equal to length of classes.\n") }
			}
		}
		else
		{ 	rf.classes = trainLabels; rf.overSampling = 0; rf.rebalancedSampling = FALSE; classCutOff =  c(0,0)  }		
		if (logX)
		{ 
		   if (is.null(whichCatVariables))  {  trainData <- generic.log(trainData) }
		   else   {  trainData[,-whichCatVariables] <- generic.log(trainData[,-whichCatVariables]) }
		}		
		if (!rf.regression)  
		{ 
			if ( (rf.featureSelectionRule[1] == "L1") | (rf.featureSelectionRule[1] == "L2") )
			{	
				rf.featureSelectionRule[1] = "entropy" 
				cat("Feature selection rule has been set to entropy.\n")
			}
		}
		else 
		{  
			if ( (rf.featureSelectionRule[1] == "entropy")  |  (rf.featureSelectionRule[1] == "gini") ) 
			{	rf.featureSelectionRule[1] = "L2" }			
			if (!is.null(depthControl) & (!is.character(depthControl)))
			{	
				if (depthControl < 1)
				{
					depthControl = 1 
					cat("Depth control option is lower than its range values. Resetting option.\n")
				}				
			}			
		}
	}
	{
		#require(parallel)	
		max_threads = detectCores()		
		if (threads == "auto")
		{	
			if (max_threads == 2) { threads = max_threads }
			else {	threads  = max(1, max_threads - 1)  }
		}
		else
		{
			if (max_threads < threads) 
			{	cat("Warning : number of threads indicated by user was higher than logical threads in this computer.\n") }
		}
		
		{
			#require(doParallel)			
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
		}
		chunkSize  <-  ceiling(ntree/getDoParWorkers())
		smpopts  <- list(chunkSize = chunkSize)
	}
	if (!rf.bootstrap & (rf.subagging == 1)) {  use.OOB = FALSE }
	if (use.OOB)
	{
		# rufObject <- foreach(icount(ntree), .export = export, .options.smp = smpopts, .inorder = FALSE, .multicombine = TRUE, 
		# .packages = "rUniformForestCppClass") %dopar%
		rufObject <- foreach(i = 1:ntree, .export = export, .options.smp = smpopts, .inorder = FALSE, .multicombine = TRUE) %dopar%
		{
			uniformDecisionTree(trainData, trainLabels, nodeMinSize = nodeMinSize, maxNodes = maxNodes, rf.features = features, 
				getSplitAt = splitAt, regression = rf.regression, bootstrap = rf.bootstrap, subagging = rf.subagging, treeDepth = depth, treeClasswt = classwt,	treeOverSampling = rf.overSampling, targetClass = rf.targetClass, OOB = use.OOB, 
				treeRebalancedSampling = rf.rebalancedSampling, x.bagging = rf.x.bagging, random.combination = rf.randomCombination, randomFeature = rf.randomFeature, treeCatVariables = whichCatVariables, outputPerturbation = rf.outputPerturbationSampling, featureSelectionRule = rf.featureSelectionRule, treeDepthControl = depthControl)
		}		
		stopCluster(Cl)
		{
			OOB.matrix = NULL 
			new.rufObject = vector("list", ntree)
			new.rufObject$Tree = lapply(rufObject, function(Z) Z$Tree) 
			for (i in 1:ntree) 	{	OOB.matrix <- rbind(OOB.matrix, cbind(rufObject[[i]]$OOB.idx, rep(i,length(rufObject[[i]]$OOB.idx))))  } 
			OOB.val = sort(unique(OOB.matrix[,1])) 
			n.OOB = nrow(trainData) 
			OOB.votes2 = matrix(Inf, n.OOB, ntree)			
			if (is.null(classwt))
			{
				OOB.votes <- randomUniformForestCore.predict(new.rufObject$Tree, trainData, pr.regression = rf.regression, 
				classes = rf.classes, OOB.idx = TRUE, pr.parallelPackage = parallelPackage[1], pr.imbalance = classCutOff, 
				pr.threads = threads )				 
				for (j in 1:ntree)
				{
					idxJ = which(OOB.matrix[,2] == j)
					if (length(idxJ) > 0) 	{	OOB.votes2[OOB.matrix[idxJ,1],j] = OOB.votes[OOB.matrix[idxJ,1],j] 	}
				}
				OOB.object <- majorityClass(OOB.votes2, rf.classes, m.imbalance = classCutOff, m.regression = rf.regression)
			}
			else
			{
				OOB.allWeightedVotes = matrix(Inf, n.OOB, ntree)
				OOB.votes <- randomUniformForestCore.predict(new.rufObject$Tree, trainData, pr.regression = rf.regression, 
				classes = rf.classes, pr.classwt = TRUE, OOB.idx = TRUE, pr.parallelPackage = parallelPackage[1], pr.imbalance = classCutOff, pr.threads = threads)						
				for (j in 1:ntree)
				{
					idxJ = which(OOB.matrix[,2] == j)
					if (length(idxJ) > 0) 	
					{	
						OOB.votes2[OOB.matrix[idxJ,1],j] = OOB.votes$all.votes[OOB.matrix[idxJ,1],j] 	
						OOB.allWeightedVotes[OOB.matrix[idxJ,1],j] = OOB.votes$allWeightedVotes[OOB.matrix[idxJ,1],j]
					}
				}				
				OOB.object <- majorityClass(OOB.votes2, rf.classes, m.classwt = OOB.allWeightedVotes, m.imbalance = classCutOff, 
					m.regression = rf.regression)			
			}			
			OOB.pred = OOB.object$majority.vote
			if (BreimanBounds)
			{
				strengthCorr.object <- strength_and_correlation(s.trainLabels = trainLabels, OOB.votes2, OOB.object, rf.classes, 
					s.regression = rf.regression, s.parallelPackage = parallelPackage[1])
			}
			if (rf.regression)
			{  
				OOB.confusion = "only prediction error for regression. See ..$pred.error"				
				MSE <- L2Dist(OOB.pred[OOB.val], trainLabels[OOB.val])/n.OOB				
				pred.error = round(MSE,4)
				percent.varExplained = max(0, 100*round(1 - MSE/var(trainLabels[OOB.val]),4))
			}
			else
			{
				OOB.pred = if ( min(trainLabels) == 0) { OOB.pred - 1 }  else { OOB.pred }				
				OOB.confusion <- confusion.matrix(OOB.pred[OOB.val], trainLabels[OOB.val])
				pred.error <- generalization.error(OOB.confusion)
			}
		}
		if (variableImportance)
		{	
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
			gen.rufObject = new.rufObject$Tree	
		}
		else
		{	
			if (rf.regression)
			{ 
				if (BreimanBounds)
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2,	pred.error = pred.error, 
					percent.varExplained = percent.varExplained, regression = rf.regression) )
				}
				else
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2,	pred.error = pred.error, percent.varExplained = percent.varExplained,  
					regression = rf.regression) )
				}
			}
			else
			{ 	
				if (BreimanBounds)
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, regression = rf.regression) )
				}
				else
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2, pred.error = pred.error, regression = rf.regression) )
				}					
			}
		}		
	} 
	else
	{
		rufObject = vector('list', ntree)			
		# rufObject <- foreach(icount(ntree), .export = export, .options.smp  = smpopts, .inorder = FALSE, .multicombine = TRUE, .packages = "rUniformForestCppClass") %dopar%
		rufObject <- foreach(i = 1:ntree, .export = export, .options.smp  = smpopts, .inorder = FALSE, .multicombine = TRUE) %dopar%
		{	
			uniformDecisionTree(trainData, trainLabels, nodeMinSize = nodeMinSize, maxNodes = maxNodes, rf.features = features, 
				getSplitAt = splitAt, regression = rf.regression, bootstrap = rf.bootstrap, subagging = rf.subagging, treeDepth = depth, treeClasswt = classwt, treeOverSampling = rf.overSampling, targetClass = rf.targetClass, x.bagging = rf.x.bagging, treeRebalancedSampling = rf.rebalancedSampling, random.combination = rf.randomCombination, randomFeature = rf.randomFeature, treeCatVariables = whichCatVariables, outputPerturbation = rf.outputPerturbationSampling, 
				featureSelectionRule = rf.featureSelectionRule, treeDepthControl = depthControl)$Tree
		}					
		if (variableImportance) {	gen.rufObject = rufObject	}
		else 
		{	
			stopCluster(Cl) 			
			return( list(object = rufObject, regression = rf.regression))	
		}
	}	
	if (variableImportance)
	{
	   	if (ntree <= 100) {  threads = 1 }						
		varImpMatrix1 <- unlist(lapply(gen.rufObject, function(Z) Z[,"split var"]))		
		if (rf.regression) 	{  varImpMatrix2 <- unlist(lapply(gen.rufObject, function(Z) Z[,"L2Dist"]))/1000 }
		else {	varImpMatrix2 <- unlist(lapply(gen.rufObject, function(Z) Z[,"Gain"])) }			
		varImpMatrix <- cbind(varImpMatrix1, varImpMatrix2)
		if (!rf.regression)
		{   
			predMatrix <- foreach(i = 1:ntree, .options.smp = smpopts, .inorder = TRUE, .combine = rbind, .multicombine = TRUE) %dopar%	
			{ 	
				predIdx = which(gen.rufObject[[i]][,"left daughter"] == 0)
				predClass = gen.rufObject[[i]][predIdx, "prediction"]				
				predVar = vector(length = length(predIdx))
				for (j in seq_along(predIdx))
				{ 
					if ((predIdx[j] %% 2) == 0) { predVar[j] = gen.rufObject[[i]][which(gen.rufObject[[i]][,"left daughter"] == predIdx[j]), "split var"] }
					else {	predVar[j] = gen.rufObject[[i]][which(gen.rufObject[[i]][,"right daughter"] == predIdx[j]), "split var"] }
				}					
				cbind(predVar, predClass)
			}	
		}
		stopCluster(Cl)
		varImpMatrix = varImpMatrix[-which(varImpMatrix[,1] == 0),]		
		na.gain = which(is.na(varImpMatrix[,2]))
		if (length(na.gain) > 0)
		{ 	varImpMatrix = varImpMatrix[-na.gain,]	}			
		rf.var = unique(sortCPP(varImpMatrix[,1]))
		n.var = length(rf.var)			
		if (!rf.regression)
		{	
			var.imp.output <- matrix(NA, n.var, 4)
			for (i in 1:n.var)
			{  	
				classTable = sort(table(predMatrix[which(rf.var[i] == predMatrix[,1]),2]),decreasing = TRUE)
				classMax = as.numeric(names(classTable))[1]
				classFreq =	(classTable/sum(classTable))[1]
				var.imp.output[i,] =  c(rf.var[i], sum(varImpMatrix[which(rf.var[i] == varImpMatrix[,1]),2]), classMax, classFreq)
			}
		}
		else
		{	
			var.imp.output <- matrix(NA, n.var, 2)
			for (i in 1:n.var)
			{  	var.imp.output[i,] = c(rf.var[i], sum(varImpMatrix[which(rf.var[i] == varImpMatrix[,1]),2]))	}
		}
		var.imp.output = var.imp.output[order(var.imp.output[,2], decreasing = TRUE),]
		var.imp.output = round(cbind( var.imp.output, 100*var.imp.output[,2]/max(var.imp.output[,2])),2)
		percent.importance = round(100*var.imp.output[,2]/sum(var.imp.output[,2]),0)
		var.imp.output[,2] = round(var.imp.output[,2],0)		
		var.imp.output = cbind(var.imp.output, percent.importance)		
		if (rf.regression) 
		{	colnames(var.imp.output) = c("variables", "score", "percent", "percent importance")	}
		else
		{	
			colnames(var.imp.output) = c("variables", "score", "class", "class frequency", "percent", "percent importance")
			row.names(var.imp.output) = NULL	
		}		
		var.imp.output = data.frame(var.imp.output)		
		if (!is.null(colnames(trainData))) 	{	var.imp.output[,1] = colnames(trainData)[var.imp.output[,1]]	}		
		if (use.OOB)
		{	
			if (rf.regression)
			{
				if (BreimanBounds)
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, 
					percent.varExplained = percent.varExplained,  variableImportance = var.imp.output, regression = rf.regression) )
				}
				else
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2, pred.error = pred.error, percent.varExplained = percent.varExplained, 
					variableImportance = var.imp.output, regression = rf.regression) )
				}
			}
			else
			{
				if (BreimanBounds)
				{
					return( list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, 
					variableImportance = var.imp.output, regression = rf.regression) )
				}
				else
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2, pred.error = pred.error, variableImportance = var.imp.output, regression = rf.regression) )	
				}
			}
		}
		else
		{	return(list(object = gen.rufObject, variableImportance = var.imp.output, regression = rf.regression) )	}
	}
}

rUniformForestPredict <- function(object, X, 
	type = c("response", "prob", "votes", "confInt", "ranking", "quantile", "truemajority", "all"),
	classcutoff = c(0,0), 
	conf = 0.95,
	whichQuantile = NULL,
	rankingIDs = NULL,
	threads = "auto", 
	parallelpackage = "doParallel", ...)
{
	object <- filter.object(object)
	X <- if (is.vector(X)) { t(X) } else  { X }
	X <- fillVariablesNames(X)	
	if (!is.null(object$formula) & (length(object$variablesNames) != ncol(X)))
	{
		mf <- model.frame(data = as.data.frame(cbind(rep(1, nrow(X)),X)))
		X <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
	}
	else
	{
		if (length(object$variablesNames) != ncol(X)) 
		{	
			cat("data to predict have not the same dimension than data that have been computed by the model.\n Relevant variables have been extracted.\n")
			newVars <- rmNA(match(object$variablesNames, colnames(X)))
			if (length(newVars) == 0) 
			{ stop("No relevant variable has been found. Please give to both train and test data the same column names.\n") }
			X = X[,newVars, drop = FALSE]
		}
	}
	classes = object$classes	
	flag1 = flag2 = 1
	for (i in 1:ncol(object$forestParams))
	{			
		classCutOffString = as.vector(object$forestParams[which(row.names(object$forestParams) == "classcutoff"), i])		
		if ((classCutOffString != "FALSE") & (as.numeric(classcutoff[2]) == 0))
		{
			classCutOffString = rm.string(classCutOffString, "%")
			classCutOffString = rm.string(classCutOffString, "Class ")
			classCutOffString = strsplit(classCutOffString, ",")
			classcutoff[1] = which(classes == classCutOffString[[1]][1])
			classcutoff[2] = 0.5/(as.numeric(classCutOffString[[1]][2])/100)
		}
		else
		{	
			if ((as.numeric(classcutoff[2]) != 0) & (i == 1)) 
			{  	
				classcutoff <- c(which(classes == as.character(classcutoff[1])), 0.5/(as.numeric(classcutoff[2]))) 
				if (length(classcutoff) == 1) { stop("Label not found. Please provide name of the label instead of its position") }
				if (i > 1) { cat("Only last cutoff will be used in incremental random uniform forest.\n") }
			}		
		}				
		classwtString = as.vector(object$forestParams[which(row.names(object$forestParams) == "classwt"),i])	
		classwt = FALSE
		if (is.na(as.logical(classwtString))) {  classwt = TRUE; flag1 = flag2 = -1  }
		else 
		{	flag2 = 1	}		
		if (flag1 != flag2) 
		{ stop("Incremental random Uniform forest. Class reweighting must remain (or miss) in all forests.\n") }	
		randomCombinationString = as.vector(object$forestParams[which(row.names(object$forestParams) == "randomcombination"),i])
		if  ((randomCombinationString != "FALSE") )
		{	
			if (i == 1)
			{
				random.combination = as.numeric(unlist(strsplit(randomCombinationString, ",")))
				nbCombination = length(random.combination)/3
				X <- randomCombination(NAfactor2matrix(X), combination = random.combination[1:(length(random.combination) - nbCombination)], 
					weights = random.combination[(length(random.combination) - nbCombination + 1):length(random.combination)])
			}
			else	{	cat("Only first random combination will be used in incremental random uniform forest.\n") 	}
		}
	}	
	r.threads = threads	
	predObject <- randomUniformForestCore.predict(object, X, rf.aggregate = TRUE, pr.imbalance = classcutoff, pr.threads = threads, 	
	pr.classwt = classwt, pr.parallelPackage = parallelpackage[1])		
	object = filter.forest(object)	
	if (type[1] == "response") 
	{  
		predictedObject = predObject$majority.vote
		if (!is.null(classes)) {  predictedObject = as.factor(predictedObject); levels(predictedObject) = classes }
	}	
	if (type[1] == "truemajority") 
	{
		if (object$regression)
		{ stop( "True majority vote can only be computed for classification.") }
		else
		{	
			nbObs = nrow(predObject$all.votes)
			trueMajorityVotes = rep(0,nbObs)
            alpha1 = (1 - conf)/2
			alpha2 = conf + alpha1			
			for (i in 1:nbObs)
			{
				expectedClasses = unlist(lapply(predObject$votes.data, function(Z) Z[i,1]))	
				votingMembers = unlist(lapply(predObject$votes.data, function(Z) Z[i,2]))
				outliers = c(quantile(votingMembers, alpha1), quantile(votingMembers, alpha2))
				idxWoOutliers = which(votingMembers <= outliers[1] | votingMembers >= outliers[2])
				votes = cbind(expectedClasses[-idxWoOutliers], votingMembers[-idxWoOutliers])
				colnames(votes) = c("expectedClasses", "votingMembers")
			    trueMajorityVotes[i] = which.max(by(votes, votes[,"expectedClasses"], sum))
			}
			predictedObject = as.factor(trueMajorityVotes); levels(predictedObject) = classes 
		}
	}	
	if (type[1] == "confInt")	
	{  
		if (!object$regression)
		{ stop( "confidence interval can only be computed for regression") }
		else
		{	
			alpha1 = (1 - conf)/2
			alpha2 = conf + alpha1
			Q_1 = apply(predObject$all.votes, 1, function(Z) quantile(Z, alpha1))
			Q_2 = apply(predObject$all.votes, 1, function(Z) quantile(Z, alpha2))
			SD = apply(predObject$all.votes, 1, function(Z) sd(Z))
			predictedObject = data.frame(cbind(predObject$majority.vote, Q_1, Q_2, SD))
			Q1Name = paste("LowerBound", "(q = ", round(alpha1,3), ")", sep ="")
			Q2Name = paste("UpperBound", "(q = ", round(alpha2,3), ")",sep ="")
			colnames(predictedObject) = c("Estimate", Q1Name, Q2Name, "standard deviation")
		}		
	}
	if (type[1] == "quantile")	
	{  
		if (!object$regression)
		{ stop( "quantile(s) can only be computed for regression") }
		else
		{	
			if (!is.numeric(whichQuantile)) { stop( "Please provide a numeric value between 1 and 99") }			
			if ( (whichQuantile < 1) | (whichQuantile > 99))  { stop( "Please provide a numeric value between 1 and 99") }			
			predictedObject = apply(predObject$all.votes, 1, function(Z) quantile(Z, whichQuantile/100))
		}		
	}	
	if (type[1] == "all")	
	{  
		predictedObject = predObject
		if (object$regression)
		{ 
			stdDev = apply(predObject$all.votes, 1, sd)
			confIntObject = data.frame(cbind(predObject$majority.vote, apply(predObject$all.votes, 1, function(Z) quantile(Z,0.025)), 
			apply(predObject$all.votes, 1, function(Z) quantile(Z,0.975)), stdDev))			
			colnames(confIntObject) = c("Estimate", "LowerBound", "UpperBound", "Standard deviation")
			predictedObject$confidenceInterval = confIntObject
		}
	}	
	if (type[1] == "votes")  {	predictedObject = predObject$all.votes }	
	if ( (type[1] == "prob") & (!object$regression) )
	{	
		numClasses = sort(unique(predObject$all.votes[,1]))
		predictedObject = round(getVotesProbability2(predObject$all.votes, numClasses), 2) 
		colnames(predictedObject) = classes
	}	
	if ( (type[1] == "prob") & (object$regression) ) { stop("Probabilities can not be computed for regression") }	
	if ( (type[1] == "ranking") & (!object$regression) )
	{
		predictedObject = predictedObject2 = predObject$majority.vote
		if (!is.null(classes)) 
		{  
			if (!is.numeric(as.numeric(classes)))
			{	stop("Class are not numeric values or factor of numeric values. For Ranking, numeric values are needed as an equivalent of each class.") }
			else
			{
				minClass = min(as.numeric(classes))
				minPred = min(predictedObject)
				maxClass = max(as.numeric(classes))
				maxPred = max(predictedObject)				
				if (minClass < (minPred - 1))
				{	
					predictedObject[predictedObject != minPred] = predictedObject[predictedObject != minPred] - 1
					predictedObject = predictedObject - 1
				}
				else
				{
					if (minClass < minPred)
					{ 
				       predictedObject = predictedObject - 1
					}
				}		
				if (maxClass > maxPred)
				{	predictedObject[predictedObject == maxPred] = maxPred  }
			}
		}
		else { stop("Class not found. Please check if model is computed as classification, by looking if Y is set as factor.") }			
		numClasses = sort(unique(predObject$all.votes[,1]))
		probabilities = round(getVotesProbability2(predObject$all.votes, numClasses), 2) 
		colnames(probabilities) = classes
		countPredObject = table(predictedObject2)		
		majorityVote = as.numeric(names(which.max(countPredObject)))				
		n = length(predictedObject)
		followIdx = 1:n		
		if (!is.null(rankingIDs))
		{	
			rankingObject = cbind(followIdx, rankingIDs, predictedObject, probabilities)	
			if (length(classes) > 2)
			{  
				minoritiesProbability = rowSums(probabilities[,-majorityVote])
				rankingObject = cbind(rankingObject, minoritiesProbability)
			}			
			colnames(rankingObject)[1] = "idx"
			if (is.vector(rankingIDs) | is.factor(rankingIDs)) 
			{ 
				lengthRankingIDS = 1 
				colnames(rankingObject)[2] = "ID"
				colnames(rankingObject)[3] = "majority vote"
				cases = sort(unique(rankingIDs))
			}
			else 
			{ 
				lengthRankingIDS = ncol(rankingIDs)
				if (is.null(colnames(rankingIDs)))
				{	colnames(rankingObject)[2:(2+lengthRankingIDS)] = "ID"	}			
				
				colnames(rankingObject)[lengthRankingIDS+2] = "majority vote"
				cases = sort(unique(rankingIDs[,1]))
			}			
			if (length(classes) > 2)
			{  minorityIdx = which(colnames(rankingObject) == "minoritiesProbability") }
			else
			{ 	minorityIdx = which(colnames(rankingObject) == classes[-majorityVote])	}			
			lengthCases = length(cases)
			subCases = vector('list', lengthCases)
			for (i in 1:lengthCases)
			{	subCases[[i]] = which(rankingObject[,2] == cases[i]) }							
			rankingOutputObject <- matrix(NA, n, ncol(rankingObject) + 1)
			for (i in 1:lengthCases)
			{
				if (length(subCases[[i]]) > 1)
				{
					rankingOutputObject[subCases[[i]],] <- as.matrix(cbind(sortMatrix(rankingObject[subCases[[i]],], minorityIdx, decrease = TRUE), 
					1:length(subCases[[i]])))
				}
				else
				{	rankingOutputObject[subCases[[i]],] <- c(rankingObject[subCases[[i]],], 1)	}
			}			
			rankingOutputObject = as.data.frame(rankingOutputObject)			
			for (j in 1:ncol(rankingOutputObject))
			{
				if (is.factor(rankingOutputObject[,j]) & is.numeric(as.numeric(as.vector(rankingOutputObject[,j]))))
				{	rankingOutputObject[,j] = as.numeric(as.vector(rankingOutputObject[,j])) }
			}			
			colnames(rankingOutputObject) = c(colnames(rankingObject), "rank")
			predictedObject = sortMatrix(rankingOutputObject,1)	
		}
		else
		{	
			rankingObject = cbind(followIdx, predictedObject, probabilities)			
			if (length(classes) > 2)
			{  
				minoritiesProbability = rowSums(probabilities[,-majorityVote])
				rankingObject = cbind(rankingObject, minoritiesProbability)
			}			
			colnames(rankingObject)[1] = "idx"
			colnames(rankingObject)[2] = "majority vote"			
			if (length(classes) > 2)
			{  minorityIdx = which(colnames(rankingObject) == "minoritiesProbability") }
			else
			{	minorityIdx = which(colnames(rankingObject) == classes[-majorityVote]) }			
			rankingOutputObject = sortMatrix(rankingObject, minorityIdx, decrease = TRUE)	
			predictedObject = cbind(rankingOutputObject, 1:n)
			colnames(predictedObject)[ncol(predictedObject)] = "rank"
			predictedObject = sortMatrix(predictedObject,1)	
		}
	}		
	if ( (type[1] == "ranking") & (object$regression) )
	{ stop("Ranking currently available for classification tasks") }	
	predictedObject
}

randomUniformForestCore.predict <- function(object, X, 
rf.aggregate = TRUE, 
OOB.idx = FALSE, 
pr.imbalance = c(0,0),
pr.regression = TRUE,
pr.classwt = FALSE,  
classes = -1,
pr.export = c("onlineClassify", "predictDecisionTree", "majorityClass", "rmNA", "mergeLists", "classifyMatrixCPP", "L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP",	"entropyInformationGainCPP", "runifMatrixCPP"), 
pr.threads = "auto", 
pr.parallelPackage = "doParallel")
{
	#require(rUniformForestCppClass)	
	if (!OOB.idx)
	{
		if (is.matrix(X))
		{
			matchNA = (length(which(is.na(X))) > 0)		
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values. Please use one of many models available if acccuracy is needed\n")
				X <- na.impute(X)
			}
		}
		else
		{ 
			X <- NAfactor2matrix(X)	
			matchNA = (length(which(is.na(X))) > 0)		
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values. Please use one of many models available if acccuracy is needed\n")
				X <- na.impute(X)
			}
		}
		if (!is.null(object$logX))
		{
			if (object$logX)
			{ 
			   if (is.null(object$categoricalvariables))  {  X <- generic.log(X) }
			   else  {  X[,-object$categoricalvariables] <- generic.log(X[,-object$categoricalvariables]) }
			}
		}		
		object = filter.forest(object)
	}				
	if (!is.null(object$OOB) & (!OOB.idx))
	{  OOB.predicts = object$OOB.predicts; pr.regression = object$regression; object = object$object }
	else
	{
		if (!OOB.idx) 	{ 	pr.regression = object$regression; object = object$object   }			
		if (!is.null(object$variableImportance) & (OOB.idx)) {	object = object$object	}		
	}			
	n = nrow(X)			
	if (!pr.regression)
	{
		if (classes[1] <  1)
		{
			classes = sort(as.numeric(rownames(table(object[[sample(1:length(object),1)]][,"prediction"]))))
			if (classes[1] == 0) {  classes = classes[-1]  }
		}			
		l.class = length(classes)
		class.occur = rep(0,l.class)
	}	
	{
		ntree = length(object)
		pred.X = vector("list", ntree)
		all.votes = nodes.length = nodes.depth = matrix(data = Inf, nrow = n, ncol = ntree)
		majority.vote = vector(length = n)		
		pred.X <- lapply(object, function(Z) predictDecisionTree(Z, X))									
		for (i in 1:ntree)
		{	
			all.votes[,i] = pred.X[[i]][,1]
			nodes.length[,i] = pred.X[[i]][,2]
			nodes.depth[,i] = pred.X[[i]][,3]	
		}		
		if (pr.classwt)
		{
			allWeightedVotes = matrix(data = Inf, nrow = n, ncol = ntree)
			for (i in 1:ntree)	{  allWeightedVotes[,i] = object[[i]][nodes.depth[,i], "avgLeafWeight"]	}
		}
		else
		{  allWeightedVotes = NULL }
			
		if (pr.regression) {  majority.vote = rowMeans(all.votes)	}
		else {	majority.vote = majorityClass(all.votes, classes, m.imbalance = pr.imbalance, m.classwt = allWeightedVotes)$majority.vote  }
	}				
	if (rf.aggregate & (!OOB.idx))
	{	
		return(list(all.votes = all.votes, majority.vote = majority.vote, nodes.depth = nodes.depth, nodes.length = nodes.length, 
		votes.data = pred.X))	
	}	
	else
	{
		if (OOB.idx) 
		{  
			if (pr.classwt)	{ return(list(all.votes = all.votes, allWeightedVotes = allWeightedVotes))	}
			else { return(all.votes)	}
		}
		else  {  return(majority.vote) 	}		
	}
}

majorityClass <- function(all.votes, classes, m.regression = FALSE, m.imbalance = c(0,0), m.classwt = NULL, m.threads = "auto")
{
	if (m.regression)
	{	
		all.votes[all.votes == Inf] = NA;  majority.vote = apply(all.votes, 1, function(all.votes) mean(rmNA(all.votes))); class.counts = NULL	
		return(list(majority.vote = majority.vote, class.counts = class.counts))
	}
	else
	{
		n = nrow(all.votes)
		majority.vote = vector(length = n)
		l.class = length(classes)
		class.occur = vector(length = l.class)
		class.countsMajorityVote = trueClass.counts = matrix(NA, n, l.class)
		for (i in 1:n)
		{ 
			if (!is.null(m.classwt))
			{	
				for (j in 1:l.class) 
				{
					idx = rmNA(which(all.votes[i,] == classes[j]))
					l.idx = length(idx)
									
					if (l.idx > 0) {   class.occur[j] = sum(m.classwt[i,idx]) }
					else	{ class.occur[j] = 0 }
				}				
			}
			else
			{
				for (j in 1:l.class) { 	class.occur[j] <- sum(all.votes[i,] == classes[j])  }
			}				
			if (m.imbalance[1] != 0)
			{	class.occur[m.imbalance[1]] = floor(class.occur[m.imbalance[1]]*m.imbalance[2]) } 
			majority.vote[i] = which.max(class.occur)
			class.countsMajorityVote[i,] = c(class.occur[majority.vote[i]], class.occur[-majority.vote[i]])
			trueClass.counts[i,]  = class.occur
		}
		return(list(majority.vote = majority.vote, class.counts = class.countsMajorityVote, trueClass.counts = trueClass.counts))
	}	
}

getVotesProbability  <- function(X, classes) (majorityClass(X, classes)$class.counts/ncol(X))

getVotesProbability2  <- function(X, classes) (majorityClass(X, classes)$trueClass.counts/ncol(X))

strength_and_correlation <- function(OOB.votes, OOB.object, 
	rf.classes, 
	s.trainLabels = NULL, 
	s.regression = FALSE, 
	output = NULL, 
	s.threads = "auto", 
	s.parallelPackage = "doParallel" )
{
	j = NULL
	dimOOBvotes = dim(OOB.votes)	
	n.OOB = dimOOBvotes[1]
	p.OOB = dimOOBvotes[2]	
	if (s.regression)
	{
		OOB.votes[OOB.votes == Inf] = NA		
		if (is.null(s.trainLabels))  { Y = rowMeans(OOB.votes, na.rm =TRUE)	}
		else {  Y = s.trainLabels  }			
		expectedSquaredErrorOverTrees = colMeans(  (OOB.votes - Y)^2, na.rm = TRUE)
		PE.tree = sum(expectedSquaredErrorOverTrees)/length(expectedSquaredErrorOverTrees)
		sd.T =  sqrt(expectedSquaredErrorOverTrees)
		mean.corr = mean(rowMeans(cor((Y - OOB.votes), use = "pairwise.complete.obs")))		
		if (is.na(mean.corr))
		{	
			cat("Not enough data to compute average correlation for trees. Error is then prediction error of a tree.\n")
			return(list(PE.forest = (mean(sd.T))^2,  PE.max = PE.tree, PE.tree = PE.tree, mean.corr = mean.corr))
		}
		else
		{	return(list(PE.forest = mean.corr*(mean(sd.T))^2,  PE.max = mean.corr*PE.tree, PE.tree = PE.tree, mean.corr = mean.corr))	}
	}
	else
	{
		#require(parallel)		
		max_threads = detectCores()		
		if (s.threads == "auto")
		{	s.threads  = max(1, max_threads - 1)  }
		else
		{
			if (max_threads < s.threads) 
			{	cat("Warning : number of threads indicated by user was higher than logical threads in this computer.\n") }
		}
		{	
			#require(doParallel)
			Cl <- makePSOCKcluster(s.threads, type ="SOCK")
			registerDoParallel(Cl)
		}
		chunkSize  <-  ceiling(p.OOB/getDoParWorkers())
		smpopts  <-  list(chunkSize = chunkSize)	
		p.new.OOB = apply(OOB.votes, 1, function(OOB.votes) sum(OOB.votes != Inf))			
		if (length(rf.classes) == 2)
		{
			Q.x.1 = OOB.object$class.counts[,1]/rowSums(OOB.object$class.counts)
			rawStrength = 2*Q.x.1 - 1			
			Tk.1 = Tk.2 = matrix (data = NA, ncol = p.OOB, nrow = n.OOB)
			OOB.votes.1 = cbind(OOB.votes, OOB.object$majority.vote)
			Tk.1 <- foreach(j = 1:p.OOB, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar%
			apply(OOB.votes.1[,-j], 1, function(Z) sum(Z[-p.OOB] == Z[p.OOB]))
			Tk.1 = Tk.1/p.new.OOB 
			Tk.2 = 1 - Tk.1
			stopCluster(Cl)
		}
		else
		{
			Q.x.j = (apply(OOB.object$class.counts[,-1], 1, max))/rowSums(OOB.object$class.counts)
			Q.x.y = OOB.object$class.counts[,1]/rowSums(OOB.object$class.counts)
			rawStrength = Q.x.y - Q.x.j			
			maj.class.j = vector(length = n.OOB)
			for (i in 1:n.OOB)
			{  
				second.class.max = which.max(OOB.object$class.counts[i,-1])
				maj.class.j[i] = rf.classes[-OOB.object$majority.vote[i]][second.class.max]
			}
			OOB.votes.1 = cbind(OOB.votes, OOB.object$majority.vote)
			OOB.votes.2 = cbind(OOB.votes,  maj.class.j)			
			ZZ <- function(j) cbind(apply(OOB.votes.1[,-j], 1, function(Z) sum(Z[-p.OOB] == Z[p.OOB])), apply(OOB.votes.2[,-j], 1, 
				function(Z) sum(Z[-p.OOB] == Z[p.OOB])))					
			Tk <- foreach(j = 1:p.OOB, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar% ZZ(j)		
			mixIdx = getOddEven(1:ncol(Tk))			
			Tk.1 = Tk[,mixIdx$odd]
			Tk.2 = Tk[,mixIdx$even]
			Tk.1 = Tk.1/p.new.OOB
			Tk.2 = Tk.2/p.new.OOB
			stopCluster(Cl) 
		}				
		p1 = colMeans(Tk.1)
		p2 = colMeans(Tk.2)		
		strength = mean(rawStrength)
		varStrength =  var(rawStrength)
		sd.T = ((p1 + p2 + (p1 - p2)^2))^0.5		
		mean.corr = varStrength / (mean(sd.T))^2			
		PE.est = mean.corr*(1 -	strength^2)/strength^2		
		return(list(PE = PE.est, avg.corr = mean.corr, strength = strength, std.strength = sqrt(varStrength)))
	}
}

monitorOOBError  <- function(OOB.votes, Y, regression = FALSE, threads = "auto")
{
	j = NULL
	n = length(Y)
	n.RS = sqrt(n)
	p = ncol(OOB.votes)	
	#require(parallel)	
	max_threads = detectCores()		
	if (threads == "auto")
	{	
		if (max_threads == 2) { threads = max_threads }
		else {	threads  = max(1, max_threads - 1)  }
	}
	else
	{
		if ((max_threads) < threads) 
		{	 cat("Warning : number of threads indicated by user was higher than logical threads in this computer\n.") }
	}	
	#require(doParallel)			
	Cl = makePSOCKcluster(threads, type = "SOCK")
	registerDoParallel(Cl)	
	chunkSize  <-  ceiling(p/getDoParWorkers())
	smpopts  <- list(chunkSize  =  chunkSize)
	if (!regression)
	{
		Y = as.numeric(Y)
		classes = sort(unique(Y))
		S1 = sample(classes,1)
		S2 = classes[-S1][1]
		
		OOBmonitor <- foreach(j = 1:(p-1), 
		.export = c("generalization.error", "confusion.matrix", "majorityClass", "rmNA"), .options.smp = smpopts, .combine = rbind) %dopar%
		{
			Estimate <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes)$majority.vote, Y))
			C1 <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes, m.imbalance =c(1, 1.5))$majority.vote, Y))
			C2 <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes, m.imbalance= c(2, 1.5))$majority.vote, Y))
			
			t(c(Estimate, C1, C2))
		}
		
		E0 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes)$majority.vote, Y))
		C1 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes, m.imbalance =c(1, 1.5))$majority.vote, Y))
		C2 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes, m.imbalance= c(2, 1.5))$majority.vote, Y))
		OOBmonitor = rbind(t(c(E0, C1, C2)),OOBmonitor)
	}
	else
	{
		OOBmonitor <- foreach(j = 1:(p-1), .export = c("generalization.error", "confusion.matrix", "majorityClass", "rmNA", "L2Dist"), 
			.options.smp = smpopts, .combine = c) %dopar%
		{
			Z = majorityClass(OOB.votes[,1:(j+1)], 0, m.regression = TRUE)$majority.vote
			NAIdx = which(is.na(Z))
			if (length(NAIdx) > 0) { L2Dist(Z[-NAIdx], Y[-NAIdx])/length(Z[-NAIdx])}
			else { L2Dist(Z, Y)/length(Z)  }
		}		
		Z = majorityClass(matrix(OOB.votes[,1]), 0, m.regression = TRUE)$majority.vote
		NAIdx = which(is.na(Z))
		E0 = if (length(NAIdx) > 0) { L2Dist(Z[-NAIdx], Y[-NAIdx])/length(Z[-NAIdx])} else { L2Dist(Z, Y)/length(Z)  }
		OOBmonitor = c(E0,OOBmonitor)
	}	
	stopCluster(Cl)
	return(OOBmonitor)
}

weightedVote <- function(all.votes, idx = 2, granularity = 2)
{
	all.votes <- round(all.votes, granularity)	
	apply(all.votes, 1, function(Z)
		{
			A = sort(table(rmInf(Z)), decreasing = TRUE)
			B = as.numeric(names(A))[1:idx]
			sum( B * (B/sum(abs(B))) )
		}
	)
}

weightedVoteModel <- function(votes, majorityVote, Y = NULL, nbModels = 1, idx = 1 , granularity = 1, train = TRUE, models.coeff = NULL)
{
	if (train)
	{
		if (is.null(Y))
		{ stop("output is neded to build model") }		
		if (nbModels == 1)
		{
			default.model = weightedVote(votes, idx = idx, granularity = granularity)
			lin.model = lm(Y ~ cbind(majorityVote, default.model))
		}
		else
		{
			models = matrix(data = NA, ncol = nbModels + 2, nrow = length(Y))			
			models[,1] = majorityVote
			models[,2] = weightedVote(votes, idx = idx, granularity = granularity)  
			for (j in 3:(nbModels+2))
			{	models[,j] = weightedVote(votes, idx = max(j,idx), granularity = j ) }
			lin.model = lm(Y ~ models)
		}		
		return(lin.model)
	}
	else
	{
		p = length(models.coeff)
		models = matrix(data = NA, ncol = p, nrow = length(majorityVote))			
		models[,1] = rep(1,length(majorityVote))
		models[,2] = majorityVote
		models[,3] = weightedVote(votes, idx = idx, granularity = granularity)		
		if (p > 3)
		{
			for (j in 4:p)
			{	models[,j] = weightedVote(votes,  idx =j  , granularity = j)	}
		}		
		newMajorityVote = apply(models,1, function(Z) sum(models.coeff*Z))		
		return(newMajorityVote )
	}
}

postProcessingVotes <- function(object, nbModels = 1, idx = 1, granularity = 1, predObject = NULL, swapPredictions = FALSE, X = NULL, imbalanced = FALSE)
{
	object <- filter.object(object)	
	if (rownames(object$forestParams)[1] == "reduceDimension") 
	{ stop("Post Processing does not work with objects coming from rUniformForest.big() function") }
	if (!object$forest$regression)
	{
		if (length(as.numeric(object$classes)) > 2)
		{	stop("Optimization curently works only for binary classification") }		
		if (is.null(X))
		{	stop("Please provide test data")	}
		else
		{
			if (is.null(predObject)) {	predObject <- predict(object, X, type = "all")	}
			else
			{
				if (is.null(predObject$all.votes)) 
				{ stop("Please provide full prediction object (option type = 'all' when calling predict() function).") }				
				if (is.null(predObject$majority.vote)) 
				{ stop("Please provide full prediction object (option type = 'all' when calling predict() function).") }
			}
			majorityVotePosition = which.max(table(predObject$majority.vote))
			numClasses = sort(unique(predObject$all.votes[,1]))
			probPred = round(getVotesProbability2(predObject$all.votes, numClasses), 2) 
			colnames(probPred) = object$classes
			if (imbalanced) 
			{ 	cutoff =  1 - ( mean(probPred[,2])/mean(probPred[,1]) )	}
			else 
			{ 	cutoff = 0.5/2*( mean(probPred[,2])/mean(probPred[,1]) + mean(probPred[,1])/mean(probPred[,2]) ) } 
			predVotes = predict(object, X, classcutoff = c(object$classes[majorityVotePosition], cutoff))			
			return(predVotes)
		}
	}	
	if (swapPredictions) 
	{ 
		object$forest$OOB.votes = predObject$forest$OOB.votes 
		object$forest$OOB.predicts = predObject$forest$OOB.predicts
	}	
	if (is.null(object$forest$OOB.votes)) 
	{ stop("No OOB data for post processing. Please enable OOB option and subsamplerate or bootstrap ('replace' option) when computing model.") }			
	if (is.null(object$predictionObject))
	{
		if (is.null(predObject)) { stop("Post processing can not be computed. Please provide prediction object") }
		else 
		{  
			if (is.null(predObject$majority.vote))
			{ stop("Post processing can not be computed. Please provide full prediction object (type = 'all') when calling predict()") }
			else
			{			
				meanEstimate = predObject$majority.vote
				allVotes = predObject$all.votes
			}
		}
	}
	else
	{   
		meanEstimate = object$predictionObject$majority.vote
		allVotes = object$predictionObject$all.votes
	}			
	Y = object$y
	OOBVotes = object$forest$OOB.votes	
	NAIdx = which.is.na(object$forest$OOB.predicts)	
	if (length(NAIdx) > 0)
	{
		Y = Y[-NAIdx]
		OOBVotes = OOBVotes[-NAIdx,]
		object$forest$OOB.predicts = object$forest$OOB.predicts[-NAIdx]
	}
	OOBMeanEstimate = object$forest$OOB.predicts
	L2DistOOBMeanEstimate <- L2Dist(OOBMeanEstimate, Y)/length(Y)	
	OOBMedianEstimate <- apply(OOBVotes, 1, function(Z) median(rmInf(Z)))
	L2DistOOBMedianEstimate <- L2Dist(OOBMedianEstimate, Y)/length(Y)	
	weightedModel <- weightedVoteModel(OOBVotes, OOBMeanEstimate, Y = Y, nbModels = nbModels, idx = idx, granularity = granularity)
	OOBWeightedVoteEstimate <- weightedVoteModel(OOBVotes, OOBMeanEstimate, train = FALSE, models.coeff = weightedModel$coefficients)	
	NAOOBWeightedVoteEstimate <- which.is.na(OOBWeightedVoteEstimate)	
	if (length(NAOOBWeightedVoteEstimate) > 0) 
	{ OOBWeightedVoteEstimate[NAOOBWeightedVoteEstimate] = OOBMeanEstimate[NAOOBWeightedVoteEstimate] }	
	L2DistOOBWeightedVoteEstimate <- L2Dist(OOBWeightedVoteEstimate, Y)/length(Y)	
	flagMedian =  ( (L2DistOOBMedianEstimate <= L2DistOOBMeanEstimate) | (L1Dist(OOBMedianEstimate, Y) <= L1Dist(OOBMeanEstimate, Y)) )
	flagWeighted = ( (L2DistOOBWeightedVoteEstimate <= L2DistOOBMeanEstimate) | (L1Dist(OOBWeightedVoteEstimate, Y) <= L1Dist(OOBMeanEstimate, Y)) )	
	if (flagMedian)
	{
		if (flagWeighted)
		{  weightedVoteEstimate <- weightedVoteModel(allVotes, meanEstimate, train = FALSE, models.coeff = weightedModel$coefficients) }
		else
		{  weightedVoteEstimate = rep(0, length(meanEstimate)) }
	}
	else
	{  	weightedVoteEstimate <- weightedVoteModel(allVotes, meanEstimate, train = FALSE, models.coeff = weightedModel$coefficients)  }	
	medianEstimate <- apply(allVotes, 1, function(Z) median(rmInf(Z)))	
	NAWeightedVoteEstimate <- which.is.na(weightedVoteEstimate)
	if (length(NAWeightedVoteEstimate) > 0) { weightedVoteEstimate[NAWeightedVoteEstimate] = meanEstimate[NAWeightedVoteEstimate] }	
	negativeBias = - (mean(OOBMeanEstimate) - Y)
	modelResiduals = Y - OOBMeanEstimate
	biasModel  = lm(modelResiduals ~ negativeBias) 
	biasSign = sign(meanEstimate - medianEstimate)   
	biasSign[which(biasSign == 0)] = 1
	residuals_hat = vector(length = ncol(OOBVotes))	
	residuals_hat <- apply(OOBVotes, 2, function(Z) sum(rmInf(Z))/length(rmInf(Z))) - mean(OOBMeanEstimate)
	MeanNegativeBias = -mean(residuals_hat^2)		
	biasCorrection = biasModel$coefficients[2]*biasSign*MeanNegativeBias + biasModel$coefficients[1]	
	if (flagMedian)
	{ 
		if (flagWeighted) { return( 0.5*(medianEstimate + weightedVoteEstimate + biasCorrection)) }
		else { return(medianEstimate) }
	}
	else
	{ return(weightedVoteEstimate + biasCorrection) }
}

localVariableImportance <- function(object, nbVariables = 2, Xtest = NULL, predObject = NULL, l.threads = "auto", l.parallelPackage = "doParallel")
{
	object = filter.object(object)
	rufObject = object$forest$object	
	if (is.null(object$predictionObject))
	{
		if (is.null(predObject)) 
		{ 
			if (is.null(Xtest))
			{ stop("Local variable importance can not be computed. please provide test data.") }
			else
			{
				predObject = predict(object, Xtest, type = "all")
				majority.vote = as.numeric(predObject$majority.vote) 
				pred.rufObject = predObject$votes.data
			}
		}
		else 
		{  
			if (is.null(predObject$majority.vote))
			{ 
				stop("Local variable importance can not be computed. Please provide full prediction object (type = 'all') when calling predict()") 
			}
			else
			{			
				majority.vote = as.numeric(predObject$majority.vote) 
				pred.rufObject = predObject$votes.data
			}
		}		
	}
	else
	{   
		majority.vote = as.numeric(object$predictionObject$majority.vote)
		pred.rufObject = object$predictionObject$votes.data
	}			
	ntree = length(rufObject)
	if (is.null(pred.rufObject))
	{ stop ("no data to evaluate importance") }
	else
	{
		n = dim(pred.rufObject[[1]])[1]
		{
			#require(parallel)
			threads = l.threads
			max_threads = detectCores()			
			if (threads == "auto")
			{	
				if (max_threads == 2) { threads = max_threads }
				else {	threads  = max(1, max_threads - 1)  }
			}
			else
			{
				if ((max_threads) < threads) 
				{	cat("Warning : number of threads indicated by user is higher than logical threads in this computer.\n") }
			}
			{
				#require(doParallel)				
				Cl = makePSOCKcluster(threads, type = "SOCK")
				registerDoParallel(Cl)
			}
			chunkSize  <- ceiling(ntree/getDoParWorkers())
			smpopts  <- list(chunkSize = chunkSize)
		}		
		varMatrix <- foreach(i = 1:ntree, .options.smp = smpopts, .combine = rbind, .multicombine = TRUE) %dopar%	
		{ 	
			predVar = vector(length = n)
			for (j in 1:n)
			{
				predIdx = pred.rufObject[[i]][j,3]
				predClass = rufObject[[i]][predIdx, "prediction"]
							
				if ((predIdx %% 2) == 0)
				{ 	predVar[j] = rufObject[[i]][which(rufObject[[i]][,"left daughter"] == predIdx),"split var"] }
				else
				{	predVar[j] = rufObject[[i]][which(rufObject[[i]][,"right daughter"] == predIdx),"split var"] }
			}						
			predVar
		}		
		varImp = varFreq = matrix(data = NA, ncol = nbVariables, nrow = n)
		varImp <- t(apply(varMatrix, 2, function(Z) as.numeric(names(sort(table(Z), decreasing = TRUE)[1:nbVariables])) ))
		varFreq <-  t(apply(varMatrix, 2, function(Z)  sort(table(Z), decreasing = TRUE)[1:nbVariables]/ntree ) )		
		NAIdx = which(is.na(varImp), arr.ind = TRUE)		
		if (length(NAIdx[,1]) > 0)
		{	
			varImp[NAIdx] = apply(NAIdx, 1, function(Z) { 
												newZ = rmNA(varImp[Z])
												if (length(newZ) == 1) { rep(newZ, length(Z)) } else  { sample(newZ,1)  }
											}
							)
			varFreq[NAIdx] = apply(NAIdx, 1, function(Z) { 
												newZ = rmNA(varFreq[Z])
												if (length(newZ) == 1) { rep(newZ, length(Z)) } else  { sample(newZ,1)  }
											}
							)				
		}		
		objectNames = objectFrequencyNames = vector()
		for (j in 1:nbVariables)
		{	
			objectNames[j] = paste("localVariable", j, sep = "")
			objectFrequencyNames[j] = paste("localVariableFrequency", j, sep = "") 
		}	
		variableImportance.object = cbind(majority.vote, varImp, varFreq)		
		stopCluster(Cl)		
		if (!object$forest$regression)
		{
			colnames(variableImportance.object) = c("class", objectNames, objectFrequencyNames)		
			classObject = as.numeric(names(table(variableImportance.object[,1])))		
			classMatrix = list()		
			for (i in 1:length(classObject))
			{		
				classMatrix[[i]] = round(sort(table(variableImportance.object[which(variableImportance.object[,1] == i),2])/
					sum(table(variableImportance.object[which(variableImportance.object[,1] == i),2])), decreasing = TRUE),2)
			}		
			orderedVar = unique(as.numeric(names(sort(unlist(classMatrix), decreasing = TRUE))))
			Variables = as.numeric(names(sort(unlist(classMatrix), decreasing = TRUE)))
			orderedImportance =  as.numeric(sort(unlist(classMatrix), decreasing = TRUE))			
			classVariableImportance = matrix(data = 0, ncol = length(classObject), nrow = length(orderedVar))
			for (i in 1:length(orderedVar))
			{	
				for (j in 1:length(classObject))
				{	
					idx = which( as.numeric(names(classMatrix[[j]])) == orderedVar[i])
					if (length(idx) > 0)
					{ classVariableImportance[i,j] = classMatrix[[j]][idx]	}
				}
			}		
			classVariableImportance = data.frame(classVariableImportance)
			rownames(classVariableImportance) = orderedVar
			colnames(classVariableImportance) = paste("Class", classObject, sep=" ")				
			return(list( obsVariableImportance = variableImportance.object, classVariableImportance = classVariableImportance ))
		}
		else
		{
			colnames(variableImportance.object) = c("estimate", objectNames, objectFrequencyNames)
			return(variableImportance.object)			
		}
	}
}
	
localTreeImportance <- function(rfObject, OOBPredicts = NULL, OOBVotes = NULL)
{
	rfObject = filter.object(rfObject)  
	if (!is.null(rfObject$OOB.votes))	{ OOBVotes = rfObject$OOB.votes }
	if (is.null( OOBVotes))	{ stop ( "OOB outputs missing. No optimization or weighted vote can be done") }
	if (!is.null(rfObject$OOB.predicts))	{	OOBPredicts = rfObject$OOB.predicts	}	
	n = length(OOBPredicts)  
	wPlus = wMinus = NULL
	for ( i in seq_along(OOBPredicts))
	{
		wPlus = c(wPlus, which(OOBVotes[i,] == OOBPredicts[i]))
		wMinus = c(wMinus, which(OOBVotes[i,] != OOBPredicts[i]))
	}	
	object = list(pweights = (table(wPlus)/n), nweights = (table(wMinus)/n))	
	return(object)
}

importance.randomUniformForest <- function(object, maxVar = 30, maxInteractions = 3, Xtest = NULL, predObject = NULL, ...)
{
	object <- filter.object(object) 
	maxInteractions = max(2, maxInteractions)	
	if (!is.null(Xtest))
	{
		if (!is.null(object$formula) & (length(object$variablesNames) != ncol(Xtest)))
		{
			mf <- model.frame(formula = object$formula, data = as.data.frame(Xtest))
			Xtest <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
			if (object$logX) { 	Xtest <- generic.log(Xtest) }	
		}
		else
		{		
			Xfactors <- which.is.factor(Xtest)
			catVarIdx <- which(rownames(object$paramsObject) == "categorical variables")
			Xtest = NAfactor2matrix(Xtest)
			matchNA = (length(which(is.na(Xtest))) > 0)				
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values\n")
				Xtest <- na.impute(Xtest)
			}
			if (object$logX)
			{ 	
				trueFactorsIdx <- which(Xfactors == 1)
				if (length(trueFactorsIdx) ==  0)  {  Xtest <- generic.log(Xtest) }
				else  
				{  
					if (!is.null(object$paramsObject[1,catVarIdx]))
					{	
						cat("Warning : all categorical variables have been ignored for the logarithm transformation.\nIf only some of them have been defined as categorical, variable importance will be altered.\n To avoid it, please define logarithm transformation outside of the model or ignore categorical variables, or set them to 'all'.\n")
						Xtest[,-trueFactorsIdx] <- generic.log(Xtest[,-trueFactorsIdx])						
					}
					else
					{ 	Xtest <- generic.log(Xtest)	}
				}
			}
		}
	}	
	if (!is.null(object$forest$variableImportance))
	{	
		par(las=1) 
		maxChar = floor(1 + max(nchar(object$variablesNames))/2)		
		par(mar=c(5, maxChar + 1,4,2)) 		
		varImportance1 = varImportance = object$forest$variableImportance
		if (!object$forest$regression)
		{
			varImportance[,"class"] = object$classes[as.numeric(varImportance[,"class"])]
			varImportance1[,"class"] = varImportance[,"class"]
		}		
		nVar = nrow(varImportance)
		if (nVar > maxVar) { varImportance = varImportance[1:maxVar,] }	
		barplot(varImportance[nrow(varImportance):1, "percent.importance"], horiz = TRUE, col = sort(heat.colors(nrow(varImportance)), decreasing = TRUE), names.arg = varImportance[nrow(varImportance):1,"variables"], xlab = "Relative Influence (%)",  
		main = "Variable importance based on information gain", border = NA)
		abline(v = 100/nVar, col='grey')		
		cat("\nGlobal Variable importance (", min(maxVar, nrow(varImportance)), "most important ) :\n")
		print(varImportance)		
		localVarImportance <- localVariableImportance(object, nbVariables = maxInteractions, Xtest = Xtest, predObject = predObject)		
		if (!object$forest$regression)
		{
			rownames(localVarImportance$classVariableImportance) = object$variablesNames[as.numeric(rownames(localVarImportance$classVariableImportance))]
			colnames(localVarImportance$classVariableImportance) = paste("Class ", 
				object$classes[as.numeric(rm.string(colnames(localVarImportance$classVariableImportance), "Class"))], sep = "")				
			obsVarImportance = data.frame(localVarImportance$obsVariableImportance)
			obsVarImportance[,1] = object$classes[obsVarImportance[,1]]			
		}
		else
		{	obsVarImportance = data.frame(localVarImportance) 	}		
		for (j in 2:(maxInteractions+1)) {	obsVarImportance[,j] = object$variablesNames[obsVarImportance[,j]] }		
		obsVarImportance2 = obsVarImportance
		fOrder = sort(table(obsVarImportance2[,2]), decreasing = TRUE)
        sOrder = sort(table(obsVarImportance2[,3]), decreasing = TRUE)				
		W1 = mean(obsVarImportance2[,grep("localVariableFrequency1", colnames(obsVarImportance2))[1]])
		W2 = mean(obsVarImportance2[,grep("localVariableFrequency2", colnames(obsVarImportance2))[1]])		
		minDim2 = min(length(fOrder), length(sOrder))
		partialDependence = matrix(NA, minDim2, minDim2)
		for (i in 1:minDim2) {  partialDependence[,i] = fOrder[i]*W1/W2 + sOrder[1:minDim2]  }		
		colnames(partialDependence) = names(fOrder)[1:minDim2]
		rownames(partialDependence) = names(sOrder)[1:minDim2]
		partialDependence = partialDependence/(2*nrow(obsVarImportance2))
		avg2ndOrder = colMeans(partialDependence)
		avg1rstOrder = c(rowMeans(partialDependence),0)
		partialDependence = rbind(partialDependence, avg2ndOrder)
		partialDependence = cbind(partialDependence, avg1rstOrder)		
		minDim = min(10, minDim2)
		varImportanceOverInteractions = vector()
		for (i in 1:minDim2)
		{
			idx = which(rownames(partialDependence)[i] == colnames(partialDependence))
			if (length(idx) > 0)  
			{	varImportanceOverInteractions[i] = 0.5*(avg1rstOrder[i] + avg2ndOrder[idx] + 2*partialDependence[i,idx])  }
			else 
			{	varImportanceOverInteractions[i] = avg1rstOrder[i]	}			
			names(varImportanceOverInteractions)[i] = rownames(partialDependence)[i]
		}				
		varImportanceOverInteractions = sort(varImportanceOverInteractions, decreasing = TRUE)
		if (!object$forest$regression)
		{
			cat("\n\n", "Variable importance over labels(", minDim, "most important variables ) :\n")
			print(localVarImportance$classVariableImportance[1:minDim,])
		}
		cat("\n\n Variables interactions over all observations and trees :", minDim, "most important variables at first (columns) and second (rows) order.") 
		cat("\nFor each variable (at each order) interaction with others is computed.\n")
		print(round(partialDependence[1:minDim,],4))	
		cat("\n\n Overall variables interaction :", minDim, "most important scores of overall interactions.\n")
		print(round(varImportanceOverInteractions[1:minDim],4))
		if (!object$forest$regression)
		{	cat("\n\n see ...$localVariableImportance$obsVariableImportance to get variable importance for each observation.")	}
		else
		{   cat("\n\n see ...$localVariableImportance to get variable importance for each observation.")	}
		cat("\n\n call partialDependenceOverResponses() function to get partial dependence over response\n for each variable among most important.\n")	
	}
	else
	{ stop("no variable importance defined in random uniform forest") }	
	importanceObject = list(globalVariableImportance = varImportance1, localVariableImportance = localVarImportance, 
		partialDependence = partialDependence, variableImportanceOverInteractions = varImportanceOverInteractions )		
	class(importanceObject) <- "importance"	
	importanceObject
}

partialImportance <- function(X, object, whichClass = NULL, threshold = NULL, thresholdDirection = c("low", "high"), border = NA, nLocalFeatures = 5)
{
	par(bg = "grey")
	if (is.list(object$localVariableImportance))
	{
		Z = object$localVariableImportance$obsVariableImportance
		whichClassNames <- rm.string(names(object$localVariableImportance$class), "Class ")
		numericClassNames = as.numeric(as.factor(whichClassNames))		
		if (is.null(whichClass) ) { stop("Please provide a class.") }
		if (is.character(whichClass)) { whichClass = numericClassNames[which(whichClassNames == whichClass)] 	}
		whichClass2 = whichClassNames[which(numericClassNames == whichClass)]	
		idx = which(Z[,1] == whichClass)
		Z = Z[idx, ,drop = FALSE]		
		if (dim(Z)[1] <= 1) { stop("not enough obervations found using this class.") }
	}
	else
	{   
		Z = object$localVariableImportance	
		if (!is.null(threshold))
		{ 
			if (thresholdDirection == "low") { 	idx = which(Z[,1] <= threshold) }
			else {  idx = which(Z[,1] > threshold)  }
			Z = Z[idx,]
			
			if (dim(Z)[1] < 1) { stop("no obervations found using this threshold.") }
		}
	}	
	Z = Z[,-1]
	idxLocalVar = grep("localVariable", colnames(X))
	idxPos = length(idxLocalVar)/2	
	countVars = sort(table(Z[,1:idxPos]), decreasing = TRUE)	
	obsObject = rmNA(countVars[1:nLocalFeatures]/sum(countVars))	
	XNames = colnames(X)
	par(las = 1)
	maxChar = floor(2 + max(nchar(XNames[as.numeric(names(sort(obsObject)))])))/2
	par(mar = c(5, maxChar + 2,4,2))
	barplot(sort(obsObject*100), horiz = TRUE, col = sort(heat.colors(length(obsObject)), decreasing = TRUE), border = border,
		names.arg = XNames[as.numeric(names(sort(obsObject)))], xlab = "Relative influence (%)", 
		main = if (is.list(object$localVariableImportance)) { paste("Partial importance based on observations over class ", whichClass2, sep ="") }
		else 
		{
			if (!is.null(threshold)) 
			{	
				if (thresholdDirection == "low") { paste("Partial importance based on observations (with Response < ", round(threshold, 4), ")", sep ="") }
				else  { paste("Partial importance based on observations (with Response > ", round(threshold, 4), ")", sep ="") }					
			}
			else {	"Partial importance based on observations" }
		}
	)	
	cat("Variance explained: ", round(rmNA(sum(obsObject))*100, 2), "%\n", sep="")
	return(obsObject)
}

partialDependenceOverResponses <- function(Xtest, importanceObject, whichFeature = NULL, whichOrder = c("first", "second", "all"), outliersFilter = FALSE, plotting = TRUE, followIdx = FALSE)
{
	FeatureValue = Response = Class = Observations = NULL
	if (!is.null(whichFeature))
	{	
		if (is.character(whichFeature)) { whichFeature = which(colnames(Xtest) == whichFeature) }		
		if (length(whichFeature) > 1)
		{ 
			whichFeature = whichFeature[1]
			cat("Only one variable can be computed at the same time\n")
		}
	}	
	if (whichOrder[1] == "first") 	{ idxOrder = 2	}	
	if (whichOrder[1] == "second") 	{ idxOrder = 3	}	
	if (whichOrder[1] == "all") 	
	{ 
		if (is.matrix(importanceObject$localVariableImportance))
		{	idxOrder = 2:length(grep("localVariableFrequency", colnames(importanceObject$localVariableImportance)))	}
		else
		{	idxOrder = 2:length(grep("localVariableFrequency", colnames(importanceObject$localVariableImportance$obsVariableImportance)))	}
	}	
	idx = list()
	if (is.matrix(importanceObject$localVariableImportance)) {	importanceObjectMatrix = importanceObject$localVariableImportance }
	else { importanceObjectMatrix = importanceObject$localVariableImportance$obsVariableImportance }	
	if (is.null(whichFeature))
	{	whichFeature = as.numeric(names(which.max(table(importanceObjectMatrix[,idxOrder[1]]))))	}			
	idx[[1]] = which(importanceObjectMatrix[,idxOrder[1]] == whichFeature)			
	if (length(idxOrder) > 1)
	{
		for (i in 1:length(idxOrder))
		{	idx[[i+1]] = which(importanceObjectMatrix[,1+idxOrder[i]]== whichFeature)	}
	}	
	partialDependenceMatrix = cbind(Xtest[unlist(idx), whichFeature], importanceObjectMatrix[unlist(idx),1], unlist(idx))
	partialDependenceMatrix = sortMatrix(partialDependenceMatrix, 1)	
	NAIdx = which(is.na(partialDependenceMatrix))	
	if (length(NAIdx) > 0) {  partialDependenceMatrix = partialDependenceMatrix[-NAIdx,] }	
	if (outliersFilter)
	{
		highOutlierIdx =  which(partialDependenceMatrix[,1] > quantile(partialDependenceMatrix[,1],0.95))
		lowOutlierIdx =  which(partialDependenceMatrix[,1] < quantile(partialDependenceMatrix[,1],0.05))
		if (length(highOutlierIdx) > 0 | length(lowOutlierIdx) > 0) 
		{	partialDependenceMatrix = partialDependenceMatrix[-c(lowOutlierIdx,highOutlierIdx),]	}
	}	
	idx = partialDependenceMatrix[,3]
	partialDependenceMatrix = partialDependenceMatrix[,-3]
	if (plotting)
	{
		if (dim(partialDependenceMatrix)[1] < 1)
		{ stop ("not enough points to plot partial dependencies. Please increase order of interaction") }		
		if (is.matrix(importanceObject$localVariableImportance))
		{
			#require(ggplot2) || install.packages("ggplot2")
			colnames(partialDependenceMatrix) = c("FeatureValue", "Response")
			partialDependenceMatrix = data.frame(partialDependenceMatrix)
			
			tt <- ggplot(partialDependenceMatrix, aes(x = FeatureValue, y = Response))
			plot(tt +  geom_point(colour = "lightblue") + stat_smooth(fill = "green", colour = "darkgreen", size = 1) + 
			labs(title = "Partial dependence over predictor", x = colnames(Xtest)[whichFeature], y = "Response"))
		}
		else
		{
			#require(ggplot2) || install.packages("ggplot2")
			colnames(partialDependenceMatrix) = c("Observations", "Class")
			partialDependenceMatrix = data.frame(partialDependenceMatrix)
			variablesNames = unique(partialDependenceMatrix$Class)
			partialDependenceMatrix$Class = factor(partialDependenceMatrix$Class)
			levels(partialDependenceMatrix$Class) = colnames(importanceObject$localVariableImportance$classVariableImportance)[sort(variablesNames)]			
			plot(qplot(Class, Observations, data = partialDependenceMatrix, geom = c("boxplot", "jitter"), outlier.colour = "green", outlier.size = 2.5, fill= Class, main = "Partial dependence over predictor", xlab = "", ylab = colnames(Xtest)[whichFeature]))
		}
	}
	else
	{
		if (!is.matrix(importanceObject$localVariableImportance))
		{	
			colnames(partialDependenceMatrix) = c("Observations", "Class")
			partialDependenceMatrix = data.frame(partialDependenceMatrix)
			variablesNames = unique(partialDependenceMatrix$Class)
			partialDependenceMatrix$Class = factor(partialDependenceMatrix$Class)
			levels(partialDependenceMatrix$Class) = colnames(importanceObject$localVariableImportance$classVariableImportance)[sort(variablesNames)]
		}
	}	
	if (followIdx)
	{	return(list(partialDependenceMatrix = partialDependenceMatrix, idx = as.numeric(idx) )) }
	else
	{  return(partialDependenceMatrix)  }
}	

partialDependenceBetweenPredictors <- function(Xtest, importanceObject, features, whichOrder = c("first", "second", "all"), perspective = FALSE, outliersFilter = FALSE)
{
	Variable1 = Variable2 = SameClass = ..level.. = Response = NULL
	if (length(features) != 2) { stop("Please provide two features.") }	
	#require(ggplot2)
	graphics.off()	
	pD1 <- partialDependenceOverResponses(Xtest, importanceObject, whichFeature = features[1], whichOrder = whichOrder, 
		outliersFilter = outliersFilter, plotting = FALSE, followIdx = TRUE)	
	pD2 <- partialDependenceOverResponses(Xtest, importanceObject, whichFeature = features[2], whichOrder = whichOrder, 
		outliersFilter = outliersFilter, plotting = FALSE, followIdx = TRUE)	
	sameIdx2 = find.idx(pD1$idx, pD2$idx)
	sameIdx1 = find.idx(pD2$idx, pD1$idx)
	minDim = length(sameIdx1)	
	if ( (minDim < 10) | (length(sameIdx2) < 10)) { stop("Not enough points. Please use option whichOrder = 'all'") }	
	pD1 = pD1$partialDependenceMatrix; pD2 = pD2$partialDependenceMatrix;
	pD11 = factor2matrix(pD1)[sameIdx1,]; pD22 = factor2matrix(pD2)[sameIdx2,]
	if (!is.matrix(importanceObject$localVariableImportance))
	{
		minN = min(nrow(pD11), nrow(pD22))
		pD11 = pD11[1:minN,]
		pD22 = pD22[1:minN,]		
		idx = ifelse(pD11[,2] == pD22[,2], 1, 0)
		Xi  = cbind(pD11[which(idx == 1), 1], pD22[which(idx == 1), 1])
		Xj  = cbind(pD11[which(idx == 0), 1], pD22[which(idx == 0), 1])		
		if (!is.character(features[1]))
		{ 
			fName1 = which(colnames(Xtest) == colnames(Xtest)[features[1]])
			fName2 = which(colnames(Xtest) == colnames(Xtest)[features[2]])
			features = colnames(Xtest)[c(fName1, fName2)]
		}		
		Xi = as.data.frame(Xi); colnames(Xi) = c("Variable1", "Variable2")
		Xj = as.data.frame(Xj); colnames(Xj) = c("Variable1", "Variable2")
		Xi = cbind(Xi, rep(1, nrow(Xi)))
		Xj = cbind(Xj, rep(0, nrow(Xj)))		
		colnames(Xj)[3] = colnames(Xi)[3] = "SameClass"
		X = rbind(Xi, Xj)
		X[,3] = ifelse(X[,3] == 1, TRUE, FALSE)		 
		dev.new() 
		tt <- ggplot(X, aes(x = Variable1, y = Variable2, colour = SameClass))	
		plot(tt + geom_point(size = 2) + labs(title = "Dependence between predictors", x = features[1], y = features[2]) 
			+  scale_colour_manual("Same class", values = c("red", "green") )
		)
		dev.new()
		cde1 <- geom_histogram(position = "fill", binwidth = diff(range(X[,1]))/4, alpha = 7/10)
		cde2 <- geom_histogram(position = "fill", binwidth = diff(range(X[,2]))/4, alpha = 7/10)
		tt1 <- ggplot(X, aes(x = Variable1, fill = SameClass))
		plot(tt1 + cde1 
			+ labs(title = "Class distribution", x = features[1], y = "Frequency")
			+ scale_fill_manual(paste("Same class as", features[2]),values = c("red", "lightgreen"))
		)		
		dev.new()
		tt1 <- ggplot(X, aes(x = Variable2, fill = SameClass))
		plot(tt1 + cde2 
			+ labs(title = "Class distribution", x = features[2], y = "Frequency")
			+ scale_fill_manual(paste("Same class as", features[1]),values = c("red", "lightgreen"))
		)				
		dev.new()
		tt2 <- ggplot(X, aes( x = Variable1, y = Variable2, z = SameClass))	
			try(plot(tt2 + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") 
			+ scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red")
			+ labs(title = "Heatmap of dependence between predictors", x = features[1], y = features[2])
		), silent = TRUE)		
		colnames(X) = c(features, "Same class")
	}
	else
	{
		intervals = cut(c(pD11[,2], pD22[,2]), minDim, labels = FALSE)
		pD11 = cbind(pD11,intervals[1:minDim])		
		if (nrow(pD22) != length(rmNA(intervals[(minDim + 1):(2*minDim)])))
		{
			sameN = min(nrow(pD22),length(rmNA(intervals[(minDim + 1):(2*minDim)])))
			pD22 = cbind(pD22[1:sameN,], rmNA(intervals[(minDim + 1):(2*minDim)])[1:sameN])
		}
		else
		{ 	pD22 = cbind(pD22, rmNA(intervals[(minDim + 1):(2*minDim)]))	}		
		minN = min(nrow(pD11), nrow(pD22))		
		Xi = sortMatrix(pD11,3)[1:minN,]
		Xj = sortMatrix(pD22,3)[1:minN,]
		Z  = (Xi[,2] + Xj[,2])/2		
		if (!is.character(features[1]))
		{ 
			fName1 = which(colnames(Xtest) == colnames(Xtest)[features[1]])
			fName2 = which(colnames(Xtest) == colnames(Xtest)[features[2]])
			features = colnames(Xtest)[c(fName1, fName2)]
		}		
		dev.new()
	    XiXj = as.data.frame(cbind(Xi[,1], Xj[,1], Z)); colnames(XiXj) = c("Variable1", "Variable2", "Response")		
		tt <- ggplot(XiXj, aes(x = Variable1, y = Variable2, z = Response))
		ttMore <- tt + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") +
			scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red") +
			labs(title = "Local heatmap of dependence (with frequency and intensity of Response)", x = features[1], y = features[2])
		try(plot(ttMore), silent = TRUE)
		dev.new()
		fourQuantilesCut = cut(Z,4)
		XiXj[,3] = fourQuantilesCut		
		dataCuts = table(fourQuantilesCut)		
		while (length(which(dataCuts < 5)) > 0)
		{
			lowDataCuts = names(dataCuts[which(dataCuts < 5)])
			rmIdx = find.idx(lowDataCuts, XiXj[,3])			
			XiXj = XiXj[,-3]
			Z = Z[-rmIdx]
			fourQuantilesCut = cut(Z, 4)
			XiXj = cbind(XiXj[-rmIdx,], fourQuantilesCut) 		
			dataCuts = table(fourQuantilesCut)
			colnames(XiXj)[3] = "Response"
		}
		X = cbind(XiXj, Z)
		colnames(X) = c(features, "Response in four quantile intervals", "Response")		
	    tt1 <- ggplot(XiXj, aes(x = Variable1, y = Variable2, z = Response)) 
		try(plot(tt1 + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") 
			+  scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red")	
			+ labs(title = "Global heatmap of dependence (with intensity of Response)", x = features[1], y = features[2]) 
		), silent = TRUE)
		dev.new()
		try(plot(tt + geom_point(aes(colour = Response, size= Response))
			+  stat_smooth(fill = "lightgrey", colour = "grey", size = 1) 
			+ labs(title = "Dependence between Predictors", x = features[1], y = features[2])		
			+ scale_colour_gradient2(low = "blue", mid = "green", high = "red") 	
		), silent = TRUE)		
		if (perspective)
		{
			gridSize = 150; gridLag = 100			
			dev.new()
			x = XiXj[,1]
			y = XiXj[,2]
			z = Z			
			n = length(x)			
			xyz = cbind(x,y,z)
			xyz.byX = sortMatrix(xyz,1)
			xyz.byY = sortMatrix(xyz,2)			
			newX = xyz.byX[,1]
			newY = xyz.byY[,2]			
			dummyForRepX = dummyForRepY = rep(0,n)
			for (i in 2:n)
			{ 
				if (newX[i] == newX[i-1])
				{	dummyForRepX[i] = 1 }
				
				if (newY[i] == newY[i-1])
				{	dummyForRepY[i] = 1 }
			}			
			newIdx = which(dummyForRepY == 0 & dummyForRepX == 0)
			newX = newX[newIdx]
			newY = newY[newIdx]			
			if ( (gridSize + gridLag) > length(newX))
			{
				interp.1 = seq(min(newX) + 0.1, max(newX) - 0.1,length = gridSize + gridLag - length(newX))
				interp.2 = seq(min(newY) + 0.1, max(newY) - 0.1,length = gridSize + gridLag- length(newY))				
				newX = sort(c(newX, interp.1))
				newY = sort(c(newY, interp.2))
			}			
			duplicatesX = duplicated(newX)
			duplicatesY = duplicated(newY)			
			if (sum(duplicatesX) > 0)
			{
				newX = newX[!duplicatesX]
				newY = newY[!duplicatesX]
			}			
			if (sum(duplicatesY) > 0)
			{
				newY = newY[!duplicatesY]
				newX = newX[!duplicatesY]
			}			
			newXYZ = cbind(newX, newY, rep(NA, length(newX)))
			proxyM = fillNA2.randomUniformForest(rbind(xyz, newXYZ))			
			xyz.dim =  dim(xyz)
			nn = length(newX)						
			newZ = matrix(NA, nn, nn)			
			for (i in 1:nn)
			{
				for (j in 1:nn)
				{	newZ[i,j] = mean(proxyM[which(newX[i] == proxyM[,1] | newY[j] == proxyM[,2]), 3])	}
			}											
			L.smoothNewZ = t(apply(newZ, 1, function(Z) lagFunction(Z, lag = gridLag, FUN  = mean, inRange = TRUE)))
			C.smoothNewZ = apply(newZ, 2, function(Z) lagFunction(Z, lag = gridLag, FUN  = mean, inRange = TRUE))			
			smoothIdx = round(seq(1, nrow(L.smoothNewZ), length = gridLag),0)
			newX = newX[-smoothIdx]
			newY = newY[-smoothIdx]
			C.smoothNewZ = C.smoothNewZ[,-smoothIdx]
			L.smoothNewZ = L.smoothNewZ[-smoothIdx,]
			newZ = 0.5*(C.smoothNewZ + L.smoothNewZ)			
			if (length(newX) > gridSize)
			{
				sampleIdx = sort(sample(nn, gridSize))
				newX = newX[sampleIdx]
				newY = newY[sampleIdx]
				newZ = newZ[sampleIdx, sampleIdx]
			}
			highOutlierIdx2 <- apply(newZ, 2, function(Z) which(Z > quantile(Z, 0.975)))
			lowOutlierIdx2 <- apply(newZ, 2, function(Z) which(Z < quantile(Z, 0.025)))
			ouliersIdx2 <- c(lowOutlierIdx2, highOutlierIdx2)
			if (length(ouliersIdx2) > 0) 
			{	
				newZ = newZ[-ouliersIdx2, -ouliersIdx2]  	
				newX = newX[-ouliersIdx2]
				newY = newY[-ouliersIdx2]
			}
			rm(lowOutlierIdx2);
			rm(highOutlierIdx2);
			highOutlierIdx2 <- apply(newZ, 1, function(Z) which(Z > quantile(Z, 0.975)))
			lowOutlierIdx2 <-  apply(newZ, 1, function(Z) which(Z < quantile(Z, 0.025)))
			ouliersIdx2 <- c(lowOutlierIdx2, highOutlierIdx2)
			if (length(ouliersIdx2) > 0) 
			{	
				newZ = newZ[-ouliersIdx2, -ouliersIdx2]  	
				newX = newX[-ouliersIdx2]
				newY = newY[-ouliersIdx2]
			}					
			par(bg = "grey")
			try(persp.withcol(newX, newY, newZ, heat.colors, nrow(newZ), theta = -40, phi = 15, xlab = features[1], ylab = features[2], zlab = "Response", main = "Dependence between Predictors and (smoothed) effect over Response", ticktype = "detailed", box = TRUE, expand = 0.5, shade = 0.15), silent = TRUE)
		}
	}	
	#if (.Platform$OS.type == "windows") { arrangeWindows("vertical", preserve = FALSE) }
	idxf1 = which(colnames(importanceObject$partialDependence) ==  features[1])
	idxf2 = which(rownames(importanceObject$partialDependence) ==  features[2])	
	cat("\nLevel of interactions between", features[1], "and" , features[2], "at first order:", 
	round(importanceObject$partialDependence[idxf2, idxf1],4), "\n", sep=" ")	
	idxf1 = which(rownames(importanceObject$partialDependence) == features[1])	
	if (length(idxf1) > 0)	
	{ 
		idxf2 = which(colnames(importanceObject$partialDependence) ==  features[2]) 
		if (length(idxf2) > 0)	
		{	
			cat("Level of interactions between", features[1], "and", features[2], "at second order:", 
			round(importanceObject$partialDependence[idxf1, idxf2],4), "\n", sep=" ")		
		}
	}
	cat("\nPlease use R menu to tile vertically windows and see all plots.\n")
	return(X)
}

twoColumnsImportance <- function(importanceObjectMatrix)
{
	idx = length(grep("localVariableFrequency", colnames(importanceObjectMatrix))) - 1
	tmpImportanceObjectMatrix = importanceObjectMatrix[,1:2]
	if (idx > 0)
	{
		for (i in 1:idx)
		{	tmpImportanceObjectMatrix = rbind(tmpImportanceObjectMatrix, importanceObjectMatrix[,c(1,2+idx[i])])	}
	}
	return(tmpImportanceObjectMatrix)
}

plot.importance <- function(x, nGlobalFeatures = 30, nLocalFeatures = 5, Xtest = NULL, whichFeature = NULL, whichOrder = "all", outliersFilter = TRUE, formulaInput = NULL, border = NA, ...)
{
	cat("\nPlease use R menu to tile vertically windows and see all plots.\n")
	object <- x
	Variable = Response = NULL
	if (!is.null(Xtest))
	{
		if (!is.null(formulaInput))
		{
			mf <- model.frame(formula = formulaInput, data = as.data.frame(Xtest))
			Xtest <- model.matrix(attr(mf, "terms"), data = mf)[,-1]	
		}
		else
		{		
			Xtest = NAfactor2matrix(Xtest)
			matchNA = (length(which(is.na(Xtest))) > 0)				
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values\n")
				Xtest <- na.impute(Xtest)
			}
		}
	}		
	maxVar = nGlobalFeatures
	maxVar2 = nLocalFeatures
	graphics.off()
	varImportance = object$globalVariableImportance
	n = nrow(varImportance)	
	if (n > maxVar) { varImportance = varImportance[1:maxVar,] }
	else { maxVar = n}	
	par(las = 1)
	maxChar = floor(2 + max(nchar(as.character(object$globalVariableImportance[,1])))/2)
	par(mar = c(5, maxChar + 1,4,2))
	barplot(varImportance[maxVar:1,"percent.importance"], horiz = TRUE, col = sort(heat.colors(maxVar), decreasing = TRUE), border = border,
	names.arg = varImportance[maxVar:1,"variables"], xlab = "Relative influence (%)", main = "Variable importance based on information gain")
	abline(v = 100/n, col = 'grey')	
	dev.new()
	par(las=1)
    par(mar = c(5,4,4,2))
	nbFeatures = ncol(object$partialDependence)
	newNbFeatures = min(maxVar2, nbFeatures -1)	
	if (newNbFeatures < (nbFeatures - 1))
	{
		OthersVariablesCol = colSums(object$partialDependence[(newNbFeatures+1):(nbFeatures -1), -nbFeatures, drop = FALSE])[1:newNbFeatures]
		OthersVariablesRow = rowSums(object$partialDependence[-nbFeatures, (newNbFeatures+1):(nbFeatures -1), drop = FALSE])[1:newNbFeatures]
		corner = mean(c(OthersVariablesCol, OthersVariablesRow))
		newPartialDependence = object$partialDependence[1:newNbFeatures,1:newNbFeatures]
		newPartialDependence =  rbind(cbind(newPartialDependence, OthersVariablesRow), c(OthersVariablesCol, corner))
		colnames(newPartialDependence)[ncol(newPartialDependence)] = "Other features"
		rownames(newPartialDependence)[nrow(newPartialDependence)] = "Other features"
		mosaicplot(newPartialDependence, color = sort(heat.colors(newNbFeatures + 1), decreasing = FALSE), 
		main = "Variables interactions over observations", ylab = "Most important variables at 2nd order", 
		xlab = "Most important variables at 1rst order", las = ifelse(maxChar > 10, 2,1), border = border)
	}
	else
	{
		mosaicplot(object$partialDependence[1:newNbFeatures,1:newNbFeatures], color = sort(heat.colors(newNbFeatures), decreasing = FALSE), 
		las = ifelse(maxChar > 10, 2, 1), main = "Variables interactions over observations", ylab = "Most important variables at 2nd order", xlab = "Most important variables at 1rst order", border = border)
	}	
	dev.new()
	par(las=1)
	par(mar=c(5,maxChar + 1,4,2))
	nbFeatures2 = min(maxVar2, length(object$variableImportanceOverInteractions))
	importanceOverInteractions = sort(object$variableImportanceOverInteractions, decreasing = TRUE)/sum(object$variableImportanceOverInteractions)*100
	barplot(importanceOverInteractions[nbFeatures2:1], horiz = TRUE, col = sort(heat.colors(nbFeatures2), decreasing = TRUE), 
		names.arg = names(importanceOverInteractions)[nbFeatures2:1], xlab = "Relative influence (%)", 
		main = "Variable importance based on interactions", border = border)
	abline(v = 100/n, col = 'grey')
	if (!is.matrix(object$localVariableImportance))
	{
		dev.new()
		par(las=1)
	    par(mar = c(5,maxChar + 1,4,2))
		nbFeatures3 = min(nbFeatures2, nrow(object$localVariableImportance$classVariableImportance))
		mosaicplot(t(object$localVariableImportance$classVariableImportance[1:nbFeatures3,]), color = sort(heat.colors(nbFeatures3), 
			decreasing = FALSE), main = "Variable importance over labels", border = border)
	}
	else
	{
		dev.new()
		#require(ggplot2) || install.packages("ggplot2")
		ggData = twoColumnsImportance(object$localVariableImportance)
		mostImportantFeatures = as.numeric(names(sort(table(ggData[,2]), decreasing = TRUE)))	
		if (length(unique(ggData[,2])) > maxVar2)
		{	
			mostImportantFeatures = mostImportantFeatures[1:maxVar2]	
			ggData = ggData[find.idx(mostImportantFeatures, ggData[,2], sorting = FALSE),]
		}
		if (is.null(Xtest)) {	textX = paste("V", mostImportantFeatures, collapse = " ", sep="")	}
		else {	textX = paste( colnames(Xtest)[mostImportantFeatures], collapse = ", ", sep="")	}
		textX = paste("Index of most important variables [", textX, "]")
		colnames(ggData)[1] = "Response"
		colnames(ggData)[2] = "Variable"
		ggData = data.frame(ggData)
		if (!is.null(Xtest)) { 	ggData[,"Variable"] = colnames(Xtest)[ggData[,"Variable"]] 	}
		ggData[,"Variable"] = as.factor(ggData[,"Variable"])
		if (nrow(ggData) > 1000) 
		{ 
			randomSample = sample(nrow(ggData), 1000) 
			gg <- qplot(Variable, Response, data = ggData[randomSample,], geom = c("boxplot", "jitter"), colour = Response, outlier.colour = "green", outlier.size = 1.5, fill = Variable)
		}
		else
		{	
			gg <- qplot(Variable, Response, data = ggData, geom = c("boxplot", "jitter"), colour = Response, outlier.colour = "green", outlier.size = 1.5, fill = Variable)	
		}		
		if (maxChar > 10)
		{
			plot(gg + labs(x ="", y = "Response", title = "Dependence on most important predictors") + 
			theme(axis.text.x = element_text(angle = 60, hjust = 1)))			
		}
		else
		{	
			plot(gg + labs(x ="", y = "Response", title = "Dependence on most important predictors"))
		}
	}	
	if (is.null(Xtest))	{	stop("partial dependence between response and predictor can not be computed without test data") }
	else
	{
		endCondition = 0
		dev.new()
		pD <- partialDependenceOverResponses(Xtest, object, whichFeature = whichFeature, whichOrder = whichOrder, outliersFilter = outliersFilter)
		#if (.Platform$OS.type == "windows") { arrangeWindows("vertical", preserve = FALSE) }
		idxMostimportant = rmNA(match(names(object$variableImportanceOverInteractions), colnames(Xtest)))[1:nbFeatures2]
		mostimportantFeatures = colnames(Xtest)[idxMostimportant]
		while (!endCondition)
		{
			ANSWER <- readline(cat("To get partial dependence of most important features\ngive a column number\namong", idxMostimportant,
			"\n(", mostimportantFeatures, ")\nPress escape to quit\n"))
			if ( is.numeric(as.numeric(ANSWER)) & !is.na(as.numeric(ANSWER)) )
			{  
				whichFeature = as.numeric(ANSWER)
				if (whichFeature %in% idxMostimportant)
				{	
					pD <- partialDependenceOverResponses(Xtest, object, whichFeature = whichFeature , whichOrder = whichOrder, 
					outliersFilter = outliersFilter)	
				}
				else
				{   stop("Please provide column index among most important. Partial Dependence can not be computed.") }				
			}
			else
			{	endCondition = 1	}
		}
	}	
}

print.importance <- function(x,...)
{
	object <- x
	minDim = min(10,length(object$variableImportanceOverInteractions))
	cat("\nGlobal Variable importance (", minDim, "most important based on information gain ) :\n")
	print(object$globalVariableImportance[1:minDim,])	
	cat("\n\nLocal Variable importance (", minDim , "most important based on interactions ) :\n")
	print(round(object$variableImportanceOverInteractions/sum(object$variableImportanceOverInteractions),2)[1:minDim])
	if (!is.matrix(object$localVariableImportance))
	{
		cat("\n\nVariable importance over labels(", minDim, "most important variables ) :\n")
		print(object$localVariableImportance$classVariableImportance[1:minDim,])
	}	
	cat("\n\nVariables interactions over all observations and trees at first and second order :\n")
	print(round(object$partialDependence[1:minDim,], 2))	
	if (!is.matrix(object$localVariableImportance))
	{ cat("\n\n see ...$localVariableImportance$obsVariableImportance to get variable importance, variable frequency\n and predicted response for each observation.") }
	else 
	{	
		cat("\n\n see ...$localVariableImportance to get variable importance, variable frequency\n and predicted response for each observation.")	
	}
	cat("\n\n call partialDependenceOverResponses() function with importance object computed and option 'Xtest' (which request test data)\n to get partial dependence over response for each variable among most important.\n")	
}

combineRUFObjects <- function(rUF1, rUF2)
{
	rUF1 <- filter.object(rUF1)
	rUF2 <- filter.object(rUF2)
	return(onlineCombineRUF(rUF1, rUF2))	
}

rankingTrainData <- function(trainData = NULL, trainLabels = NULL,  testData = NULL, testLabels = NULL, ntree = 100,  thresholdScore = 2/3, nTimes = 2, ...)
{	
	score = tmpScore = rep(0, nrow(testData))	
	classes = unique(trainLabels)
	majorityClass = modX(trainLabels)	
	i = 1;	rmIdx = NULL; idx = 1:nrow(testData)
	while  ( (i <= nTimes)  & (length(testData[,1]) >= 1) )
	{
		rUF <- randomUniformForestCore(trainData, trainLabels = as.factor(trainLabels), ntree = ntree, use.OOB = FALSE, rf.overSampling = -0.75, rf.targetClass = majorityClass, rf.subagging = 2/3)
		predictRUF <- randomUniformForestCore.predict(rUF, testData)	
		tmpScore = tmpScore + apply(predictRUF$all.votes, 1, function(Z)  length(which(Z != majorityClass)))
		score[idx] = score[idx] + tmpScore								
		rmIdx = which(tmpScore < (thresholdScore*ntree))		
		if (length(rmIdx) > 0)  {  testData = testData[-rmIdx,];  idx = idx[-rmIdx]; tmpScore = tmpScore[-rmIdx];   }
		i = i + 1
	}	
	return(score)
}
		
plotTreeCore <- function(treeStruct, rowNum = 1, height.increment = 1)
{
  if ( (treeStruct[rowNum, "status"] == -1) )
  {
    treeGraphStruct <- list()
    attr(treeGraphStruct, "members") <- 1
    attr(treeGraphStruct, "height") <- 0
    attr(treeGraphStruct, "label") <- if (treeStruct[rowNum,"prediction"] == 0) { "next node" } else 
		{ 	
			if (is.numeric(treeStruct[rowNum,"prediction"]))	
			{ round(treeStruct[rowNum,"prediction"],2) } else { treeStruct[rowNum,"prediction"]	}	
		}
	attr(treeGraphStruct, "leaf") <- TRUE
  }
  else
  {
	left <- plotTreeCore(treeStruct, treeStruct[rowNum, "left.daughter"], height.increment) #,  incDepth+1
    right <- plotTreeCore(treeStruct, treeStruct[rowNum, "right.daughter"], height.increment)#,  incDepth+1
    treeGraphStruct <- list(left,right)
    attr(treeGraphStruct, "members") <- attr(left, "members") + attr(right,"members")
    attr(treeGraphStruct,"height") <- max(attr(left, "height"),attr(right, "height")) + height.increment
    attr(treeGraphStruct, "leaf") <- FALSE	
	if (rowNum != 1)
	{     attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), " ?", sep ="")  }
	else
	{	  
		attr(treeGraphStruct, "edgetext") <- paste(".[no].  .", 
		treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), ".                                                   .[yes]", sep="") 
	}
  }
  class(treeGraphStruct) <- "dendrogram"
  return(treeGraphStruct)
}

plotTreeCore2 <- function(treeStruct, rowNum = 1, height.increment = 1,  maxDepth = 100 )
{
  if ((treeStruct[rowNum, "status"] == -1) | (rowNum > maxDepth))
  {
    treeGraphStruct <- list()
    attr(treeGraphStruct, "members") <- 1
    attr(treeGraphStruct, "height") <- 0
    attr(treeGraphStruct, "label") <- if (treeStruct[rowNum,"status"] == 1) { "next node" } else 
		{ 	
			if (is.numeric(treeStruct[rowNum,"prediction"]))	
			{ round(treeStruct[rowNum,"prediction"],2) } else { treeStruct[rowNum,"prediction"]	}	
		}
	attr(treeGraphStruct, "leaf") <- TRUE
  }
  else
  {
    left <- plotTreeCore2(treeStruct, treeStruct[rowNum, "left.daughter"], height.increment) #,  incDepth+1
    right <- plotTreeCore2(treeStruct, treeStruct[rowNum, "right.daughter"], height.increment)#,  incDepth+1
    treeGraphStruct <- list(left,right)
    attr(treeGraphStruct, "members") <- attr(left, "members") + attr(right,"members")
    attr(treeGraphStruct,"height") <- max(attr(left, "height"),attr(right, "height")) + height.increment
    attr(treeGraphStruct, "leaf") <- FALSE	
	if (rowNum != 1)
	{     attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), "?", sep ="")  }
	else
	{	  
		attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), "?", ".                                     .", sep="") 
	}
  }
  class(treeGraphStruct) <- "dendrogram"
  return(treeGraphStruct)
}

plotTree <- function(treeStruct, rowNum = 1, height.increment = 1, maxDepth = 100, fullTree = FALSE, xlim = NULL, ylim= NULL, center = TRUE)
{
	if (fullTree)
	{	drawTreeStruct <- plotTreeCore(treeStruct, rowNum = rowNum , height.increment = height.increment)	}
	else
	{    drawTreeStruct <- plotTreeCore2(treeStruct, rowNum = rowNum , height.increment = height.increment,  maxDepth = maxDepth)	}	
	nP <- list(col = 3:2, cex = c(2.0, 0.75), pch = 21:22, bg =  c("light blue", "pink"), lab.cex = 0.75, lab.col = "tomato")	
	if (is.null(xlim))
	{	
		if (is.null(ylim)) { ylim = c(0,8) }
		plot(drawTreeStruct, center = center, leaflab ='perpendicular', edgePar = list(t.cex = 2/3, p.col = NA, p.lty = 0, lty =c( 2,5), 
			col = c("purple", "red"), lwd = 1.5), nodePar = nP, ylab = "Tree depth", xlab = "Predictions", 
			xlim = c(0, min(30,floor(nrow(treeStruct)/2))), ylim = ylim)
	}
	else
	{
		if (is.null(ylim)) { ylim = c(0,8) }
		plot(drawTreeStruct, center = center, leaflab ='perpendicular', edgePar = list(t.cex = 2/3, p.col = NA, p.lty = 0, lty =c( 2,5), 
			col = c("purple", "red"), lwd = 1.5), nodePar = nP, ylab = "Tree depth", xlab = "Predictions", xlim = xlim, ylim = ylim )
	}
}

fillNA2.randomUniformForest <- function(X, Y = NULL, ntree = 100, mtry = 1, nodesize = 10, NAgrep = "", threads = "auto")
{	
    i = NULL 
	if (ntree > 100 & nodesize == 1)
	{ 
		cat("For nodesize = 1, ntree cannot currently set higher to 100. Resetting to default value.\n")
		ntree = 100
	}	
	n <- nrow(X)
	X <- fillVariablesNames(X)	
	if (!is.null(Y)) { trueY = Y }
	if (!is.matrix(X))
	{	
		trueX = X
		flag = TRUE
		X.factors <- which.is.factor(X)
		X <- NAfactor2matrix(X, toGrep = NAgrep)
	}
	else
	{
		flag = FALSE
		X.factors = rep(0,ncol(X))
	}	
	NAIdx = which(is.na(X), arr.ind = TRUE)	
	if (sum(NAIdx) == 0) { stop("No missing values in data.") }	
	processedFeatures = unique(NAIdx[,2])	
	nFeatures = length(processedFeatures)	
	idx <- lapply(processedFeatures, function(Z) NAIdx[which(NAIdx[,2] == Z),1]) 
	validIdx <- sapply(idx, function(Z) (length(Z) > nodesize))
	processedFeatures <- processedFeatures[which(validIdx == TRUE)]
	nFeatures = length(processedFeatures)
	invalidIdx <- which(validIdx == FALSE)
	if (length(invalidIdx) > 0)
	{	idx <- rm.InAList(idx, invalidIdx)	}		
	if (nFeatures == 0)
	{	
		cat("Not enough missing values. Rough imputation is done. Please lower 'nodesize' value to increase accuracy.\n")
		return(na.replace(trueX, fast = TRUE))	
	}
	else
	{
		if (!is.null(Y)){  X = cbind(X,Y)	}
		X <- fillVariablesNames(X)
		X <- na.replace(X, fast = TRUE)			
		{
			#require(parallel)	
			max_threads = min(detectCores(),4)		
			if (threads == "auto")
			{	
				if (max_threads == 2) { threads = min(max_threads, nFeatures) }
				else {	threads  = max(1, max_threads)   }
			}
			else
			{
				if (max_threads < threads) 
				{	cat("Warning : number of threads indicated by user was higher than logical threads in this computer.\n") }
			}			
			threads =  min(nFeatures, max_threads)
			#require(doParallel)			
			Cl <- makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)		
			chunkSize  <-  ceiling(nFeatures/getDoParWorkers())
			smpopts  <- list(chunkSize = chunkSize)
		}				
		export = c("randomUniformForest.default", "rUniformForest.big", "randomUniformForestCore.big", "randomUniformForestCore", 		"predict.randomUniformForest", "rUniformForestPredict", "uniformDecisionTree", "CheckSameValuesInAllAttributes", "CheckSameValuesInLabels", "fullNode", "genericNode", "leafNode", "filter.object", "filter.forest", 
		"randomUniformForestCore.predict", "onlineClassify", "overSampling", "predictDecisionTree", "options.filter", "majorityClass", "randomCombination", "randomWhichMax", "vector2matrix", "which.is.na", "which.is.factor", "factor2vector", "outputPerturbationSampling", "rmNA", "count.factor", "find.idx", "genericOutput", "fillVariablesNames", "is.wholenumber", "rm.tempdir", "setManyDatasets", "onlineCombineRUF", "mergeLists", "classifyMatrixCPP", "L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP", "entropyInformationGainCPP", "runifMatrixCPP")
		if (nrow(X) < 2001)
		{
			newX <- foreach( i= 1:nFeatures, .options.smp = smpopts, .inorder = FALSE, .combine = cbind, .multicombine = TRUE, 
			.export = export) %dopar%
			{
				if (X.factors[processedFeatures[i]] == 1)
				{
					rufObject <- randomUniformForest.default(X[-idx[[i]],-processedFeatures[i]], Y = as.factor(X[-idx[[i]], 
						processedFeatures[i]]), OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, 
							threads = 1)
							
					X[idx[[i]], processedFeatures[i]] <- as.numeric(as.vector(predict.randomUniformForest(rufObject, 
						X[idx[[i]], -processedFeatures[i]])))
				}
				else
				{
					rufObject <- randomUniformForest.default(X[-idx[[i]],-processedFeatures[i]], Y = X[-idx[[i]], processedFeatures[i]], 
						OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1)
						
					X[idx[[i]], processedFeatures[i]] <- predict.randomUniformForest(rufObject, X[idx[[i]], -processedFeatures[i]])	
				}							
				if (mean(is.wholenumber(rmNA(X[-idx[[i]], processedFeatures[i]]))) == 1) {  round(X[,processedFeatures[i]]) }
				else {  X[,processedFeatures[i]]  }	
			}			
			X[,processedFeatures] = newX
			stopCluster(Cl)				
		}
		else
		{
			newX <- foreach( i = 1:nFeatures, .options.smp = smpopts, .inorder = FALSE, .combine = cbind, .multicombine = TRUE, 
			.export = export) %dopar%
			{
				if (X.factors[processedFeatures[i]] == 1)
				{
					rufObject <- rUniformForest.big(X[-idx[[i]],-processedFeatures[i]], Y = as.factor(X[-idx[[i]], 
						processedFeatures[i]]), nforest = floor(nrow(X)/1000), replacement = TRUE, randomCut = TRUE, OOB = FALSE, 
						importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1)
							
					X[idx[[i]], processedFeatures[i]] <- as.numeric(as.vector(predict.randomUniformForest(rufObject, 
						X[idx[[i]], -processedFeatures[i]])))
				}
				else
				{
					rufObject <- rUniformForest.big(X[-idx[[i]],-processedFeatures[i]], Y = X[-idx[[i]], processedFeatures[i]], 
						nforest = floor(nrow(X)/1000), replacement = TRUE, randomCut = TRUE, OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1)
						
					X[idx[[i]], processedFeatures[i]] <- predict.randomUniformForest(rufObject, X[idx[[i]], -processedFeatures[i]])	
				}
							
				if (mean(is.wholenumber(rmNA(X[-idx[[i]], processedFeatures[i]]))) == 1) {  round(X[,processedFeatures[i]]) }
				else {  X[,processedFeatures[i]]  }	
			}
			X[,processedFeatures] = newX
			stopCluster(Cl)
		}			
		if (sum(X.factors) != 0)
		{
			factorsIdx = which(X.factors != 0)
			X = as.true.matrix(X)
			X = as.data.frame(X)
			for (j in 1:length(factorsIdx))
			{	
				k = factorsIdx[j]
				X[,k] = as.factor(X[,k])
				Xlevels = as.numeric(names(table(X[,k])))
				levels(X[,k]) = levels(as.factor(trueX[,k]))[1:length(Xlevels)]
			}				
		}
		return(X)
	}
}
# END OF FILE