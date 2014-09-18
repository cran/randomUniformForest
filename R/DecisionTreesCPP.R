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

CheckSameValuesInLabels <- function(Y, n, m, nodeMinSize)
{
	flagZeros = TRUE
	if (n > nodeMinSize)
	{	if (sum(Y) !=  n*Y[1])	{  flagZeros = FALSE }	}	
	return(flagZeros)
}

CheckSameValuesInAllAttributes <- function(X, n, m, nodeMinSize)
{
	flagZeros = TRUE
	if (n > nodeMinSize)
	{ 	
		i = 1
		Cte <- sum(X[i,])/m		
		while(flagZeros & (i < n))
		{
			if (Cte != X[i,1])
			{  flagZeros = FALSE  }
			else
			{
				if (Cte != X[i+1,1])
				{	flagZeros = FALSE	}
			}
			i = i + 1
		}			
	}	
	return(flagZeros)
}

options.filter <- function(X, Y, nodeMinSize, o.bootstrap = FALSE, o.subagging = FALSE, o.treeOverSampling = 0, o.targetClass = -1, 
	o.OOB = FALSE, o.treeRebalancedSampling = FALSE) 
{
	dimX = dim(X)
	n = dimX[1]
	p = dimX[2]
	follow.idx = 1:n	
	if (!is.numeric(o.treeRebalancedSampling))
	{
		if (o.treeRebalancedSampling)
		{
			classesObject = table(Y)
			minLengthClasses =  which.min(classesObject)
			classes = as.numeric(names(classesObject))
			sampleSize= classesObject[minLengthClasses]
			classesSample = NULL			
			for (i in 1:length(classes))
			{ 	classesSample = c(classesSample, sample(which(Y == classes[i]), sampleSize)) 	}			
			follow.idx = sort(classesSample)
			n = length(follow.idx)
		}
	}
	else
	{
		if (max(o.treeRebalancedSampling) < 1)
		{	o.treeRebalancedSampling = round(o.treeRebalancedSampling*n,0) }
		if (o.subagging != 1) { stop("rebalanced sampling is not currently compatible with sub-sampling")  }
		classesObject = table(Y)
		minLengthClasses =  which.min(classesObject)
		classes = as.numeric(names(classesObject))
		classesSample = NULL		
		for (i in 1:length(classes))
		{ 	
			classesSampleSize = which(Y == classes[i])
			classesSample = c(classesSample, sample(classesSampleSize, o.treeRebalancedSampling[i], 
			replace = if (classesSampleSize >  o.treeRebalancedSampling[i]) { FALSE } else { TRUE })) 	
		}	
		follow.idx = sort(classesSample)
		n = length(follow.idx)
	}	
	if (o.treeOverSampling != 0)
	{ 
		if (o.targetClass == -1) { stop ("no class to target") }			
		oversample <- overSampling(Y[follow.idx], which.class = o.targetClass, proportion = o.treeOverSampling)	
		follow.idx = sort(c(oversample$C1, oversample$C2))
		n = length(follow.idx)
		
		if (o.treeOverSampling == -1) { o.subagging = 1 }
	}	
	OOBIdx = 0
	if (o.subagging != 1)
	{	
		follow.idx = sortCPP(follow.idx)
		nn = sort(sample(follow.idx, floor(o.subagging*length(follow.idx)), replace = o.bootstrap))
		if (o.OOB)	{  OOBIdx = follow.idx[-nn]	}
		follow.idx = rmNA(follow.idx[nn])		
		if (o.OOB)
		{
			idx = vector()
			for (i in 1:length(follow.idx))
			{	idx = c(idx, which(follow.idx[i] == OOBIdx)) }		
			if (length(idx) > 0) { OOBIdx = OOBIdx[-idx]	}	
		}
		n = length(follow.idx)
	}	
	if ( (o.bootstrap) & (o.subagging == 1) )
	{  
	    bootstrap.idx = sortCPP(sample(n, n, replace = TRUE))				
		follow2.idx = follow.idx[bootstrap.idx]
		
		while(CheckSameValuesInAllAttributes(X[follow2.idx,], n, p, nodeMinSize) | CheckSameValuesInLabels(Y[follow2.idx], n, p, nodeMinSize) )
		{	bootstrap.idx = sample(n, n, replace = TRUE); follow2.idx = follow.idx[bootstrap.idx]	}
		if (o.OOB)
		{ 	 
			OOBIdx = follow.idx[1:n][-bootstrap.idx]
			rm.idx = find.idx(follow2.idx, OOBIdx, sorting =FALSE)
			if (length(rm.idx) > 0)  { OOBIdx = OOBIdx[-rm.idx] }			
		}		
		follow.idx = follow2.idx
	}	
	return(list(X = X[follow.idx,], Y = Y[follow.idx], n = n, p = p, OOBIdx = OOBIdx))
}

uniformDecisionTree <- function(X, Y, nodeMinSize = 1, maxNodes = Inf, rf.features = floor(4/3*ncol(X)), 
treeDepth = Inf, 
treeDepthControl = NULL, 
getSplitAt = "random", 
featureSelectionRule = c("entropy", "gini", "random", "L2", "L1"),
regression = FALSE, 
bootstrap = FALSE, 
subagging = 1, 
treeClasswt = NULL, 
treeOverSampling = 0,
targetClass = -1,  
treeRebalancedSampling = FALSE, 
OOB = FALSE, 
x.bagging = FALSE, 
random.combination = 0,
randomFeature = FALSE, 
treeCatVariables = NULL,
outputPerturbation = FALSE) 
{
	{
		set.seed(sample(100000,1))
		gainFunction = featureSelectionRule[1] 
		randomSubspace = FALSE
		if (is.character(rf.features))
		{ 
			rf.features = sample(2:(2*ncol(X)), 1) 
			randomSubspace = TRUE
		}		
		randomNodesize = FALSE
		if (is.character(nodeMinSize))
		{ 
			nodeMinSize = sample(max(floor(nrow(X)/2),1),1) 
			randomNodesize = TRUE
		}		
		if (is.character(maxNodes))	{ maxNodes = sample(3:floor(nrow(X)),1) }				
		if (rf.features == 1) { randomFeature = TRUE }				
		if (outputPerturbation) 
		{  
			minValue = if (treeOverSampling == 0) { 0.05 } else { 1/(length(Y)-1) }
			Y <- outputPerturbationSampling(Y, whichClass = targetClass, sampleSize = max(min(treeOverSampling, 1), minValue), regression = regression)
			treeOverSampling = 0 
		}
		object.filtered <- options.filter(X, Y, nodeMinSize, o.bootstrap = bootstrap, o.subagging = subagging, o.treeOverSampling = treeOverSampling, o.targetClass = targetClass, o.OOB = OOB, o.treeRebalancedSampling = treeRebalancedSampling)			
		init.Y = rep(0, nrow(X))
		init.X = matrix(0, nrow(X), ncol(X)) 
		init.Y = object.filtered$Y
		init.X = object.filtered$X
		n = object.filtered$n
		p = object.filtered$p
		OOBIdx = object.filtered$OOBIdx
		if (nodeMinSize >= n) { stop("minimal number of observations (nodesize) can not be greater than sample size") }
		m = rf.features
		if (x.bagging) { m = min(p, rf.features) }		
		mCoeff = m/p
		iGain = iGain.tmp = rep(0, rf.features)
		minNodes = 3
		if (!regression) 
		{ 
			if (gainFunction == "random") {	gainFunction = sample(  c("entropy", "gini"), 1) }
			classes = sortCPP(unique(init.Y)) 
			nClasses = length(classes) 
			min_class = min(classes)			
			if (is.character(treeDepthControl)) {  thresholdLimit <- runif(1,0,0.01) }
			else  
			{  
				if (is.null(treeDepthControl)) 	{ 	thresholdLimit = 0  }
				else { thresholdLimit = treeDepthControl }
			}	
		}
		else
		{  
			if (gainFunction == "random") {	gainFunction = sample( c("L1", "L2"), 1) }
			classes = init.Y 
			nClasses = 1		   
			if (!is.null(treeDepthControl))
			{	
				if (gainFunction == "L2") 
				{ 
					L2Limit <- L2DistCPP(sum(init.Y)/n, init.Y)
					if (is.character(treeDepthControl)) {  limitCoeff = sample(4:16, 1) }
					else  {  limitCoeff = treeDepthControl	}
				} 
				else 
				{ 					
					L2Limit <- L1DistCPP(sum(init.Y)/n, init.Y) 
					if (is.character(treeDepthControl)) {  limitCoeff = sample(2:4, 1) }
					else  {  limitCoeff = treeDepthControl	}
				}				
			}
			else
			{ thresholdLimit = 0 }
		}
		nrowTree = floor(log(n)/log(2))
		nodes.nb = averageLeafWeight = vector(length = nrowTree)
		Tree = new.Tree = matrix(data = Inf, ncol = 7, nrow = nrowTree)
		rmVarList = list()
		dataSpaceList = vector("list", floor(log(n)/log(2)))
		dataSpaceList = lapply(dataSpaceList, function(Z) 1:n)
		dataSpaceList[[1]] = (1:n)
		rmVarList[[1]] = updated.nodes = 0L	
		endCondition = k = idxBestFeature = allNodes = 1L 
		proba.Y = vector(length = nClasses)		
		n.dup = n
		prev.iGain = vector()
		rmVar = NULL
	}
	rm(object.filtered)
	allDim = 1:p
	flagOptimization_1 = 0
	if (x.bagging)  { flagOptimization_1 = 1 }	
	if (randomFeature | outputPerturbation | x.bagging) { lowCorr <- FALSE }
	else
	{
	   lowCorr <- if (sample(0:1,1) == 1) { TRUE } else { FALSE }
	}
	while (endCondition)
	{		
		set.seed(sample(100000,1))
		{
			n = length(dataSpaceList[[k]])
			nodeMinSizeCandidate = max(floor(n/2),1)
			if ( randomNodesize & (k > nodeMinSizeCandidate) ) { nodeMinSize <- sample(nodeMinSizeCandidate,1)	}
			Y <- init.Y[dataSpaceList[[k]]] 
			X <- init.X[dataSpaceList[[k]], allDim, drop = FALSE] 
		}		
		pureNode = FALSE
		flag = 0		
		if ( (k >= minNodes) & (CheckSameValuesInAllAttributes(X, n, p, nodeMinSize) | (k >= (2^treeDepth - 1)) |  (k >= maxNodes)) )
		{	
			if (n != 0)
			{  
				Leaf <- leafNode(X, Y, k, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)				
				if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
				else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }
				nodes.nb[k] = n
			}
			else
			{	
				previous.k <- which(Tree == k, arr.ind = TRUE)[1] 
				Y <- init.Y[dataSpaceList[[previous.k]]]
				X <- matrix(init.X[dataSpaceList[[previous.k]], , drop = FALSE], length(previous.k), p)				
				Leaf <- leafNode(X, Y, k, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)				
				if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
				else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }
				nodes.nb[k] = length(Y)
			}
			updated.nodes = updated.nodes + 2			
			if (!is.null(treeClasswt) & (!regression)) 	{  averageLeafWeight[k] <- sum(treeClasswt*Leaf$nb)/nodes.nb[k]	}
		}
		else
		{
			if ((k >= minNodes) & CheckSameValuesInLabels(Y, n, p, nodeMinSize) ) 
			{ 
				Tree[k,] <- c( genericNode(k, -1), leafNode(X, Y, k, n, nClasses, classes, l.classwt = treeClasswt, 
					which.values = "output", l.regression = regression) )				
				nodes.nb[k] = n 
				pureNode = TRUE				
				if (!is.null(treeClasswt) & (!regression)) 	{  	averageLeafWeight[k] = treeClasswt[Y[1]]	}
			}
			else
			{
				if ( rmVarList[[k]][1] != 0 )	
				{	
					featuresIdx = allDim[-rmVarList[[k]]]  
					m = max(1, floor(mCoeff*length(featuresIdx)))  
				}
				else	
				{ 	
					featuresIdx = allDim
					m = rf.features	
				}
				candidateVariables <- sample(featuresIdx, m, replace = !x.bagging)
				if (!flagOptimization_1 & !randomSubspace & lowCorr)
				{
					
					randomLimit <- if (regression)  { sample(3:sample(3:7,1),1) }  else { 3 }								
					if ( (k <= randomLimit) & (m >= 2) )
					{   U.coeff = runif(1); candidateVariables = sample(featuresIdx, max(2,floor(U.coeff*2*p)), replace = TRUE)	}
					else
					{ 	flagOptimization_1 = 1 }					
				}				
				if (randomSubspace)
				{	
					U.coeff <- 2*runif(1)
					candidateVariables <- sample(featuresIdx, max(2,floor(U.coeff*p)), replace = TRUE)	
				}
				candidateVariables <- sortCPP(candidateVariables)
				colXTmp <- length(candidateVariables)
				catVar = 0 
				if (!is.null(treeCatVariables))
				{
					catVar = NULL 
					for (i in seq_along(treeCatVariables)) 
					{ 
						matchCatVar <- which(candidateVariables == treeCatVariables[i])
						if (length(matchCatVar) > 0) {  catVar = c(catVar, matchCatVar)  }
					}					
					if (!is.null(catVar)) 
					{ 
						catVar <- sortCPP(unique(candidateVariables[catVar]))												
						for (j in 1:length(catVar))
						{
							if (length(unique(X[,catVar[j]])) > 1)
							{
								dummy <- sample(unique(X[,catVar[j]]), 2)
								dummyIdx = which(X[,catVar[j]] == dummy[1])
								X[dummyIdx, catVar[j]] = dummy[1]
								X[-dummyIdx, catVar[j]] = dummy[2]
							}
						}
					}
					else
					{ catVar = 0 }
				}
				XTmp <- X[ ,candidateVariables, drop = FALSE]
				variablesLength <- checkUniqueObsCPP(XTmp)
				idxVariableLength <- which(variablesLength > 1)
				uniqueObs <- which(variablesLength == 1)
				lengthUniqueObs <- length(uniqueObs)
				if (!regression)
				{
					if (nClasses == 2) 
					{  
						proba.Y[1] = sum(Y == classes[1])
						proba.Y[2] = n - proba.Y[1]  
					}
					else  
					{	
						for (i in 1:nClasses) 
						{  proba.Y[i] = sum(Y == classes[i]) }	
					}						
					if (is.null(treeClasswt))
					{
						if (gainFunction == "entropy") { YGain <- crossEntropyCPP(proba.Y/n) } 
						else {  YGain <- giniCPP(proba.Y/n) }
					}
					else
					{
						if (gainFunction == "entropy") { YGain <- crossEntropyCPP(treeClasswt*proba.Y/n) } 
						else {  YGain <- giniCPP(treeClasswt*proba.Y/n) }
					}
				}
				else
				{  YGain =  sample(c(0,0.01),1);  classes = 0 }  								
				if (lengthUniqueObs > 0) { 	rmVar <- sortCPP(unique(candidateVariables[uniqueObs])) }				 				
				nCurrentFeatures <- length(variablesLength)								
				if (!is.na(idxVariableLength[1])) 				
				{		
					candidateVariables = candidateVariables[idxVariableLength]
					nIdxVariableLength <- length(idxVariableLength)
					if (nIdxVariableLength == 1)
					{
						XTmp <- XTmp[ ,idxVariableLength, drop = FALSE]
						nCurrentFeatures = 1
						thresholds <- runif(1, min(XTmp), max(XTmp))								
						if (!randomFeature)							
						{ 
							if (regression) 
							{ 
								if (gainFunction == "L2") {	iGain <- YGain - L2InformationGainCPP(Y, XTmp, thresholds) }
								else { iGain <- YGain - L1InformationGainCPP(Y, XTmp, thresholds) }
							}
							else  
							{  
								if ( is.null(treeClasswt) ) { newTreeClasswt <- rep(0, nClasses) }
								else  {  newTreeClasswt = treeClasswt }								
								flagIGFunction = 1
								if (gainFunction != "entropy") { flagIGFunction = -1 }								
								iGain <- YGain - entropyInformationGainCPP(Y, XTmp, thresholds, classes, nClasses, newTreeClasswt, flagIGFunction)
							}
						}										
					}					
					else 
					{
						XTmp <- XTmp[,idxVariableLength, drop = FALSE]						
						p.XTmp = nCurrentFeatures = nIdxVariableLength																		
						thresholds <- runifMatrixCPP(XTmp)
						if (!randomFeature)
						{
							if (regression)  
							{ 	
								if (gainFunction == "L2") {	iGain <- YGain - L2InformationGainCPP(Y, XTmp, thresholds) }
								else { iGain <- YGain - L1InformationGainCPP(Y, XTmp, thresholds) }								
							}
							else 
							{	
								if (is.null(treeClasswt)) 	{ newTreeClasswt = rep(0, nClasses) }
								else  {  newTreeClasswt = treeClasswt }
								
								flagIGFunction = 1
								if (gainFunction != "entropy") { flagIGFunction = -1 }
								
								iGain <- YGain - entropyInformationGainCPP(Y, XTmp, thresholds, classes, nClasses, newTreeClasswt, flagIGFunction)
							}
						}
					}
				}
				na.gain = which(is.na(iGain) | (iGain == Inf) | (iGain == -Inf))
				na.gain.length = length(na.gain)											
				if (!randomFeature)
				{
					if (na.gain.length == nCurrentFeatures) {  iGain = rep(0, nCurrentFeatures) }
					else {	if (na.gain.length > 0) {  iGain[na.gain] = 0 }  }
				}				
				if (randomFeature & !is.na(idxVariableLength[1]))
				{
				    idxBestFeature <-  if (length(thresholds) > 1) { sample(seq_along(thresholds), 1) } else  { 1 }
					if (regression)
					{	
						if (gainFunction == "L2") 
						{ 	
							iGain[idxBestFeature] <- YGain - 
							L2InformationGainCPP(Y, XTmp[,idxBestFeature, drop = FALSE], thresholds[idxBestFeature])
						}
						else
						{
							iGain[idxBestFeature] <- YGain - 
							L1InformationGainCPP(Y, XTmp[,idxBestFeature, drop = FALSE], thresholds[idxBestFeature])
						}
					}
					else  
					{  
						if (is.null(treeClasswt)) { newTreeClasswt = rep(0, nClasses) }
						else  { newTreeClasswt = treeClasswt }								
						flagIGFunction = 1
						if (gainFunction != "entropy") { flagIGFunction = -1 }								
						iGain[idxBestFeature] <- YGain - entropyInformationGainCPP(Y, XTmp[, idxBestFeature, drop = FALSE], thresholds[idxBestFeature], classes, nClasses, newTreeClasswt, flagIGFunction)
					}
				}
				else
				{	
					if (randomFeature) { idxBestFeature <- sample(candidateVariables,1) }
					else { idxBestFeature <- randomWhichMax(iGain) }
				}					
				iGain <- abs(iGain)
				if (regression & !is.null(treeDepthControl)) {  thresholdLimit = L2Limit/(limitCoeff*2*k)  }								
				if  (max(iGain) <= thresholdLimit)
				{	
					Leaf <- leafNode(X, Y, k, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)					
					if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
					else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }					
					if (!is.null(treeClasswt) & (!regression)) 	{  averageLeafWeight[k] = sum(treeClasswt*Leaf$nb)/n }					
					nodes.nb[k] = n
					pureNode = TRUE
				}
			}
			if (pureNode)
			{	updated.nodes = updated.nodes + 2; iGain = iGain.tmp	}
			else
			{
				belowK = 2*k - updated.nodes
				aboveK = 2*k + 1 - updated.nodes				
				if (lengthUniqueObs > 0) {  rmVarList[[belowK]] = rmVarList[[aboveK]] = rmVar }
				else { 	rmVarList[[belowK]] = rmVarList[[aboveK]] = 0	}								
				bestFeature = candidateVariables[idxBestFeature]
				bestThreshold = thresholds[idxBestFeature]
				XObsOfBestFeature <- X[,bestFeature]				
				subIdx = idxLow = which(XObsOfBestFeature <= bestThreshold)					
				nLowIdx = length(idxLow)  
				nHighIdx = n - nLowIdx				
				Tree[k,] <- fullNode(k, belowK, aboveK, bestFeature, bestThreshold)				
				nodes.nb[k] = n																	
				if ( (nLowIdx == n) | (nLowIdx == 0) )
				{ dataSpaceList[[belowK]] = 0L}
				else 
				{						
					dataSpaceList[[belowK]] = dataSpaceList[[k]][subIdx]
					flag = 1					
					if ( sum(XObsOfBestFeature[idxLow]) == (nLowIdx*XObsOfBestFeature[idxLow][1]) )  
					{  rmVarList[[belowK]] = bestFeature }
					allNodes = allNodes + 1L
				}							
				if ( (nHighIdx == 0) | (nHighIdx == n) )  
				{ dataSpaceList[[aboveK]] = 0L; allNodes = allNodes + 1L }
				else 
				{	
					dataSpaceList[[aboveK]] = (dataSpaceList[[k]][-subIdx])
					if ( sum(XObsOfBestFeature[-idxLow]) == (nHighIdx*XObsOfBestFeature[-idxLow][1])) 
					{  rmVarList[[aboveK]] = bestFeature }
					allNodes = allNodes + 1L
				}
			}
			prev.iGain[k] = iGain[idxBestFeature] 
			iGain = iGain.tmp
		}
		dataSpaceList[[k]] = 0L		
		k = k + 1L
		if (k > allNodes) {  endCondition = 0	}
		if ( (k < minNodes) & (endCondition == 0) )
		{ 
			set.seed(sample(100000,1))
			iGain = iGain.tmp = rep(0, rf.features)
			prev.iGain = vector()
			proba.Y = vector(length = nClasses)
			k = endCondition = idxBestFeature = allNodes = 1L			
			# dataSpaceList = rmVarList = list()
			rmVarList = list()
			dataSpaceList = vector("list", floor(log(n.dup)/log(2)))
			dataSpaceList = lapply(dataSpaceList, function(Z) 1:n.dup)			
			rmVarList[[1]] = updated.nodes = flagOptimization_1 = 00
			dataSpaceList[[1]] = 1:n.dup			
			Tree = new.Tree = matrix(data = Inf, ncol = 7, nrow = nrowTree)
			nodes.nb = averageLeafWeight = vector(length = nrowTree)			
			rmVar = NULL
			allDim = 1:p
			if (x.bagging)  { flagOptimization_1 = 1 }	
			if (!regression) 
			{ 
				if (is.character(treeDepthControl)) {  thresholdLimit <- runif(1,0,0.01) }
				else  
				{  
					if (is.null(treeDepthControl)) 	{ 	thresholdLimit = 0  }
					else { thresholdLimit = treeDepthControl  }
				}	
			}
			else
			{
				if (!is.null(treeDepthControl))
				{	
					if (gainFunction == "L2") 
					{ 
						L2Limit <- L2DistCPP(sum(init.Y)/n, init.Y)
						if (is.character(treeDepthControl)) {  limitCoeff = sample(4:16, 1) }
						else  {  limitCoeff = treeDepthControl	}
					} 
					else 
					{ 					
						L2Limit <- L1DistCPP(sum(init.Y)/n, init.Y) 
						if (is.character(treeDepthControl)) {  limitCoeff = sample(2:4, 1) }
						else  {  limitCoeff = treeDepthControl	}
					}				
				}
				else
				{ thresholdLimit = 0 }
			}					
		}		
		if ( (dim(Tree)[1] - k) < 2) 
		{ 
			Tree = rbind(Tree, new.Tree)
			new.nodes.nb = vector(length = nrowTree); nodes.nb = c(nodes.nb,new.nodes.nb)
		}
	}
	Tree.size = which((Tree == Inf), arr.ind = TRUE)[1] - 1
	Tree = Tree[1:Tree.size,]		
	if (Tree.size == 1) {  Tree = matrix(Tree, nrow = 1) }
	if (is.null(treeClasswt)) 
	{ 
		Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size])
		rownames(Tree) = Tree[,1]
		Tree = Tree[,-1]
		if (regression) 
		{ colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "L2Dist") }
		else 
		{	colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "Gain")	}
	}
	else
	{
		if (regression) { Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size]) }
		else { 	Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size], averageLeafWeight[1:Tree.size]) }
		
		rownames(Tree) = Tree[,1]
		Tree = Tree[,-1]
		if (regression) 
		{ colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "L2Dist") }
		else 
		{	
			colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "Gain", "avgLeafWeight" ) 
		}
	}	
	if (OOB) {	return(list(Tree = Tree, OOB.idx = OOBIdx))	}
	else {	return(list(Tree = Tree)) }
}

genericNode <- function(k, status) c(k, 0, 0, 0, 0, status)
fullNode <- function(k, left.node, right.node, attribute, split.point)  c(k, left.node, right.node, attribute, split.point, 1, 0) 
leafNode <- function(X, Y, k, n, nClasses, classes, which.values = "input", l.regression = FALSE, l.classwt = NULL)
{
	if (which.values == "input")
	{
		if (l.regression) 
		{ 
			if (n == 1) { return(Y) }
			else  
			{ 	return(sum(Y)/n) }		
		}
		else
		{ 	
			classNb = rep(0, nClasses)
			classNb[1] = sum(Y == classes[1])
			if (nClasses == 2) { classNb[2] =  n - classNb[1] }
			else
			{
				for (i in 2:nClasses) 
				{ classNb[i] = sum(Y == classes[i]) }
			}		
			if (!is.null(l.classwt)) {	classProb = l.classwt*classNb	}
			else  {  classProb = classNb  }	
			return (list(Class = randomWhichMax(classProb), nb = classNb))
		}			
	}
	else
	{ return(Y[1]) }
}

onlineClassify <- function(treeObject, X, c.regression = 0)
{
	p = length(X)
	YonlineLearning = X[p]; X = X[-p]
	classifyObject <- classifyMatrixCPP(treeObject, matrix(X, 1, p))
	nodes.idx = 1
	if (c.regression != 0)  
	{ 	
		if ( (YonlineLearning - classifyObject[1])^2  > c.regression )  {  nodes.idx = 0  } 	
	}
	else 
	{ 	
		if (YonlineLearning == classifyObject[1]) 	{	nodes.idx = 0 	} 	
	}
	return(nodes.idx)
}

predictDecisionTree <- function(treeObject, X, p.onlineLearning = FALSE, Y.p.onlineLearning = 0, p.regression = 0) 
{
	n = nrow(X)
	if (!is.matrix(treeObject)) { 	treeObject = treeObject$Tree	}
	if (p.onlineLearning) 
	{  
		X = cbind(X, Y.p.onlineLearning)
		predictedObject = vector(length = n) 
		for (i in 1:n) 
		{  predictedObject[i] <- onlineClassify(treeObject, X[i,], c.regression = p.regression)	}
	}
	else 
	{  
		predictedObject = matrix(data = 0, ncol = 3, nrow = n)
		colnames(predictedObject) = c("value", "nodes.data", "depth") 
		predictedObject <- classifyMatrixCPP(treeObject, X)	
	}	
	return(predictedObject)
}
# END OF FILE