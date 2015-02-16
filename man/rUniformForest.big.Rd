\name{rUniformForest.big}
\alias{rUniformForest.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random Uniform Forests for Classification and Regression with large data sets}

\description{
  Implements random uniform forests for data sets that are too large too fit in physical memory but enough too fit in virtual memory. data set is randomly (or not) cut in many sub-samples and each one is processed, getting many base (but ensemble) models per sub-sample. At the end, all base models are combined to obtain one ensemble of ensembles model. If data can not reside in physical memory, but can reside in virtual memory (physical memory + swap file) then consider R packages 'bigmemory', 'data.table' of 'ff' to load data. To save memory (and computing time), subsamples of data (e.g. using an Hadoop environment like) are suited, computing then one forest per sub-sample and combining all trees from all forests using \code{\link{rUniformForest.combine}}. Note that rUniformForest.big( ) is first designed to compute large files on a small computer, at the expense of accuracy. But, in case of a shifting distribution, model may be more robust than the standard one (at least for regression).
}

\usage{
rUniformForest.big(X, Y = NULL, 
	xtest = NULL,
	ytest = NULL,
	nforest = 2,
	randomCut = FALSE,
	reduceDimension = FALSE,
	reduceAll = FALSE,
	replacement = FALSE,
	subsample = FALSE,
	ntree = 100,
    nodesize = 1, 
	maxnodes = Inf,
	mtry = ifelse(bagging,ncol(X),floor(4/3*ncol(X))),
	regression = ifelse(is.factor(Y),FALSE, TRUE),
    subsamplerate = ifelse(regression, 0.7, 1),
	replace = ifelse(regression,FALSE,TRUE),
	OOB = TRUE,
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
    depth = Inf,
    depthcontrol = NULL,
    importance = TRUE,	
	bagging = FALSE,
	unsupervised = FALSE,
	unsupervisedMethod = c("uniform univariate sampling", 
	"uniform multivariate sampling", "with bootstrap"),
	classwt = NULL,	
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
    rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L1", "L2"),	
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	subset = NULL,
	usesubtrees = FALSE,
	threads = "auto",
	parallelpackage = "doParallel")
}

\arguments{
  \item{X}{
	    a large data frame, or matrix, of predictors describing the model to be fitted. 
	}
	\item{Y}{
	   a responses vector. If it is a factor, classification is assumed, otherwise regression is computed.
	}
	\item{xtest}{
	     an optional data frame or matrix (like X) containing predictors for the test set.
	}
	\item{ytest}{
	     optional responses for the test set.
	}
	  \item{nforest}{
	     number of forests to compute. Size of each subsample will be 'number of observations in the original data set / nforest'. if 'nforest' is too high, accuracy will be lost. If too low, computing time will be long. 
	}
	  \item{randomCut}{
	     should original data be cut randomly ?
	}
	  \item{reduceDimension }{
	     should dimension be reduced in original data (or subsamples, if nforest > 1) ? useful for speed, but can reduce dramatically accuracy.
	}
	  \item{reduceAll}{
	     should subsamples have lower dimension ? It is recommended to use it only to quickly get a base result.
	}
	 \item{replacement}{
		   should sample of data be done with replacement, e.g. sample n observations between n with replacement then divide them in 'nforest' subsamples. 
	}
	 \item{subsample}{
		   value of subsample rate (m/n), e.g. sample m observations (m < n) between n, then divide them in 'nforest' subsamples. Note than 'subsample' can be combine with 'replacement' option.
	}
	\item{ntree}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{mtry}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{nodesize}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{maxnodes}{
	 	  see \code{\link{randomUniformForest}}.
	}
	\item{depth}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{depthcontrol}{
	   see \code{\link{randomUniformForest}}.
	}
	\item{regression}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{replace}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{OOB}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{BreimanBounds}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{subsamplerate}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{importance}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{bagging}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{unsupervised}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{unsupervisedMethod}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{classwt}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{oversampling}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{targetclass}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{outputperturbationsampling}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{rebalancedsampling}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{featureselectionrule}{
	   see \code{\link{randomUniformForest}}.
	}
	\item{randomcombination}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{randomfeature}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{categoricalvariablesidx}{
	     see \code{\link{randomUniformForest}}.
	}
	\item{na.action}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{logX}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{classcutoff}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{subset}{
		see \code{\link{randomUniformForest}}.
	}
	\item{usesubtrees}{
		see \code{\link{randomUniformForest}}.
	}
	\item{threads}{
	    see \code{\link{randomUniformForest}}.
	}
	\item{parallelpackage}{
	     see \code{\link{randomUniformForest}}.
	}
}

\value{
  an object of class randomUniformForest.
}

\author{
  Saip CISS \email{saip.ciss@wanadoo.fr}
}

\seealso{
  \code{\link{rUniformForest.combine}}, \code{\link{rUniformForest.grow}}
}

\examples{
## not run
## Classification: synthetic data
# n = 100;  p = 10 ## for ease of use, we consider small 'n'
## Simulate 'p' gaussian vectors with random parameters between -10 and 10.
#X <- simulationData(n,p)

## Make a rule to create response vector
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]

# Y <- as.factor(ifelse(rule > mean(rule), 1, 0))

# big.ruf <- timer(rUniformForest.big(X, Y, nforest = 2, 
# threads = 1, BreimanBounds = FALSE, replacement = TRUE, importance = FALSE))

## elapsing time
# big.ruf$time 

## OOB accuracy
# big.ruf$object

## standard model
# std.ruf <- timer(randomUniformForest(X, Y, threads = 1, ntree = 20, BreimanBounds = FALSE))

## elapsing time. Note that for small 'n' standard case will be faster.
# std.ruf$time  

## OOB accuracy
#std.ruf$object

## not run
##  regression
# Y = rule
# big.ruf <- timer(rUniformForest.big(X, Y, nforest = 2, 
# threads = 2, BreimanBounds = FALSE, subsample = 0.7))
# big.ruf  

## classic random uniform forest 
# std.ruf <- timer(randomUniformForest(X, Y, threads = 2, BreimanBounds = FALSE))
# std.ruf  # accuracy gap is much larger in case of regression

## but, one can consider a new case, e.g. shifting distribution, to see how it works
# newX <- simulationData(n,p)
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# newRule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]
# newY = newRule 

## predict using standard model
# pred.std.ruf <- predict(std.ruf$object, newX)

## get mean squared error
# sum( (pred.std.ruf - newY)^2 )/length(newY)

## predict using rUniformForest.big
# pred.big.ruf <- predict(big.ruf$object, newX)

## get mean squared error : both errors will be more closer, and for large 'n' (and more trees), 
## rUniformForest.big might have lower error
# sum( (pred.big.ruf - newY)^2 )/length(newY)
}