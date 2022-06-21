\name{postProcessingVotes}
\alias{postProcessingVotes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-processing for Regression}

\description{
  Post-processing use OOB votes and predicted values to build more accurate estimates of Response values. Note that post-processing can not ensure that new estimate will have a lower error. It works for many cases but not all.
}
\usage{
postProcessingVotes(object,
	nbModels = 1,
	idx = 1,
	granularity = 1,
	predObject = NULL,
	swapPredictions = FALSE,
	X = NULL,
	Xtest = NULL,
	imbalanced = FALSE,
	OOB = FALSE,
	method = c("cutoff", "bias", "residuals"),
	keep2ndModel = FALSE, 
	largeData = FALSE,
	\dots)
}

\arguments{
	\item{object}{
 		a randomUniformForest object with OOB data.
}
	\item{nbModels}{
 		how many models to build for new estimates. Usually one is enough.
}
	\item{idx}{
 		how many values to choose in OOB model for each new predicted value. Usually one is enough.
}
	\item{granularity}{
 		degree of precision needed for each old estimate value. Usually one is enough.
}
	\item{predObject}{
 		if current model is built with full sample, then using an old model 'predObject' (a randomUniformForest object) that have OOB data can help to reduce error. Must be used with 'swapPredictions = TRUE'
}
	\item{swapPredictions }{
 		set it to TRUE if two models, current one without OOB data and old one with OOB data, have to be used for trying to reduce prediction error.
}
	\item{X}{
		not currently used.
}
	\item{Xtest}{
		test data in the case of regression, for a more friendly output of the model.
	}
	\item{imbalanced}{
		if TRUE, may improve metrics in the case of imbalanced datasets.
}
\item{OOB}{
	if FALSE, does not use OOB informations.
}
	\item{method}{
    for classification, if expected bias is enough high, one may use it as a method to improve AUC. Otherwise, use the default one, 'cutoff'. Both tend to get the same results, despite a few tests for now, but 'bias' method seems more robust. 'residuals' is used in regression only as a powerful but computationally intensive method and replaces the default internal one.
}
\item{keep2ndModel}{
	if TRUE, and for regression, keep the model based on residuals for further modelling and predictions.
}
\item{largeData}{
	if TRUE, and for regression, use \code{\link{rUniformForest.big}} to compute model for the residuals.
}
\item{\dots}{
	arguments to use for the computation of the model from the residuals.
}
}

\references{
	Xu, Ruo, Improvements to random forest methodology (2013). Graduate Theses and Dissertations. Paper 13052.
}

\value{
  a vector of predicted values.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
# Note that post-processing works better with enough trees, at least 100, and enough data
n = 200;  p = 20
# Simulate 'p' gaussian vectors with random parameters between -10 and 10.
X <- simulationData(n,p)

# give names to features
X <- fillVariablesNames(X)

# Make a rule to create response vector
epsilon1 = runif(n,-1,1)
epsilon2 = runif(n,-1,1)

# a rule with many noise (only four variables are significant)
Y = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]

# randomize then  make train and test sample
twoSamples <- cut(sample(1:n,n), 2, labels = FALSE)

Xtrain = X[which(twoSamples == 1),]
Ytrain = Y[which(twoSamples == 1)]
Xtest = X[which(twoSamples == 2),]
Ytest = Y[which(twoSamples == 2)]

# compute an accurate model (in this case bagging and log data works best) and predict
rUF.model <- randomUniformForest(Xtrain, Ytrain, xtest = Xtest, ytest = Ytest, 
bagging = TRUE, logX = TRUE, ntree = 60, threads = 2)

# get mean squared error
rUF.model

# post process
newEstimate <- postProcessingVotes(rUF.model)

# get mean squared error
sum((newEstimate - Ytest)^2)/length(Ytest)

## regression do not use all data but sub-samples.
## Comparing, when using full sample (but then, we do not have OOB data)

# rUF.model.fullsample <- randomUniformForest(Xtrain, Ytrain, xtest = Xtest, ytest = Ytest,
# subsamplerate = 1, bagging = TRUE, logX = TRUE, ntree = 60, threads = 2)
# rUF.model.fullsample

## Nevertheless we can use old model with OOB data to fit a new estimate
## newEstimate.fullsample <- postProcessingVotes(rUF.model.fullsample,
# predObject = rUF.model, swapPredictions = TRUE)

## get mean squared error
# sum((newEstimate.fullsample - Ytest)^2)/length(Ytest)
}