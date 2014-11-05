\name{generic.cv} 
\alias{generic.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generic k-fold cross-validation}

\description{
 Performs k-fold cross-validation 'n' times for any specified algorithm, using two of many metrics(test error, AUC, precision,...)
}

\usage{
generic.cv(X, Y, 
nTimes = 1, 
k = 10, 
seed = 2014, 
regression = TRUE, 
genericAlgo = NULL, 
specificPredictFunction = NULL, 
metrics = c("none", "AUC", "precision", "F-score", "L1", "geometric mean", 
"geometric mean (precision)"))
}

\arguments{
   \item{X}{
    a matrix or dataframe of observations
}
 \item{Y}{
    a response vector for the observed data.
}
 \item{nTimes}{
   number of times that k-fold cross-validation need to be performed.
}
 \item{k}{
    how many folds ?
}
\item{seed}{
    the seed for reproducibility.
}
\item{regression}{
	if TRUE, performs regression.
}
 \item{genericAlgo}{
    wrapper function to embed the algorithm needed to be assesses. One can eventually add options. NULL is only for convenience. Wrapper function is needed to assess cross-validation.
}
 \item{specificPredictFunction}{
    if the assessed model does not support the R generic method 'predict', one has to define here, with a function, how predictions have to be generated.
}
\item{metrics}{
    One of many other metrics one can call with the standard one, test error (or MSE for regression). 
}
}

\value{
a list with the following components :
   \item{testError}{the values of test error.}
  \item{avgError}{mean of test error.}
  \item{stdDev}{standard deviation of test error.}
  \item{metric}{values of the metric chosen.}
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
## not run
# data(iris)
# Y <- iris$Species
# X <- iris[,-which(colnames(iris) == "Species")]

## 10-fold cross-validation for randomForest algorithm
## create the wrapper function
# genericAlgo <- function(X, Y) randomForest(X, Y)

## run 
# RF.10cv.iris <- generic.cv(X, Y, genericAlgo = genericAlgo, regression = FALSE)

## 10-fold cross-validation for Gradient Boosting Machines algorithm (gbm package)
## create the wrapper function

# require(gbm) || install.packages("gbm")
# genericAlgo <- function(X, Y) gbm.fit(X, Y, distribution = "multinomial", n.trees = 500,
# shrinkage = 0.05, interaction.depth = 24, n.minobsinnode = 1) 

## create a wrapper for the prediction function of gbm
# nClasses = length(unique(Y))
# specificPredictFunction <- function(model, newdata)
# {
#	modelPrediction = predict(model, newdata, 500) 
#	predictions = matrix(modelPrediction, ncol = nClasses )
#	colnames(predictions) = colnames(modelPrediction)
#	return(as.factor(apply(predictions, 1, function(Z) names(which.max(Z)))))
# }

## run
# gbm.10cv.iris <- generic.cv(X, Y, genericAlgo = genericAlgo, 
# specificPredictFunction = specificPredictFunction, regression = FALSE)
}