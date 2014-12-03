\name{fillNA2.randomUniformForest} 
\alias{fillNA2.randomUniformForest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Missing values imputation by randomUniformForest}

\description{
 Impute missing values using randomUniformForest. Each variable containing missing values is, in turn, considered as a responses vector, where non-missing values are training responses and missing values, responses to predict.
}

\usage{
fillNA2.randomUniformForest(X, Y = NULL, ntree = 100, 
	mtry = 1, 
	nodesize = 10,
	categoricalvariablesidx = NULL,
	NAgrep = "", 
	maxClasses = floor(0.01*min(3000, nrow(X))+2),
	threads = "auto")
}

\arguments{
   \item{X}{
     A data frame or matrix of predictors.
}
 \item{Y}{
     response vector, e.g. training labels. Note that it is strongly recommended to use the default values in order 
	 to not creating bias, for example if test set needs also to be imputed.
}
 \item{ntree}{
     number of trees to grow for each predictor with missing values. Default value usually works, unless one wants a better level of accuracy.
}
 \item{mtry}{
     number of predictors to try for growing the forest. Default value means full randomized trees will be grown, relying more on convergence of random forests (to achieve a good estimation) and less to the data.
}
 \item{nodesize}{
     smallest number of observations in each leaf node.
}
 \item{categoricalvariablesidx}{
    see \code{\link{randomUniformForest}}.
}
 \item{NAgrep}{
     which symbol(s) (in case of a data frame), e.g. "?", in addition to "NA" (that is automatically matched) have to be considered as missing values. One can put many symbols, for example 'NAgrep = c(" ", "?")'. 
}
\item{maxClasses}{
	maximal number of classes for a categorical variable. Note that the function will look first the structure of the dataset to detect categorical variables. One  has to take care of categories with many classes but not defined as factors. Hence in case of doubt, use a large value of 'maxClasses' is preferable. 
}
\item{threads}{
    how many logical cores to complete data. Default values will let algorithm use all available cores.  
}
}

\details{
	Algorithm uses randomUniformForest to complete matrix. At the first step, all missing values are identified and rough imputation is done using  most frequent (or median) values for each predictor with missing values. Then, these predictors are, one by one, considered as responses, where training data are the ones with non missing values (or roughly fixed) and test sample are data whose values are really missing. There is only one iteration (over all predictors) since this version already requires computing time. Hence, it is strongly  recommended to use it only if others imputation models do not work or if speed is not mandatory. Note that the function will try to render exactly the same structure than the one in the original dataset.
}

\value{
 a matrix or a data frame containing completed values.
}

\author{
 Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
## not run
## get same example as rfImpute() function from randomForest package
# data(iris)
# iris.na <- iris
# set.seed(111)
## artificially drop some data values.
# for (i in 1:4) iris.na[sample(150, sample(20)), i] <- NA

## imputation 
# iris.imputed <- fillNA2.randomUniformForest(iris.na, threads = 1)

## model with imputation
# iris.NAfixed.ruf <- randomUniformForest(Species ~ ., iris.imputed, 
# BreimanBounds = FALSE, threads = 1)
# iris.NAfixed.ruf

## Compare with true data (OOB evaluation)
# iris.ruf <- randomUniformForest(Species ~ ., iris, BreimanBounds = FALSE, threads = 1)
# iris.ruf
}