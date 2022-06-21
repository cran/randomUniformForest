\name{fillNA2.randomUniformForest}
\alias{fillNA2.randomUniformForest}
\alias{rufImpute}
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
	maxClasses = max(10, which.is.factor(X[, , drop = FALSE], count = TRUE)),
	na.action = c("accurateImpute", "fastImpute", "veryFastImpute"),
	threads = "auto",
	\dots)
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
     number of predictors to try when growing the forest. Default value means full randomized trees will be grown, relying more on convergence of random forests (to achieve a good estimation) and less to the data.
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
\item{na.action}{
	Type of imputation to make.
}
\item{threads}{
    how many logical cores to complete data. Default values will let algorithm use all available cores.  
}
\item{\dots}{
     not currently used.
}
}

\details{
	Algorithm uses randomUniformForest to complete matrix. At the first step, all missing values are identified and rough imputation is done using  most frequent (or median) values for each predictor with missing values. Then, these predictors are, one by one, considered as response vectors, where training data are the ones with non missing values (or roughly fixed) and test sample are data whose values are really missing. There is only one iteration (over all predictors) since this version already requires computing time. Hence, it is strongly  recommended to use it only if others imputation models do not work or if speed is not mandatory. 
}

\value{
 a matrix or a data frame containing completed values.
}

\note{
The function will try to render exactly the same structure than the one in the original dataset. Note also that '""' (void character) will be automatically considered as a true category; hence if it is really a missing value 
one should replace it, for example, by '"?"' and type 'NAgrep = "?"' in the options of the function. See examples.

Note also that rufImpute( ) is an alias of the function.
}

\author{
 Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
## not run

## A - usual case
# get same example as rfImpute() function from randomForest package
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

## B - hard case : titanic dataset
## see http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic.html
## for more informations

# URL = "http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/"
# dataset = "titanic3.csv"
# titanic3 = read.csv(paste(URL, dataset, sep =""))

## 1309 observations, 14 variables, 2 classes
# XY = titanic3
# Y = XY[,"survived"]
# X = XY[,-which(colnames(XY) == "survived")]

## remove name :
# XX = X = X[,-2]

## first imputation : not working good because missing values are "" and NA
# X.imputed.1 = rufImpute(X)

# head(X.imputed.1)
# head(X)

## 1 - one first has to replace all factors by characters in order to handle ""
# for (j in 1:ncol(X))
# {	if (is.factor(X[,j])) { XX[,j] = as.character(X[,j]) } }

## 2 - replace "" by "?"
# XX[which(XX == "", arr.ind  = TRUE)] = "?"

## 3 - impute by including "?" in missing values, 
## increasing 'maxClasses' to be sure to match all categorical variables
## setting categorical variables explicitly (they may be integers, factors or characters)
## Note : integers may also be viewed as numerical values

# str(XX) ##gives the type of all variables

# categoricalVariables = vector();  i = 1
# for (j in 1:ncol(X)) 
# { 
#	if (class(X[,j]) != "numeric") { categoricalVariables[i] = j; i = i + 1 }
# }

# X.imputed.1 = rufImpute(XX, NAgrep = "?", maxClasses = 1200, 
# categorical = categoricalVariables)

# Take a random sample and compare
# idx = sample(nrow(X), 20)
# X[idx,]
# X.imputed.1[idx,]

## modify eventually some numeric values like 'age' to match cases
## use more trees and less randomization to possibly increase accuracy 
## at the risk of less consistency.

# X.imputed.2 = rufImpute(XX, NAgrep = "?", maxClasses = 1200, mtry = 4, nodesize = 5,
# ntree = 200, categorical = categoricalVariables)

## 4- assess the imputed matrix : OOB evaluation
## - base model : omit missing values is not possible (too many ones). Roughly impute instead.
# titanic.baseModel.ruf = randomUniformForest(X, as.factor(Y), na.action = "fastImpute",
# categorical = categoricalVariables)

# titanic.baseModel.ruf

## - imputed model
# titanic.imputedModel.ruf = randomUniformForest(X.imputed.1, as.factor(Y), 
# categorical = categoricalVariables)

# titanic.imputedModel.ruf

## roughly (and internal) imputation works better in this case and one should investigate
## models and data to understand possible reasons (one influential feature, many categories, ...)
}