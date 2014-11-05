\name{predict.randomUniformForest}
\alias{predict.randomUniformForest}
\alias{predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict method for random Uniform Forests objects
}

\description{Prediction of test data using random Uniform Forests.
}

\usage{
\method{predict}{randomUniformForest}(object, X, 
	type = c("response", "prob", "votes", "confInt", 
	"ranking", "quantile", "truemajority", "all"),
	classcutoff = c(0,0), 
	conf = 0.95,
	whichQuantile = NULL,
	rankingIDs = NULL,
	threads = "auto", 
	parallelpackage = "doParallel",
	\dots)
}

\arguments{
	\item{object}{
    an object of class randomUniformForest, as one created by the randomUniformForest( ) function.
}
	\item{X}{
     a data frame or matrix containing new data (without response values).
}
	\item{type}{
    one of response, prob, votes, prediction intervals, ranking, quantile, true majority or raw outputs, indicating respectively the type of output: 'predicted values', 'matrix of class probabilities', 'matrix of vote counts', 'prediction interval for each response', 'ranked responses' in case of recommendation scheme (note that it is currently experimental and it, first, requires classification modelling), 'quantile regression', 'predicted values based on nodes votes instead of tree votes', or lastly, 'raw outputs' of a random uniform forest. Note that "confInt" and "quantile" use an estimate of the conditional expectation for each tree parameter (the randomness), which contains many duplicates, due to bootstrap and the nature of ensemble learning, to produce estimates. These ones are usually loose since forest replicates distribution of Y (knowing X) to take decisions. A more accurate option is to use \code{\link{bCI}} function.
}
	\item{classcutoff}{
    see \code{\link{randomUniformForest}}.
}
	\item{conf}{
     if type == "confInt", value of 'conf' (greater than 0 and lesser than 1) is the desired level of confidence for prediction interval.
}
	\item{whichQuantile}{
    if type == "quantile", value of 'whichQuantile' is the desired quantile (greater than 0 and lesser than 1).
}
	\item{rankingIDs}{
     experimental. If type == "ranking", vector of ID if one wants to rank responses instead of classify them. Ranking usually involves many same users (or objects, or items) whose choices have to be ranked, e.g., by a recommendation engine. For example, if 'user1' is present four times in test sample, its ID gives to the model a way to identify him (or it) and to order response values, beginning by the most relevant. Then, for many users, the process is repeated giving, finally, predictions from most to least relevant, taking into account all users. random Uniform Forests approach currently uses class probabilities to order predictions values. Note that optimizing AUC is one of the methods leading to good ranking methods if choices are binary, and NDCG is a way to assess model. 
}
	\item{threads}{
    see \code{\link{randomUniformForest}}.
}
	\item{parallelpackage}{
    see \code{\link{randomUniformForest}}.
}
\item{\dots}{
    not used currently.
}
}

\value{
  \item{response }{predicted values. Default option that returns values in the same way than original training responses.}
  \item{prob }{for classification only. Matrix of class probabilities.}
  \item{votes }{matrix of vote counts. Each row is an observation and each columns is tree output.}
  \item{confInt }{for regression only. Matrix where each row is an observation and each column one of the prediction interval bounds. }
  \item{quantile }{for regression only. Vector of predicted quantiles for 'conf' (value of the option) level of confidence.}
  \item{ranking }{a matrix or data frame. Description will be updated soon.}
  \item{truemajority }{predicted values, using raw outputs of trees (not majority vote). This option makes sense if one set 'nodesize' option greater than 1. Hence, aggregation is participative (at the leaf level) and not representative (at the tree level).}
  \item{all }{raw outputs of the model. Not useful, unless further computation is needed, for example in case of Post-processing.}
}

\author{
  Saip Ciss  \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{postProcessingVotes}}, \code{\link{bCI}}, \code{\link{model.stats}}  
}

\examples{
## same as randomForest example
#  data(iris)
# set.seed(111)
# ind <- sample(2, nrow(iris), replace = TRUE, prob = c(0.8, 0.2))

# iris.ruf <- randomUniformForest(Species ~ ., data = iris[ind == 1,], OOB = FALSE, 
# importance = FALSE, threads = 1)
# iris.pred <- predict(iris.ruf, iris[ind == 2,])

# table(observed = iris[ind == 2, "Species"], predicted = iris.pred)

## get all votes : note that aliases of classes are used internally and, for intermediate
## results, are not converted to their true values
# iris.all.votes <- predict(iris.ruf, iris[ind == 2,], type = "votes")

## get class probabilities
# iris.class.prob <- predict(iris.ruf, iris[ind == 2,], type = "prob")
# iris.class.prob
}