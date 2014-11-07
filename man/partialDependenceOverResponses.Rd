\name{partialDependenceOverResponses} 
\alias{partialDependenceOverResponses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial Dependence Plots and Models}

\description{
  Computes partial dependence between expected conditional response and all values of its target feature, knowing
the distribution of all others features in the data (e.g. marginal effect of the target feature over the response)
}

\usage{
partialDependenceOverResponses(Xtest, importanceObject, 
	whichFeature = NULL, 
	whichOrder = c("first", "second", "all"), 
	outliersFilter = FALSE, 
	plotting = TRUE, 
	followIdx = FALSE, 
	maxClasses = if (is.null(whichFeature)) { 10 } else { max(10, 
	which.is.factor(Xtest[, whichFeature, drop = FALSE], count = TRUE)) })
}

\arguments{
  \item{Xtest}{
     a matrix or data frame specifying test (or train) data.
}
 \item{importanceObject}{
     an object of class importance.
}
 \item{whichFeature}{
    feature that one needs to see dependence with responses (either train responses or predicted values).
}
  \item{whichOrder}{
     at which order, partial dependence does it need to be computed ? 
}
  \item{outliersFilter}{
  filter outliers ?
}
  \item{plotting}{
  plot partial dependence ?
}
\item{followIdx}{
  not currently used.
}
\item{maxClasses}{
  for variables with discrete values that need to be treated as categorical for a better visualization and
  that have more than 'maxClasses' unique values. Or for categorical variables that one knows to be categorical but
  whose data are stored as a R matrix.
}
}

\details{
  Partial dependence shows how expected conditional response is evolving with its target feature, knowing the distribution of all others features. For example, if one wants to know how response is moving, on average, depending on one target feature and for all possible values of the others features. For Classification, partial dependence shows the dependence of each important variable over classes. One key advantage of partial dependence is to allow extrapolation for ensemble models. Note that this version is inspired by Friedman (2001) ideas, but uses a different implementation.
}

\value{
 a matrix containing feature values and expected conditional responses.
}

\references{
Friedman, J.H., 2001. \emph{Greedy function approximation: A gradient boosting machine}. Ann. Statist. 29, 1189-1232.
}

\author{
 Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\seealso{
\code{\link{partialDependenceBetweenPredictors}}, \code{\link{partialImportance}}
}

\examples{
## not run
## NOTE:  please remove comments to run
#### Classification: "car evaluation" data (http://archive.ics.uci.edu/ml/datasets/Car+Evaluation)
# data(carEvaluation)
# car.data <- carEvaluation

# n <- nrow(car.data)
# p <- ncol(car.data)

# trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

## train examples
# car.data.train <- car.data[trainTestIdx == 1, -p]
# car.class.train <- as.factor(car.data[trainTestIdx == 1, p])

## test data
# car.data.test <- car.data[trainTestIdx == 2, -p]
# car.class.test <- as.factor(car.data[trainTestIdx == 2, p])

## compute model : train then test in the same function
# car.ruf <- randomUniformForest(car.data.train, car.class.train, 
# xtest = car.data.test, ytest = car.class.test, threads = 2)
# car.ruf

## compute importance
# car.ruf.importance <- importance.randomUniformForest(car.ruf, 
# Xtest = car.data.train, threads = 2)

## plot partial dependence, at all orders, between classes and 'safety' feature on train data.
## Note that data are mainly categorical, ordinal and characters (e.g. v-high, high, med, low 
## for 'safety'). Model deals internally with characters and compute them as numerical values 
## (e.g.  "v-high, high, med, low" leads to 4,1,3,2.)

## how safety affects labels, knowing all others features ?
# pD.car.safety <- partialDependenceOverResponses(car.data.train, car.ruf.importance, 
# whichFeature = "safety", whichOrder = "all")

## What's happen at first order ?
# pD.1rstOrder.car.safety <- partialDependenceOverResponses(car.data.train, car.ruf.importance,
# whichFeature = "safety", whichOrder = "first")

## plot partial dependence, at first order, default, 
## between classes and 'buying' feature on train data
# pD.1rstOrder.car.buying <- partialDependenceOverResponses(car.data.train, car.ruf.importance, 
# whichFeature = "buying")

## plot partial dependence, at second order, between classes and 'buying' feature on train data.
## Second order means 'buying' feature is supposed to be the second most important feature
## and first one is unknown.
# pD.2ndOrder.car.buying <- partialDependenceOverResponses(car.data.train, car.ruf.importance, 
# whichFeature = "buying", whichOrder = "second")

## if one wants to assess test set, e.g. test responses are unkown
# car.ruf <- randomUniformForest(car.data.train, car.class.train, 
# xtest = car.data.test, threads = 2)
# car.ruf

## compute importance object deeper (increasing level of interactions) and on test data
# car.ruf.importance <- importance.randomUniformForest(car.ruf, Xtest = car.data.test, 
# maxInteractions = 6, threads = 2)

# pD.1rstOrder.car.safety <- partialDependenceOverResponses(car.data.test, car.ruf.importance, 
# whichFeature = "safety")

# pD.2ndOrder.car.buying <- partialDependenceOverResponses(car.data.test, car.ruf.importance, 
# whichFeature = "buying", whichOrder = "second")

# pD.allOrders.car.priceOfMaintenance <- partialDependenceOverResponses(car.data.test,
# car.ruf.importance, whichFeature = "priceOfMaintenance", whichOrder = "all")

# pD.2ndOrder.car.priceOfMaintenance <- partialDependenceOverResponses(car.data.test,
# car.ruf.importance, whichFeature = "priceOfMaintenance", whichOrder = "second")

#### Regression : "Concrete Compressive Strength" data 
## (http://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength)

# data(ConcreteCompressiveStrength)
# ConcreteCompressiveStrength.data = ConcreteCompressiveStrength

# n <- nrow(ConcreteCompressiveStrength.data)
# p <- ncol(ConcreteCompressiveStrength.data)

# trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

## train examples
# Concrete.data.train <- ConcreteCompressiveStrength.data[trainTestIdx == 1, -p]
# Concrete.responses.train <- ConcreteCompressiveStrength.data[trainTestIdx == 1, p]

## model
# Concrete.ruf <- randomUniformForest(Concrete.data.train, Concrete.responses.train,
# featureselectionrule = "L1", threads = 2)
# Concrete.ruf

## importance object only for train examples
# Concrete.ruf.importance <- importance.randomUniformForest(Concrete.ruf,
# Xtest = Concrete.data.train, maxInteractions = 4, threads = 2)

## partial dependence only for train examples :
## at all orders
# pD.Concrete.Cement <- partialDependenceOverResponses(Concrete.data.train, 
# Concrete.ruf.importance, whichFeature = "Cement", whichOrder = "all")

# pD.Concrete.Age <- partialDependenceOverResponses(Concrete.data.train, 
# Concrete.ruf.importance, whichFeature = "Age", whichOrder = "all")

## at first order
# pD.1rstOrder.Concrete.Water <- partialDependenceOverResponses(Concrete.data.train, 
# Concrete.ruf.importance, whichFeature = "Water")
}

\concept{partial dependence}