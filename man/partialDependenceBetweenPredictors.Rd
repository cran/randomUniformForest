\name{partialDependenceBetweenPredictors} 
\alias{partialDependenceBetweenPredictors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Partial Dependence between Predictors and effect over Response}

\description{
  Computes partial dependence between two predictors, and their effects on response values.
}

\usage{
partialDependenceBetweenPredictors(Xtest, importanceObject, features, 
	whichOrder = c("first", "second", "all"),
	perspective = FALSE,
	outliersFilter = FALSE)
}

\arguments{
  \item{Xtest}{
     a matrix or data frame specifying test (or train) data.
}
  \item{importanceObject }{
     an object of class importance.
}
  \item{features}{
    features that one needs to see dependence with responses (either train responses or predicted values).
}
  \item{whichOrder}{
     at which order, partial dependence does it need to be computed ? 
}
  \item{perspective}{
     plot dependence in 3D ? 
}
  \item{outliersFilter}{
  filter outliers ?
}
}

\details{
  partial dependence shows how response values are evolving depending to two target features and knowing the distribution of all others features. Note that it is essential to first have a view on target features importance.
  Steps can be given as this :
  1- get importance object (and plot it) that shows almost all objects that could explain the link between features and reponse.
  2- compute dependence between the two target features and response values. Link that with step 1 to obtain precise effects on response feature.
}

\value{
 A matrix containing values of the two features and values of expected conditional response, for regression. For classification, responses column is whether or not the values of the two features share the same class. partialDependenceBetweenPredictors() function also returns a set of figures representing how the dependency between the two features is affecting classes of the problem (or responses values in case of regression). It also retruns a measure of dependence between the two predictors at first (one of the two feature is supposed to be the most important one in the data) and second order (one of the two feature is supposed to be the second most important one). 
}

\note{
Please note that many plots (4 or 5, depending on task and option) will be produced. One has to first tile windows in R, to get the big picture.
There is currently not a generic way (at least, I did not find it) to automatically tile windows in R in all platforms.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{partialImportance}}
}

\examples{
## not run
## NOTER: please remove comments to run
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

## compute model : train then test in the same function. 
## Note that OOB error hardly works due to imbalanced classes.
## Handle it, using second upper Breiman's bound...
# car.ruf <- randomUniformForest(car.data.train, car.class.train, 
# xtest = car.data.test, ytest = car.class.test)
# car.ruf

## get for example two most important features, using table of features
## we choose "buying" and "safety"
# summary(car.ruf) 

## compute importance object (getting almost all objects wich could lead to a better explanation) 
## with deepest level of interactions to get enough points
# car.ruf.importance <- importance.randomUniformForest(car.ruf,
# Xtest = car.data.train, maxInteractions = 8)

## compute and plot partial dependence between "buying" and "safety" on train data
## heatmap leads to clusters, dependence tells what happens there, class distributions 
## tell how it happens and one has to connect these with 'importanceObject' figures
## trick : enlarge then reduce window to see figures
# pDbetweenPredictors.car.buyingAndSafety <- partialDependenceBetweenPredictors(car.data.train,
# car.ruf.importance, c("buying", "safety"), whichOrder = "all") 


#### Regression : "Concrete Compressive Strength" data 
## (http://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength)

# data(ConcreteCompressiveStrength)
# ConcreteCompressiveStrength.data = ConcreteCompressiveStrength

# n <- nrow(ConcreteCompressiveStrength.data)
# p <- ncol(ConcreteCompressiveStrength.data)

# trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

## train examples
# concrete.data.train <- ConcreteCompressiveStrength.data[trainTestIdx == 1, -p]
# concrete.responses.train <- ConcreteCompressiveStrength.data[trainTestIdx == 1, p]

## test data
# concrete.data.test <- ConcreteCompressiveStrength.data[trainTestIdx == 2, -p]
# concrete.responses.test <- ConcreteCompressiveStrength.data[trainTestIdx == 2, p]

## model
# concrete.ruf <- randomUniformForest(concrete.data.train, concrete.responses.train,
# featureselectionrule = "L1")
# concrete.ruf

## importance at the deepest level of interactions
# concrete.ruf.importance <- importance.randomUniformForest(concrete.ruf,
# Xtest = concrete.data.train, maxInteractions = 8)

## compute and plot partial dependence between "Age" and "Cement" on train data, 
## without 3D representation and with filter upon outliers
# pDbetweenPredictors.concrete.cementAndwater <- 
# partialDependenceBetweenPredictors(concrete.data.train,
# concrete.ruf.importance, c("Age", "Cement"), whichOrder = "all", outliersFilter = TRUE)

## compute and plot partial dependence between "Age" and "Cement" on train data, 
## with 3D representation (slower)
# pDbetweenPredictors.concrete.cementAndwater <- 
# partialDependenceBetweenPredictors(concrete.data.train,
# concrete.ruf.importance, c("Age", "Cement"), whichOrder = "all", perspective = TRUE)
}