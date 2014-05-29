\name{partialImportance} 
\alias{partialImportance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial Importance for random Uniform Forests
}

\description{
  Describe what features are the most important for one specifically class (in case of classification) or explain features that are affecting, the most, variability of the response (for regression), either for train or test sample.
}

\usage{
partialImportance(X, object, 
	whichClass = NULL,
	threshold = NULL,
	thresholdDirection = c("low", "high"),
	border = NA,
	nLocalFeatures = 5)
}

\arguments{
  \item{X}{
     a matrix or data frame specifying test (or train) data.
}
  \item{object}{
     an object of class importance.
}
  \item{whichClass}{
    for classification only. The index of the class that needs to be computed.
}
  \item{threshold}{
    for regression only. The threshold point. Partial importance will compute importance below or above this point (see below).
}
  \item{thresholdDirection }{
     where importance does it need to be fitted ? "low" will be lower than 'threshold' and "high" will be above.
}
 \item{border}{
  visualization option. Draw border around barplots ? 
}
  \item{nLocalFeatures}{
  how many features does it need to be assessed ?
}
}

\details{
 Partial importance must be used in conjunction with importance object, since it explains what features have influence on a class (or on
the variability of the response) but not how (that can be lower, higher or all values of one or many features). Partial Importance produces final (and local) rules that lead to the decision (observed or predicted class, observed or predicted value in cas of regression).
}

\value{
 A vector containing relative influence of the  most important features in a decreasing order. 
}

\author{
 Saip Ciss \email{saip.ciss@wanadoo.fr}
}


\examples{
## not run 
## NOTE : please remove comments to run 
#### Classification: "car evaluation" data (http://archive.ics.uci.edu/ml/datasets/Car+Evaluation)

# data(carEvaluation)
# car.data = carEvaluation 

# n <- nrow(car.data)
# p <- ncol(car.data)

# trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

## train examples
# car.data.train <- car.data[trainTestIdx == 1, -p]
# car.class.train <- as.factor(car.data[trainTestIdx == 1, p])

## test data
# car.data.test <- car.data[trainTestIdx == 2, -p]
# car.class.test <- as.factor(car.data[trainTestIdx == 2, p])

## compute model : train and predict in the same function
# car.ruf <- randomUniformForest(car.data.train, car.class.train, xtest = car.data.test)
# car.ruf

## compute importance object deeper (increasing level of interactions) on test data
# car.ruf.importance <- importance(car.ruf, Xtest = car.data.test,
# maxinteractions = 6, threads = 2)

## get partial importance for the three most important features 
## in test set for "unacceptable" cars
# car.ruf.partialImportance.unacc <- partialImportance(car.data.test, 
# car.ruf.importance, whichClass = "unacc")
									
## get partial importance for the three most important features in test set for "acceptable" cars
# car.ruf.partialImportance.acc <- partialImportance(car.data.test, car.ruf.importance, 
# whichClass = "acc")									

#### Regression : "Concrete Compressive Strength" data
## (http://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength)

# data(ConcreteCompressiveStrength.data)
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

## importance on train data
# concrete.ruf.importance <- importance.randomUniformForest(concrete.ruf, 
# Xtest = concrete.data.train, maxInteractions = 8, threads = 2)

## partial importance for concrete compressive strength higher than 40.
# concrete.ruf.partialImportance <- partialImportance(concrete.data.train, 
# concrete.ruf.importance, threshold = 40, thresholdDirection = "high")

## partial importance for concrete compressive strength lower than 30.
# concrete.ruf.partialImportance <- partialImportance(concrete.data.train, 
# concrete.ruf.importance, threshold = 30, thresholdDirection = "low")
}