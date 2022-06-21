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
	outliersFilter = FALSE,
	maxClasses = max(10, which.is.factor(Xtest[,, drop = FALSE], 
	count = TRUE)),
	bg = "grey")
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
    if TRUE, plot dependence in 3D with animation or not. Note that model will lead to interpolation. If one needs extrapolation, 'outputperturbationsampling' option can help?
}
  \item{outliersFilter}{
  filter outliers ?
}
 \item{maxClasses}{
  for variables with discrete values that need to be treated as categorical for a better visualization and
  that have more than 'maxClasses' unique values. Or for categorical variables that one knows to be categorical but
  whose data are stored as a R matrix.
}
  \item{bg}{
  background color for the plot. Type 'bg = "none"' to get a white background.
}
}

\details{
  partial dependence shows how response values are evolving depending to a pair of variables and knowing the distribution of all others covariates. Note that it is essential to first have a view on variable importance object.
  Steps can be given as this :\cr
  1- get importance object (and plot it), that shows almost all objects that could explain the link between features and response.\cr
  2- compute dependence between two target features and response values. Link this point with step 1 to obtain precise effects on response feature.
  
  Note that the function shows both the dependence between the pair of variables and the effect over the whole response values.
}

\value{
 A matrix containing values of the two features and values of expected conditional response, for regression. For classification, responses column is whether or not the values of the two features share the same class. partialDependenceBetweenPredictors() function also returns a set of figures representing how the dependency between the two features is affecting classes of the problem (or responses values in case of regression). It also returns a measure of dependence between the two predictors at first (one of the two feature is supposed to be the most important one in the data) and second order (one of the two features is supposed to be the second most important one). 
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

#### Classification: "car evaluation" data 
## (http://archive.ics.uci.edu/ml/datasets/Car+Evaluation)
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
## use 'categorical' option to better handle  variable importance
# car.ruf <- randomUniformForest(car.data.train, car.class.train, 
# xtest = car.data.test, ytest = car.class.test, categorical = "all")
# car.ruf

## get for example two most important features, using the table
## we choose "buying" and "safety"
# summary(car.ruf) 

## compute importance with deepest level of interactions to get enough points
# car.ruf.importance <- importance.randomUniformForest(car.ruf,
# Xtest = car.data.train, maxInteractions = 6)

## compute and plot partial dependence between "buying" and "safety" on train data
## is interaction leading to the same class (underlying structure)?
## is dependence linear or, for categorical variables, how is the effect of cross-tabulating 
## the variables over each class ? 
## for which values of the pair is the dependence most effective ?

# pDbetweenPredictors.car.buyingAndSafety <- partialDependenceBetweenPredictors(car.data.train,
# car.ruf.importance, c("buying", "safety"), whichOrder = "all") 

## Interpretation : 
## 1 - if "safety" is 'low', with average frequency of (around) 63 percent the 
## label of the evaluated car would be the same than the one coming from any buying category
## meaning that, the buying price will not have influence.
## If safety is 'high', the evaluated car will depend on the buying price...
## In the same sense, if the buying price is 'very high', the label of the evaluated car 
## will also depend to the safety category most of the time. 
## For a high (or lower) buying price, it will be much less dependent.
## But:
## 2 - the "heatmap..." states that confidence (more data) will be greater for cases 
## where buying price is very high (with safety low) or high (with safety high)
## 3 - Hence looking to the "dependence between predictors",  shows that if buying price is high, 
## it means that safety must be also high before one draws a conclusion.
## 4 - Looking "Variable Importance based on interactions" tells which variable dominates 
## the other when a car is evaluated.

#### Regression : "Concrete Compressive Strength" data 
## (http://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength)

# data(ConcreteCompressiveStrength)
# ConcreteCompressiveStrength.data = ConcreteCompressiveStrength

# n <- nrow(ConcreteCompressiveStrength.data)
# p <- ncol(ConcreteCompressiveStrength.data)

# set.seed(2015)
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

## Assessing test set only 
## importance at the deepest level of interactions
# concrete.ruf.importance <- importance.randomUniformForest(concrete.ruf,
# Xtest = concrete.data.test, maxInteractions = 8)

## compute and plot partial dependence between "Age" and "Cement", 
## without 3D representation and with filter upon outliers

# pDbetweenPredictors.concrete.cementAndAge <- 
# partialDependenceBetweenPredictors(concrete.data.test,
# concrete.ruf.importance, c("Age", "Cement"), whichOrder = "all", outliersFilter = TRUE)

## compute and plot partial dependence between "Age" and "Cement", 
## with 3D representation (slower)

# pDbetweenPredictors.concrete.cementAndAge <- 
# partialDependenceBetweenPredictors(concrete.data.test,
# concrete.ruf.importance, c("Age", "Cement"), whichOrder = "all", perspective = TRUE)
}

\concept{partial dependence}