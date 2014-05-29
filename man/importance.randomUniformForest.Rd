\name{importance.randomUniformForest} 
\alias{importance.randomUniformForest}
\alias{importance}
\alias{plot.importance}
\alias{print.importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variables Importance for random Uniform Forests}

\description{
  Compute object thats leads to a full analysis of features (importance, dependency, interactions, selection, ...).
}

\usage{
	\method{importance}{randomUniformForest}(object,
	maxVar = 30,
	maxInteractions = 3,
	Xtest = NULL,
	predObject = NULL,
	\dots)
	\method{plot}{importance}(x,
	nGlobalFeatures = 30,
	nLocalFeatures = 5,
	Xtest =	NULL,
	whichFeature = NULL,
	whichOrder = "all",
    outliersFilter = TRUE,
	formulaInput = NULL,
	border = NA,
	\dots)
	\method{print}{importance}(x, \dots)
}

\arguments{
	\item{x, object}{
     an object of class randomUniformForest (for the 'print' and  'plot' method, an object of class 'Importance').
}
	\item{nGlobalFeatures, nLocalFeatures}{ for the 'plot' method, number of global and local features to show.
}
	\item{maxVar}{
     maximum  number of features to display.
}
	\item{maxInteractions }{
     maximum  order of interactions. Default value is 3, meaning function will compute interactions for each variable at first order (current variable is supposed to be the most important), second order (current variable is supposed to be the second most important, but the first is unknown) and third order (current variable is supposed to be the third most important, but the first and the second are unknown) and so on, depending on the value of 'maxInteractions'.
}
	\item{Xtest}{
     current matrix used to compute 'object' model. Can be either training or test matrix. If it is the latter, please read below.
}
	\item{whichFeature}{
    for the 'plot' method, the feature (by its name or position) that one need to be asessed. It will be used in partial dependence. Useful only if the feature is not an important one for the model.
}
	\item{whichOrder}{ 
	for the 'plot' method, the order(s) at which some of the computation, e.g. partial dependence, has to be done. At "first" order, computation is done considering only each feature as the most important. At "second", it is considered at the second most important (and the first is unknown). At "all", each feature is assessed considering it is the most important, then the second most (and the first is unknown), then the third most (and the first and the second are unknown), and so on until 'maxInteractions' value is reached. Then one has to first look to this latter value to know at which maximum order a feature will be assessed. Influence of the variable will be a combination of all orders of local importance, depending on the level of interactions specified by 'maxInteractions'.
}
	\item{outliersFilter}{
	for the 'plot' method, do outliers of a feature need to be removed ? if TRUE, observations above 0.95 quantile and below 0.05 quantile will be removed.
}
	\item{formulaInput}{
	for the 'plot' method, if one uses formula, it has to be copy there in order to match figures. Not recommended, since formula can lead to unexpected effects.
}
	\item{border}{
	for the 'plot' method, if positive value, draw borders around barplots.
}
  \item{predObject}{
  	if 'object' was computed without an evaluation of test data, one must provide test data in 'Xtest' option and a full prediction object (calling option type = "all", when using predict() function. See examples).
}
 \item{\dots}{
  	others options currently not used.
}
}

\details{
  'Importance' retrieves global importance from 'object' which is "variable importance from the model and for training data". Then it computes, either for train or test data and for at least two orders :
  1- interactions, which show a kind of dependence between two features, for all couples of variables. Interactions shows links between variables by using informations at the leaf level. Level of interactions is standardized between 0 and 1. Most important (highly predictive) features usually have more interactions, but a low predictive feature can have many interactions that gives insights about the big picture.
  2- overall interactions, that lead to local variable importance which is "variable importance from the model and interactions between features". Local variable importance shows how response values can be explained, using more than just predictive features. 
  3 - partial dependence, that leads to show how response values are evolving according to one feature, and knowing distribution of all others.
  4- variable importance over labels (in classification), that shows, on the same figure, how each important feature is matching labels from the view point of the model. Features with influence are the most discriminant over labels. In regression, object is called "Dependence on most important predictors". 
  
  Summarized, importance object tells (or leads to, using others functions) "what, how and when" features affect the response. Note that importance depends of hyper-parameters of the model (mtry, depth, nodesize, ...) and one has to take care. Categorical features (and all others variables or types of variable) are treated by randomUniformForest() function as numeric values (default option). Accuracy is usually not affected, but their importance has to be assessed more carefully, using all objects provided by importance() and dependence functions (see examples, for a case study summary).
}

\value{
 an object of class importance.
}

\note{
Please note that 'plot' method of importance object will produce many plots (5) and a prompt
entering a loop which will wait user input to plot one by one partial dependencies for most
important features. One has to first tile windows in R, to get the big picture.
There is currently not a generic way (at least, I did not find it) to automatically tile windows in R 
in all platforms.
}

\author{
  Saip  CISS \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{partialDependenceOverResponses}}, \code{\link{partialDependenceBetweenPredictors}},
 \code{\link{partialImportance}}
}

\examples{
## not run 
## NOTES: please remove comments to run, since importance and plot method will draw many graphics,
## use of option 'threads = 1' (disabling parallel processing) to speed up computing, 
## since parallel processing is slower for small sample.
## 1 - Importance for Classification and Regression with formula
####  Classification

# data(iris)
# iris.ruf <- randomUniformForest(Species ~ ., data = iris, threads = 1)

## global importance: giving 2 ways to explain importance (table and visualization)
# summary(iris.ruf)

## much more about importance : giving 4 ways to explain importance in a complementary fashion 
## of global importance
# iris.ruf.importance <- importance(iris.ruf, Xtest = iris, threads = 1)

## get importance summary
# iris.ruf.importance

## visualizing all in one
# plot(iris.ruf.importance, Xtest = iris)

## get details about observations (link between observations and features) : 
## note that class and features are replaced by their internal numeric values.
## Retrieve model of observations importance
# iris.observationsImportance <- 
# data.frame(iris.ruf.importance$localVariableImportance$obsVariableImportance)

## first impute true labels
## iris.observationsImportance$class <- as.factor(iris.observationsImportance$class)
# levels(iris.observationsImportance$class) <- iris.ruf$variables

## draw a sample of 10 observations to see which features affect the most every observation 
## and the label associated (either from train or test sample)
# iris.observationsImportance[sample(nrow(iris), 10), ]

#### Regression

# data(airquality)
# airquality.data = airquality

## impute NA
# airquality.NAimputed <- fillNA2.randomUniformForest(airquality.data)

## compute model
# ozone.ruf <- randomUniformForest(Ozone ~ ., data = airquality.NAimputed, threads = 1)
# ozone.ruf

## print and plot (global) variable importance
# summary(ozone.ruf)

## much more about importance
# ozone.ruf.importance <- importance(ozone.ruf, Xtest = airquality.NAimputed, threads = 1)

## visualizing all in one : in case of formula, 'formulaInput' is needed for the 'plot' method
# plot(ozone.ruf.importance, Xtest = airquality.NAimputed, formulaInput = ozone.ruf$formula)

## 2- Importance for Classification and Regression without formula (more usual and recommended 
# for random Uniform Forests)

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
# xtest = car.data.test, ytest = car.class.test, threads = 1)
# car.ruf

## global importance: note that 'safety' does not appear to be an important feature in the barplot
## but in the table, it is, by far, the most important feature of unacceptable (unacc) cars.
# summary(car.ruf)

## interactions and local importance tell most of the story...
## (not run)
# car.ruf.importance <- importance(car.ruf, Xtest = car.data.train, threads = 1)

## ...that can be used to explain train data
# plot(car.ruf.importance, Xtest = car.data.train)

## or explain test data
# car.ruf.importance.test <- importance(car.ruf, Xtest = car.data.test, threads = 1)
# plot(car.ruf.importance.test, Xtest = car.data.test)

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

## test data
# Concrete.data.test <- ConcreteCompressiveStrength.data[trainTestIdx == 2, -p]
# Concrete.responses.test <- ConcreteCompressiveStrength.data[trainTestIdx == 2, p]

## model
# Concrete.ruf <- randomUniformForest(Concrete.data.train, Concrete.responses.train,
# featureselectionrule = "L1", threads = 1)
# Concrete.ruf

## predictions : option ' type = "all" ' is needed to manually assess importance of a test set 
# Concrete.ruf.pred <- predict(Concrete.ruf, Concrete.data.test, type = "all")

## more options about importance : more interactions
# Concrete.ruf.importance <- importance(Concrete.ruf, Xtest = Concrete.data.test,
# maxInteractions = 6, predObject = Concrete.ruf.pred, threads = 1)

## and more features
# plot(Concrete.ruf.importance, nLocalFeatures = 7, Xtest = Concrete.data.test)
}

\concept{variable importance}
\concept{variable selection}