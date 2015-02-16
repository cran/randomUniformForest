\name{importance.randomUniformForest} 
\alias{importance.randomUniformForest}
\alias{importance}
\alias{plot.importance}
\alias{print.importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variable Importance for random Uniform Forests}

\description{
  Compute object that leads to a full analysis of features (importance, dependency, interactions, selection, ...).
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
	whichOrder = if (ncol(x$globalVariableImportance) > 5) 
	{  
		if (nrow(x$localVariableImportance$obsVariableImportance) > 1000) 
		"first" else "all" 
	} 
	else
	{  if (nrow(x$localVariableImportance) > 1000) "first" else "all" },
    outliersFilter = FALSE,
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
    for the 'plot' method, the feature (by its name or position) that one need to be assessed. It will be used in partial dependence. Useful only if the feature is not an important one for the model.
}
	\item{whichOrder}{ 
	for the 'plot' method, the order(s) at which some of the computation, e.g. partial dependence, has to be done. At "first" order, computation is done considering only each feature as the most important. At "second", the feature is considered at the second most important (and the first is unknown). At "all", each feature is assessed considering it is the most important, then the second most (and the first is unknown), then the third most (and the first and the second are unknown), and so on until 'maxInteractions' value is reached. Then one has to first look to this latter value to know at which maximum order a feature will be assessed. Influence of the variable will be a combination of all orders of local importance, depending on the level of interactions specified by 'maxInteractions'.
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
  'Importance' retrieves (or computes) global importance from 'object' which is "variable importance from the model and for training data based on most predictive (by score) and discriminant(by class and class.frequency) features". 
  
  Then it computes, either for train or test data and for at least two orders :\cr
  1 - 'interactions', which show a kind of dependence between two features, for all pairs of variables. Interactions shows links between variables by using informations at the leaf level. Most important (highly predictive) features usually have more interactions, but a low predictive feature can have many interactions that provide insights about the big picture. Interactions are currently computed and plotted at the first and second order.\cr
  
  2 - 'overall interactions' that leads to "variable importance based on interaction". It shows how response values can be explained, using more than just predictive features and it is an aggregated function of interactions.\cr
  
  3 - 'partial dependence', that leads to show how response values are evolving according to one feature, and knowing the distribution of all others. See \code{\link{partialDependenceOverResponses}}.\cr
  
  4 - 'variable importance over labels' (in classification), that shows, on the same figure, how each important feature is matching labels from the point of view of the model. More precisely, this object acts a conditional importance. First a label is selected. Then importance is assessed conditionally. Hence features take influence, once the class is known (or predicted), meaning that there possibly exists more influential features (that lead to choose a class rather than another). In regression, the plotted object is called 'Dependence on most important predictors'.\cr
  
  5 - 'Importance( )' function allows a deeper analysis, depending on options and functions like \code{\link{partialDependenceBetweenPredictors}}, \code{\link{partialImportance}} or \code{\link{clusterAnalysis}}.\cr
  
  6 - To neutralize bias induced by correlated variables, it is strongly recommended :\cr
	a - to use \code{\link{randomUniformForest}} with 'mtry = 1', e.g. as a purely random forest (nodes are built and chosen completely at random)\cr
	b - then, to compute importance and use the resulted object as a benchmark of Variable Importance. If features have influence in the purely random case, then they also will have one in the weaker random case... The reciprocal is not valid.\cr
 
  Summarized, importance object tells (or leads to, using others functions) "which, how, when and where" features affect the response and interact each other. Note that importance strongly depends of hyper-parameters of the model (mtry, depth, nodesize, ...) and one has to take care. Categorical features are treated like numeric ones (default option) except if 'categorical' option is specified in \code{\link{randomUniformForest}}. 'categorical' option leads to a better assessment but might hurt accuracy. One possible reason is that the cross-entropy function (optimization criterion) would be more efficient for variables with high variability.
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

\references{
 Ciss, S., 2015b. Variable Importance in Random Uniform Forests. hal-1104751
 https://hal.archives-ouvertes.fr/hal-01104751
}

\author{
  Saip  CISS \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{partialDependenceOverResponses}}, \code{\link{partialDependenceBetweenPredictors}},
 \code{\link{partialImportance}},  \code{\link{clusterAnalysis}}
}

\examples{
## not run 
## 1 - Synthetic data:
## generate data

# set.seed(2014)
# n = 1000; p = 100 
# X = simulationData(n, p)
# X = fillVariablesNames(X)
# epsilon1 = runif(n,-1,1); epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]
# Y = as.factor(ifelse(rule > mean(rule), "+","-"))

## X1, X2, X3, X4 are the most important ones since they generate the labels
## Other ones are noise.

## run model
# synth.model.ruf = randomUniformForest(X, as.factor(Y))

## a - summary of the model provides (by default) global variable importance 
## (most predictive and most discriminant ones (in the case of classification)
# summary(synth.model.ruf)

## b - get details of importance : local variable importance and partial dependencies
## we choose 'maxInteractions' between covariates to be equal to 2
# synth.importance.ruf = importance(synth.model.ruf, Xtest = X, maxInteractions = 2)

## show interactions, variable importance based on interactions and variable importance
## conditionally to each label
# synth.importance.ruf

## c - plot details : 
## - adds partial dependence, showing for which values of a variable, the class is changing,
## - displays interactions, showing how influential variables are covering all the possible
## interactions
## - displays variable importance based on interactions, showing their relative influence 
## - displays variable importance over labels, showing how each variable is influencing a class
## and how it is discriminant in the separation between two or more classes.

# plot(synth.importance.ruf, Xtest = X, nLocalFeatures = 6)

## d - complete analysis with other tools, for example clusterAnalysis() :
## synth.Analysis.ruf = clusterAnalysis(synth.importance.ruf, X, components = 3, 
## clusteredObject = synth.model.ruf, OOB = TRUE)

## or clusteringObservations() : table and plot
# synth.Analysis2.ruf = clusteringObservations(synth.model.ruf, X, OOB = TRUE,  
# importanceObject = synth.importance.ruf) 

## 2 - Importance for Classification and Regression (with formula)
####  Classification

# data(iris)
# iris.ruf <- randomUniformForest(Species ~ ., data = iris, threads = 1)

## global importance : 
# summary(iris.ruf)

# iris.ruf.importance <- importance(iris.ruf, Xtest = iris, threads = 1)

## get importance summary
# iris.ruf.importance

## visualizing all in one
# plot(iris.ruf.importance, Xtest = iris)

#### Regression

# data(airquality)
# airquality.data = airquality

## impute NA
# airquality.NAimputed <- fillNA2.randomUniformForest(airquality.data)

## compute model
# ozone.ruf <- randomUniformForest(Ozone ~ ., data = airquality.NAimputed, threads = 1)

# ozone.ruf
# summary(ozone.ruf)
# ozone.ruf.importance <- importance(ozone.ruf, Xtest = airquality.NAimputed, threads = 1)

## visualization: in case of formula, 'formulaInput' is needed for the 'plot' method
# plot(ozone.ruf.importance, Xtest = airquality.NAimputed, formulaInput = ozone.ruf$formula)

## 3 - Importance for Classification and Regression without formula (more usual and recommended) 

#### Classification: "car evaluation" data 
## (http://archive.ics.uci.edu/ml/datasets/Car+Evaluation)

# data(carEvaluation)
# car.data <- carEvaluation
# n <- nrow(car.data)
# p <- ncol(car.data)

# trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

## training examples
# car.data.train <- car.data[trainTestIdx == 1, -p]
# car.class.train <- as.factor(car.data[trainTestIdx == 1, p])

## test data
# car.data.test <- car.data[trainTestIdx == 2, -p]
# car.class.test <- as.factor(car.data[trainTestIdx == 2, p])

## compute model : train then test in the same function
## option 'categorical' may be used for categorical variables and should be consistent
## with Variable Importance but may lead to less accuracy.

# car.ruf <- randomUniformForest(car.data.train, car.class.train,
# xtest = car.data.test, ytest = car.class.test, categorical = "all")
# car.ruf

## global importance: note that 'safety' does not appear to be an important feature 
## in the barplot but in the table, it is the most important of unacceptable (unacc) cars.
# summary(car.ruf)

## interactions and local importance tell most of the story...
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

## more interactions
# Concrete.ruf.importance <- importance(Concrete.ruf, Xtest = Concrete.data.test,
# maxInteractions = 6, predObject = Concrete.ruf.pred, threads = 1)

## or more features to plot
# plot(Concrete.ruf.importance, nLocalFeatures = 7, Xtest = Concrete.data.test)
}

\concept{variable importance}
\concept{variable selection}