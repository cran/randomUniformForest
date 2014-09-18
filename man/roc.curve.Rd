\name{roc.curve} 
\alias{roc.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ROC and precision-recall curves for random Uniform Forests
}

\description{
 plot ROC and precision-recall curves for objects of class randomUniformForest and compute F-beta score. It also works for any other model that provides predicted labels (but only for ROC curve).
}

\usage{
roc.curve(X, Y, classes,
	positive = classes[2],
	ranking.threshold = 0,
	ranking.values = 0,
	falseDiscoveryRate = FALSE,
	plotting = TRUE,
	printValues = TRUE,
	Beta = 1)
}

\arguments{
   \item{X}{
     a vector (or a factor) of predictions (with two classes) or an object of class randomUniformForest (with OOB option enabled).
}
 \item{Y}{
     a vector of numeric (integer) responses, or a factor, with two classes.
}
 \item{classes}{
    a vector (or a factor) of values that designate the class.
}
 \item{positive}{
    convention for the positive class (e.g. the minority class).
}
\item{ranking.threshold}{
    option currently implemented but not fully tested.
}
\item{ranking.values}{
    option currently implemented but not fully tested.
}
 \item{falseDiscoveryRate}{
    if TRUE, precision-recall curve is plotted. if FALSE, default value, ROC curve is plotted.
}
\item{plotting}{
    plotting the ROC curve ?
}
\item{printValues}{
   display values to screen ?
}
\item{Beta}{
    'beta' value for F-beta score.
}
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\seealso{
\code{\link{importance.randomUniformForest}}
}

\examples{
## Classification : "breast cancer" data 
# http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)

data(breastCancer)
breastCancer.data <- breastCancer

# remove ID (first column) and divide data in train and test set
breastCancer.data = breastCancer.data[,-1]

n <- nrow(breastCancer.data)
p <- ncol(breastCancer.data)

trainTestIdx <- cut(sample(1:n, n), 2, labels= FALSE)

# train examples
breastCancer.data.train <- breastCancer.data[trainTestIdx == 1, -p]
breastCancer.class.train <- as.factor(breastCancer.data[trainTestIdx == 1, p])

# rename class in benign (class 2) and malignant (class 4) to have a better view
levels(breastCancer.class.train) = c("benign", "malignant")

# test data
breastCancer.data.test <- breastCancer.data[trainTestIdx == 2, -p]
breastCancer.class.test <- as.factor(breastCancer.data[trainTestIdx == 2, p])

levels(breastCancer.class.test) = c("benign", "malignant")

# compute model : train then test in the same function and assign class weights 
# to better match the distribution (OOB errors and Breiman's bounds can help to choose weights)
# Note that in this case 'recall' (or sensitivity) is the objective, 
# e.g. match all possible cases of malignant tumours even if false positive rate increase 
#(in this latter case, further steps will reveal the truth). If malignant tumour is not detected, 
# then diagnosis error is, by far, more critical.
breastCancer.ruf  <- randomUniformForest(breastCancer.data.train, breastCancer.class.train, 
xtest = breastCancer.data.test, ytest = breastCancer.class.test, 
classwt = c(1, 3.5), threads = 2, ntree = 40, BreimanBounds = FALSE)

# get a summary of model
breastCancer.ruf

## plot ROC Curve for test data
# roc.curve(breastCancer.ruf, breastCancer.class.test, levels(breastCancer.class.test)) 

## plot precision-recall curve for test data
# roc.curve(breastCancer.ruf, breastCancer.class.test, levels(breastCancer.class.test),
# falseDiscoveryRate = TRUE)

## associate cut-off and purely random forest as an alternative to find maximum malignant cases 
## with a low false positive rate. 
## 'classcutoff' option is a bit tricky. Let's take the example we will use below.
## classcutoff = c("benign", 1.25) means that number of votes for class "benign" 
## will be weighted by 0.4 (= Cte/1.25, where Cte = 0.5) for each response.
## Hence "benign" will never have majority unless it has 2.5 (1.25/0.5) times more votes 
## than "malignant" class and all votes sum to total number of trees.
# breastCancer.cutOff.ruf <- randomUniformForest(breastCancer.data.train, breastCancer.class.train, 
# xtest = breastCancer.data.test, ytest = breastCancer.class.test, classcutoff = c("benign", 1.25), 
# randomfeature = TRUE, ntree = 50, threads = 2, BreimanBounds = FALSE)

# roc.curve(breastCancer.cutOff.ruf, breastCancer.class.test, levels(breastCancer.class.test)) 
# roc.curve(breastCancer.cutOff.ruf, breastCancer.class.test, levels(breastCancer.class.test),
# falseDiscoveryRate = TRUE)

## evaluate OOB data, when there is no test set
# breastCancer.ruf <- randomUniformForest(breastCancer.data.train, breastCancer.class.train,
# classwt = c(1, 3.5), threads = 2)
# roc.curve(breastCancer.ruf, breastCancer.class.train, levels(breastCancer.class.train)) 
}