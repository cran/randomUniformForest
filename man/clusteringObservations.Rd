\name{clusteringObservations}
\alias{clusteringObservations}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cluster observations of a (supervised) randomUniformForest object
}
\description{
   Provide a clustering scheme for observations (in training or test sample) in a (supervised) randomUniformForest object. Observations are clustered using their labels or predictions. The rest of the process is the same than in the unsupervised case. clusteringObservations() may be needed when one wants to know how the classifier can have some troubles in the prediction task, especially when a link between features
   and decision rule has to be made.
}
\usage{
clusteringObservations(object, X,
OOB = FALSE,  
predObject = NULL, 
importanceObject = NULL, 
baseModel = c("proximity", "proximityThenDistance",  "importanceThenDistance"),
MDSmetric = c("metricMDS", "nonMetricMDS"),
endModel = c("MDS", "MDSkMeans", "MDShClust"),
\dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
		an object of class randomUniformForest.
	}
	\item{X}{
		the dataset to predict.
	}
	\item{OOB}{
		set it to TRUE, if one wants to assess the OOB classifier
	}
	\item{predObject}{
	    full prediction object expected for test set, if no predictions in 'object'.
	}
	\item{importanceObject}{
	    object of class importance, if 'baseModel' equals 'importanceThenDistance'.
	}
	\item{baseModel}{
	   see \code{\link{unsupervised.randomUniformForest}}.
	}
	\item{MDSmetric}{
	   see \code{\link{unsupervised.randomUniformForest}}.
	}
	\item{endModel}{
	   see \code{\link{unsupervised.randomUniformForest}}.
	}
	\item{\dots}{
	   others parameters of \code{\link{unsupervised.randomUniformForest}}.
	}
}

\details{
clusteringObservations() function leads to many modelling paradigms, depending on what is called.\cr
- To summarize, the first level is a clustering scheme of the learning task, e.g. 'OOB = FALSE'. In this case, the modelling links the learning with the true labels, providing a picture on how the classifier is "seeing" the training sample. In particular, this representation provides insights on how the classification task might be difficult.\cr
- The second level calls the OOB point of view, e.g. how the classifier is generalizing its current knowledge. In this case, one will simply visualize where errors will probably appear the most.\cr
- The third level involves the prediction task if 'predObject' is provided or if a test set is provided when doing the supervised task (option 'xtest' in \code{\link{randomUniformForest}}). In this case, the classifier will map its own predictions on the low dimensional space. 

In the three levels, there is at least one layer (the dissimilarity matrix), in addition, involved. If one overrides the default 'endModel = MDS' option, then more layers will be added, leading to a more complex object.
}

\value{
  An object of class unsupervised, embedding the supervised model.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}


\examples{
## not run

## Vehicle dataset : 
# data(Vehicle, package = "mlbench")
# XY = Vehicle
# Y = XY[,"Class"]
# X = XY[, -ncol(XY)]

## 1- learn supervised model and visualize only true labels 
# vehicle.ruf = randomUniformForest(X, as.factor(Y), ntree = 200)

## compute, for more details, importance
# vehicle.importance = importance(vehicle.ruf, Xtest = X)

## clustering Object  
## 'OOB = FALSE' ensures that OOB classifier will not be called
## One may note that class 'opel' and 'saab' are hard to separate
# vehicle.clusterObject = clusteringObservations(vehicle.ruf, X, OOB = FALSE,
# importanceObject = vehicle.importance, clusters = 4)

# vehicle.clusterObject 

## 2 - visualize the OOB point of view
# vehicle.clusterObjectOOB = clusteringObservations(vehicle.ruf, X, OOB = TRUE,
# importanceObject = vehicle.importance, clusters = 4)

## note that true labels are left to their numeric values in order to not modify
## the internal structure of the supervised model
# vehicle.clusterObjectOOB

## 3 - learn supervised model, predict and visualize predictions
## set training and test sets

# set.seed(2014)
# train_test = init_values(X, Y, sample.size = 1/2)
# X1 = train_test$xtrain
# Y1 = train_test$ytrain
# X2 = train_test$xtest
# Y2 = train_test$ytest

## learn and predict
# vehicle.new.ruf = randomUniformForest(X1, as.factor(Y1), xtest = X2, ntree = 200)

## compute importance for the test set
# vehicle.new.importance = importance(vehicle.new.ruf, Xtest = X2)

## clustering predictions : no satisfaction. Clustering seems to the reside in a higher dimension
# vehicle.new.clusterObject = clusteringObservations(vehicle.new.ruf, X2, 
# importanceObject = vehicle.new.importance, clusters = 4)

## 3.1 - ask help : let's kMeans talk, i.e. can the problem have a dimension as small as 2 ?
## note that the clustering process generates its own predictions and the labels of clusters 
## generated have, usually, no direct links with the true labels. 
## More precisely cluster 1 is usually not label 1 of the supervised case.
## For the same reasons, importanceObject must be computed in another manner

# vehicle.more.clusterObject = clusteringObservations(vehicle.new.ruf, X2, 
# endModel = "MDSkMeans", clusters = 4)

## 3.2 - Clustering is achieved, but what is the cost for the predictions ?
## compare with true test labels : first let's look to the supervised model

# stats.new = model.stats(vehicle.new.ruf$predictionObject$majority.vote, as.numeric(Y2))

## for the clustering model, pairing true test labels and cluster labels is done as a shortcut.
## Finally, even if clustering is achievable it does not necessarily correspond to the true labels

# clusterPredictions = mergeOutliers(vehicle.more.clusterObject)

# 'bus'
# table(clusterPredictions[which(as.numeric(Y2) == 1)])

# 'opel'
# table(clusterPredictions[which(as.numeric(Y2) == 2)])

# 'saab'
# table(clusterPredictions[which(as.numeric(Y2) == 3)])

# 'van'
# table(clusterPredictions[which(as.numeric(Y2) == 4)])
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}