\name{update.unsupervised}
\alias{update.unsupervised}
\alias{update}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Update Unsupervised Learning object
}
\description{
   Update unsupervised learning object with new data in order to achieve incremental learning. New MDS (or spectral) points are predicted with new data and learning of MDS (spectral) points of the former unsupervised object. Note that
   new data are expected to have the same distribution than previous ones.
}
\usage{
\method{update}{unsupervised}(object, X, 
	oldData = NULL, 
	mapAndReduce = FALSE, 
	updateModel = FALSE, 
	\dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	  (vector of) objects of class unsupervised, coming from \code{\link{unsupervised.randomUniformForest}}, that needs to be updated.
	}
	\item{X}{
	  new data frame or matrix with same variables than the ones used to build 'object'.
	}
	\item{oldData}{
	  former data frame or matrix needed, if 'object' was built from a small dataset (less than 10000 rows).
	}
	\item{mapAndReduce}{
	  if TRUE, 'X' will be learned by chunks that will be combined to build the forest classifier for the MDS points.
	}
	\item{updateModel}{
	   if TRUE, the resulted object will also embed a learned model of all the MDS points (hence learning its own predictions). Note that this option is recommended to be enabled since it will allow to assess the resulted 
	   model against, for example, any unsupervised algorithm, or itself, based on a subsample of 'X'.
	}
	\item{\dots}{
	  all hyper parameters of \code{\link{randomUniformForest}}.
	}
}


\value{
  An object of class unsupervised, which is a list with the following components:
  
  \item{proximityMatrix}{the resulted dissimilarity matrix.}
  \item{MDSModel}{the resulted Multidimensional scaling model.}
  \item{unsupervisedModel}{the resulted unsupervised model with clustered observations in unsupervisedModel$cluster.}
  \item{largeDataLearningModel}{ if the dataset is large, the resulted model that learned a sample of the MDS points, and predicted others points.}
  \item{gapStatistics}{if K-means algorithm has been called, the results of the gap statistic. Otherwise NULL.}
  \item{rUFObject}{Random Uniform Forests object.}
  \item{nbClusters}{Number of clusters found.}
  \item{params}{options of the model.}
}


\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{combineUnsupervised}}, \code{\link{modifyClusters}}, \code{\link{mergeClusters}}, 
 \code{\link{clusteringObservations}}, \code{\link{as.supervised}} 
}

\examples{
## not run
## Water Treatment Plant Data Set
## Data can be download at https://archive.ics.uci.edu/ml/datasets/Water+Treatment+Plant

# URL = "http://archive.ics.uci.edu/ml/machine-learning-databases/water-treatment/"
# dataset = "water-treatment.data"

# X = read.table(paste(URL, dataset, sep= ""), sep = ",")

## 1- Preprocessing
## first, look at the first column and format date
#  Dates = rm.string(as.character(X[,1]), "D-")
#  DatesAsStringTable = do.call(rbind, strsplit(Dates, "/"))
#  DatesasNumericTable = t(apply(DatesAsStringTable, 1, as.numeric))

##  Then, transform data as a pure R matrix and add new dates
#  XX = as.true.matrix(X)[,-1]
#  XX = cbind(DatesasNumericTable, XX)
#  colnames(XX)[1:3] = c("day", "month", "year")

# Look the new data
# head(XX)
# str(XX)

## and fill missing values,
# X.imputed = fillNA2.randomUniformForest(XX)

## 2 - run unsupervised analysis on the first half of dataset 
# subset.1 = 1:floor(nrow(X.imputed)/2)
# WaterTreatment.model.1 = unsupervised.randomUniformForest(X.imputed, subset = subset.1, 
# baseModel = "proximityThenDistance",  seed = 2014)

## assess roughly the model and visualize
#  WaterTreatment.model.1

## 3 - update model with the second half of dataset
# WaterTreatment.updated = update.unsupervised(WaterTreatment.model.1, 
# X.imputed[-subset.1,], oldData = X.imputed[subset.1,])

# WaterTreatment.updated

## view how MDS points have been learned :
## first component
# WaterTreatment.updated$largeDataLearningModel[[1]]

## second component
# WaterTreatment.updated$largeDataLearningModel[[2]]
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}