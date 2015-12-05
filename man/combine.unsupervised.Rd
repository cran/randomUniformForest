\name{combine.unsupervised}
\alias{combine.unsupervised}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Combine Unsupervised Learning objects
}
\description{
   Combine unsupervised learning objects in order to achieve incremental learning. Only the MDS points are used then standardized before calling a clustering algorithm on all. Note that the function is currently highly experimental
   with, currently, a lack of applications.
}
\usage{
combine.unsupervised(\dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{\dots}{
	  (vector of) objects of class unsupervised, coming from \code{\link{unsupervised.randomUniformForest}}, that needs to be combined.
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
 \code{\link{update.unsupervised}}, \code{\link{modifyClusters}}, \code{\link{mergeClusters}}, 
 \code{\link{clusteringObservations}}, \code{\link{as.supervised}} 
}

\examples{
## not run
## Wine Quality Data Set
## http://archive.ics.uci.edu/ml/datasets/Wine+Quality

# data(wineQualityRed)
# X = wineQualityRed[, -ncol(wineQualityRed)]

## 1 - run unsupervised analysis on the first half of dataset 

# subset.1 = 1:floor(nrow(X)/2)
# wineQualityRed.model.1 = unsupervised.randomUniformForest(X, subset = subset.1, 
# baseModel = "proximityThenDistance")

## assess roughly the model and visualize
#  wineQualityRed.model.1

# plot(wineQualityRed.model.1)

## 2 - run unsupervised analysis on the second half of dataset 
# wineQualityRed.model.2 = unsupervised.randomUniformForest(X, subset = -subset.1)

## 3 - combine
# wineQualityRed.combinedModel = 
# combine.unsupervised(wineQualityRed.model.1, wineQualityRed.model.2)

## visualize and plot
# wineQualityRed.combinedModel
# plot(wineQualityRed.combinedModel)

# compare with the full data
# wineQualityRed.model = unsupervised.randomUniformForest(X, 
# baseModel = "proximityThenDistance")
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}