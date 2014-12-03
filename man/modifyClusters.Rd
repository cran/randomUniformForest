\name{modifyClusters}
\alias{modifyClusters}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Change number of clusters (and clusters shape) on the fly
}
\description{
   Modify on the fly the number of clusters in the unsupervised mode of Random Uniform Forests. Once the unsupervised model is built, one can change the clustering scheme, adding or merging clusters without computing again the model. Average silhouette coefficient or inter-classes variance are automatically adjusted in order to assess the new scheme.
}
\usage{
modifyClusters(object, decreaseBy = NULL, increaseBy = NULL, seed = 2014)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	    an object of class unsupervised.
	}
	\item{decreaseBy}{
	   decrease the current number of clusters by the value of 'decreaseBy'.
	}
	\item{increaseBy}{
	  decrease the current number of clusters by the value of 'increaseBy'.
	}
	\item{seed}{
	  see \code{\link{unsupervised.randomUniformForest}}.
	}
}

\value{
  An object of class unsupervised. 
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{mergeClusters}} 
}

\examples{
## not run
## load iris data
# data(iris)

## run unsupervised modelling, removing labels and committing 4 clusters
# iris.rufUnsupervised = unsupervised.randomUniformForest(iris[,-5], threads = 1, clusters = 4)

## view a summary
# iris.rufUnsupervised

## plot clusters 
# plot(iris.rufUnsupervised)

## modify clusters, decreasing them by  one
# iris.rufUnsupervisedNew = modifyClusters(iris.rufUnsupervised, decreaseBy = 1)

## assess fitting comparing average Silhouette before and after
# iris.rufUnsupervisedNew

## plot to see the new clusters
# plot(iris.rufUnsupervisedNew) 
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}