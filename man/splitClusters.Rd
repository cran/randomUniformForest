\name{splitClusters}
\alias{splitClusters}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Split a cluster on the fly
}
\description{
   Given one (or many) cluster(s), the function splits it in two new clusters.
}
\usage{
splitClusters(object, whichOnes, seed = 2014, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	    an object of class unsupervised.
	}
	\item{whichOnes}{
	   a vector defining which cluster(s), given by its label, has to be split ?
	}
	\item{seed}{
	  see \code{\link{unsupervised.randomUniformForest}}.
	}
	\item{\dots}{
	  not currently used.
	}
}

\value{
  An object of class unsupervised. 
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{mergeClusters}}, \code{\link{modifyClusters}}  
}

\examples{
## not run
## load iris data
# data(iris)

## run unsupervised modelling, removing labels and committing 2 clusters
# iris.uruf = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2, 
# threads = 1, clusters = 2)

## view a summary
# iris.uruf

## plot clusters 
# plot(iris.uruf)

## split the cluster which has the highest count (cluster 1, in our case)
# iris.urufSplit = splitClusters(iris.uruf, 1)

## assess fitting comparing average Silhouette before and after
# iris.urufSplit

## plot to see the new clusters
# plot(iris.urufSplit) 
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}