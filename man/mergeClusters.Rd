\name{mergeClusters}
\alias{mergeClusters}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Merge two arbitrary, but adjacent, clusters
}
\description{
   merge, on the flight, two adjacent clusters in order to allow better clustering scheme, if needed, and to avoid new computation of the unsupervised mode.
}
\usage{
mergeClusters(object, whichOnes)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	    an object of class unsupervised.
	}
	\item{whichOnes}{
	   which clusters are needed to be merged ? Should be two adjacent ones.
	}
}

\value{
  An object of class unsupervised. 
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{modifyClusters}} 
}

\examples{
## not run

## Crabs data 
# data(crabs, package = "MASS")

## and more about
# ?crabs

## model : commit 4 clusters
# crabs.rufUnsupervised = unsupervised.randomUniformForest(crabs, 
# categoricalvariablesidx = "all", nodesize = 5, threads = 1, clusters = 4)

## visualize clusters and merge adjacent clusters
# plot(crabs.rufUnsupervised)

## we can first merge clusters 1 and 4
## note that clusters may change if run again
# crabs.rufUnsupervisedNew = mergeClusters(crabs.rufUnsupervised, c(1,4))

## one can assess the fitting, comparing old and new model
# crabs.rufUnsupervised
# crabs.rufUnsupervisedNew

## visualize new model
# plot(crabs.rufUnsupervisedNew)

## merge new clusters 1 and 2 and look if it will be better
# crabs.rufUnsupervisedNewest = mergeClusters(crabs.rufUnsupervisedNew, c(1,2))

# crabs.rufUnsupervisedNewest
# plot(crabs.rufUnsupervisedNewest)

## NOTE : mergeClusters() provides choice on how to choose and assess clusters
## using simply visualization. 
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}