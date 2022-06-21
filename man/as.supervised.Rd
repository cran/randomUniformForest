\name{as.supervised}
\alias{as.supervised}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Conversion of an unsupervised model into a supervised one
}
\description{
   Turn an unsupervised object of class unsupervised into a supervised one of class RandomUniformForest, allowing prediction of next unlabelled datasets, full analysis of variable importance in the unsupervised case and incremental unsupervised learning.
}
\usage{
as.supervised(object, X, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	    an object of class unsupervised 
	}
	\item{X}{
	   the dataset previously learned in the unsupervised mode.
	}
	\item{\dots}{
	  allow all options of \code{\link{randomUniformForest}} to be passed to the as.supervised() function, e.g. ntree, mtry,  nodesize, ...
	}
}

\details{
 The function get clusters labels and send them to the randomUniformForest classifier in order to be learnt with the data. The resulting object can be used for any coming dataset for fast prediction or clustering. Note that main argument of the method is that the model which generates the unsupervised model also generates the supervised one. The process to clustering, in both unsupervised and supervised cases, is filtered by the dissimilarity matrix in conjunction with the MDS function. 
}

\value{
  An object of class randomUniformForest. See \code{\link{randomUniformForest}}.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{modifyClusters}}, \code{\link{mergeClusters}} \code{\link{clusteringObservations}}
}

\examples{
# see unsupervised.randomUniformForest() function.
}

\keyword{unsupervised}
\keyword{learning}
\keyword{supervised}