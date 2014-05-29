\name{simulationData} 
\alias{simulationData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of gaussian vector
}

\description{
 Simulate a gaussian vector with 'p' independent components of length 'n'. Parameters of each component are uniformly random and are taken between -10 and 10, with (absolute) standard deviation equals mean.
}

\usage{
simulationData(n, p, distrib = rnorm, colinearity = FALSE)
}

\arguments{
   \item{n}{
    number of observations to draw.
}
 \item{p}{
    number of variables to draw.
}
 \item{distrib}{
   distribution to use. Currently only gaussian one is supported.
}
 \item{colinearity}{
    not currently used.
}
}

\value{
a matrix with 'n' observations and 'p' variables.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
X <- simulationData(100,10)
summary(X)
}