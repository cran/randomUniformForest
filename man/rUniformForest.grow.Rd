\name{rUniformForest.grow}
\alias{rUniformForest.grow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add trees to a random Uniform Forest 
}
\description{
  add more trees to the ensemble model
}
\usage{
rUniformForest.grow(object, X, Y = NULL, ntree = 100, threads = "auto")
}

\arguments{
  \item{object}{
     an object of class randomUniformForest.
}
  \item{X}{
     current matrix or data frame used to compute 'object' model. One can also choose another matrix with same features and appropriate training responses 	(see below). All parameters will still remain. To combine heterogeneous models, use randomUniformForest.combine() function.
}
  \item{Y}{
     vector of training reponses, if one wants to add trees with another matrix (with same features, but different observations). Otherwise, let it to NULL if same model has to be computed on same data.
}
  \item{ntree}{
     number of trees to add to the model.
}
  \item{threads}{
  compute model in parallel for computers with many cores. Default value is "auto", letting model running
	on all logical cores minus 1. User can set 'threads' to any values >= 1, depending on the number of cores (including logical).
}
}

\details{
  rUniformForest.grow allows both to add new trees or new model (by adding trees on a new matrix and training responses) built on
the same bases than the former. Note that with formula, only new trees can be added, not new model.
}
\value{
 an object of class randomUniformForest, containing new and old trees.
}

\author{
 Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{rUniformForest.big}}, \code{\link{rUniformForest.combine}}
}

\examples{
# data(iris)
# iris.rf <- randomUniformForest(Species ~ ., iris, ntree = 50, threads = 1, BreimanBounds = FALSE)
# iris.rf <- rUniformForest.grow(iris.rf, iris, ntree = 20, threads = 1)
# iris.rf
}