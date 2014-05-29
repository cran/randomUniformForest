\name{getTree.randomUniformForest} 
\alias{getTree.randomUniformForest}
\alias{getTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract a tree from a forest}

\description{
 get the structure of a tree from a randomUniformForest object.
}

\usage{
getTree.randomUniformForest(object, whichTree, labelVar = TRUE)
}

\arguments{
   \item{object}{
    an object of class randomUniformForest.
}
 \item{whichTree}{
    which tree have to be printed ?
}
 \item{labelVar}{
   if TRUE, names of the variables will be printed, rather than their positions.
}
}

\value{
a data frame representing the raw tree structure
}

\note{
Trees in randomUniformForest package have almost the same structure than trees in the randomForest package.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
data(iris)
iris.ruf <- randomUniformForest(Species ~ ., data = iris, 
threads = 1, ntree = 20, BreimanBounds = FALSE) 

# get the 10th tree
OneTree <- getTree.randomUniformForest(iris.ruf, 10)
}