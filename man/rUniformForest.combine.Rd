\name{rUniformForest.combine}
\alias{rUniformForest.combine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Incremental learning for random Uniform Forests
}
\description{
  Combine many random uniform forests on same or different data and/or same and different parameters to obtain one ensemble model. Function acts as the "Reduce" task of MapReduce paradigm. "Map" task has here no defined length, since rUniformForest.combine() is an ensemble of ensembles of any size, parameters or data that can be continuously merged and evaluated as a single model. Native incremental learning is then achieved adding and removing, but never modifying, trees.
}
\usage{
rUniformForest.combine(...)
}

\arguments{
  \item{\dots}{
     vector enumerating objects of class randomUniformForest
}
}

\details{
  Incremental learning is assumed as a sub-sampling process if distribution is not switching. Data arrive continuously and random Uniform Forests
  are trained and updated using chunks of the data. Computing time and memory are constant and only prediction time will increase depending on the size (number of trees) of the model. Hence "split-apply-combine" ideas will matter and chunks just have to be processed on all available workstations. For a shifting distribution, changes can happen in the current processed chunk or between one chunk and another. Then, one
  needs to choose how the model has to be updated. Note that incremental learning is effective because cut-points are chosen randomly, according to the Uniform distribution on the support of each candidate variable. rUniformForest.combine() acts as "compute everywhere, combine in one place" and there is no need to split data here nor to apply the same model. Data are assumed to come from any source, continuously and by chunks, sharing only the same or some common features. Parameters in any model are allowed to change.
}


\value{
  an object of class randomUniformForest
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\seealso{
 \code{\link{rm.trees}}
}

\examples{
## not run
## Classification : synthetic data 
## get many forests and combine them

# n = 200;  p = 10
## Simulate 'p' gaussian vectors with random parameters between -10 and 10.
# X <- simulationData(n,p)

## Make a rule to create response vector
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]

# Y <- as.factor(ifelse(rule > mean(rule), 1, 0))

## create many subsamples
# manyCuts <- cut(1:n, 5, labels = FALSE)

## compute many different models 
# ruf1 <- randomUniformForest(X[which(manyCuts == 1),], Y[which(manyCuts == 1)], 
# ntree = 30, mtry = 2, BreimanBounds = FALSE, importance = FALSE, threads = 1)

# ruf2 <- randomUniformForest(X[which(manyCuts == 2),], Y[which(manyCuts == 2)],
# ntree = 40, nodesize = 10, BreimanBounds = FALSE, importance = FALSE, threads = 1)

# ruf3 <- randomUniformForest(X[which(manyCuts == 3),], Y[which(manyCuts == 3)],
# ntree = 50, bagging = TRUE, BreimanBounds = FALSE, importance = FALSE, threads = 1)

## combine them in one ensemble 
# ruf.combined <- rUniformForest.combine(ruf1, ruf2, ruf3)
# ruf.combined

## or 
# rufObjects <- list(ruf1, ruf2, ruf3)
# ruf.combined <- rUniformForest.combine(rufObjects)

## remove 10 older trees
# ruf.combined <- rm.trees(ruf.combined, method = "oldest", howMany = 10)
# ruf.combined

## compute a new model and update
# ruf4 <- randomUniformForest(X[which(manyCuts == 4),], Y[which(manyCuts == 4)],
# ntree = 40, rebalancedsampling = TRUE, BreimanBounds = FALSE, threads = 1)
# ruf.updateCombined <- rUniformForest.combine(ruf.combined, ruf4)
# ruf.updateCombined

# ruf.pred <- predict(ruf.updateCombined, X[which(manyCuts == 5),])

## confusion matrix
# table(ruf.pred, Y[which(manyCuts == 5)])

## comparison with offline learning (that will always be better, 
## except in the case of a shifting distribution)
# ruf.offline <- randomUniformForest(X[which(manyCuts != 5),], Y[which(manyCuts !=5)], 
# threads = 1, ntree = 150, BreimanBounds = FALSE)
# ruf.offline.pred <- predict(ruf.offline, X[which(manyCuts == 5),])

## confusion matrix
# table(ruf.offline.pred, Y[which(manyCuts == 5)])
}

\keyword{incremental}
\keyword{learning}