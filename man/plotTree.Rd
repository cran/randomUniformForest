\name{plotTree} 
\alias{plotTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a Random Uniform Decision Tree
}

\description{
 plot the tree structure, showing nodes, variables, cut-points and predictions.
}

\usage{
plotTree(treeStruct, 
rowNum = 1, 
height.increment = 1, 
maxDepth = 100, 
fullTree = FALSE, 
xlim = NULL, 
ylim= NULL, 
center = TRUE)
}

\arguments{
   \item{treeStruct}{
    a data frame same at the output of getTree.randomUniformForest() function.
}
 \item{rowNum}{
    internal option not currently used.
}
 \item{height.increment}{
   internal option not currently used.
}
 \item{maxDepth}{
   internal option not currently used.
}
\item{fullTree}{
    if TRUE, draw full tree, but if the tree is deep or large, one will not see all. if FALSE, only parts of the tree
	will be shown but visualization will be clearer. Note that using fullTree = TRUE, simultaneously with 'xlim' and 'ylim', one can see the whole tree (but not the details).
}

\item{xlim}{
	a vector of two values, giving the width of the tree to show. One can, then, explore tree structure more in details.
}
 \item{ylim}{
    a vector of two values, giving the depth of the tree to show. One can, then, explore tree structure more in details.
}
\item{center}{
    if TRUE, the tree will be centered.
}
}


\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
# not run
data(airquality)
ozone.ruf <- randomUniformForest(Ozone ~ ., data = airquality, 
ntree = 20, BreimanBounds = FALSE, threads = 1)

OneTree <- getTree.randomUniformForest(ozone.ruf, 10)

# plotTree(OneTree) ##only a part is visible

## full tree : 
# plotTree(OneTree, fullTree = TRUE, xlim = c(1,55), ylim = c(0, 11))
}