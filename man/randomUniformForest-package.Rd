\name{randomUniformForest-package}
\alias{randomUniformForest-package}
\docType{package}
\title{
Random Uniform Forests for Classification and Regression
}
\description{
Ensemble model, for classification and regression, based on a forest of of unpruned and randomized binary trees. Each tree is grown by sampling, with replacement, a set of variables at each node. Each cut-point is generated randomly, according to the Uniform law on the support of each candidate variable. Optimal random node is, then, choosen by maximizing information gain (classification) or minimizing 'L2' (or 'L1') distance (regression). Data are either bootstrapped or subsampled for each tree. Random Uniform Forests are aimed to lower correlation between trees, to offer more details about variable importance and selection and to allow native incremental learning.
}
\details{
\tabular{ll}{
Package: \tab randomUniformForest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-12\cr
License: \tab BSD_3_clause\cr
}
}
\references{
 Biau, G., Devroye, L., Lugosi, G., 2008. \emph{Consistency of random forests and other averaging classifiers}. The Journal of Machine Learning Research 9, 2015-2033.
 
 Bousquet, O., Boucheron, S., Lugosi, G., 2004. \emph{Introduction to Statistical Learning Theory}, in: Bousquet, O., Luxburg, U. von, Ratsch, G. (Eds.), Advanced Lectures on Machine Learning, Lecture Notes in Computer Science. Springer Berlin Heidelberg, pp. 169-207.
 
 Breiman, L, 1996. \emph{Heuristics of instability and stabilization in model selection}. The Annals of Statistics 24, no. 6, 2350--2383.
 
 Breiman, L., 1996. \emph{Bagging predictors}. Machine learning 24, 123-140.
 
 Breiman, L., 2001. \emph{Random Forests}, Machine Learning 45(1), 5-32.
 
 Breiman, L., Friedman, J.H., Olshen, R.A., Stone, C., 1984. \emph{Classification and Regression Trees}. New York: Chapman and Hall. 
 
 Ciss, S., 2014. PhD thesis: \emph{Random uniform forests and irregularity detection in social security contributions}. Submitted. 
 
 Devroye, L., Gyorfi, L., Lugosi, G., 1996. \emph{A probabilistic theory of pattern recognition}. New York: Springer. 
 
 Efron, B., 1979. \emph{Bootstrap Methods: Another Look at the Jackknife}. The Annals of Statistics 7, 1-26. 
 
 Hastie, T., Tibshirani, R., Friedman, J.J.H., 2001. \emph{The elements of statistical learning}. New York: Springer. 
 
 Ho, T.K., 1998. \emph{The random subspace method for constructing decision forests}. IEEE Transactions on Pattern Analysis and Machine Intelligence 20, 832-844. 
 
 Lin, Y., Jeon, Y., 2002. \emph{Random Forests and Adaptive Nearest Neighbors}. Journal of the American Statistical Association 101-474.
 
 Vapnik, V.N., 1995. \emph{The nature of statistical learning theory}. Springer-Verlag New York, Inc., New York, NY, USA.
}
\author{
Saip Ciss

Maintainer: Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\keyword{ 
randomUniformForest
rUniformForest.grow
rUniformForest.big
rUniformForest.combine
predict.randomUniformForest
importance.randomUniformForest
partialDependenceBetweenPredictors
partialDependenceOverResponses
partialImportance
postProcessingVotes
bCI
roc.curve
rm.trees
fillNA2.randomUniformForest
}