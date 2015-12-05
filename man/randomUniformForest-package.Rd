\name{randomUniformForest-package}
\alias{randomUniformForest-package}
\docType{package}
\title{
Random Uniform Forests for Classification, Regression and Unsupervised Learning
}
\description{
Ensemble model for classification, regression and unsupervised learning, based on a forest of unpruned and randomized binary decision trees. Unlike \emph{Breiman's Random Forests}, each tree is grown by sampling, \emph{with replacement}, a set of variables at each node. Each cut-point is generated randomly, according to the \emph{continuous Uniform distribution between two random points of each candidate variable or using its whole current support}. Optimal random node is, then, selected among many full random nodes by maximizing information gain (classification) or minimizing 'L2' (or 'L1') distance (regression). Unlike \emph{Extremely Randomized Trees}, data are either \emph{bootstrapped or sub-sampled for each tree}. Random Uniform Forests are aimed to lower correlation between trees, to offer a deep analysis of variable importance and to allow native distributed and incremental learning. The unsupervised mode introduces clustering and dimension reduction, using a three-layer engine (dissimilarity matrix, Multidimensional Scaling and k-means or hierarchical clustering). 
}
\details{
\tabular{ll}{
Package: \tab randomUniformForest\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2015-01-05\cr
License: \tab BSD_3_clause\cr
}

Installation:  install.packages("randomUniformForest")\cr
Usage:         library(randomUniformForest)\cr


}
\references{
 Amit, Y., Geman, D., 1997. Shape Quantization and Recognition with Randomized Trees. \emph{Neural Computation} 9, 1545-1588. 
 
 Biau, G., Devroye, L., Lugosi, G., 2008. Consistency of random forests and other averaging classifiers. \emph{The Journal of Machine Learning Research} 9, 2015-2033.
 
 Bousquet, O., Boucheron, S., Lugosi, G., 2004. \emph{Introduction to Statistical Learning Theory}, in: Bousquet, O., Luxburg, U. von, Ratsch, G. (Eds.), Advanced Lectures on Machine Learning, Lecture Notes in Computer Science. Springer Berlin Heidelberg, pp. 169-207.
 
 Breiman, L, 1996. Heuristics of instability and stabilization in model selection. \emph{The Annals of Statistics} 24, no. 6, 2350-2383.
 
 Breiman, L., 1996. Bagging predictors. \emph{Machine learning} 24, 123-140.
 
 Breiman, L., 2001. Statistical Modeling: The Two Cultures (with comments and a rejoinder by the author). \emph{Statistical Science} 16, no. 3, 199-231.
 
 Breiman, L., 2001. Random Forests, \emph{Machine Learning} 45(1), 5-32.
 
 Breiman, L., Friedman, J.H., Olshen, R.A., Stone, C., 1984. \emph{Classification and Regression Trees}. New York: Chapman and Hall. 
 
 Ciss, S., 2014. PhD thesis: \emph{Forets uniformement aleatoires et detection des irregularites aux cotisations sociales}. Universite Paris Ouest Nanterre, France. In french.\cr
 English title : \emph{Random Uniform Forests and Irregularity Detection in social Security contributions}.\cr
 Link : https://www.dropbox.com/s/q7hbgeafrdd8qtc/Saip_Ciss_These.pdf?dl=0
 
 Ciss, S., 2014a. Random Uniform Forests. Pre-print. 
 
 Ciss, S., 2014b. Variable Importance in Random Uniform Forests. Pre-print. 
 
 Cox, T. F., Cox, M. A. A., 2001. \emph{Multidimensional Scaling. Second edition}. Chapman and Hall. 
 
 Devroye, L., Gyorfi, L., Lugosi, G., 1996. \emph{A probabilistic theory of pattern recognition}. New York: Springer.
 
 Dietterich, T.G., 2000. \emph{Ensemble Methods in Machine Learning}, in : Multiple Classifier Systems, Lecture Notes in Computer Science.
 Springer Berlin Heidelberg, pp. 1-15. 
 
 Efron, B., 1979. Bootstrap Methods: Another Look at the Jackknife. \emph{The Annals of Statistics} 7, 1-26.
 
 Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika} 53, 325-328.
 
 Hastie, T., Tibshirani, R., Friedman, J.J.H., 2001. \emph{The elements of statistical learning}. New York: Springer. 
 
 Ho, T.K., 1998. The random subspace method for constructing decision forests. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence} 20, 832-844. 
 
 Lin, Y., Jeon, Y., 2002. Random Forests and Adaptive Nearest Neighbors. \emph{Journal of the American Statistical Association} 101-474.
 
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
unsupervised.randomUniformForest
as.supervised
clusteringObservations
}