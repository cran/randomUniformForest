\name{randomUniformForest-package}
\alias{randomUniformForest-package}
\docType{package}
\title{
Random Uniform Forests for Classification, Regression and Unsupervised Learning
}
\description{
Ensemble model for classification, regression and unsupervised learning, based on a forest of unpruned and randomized binary decision trees. Unlike \emph{Breiman's Random Forests}, each tree is grown by sampling, \emph{with replacement}, a set of variables before splitting each node. Each cut-point is generated randomly, according to the \emph{continuous Uniform distribution between two random points of each candidate variable or using its whole current support}. Optimal random node is, then, selected among many full random ones by maximizing Information Gain (classification) or minimizing a distance (regression), 'L2' (or 'L1'). Unlike \emph{Extremely Randomized Trees}, data are either \emph{bootstrapped or sub-sampled for each tree}. From the theoretical side, Random Uniform Forests are aimed to lower correlation between trees and to offer a deep analysis of variable importance. The unsupervised mode introduces clustering and dimension reduction, using a three-layer engine: dissimilarity matrix, Multidimensional Scaling (or spectral decomposition) and k-means (or hierarchical clustering). From the practical side, Random Uniform Forests are designed to provide a complete analysis of (un)supervised problems and to allow native distributed and incremental learning.
}
\details{
\tabular{ll}{
Package: \tab randomUniformForest\cr
Type: \tab Package\cr
Version: \tab 1.1.5\cr
Date: \tab 2015-02-16\cr
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
 
 Ciss, S., 2015a. Random Uniform Forests. Preprint. hal-01104340.
 
 Ciss, S., 2015b. Variable Importance in Random Uniform Forests. Preprint. hal-01104751. 
 
 Ciss, S., 2015c. Generalization Error and Out-of-bag Bounds in Random (Uniform) Forests. Preprint. hal-01110524.  
 
 Cox, T. F., Cox, M. A. A., 2001. \emph{Multidimensional Scaling. Second edition}. Chapman and Hall. 
 
 Devroye, L., Gyorfi, L., Lugosi, G., 1996. \emph{A probabilistic theory of pattern recognition}. New York: Springer.
 
 Dietterich, T.G., 2000. \emph{Ensemble Methods in Machine Learning}, in : Multiple Classifier Systems, Lecture Notes in Computer Science. Springer Berlin Heidelberg, pp. 1-15. 
 
 Efron, B., 1979. Bootstrap Methods: Another Look at the Jackknife. \emph{The Annals of Statistics} 7, 1-26.
 
 Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika} 53, 325-328.
 
 Gyorfi, L., 2002. \emph{A distribution-free theory of nonparametric regression}. Springer Science & Business Media.
 
 Hastie, T., Tibshirani, R., Friedman, J.J.H., 2001. \emph{The elements of statistical learning}. New York: Springer. 
 
 Ho, T.K., 1998. The random subspace method for constructing decision forests. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence} 20, 832-844. 
 
 Lin, Y., Jeon, Y., 2002. Random Forests and Adaptive Nearest Neighbors. \emph{Journal of the American Statistical Association} 101-474.
 
 Ng, A. Y., Jordan, M. I., Weiss, Y., 2002. On spectral clustering: Analysis and an algorithm. \emph{Advances in neural information processing systems}, 2, 849-856.
 
 Vapnik, V.N., 1995. \emph{The nature of statistical learning theory}. Springer-Verlag New York, Inc., New York, NY, USA.
}
\author{
Saip Ciss

Maintainer: Saip Ciss \email{saip.ciss@wanadoo.fr}
}
\examples{
## Presenting some of the core functions
## Not run

## 1 - Classification: iris data set (assess whole data)

## load data (included in R):
# data(iris)
# XY = iris
# p = ncol(XY)
# X = XY[,-p]
# Y = XY[,p]

## Train a model : using formula
# iris.ruf = randomUniformForest(Species ~., XY, threads = 1)

## or using matrix
## iris.ruf = randomUniformForest(X, as.factor(Y), threads = 1)

## Assess model : Out-of-bag (OOB) evaluation
# iris.ruf

## Variable Importance : base assessment
# summary(iris.ruf)

## Variable Importance : deeper assessment (explain the modelling)
# iris.importance = importance(iris.ruf, Xtest = X, maxInteractions = p - 1)

## Visualize : details of Variable Importance 
## (tile windows vertically, using the R menu, to see all plots)
# plot(iris.importance, Xtest = X)

## Analyse : get an interpretation of the model results
# iris.ruf.analysis = clusterAnalysis(iris.importance, X, components = 3, 
# clusteredObject = iris.ruf, OOB = TRUE)

## Dimension reduction, clustering and visualization : OOB evaluation
# iris.clust.ruf = clusteringObservations(iris.ruf, X, importanceObject = iris.importance)

## 2 - Regression: Boston Housing (assess a test set)

## load data :
# install.packages("mlbench") ##if not installed
# data(BostonHousing, package = "mlbench")

# XY = BostonHousing
# p = ncol(XY)
# X = XY[,-p]
# Y = XY[,p]

## get random training and test sets :
## reproducibility :
# set.seed(2015)

# train_test = init_values(X, Y, sample.size = 1/2)
# Xtrain = train_test$xtrain
# Ytrain = train_test$ytrain
# Xtest = train_test$xtest
# Ytest = train_test$ytest

## Train a model : 
# boston.ruf = randomUniformForest(Xtrain, Ytrain)

## Assess (quickly) the model : 
# boston.ruf
# plot(boston.ruf)
# summary(boston.ruf)

## Predict the test set :
# boston.pred.ruf = predict(boston.ruf, Xtest)

## or predict quantiles
# boston.predQunatile_97.5.ruf = predict(boston.ruf, Xtest, type = "quantile", 
# whichQuantile = 0.975)

## or prediction intervals
# boston.predConfInt_95.ruf = predict(boston.ruf, Xtest, type = "confInt", conf = 0.95)

## Assess predictions :
# statsModel = model.stats(boston.pred.ruf, Ytest, regression = TRUE)

## Avoiding overfitting : under the i.i.d. assumption, OOB error
## is expected to be an upper bound of MSE. Convergence is first needed.
## Convergence needs low correlation between trees residuals.

# boston.ruf

## The easy way; reduce correlation(decreasing 'mtry' value) + post-processing 
# bostonNew.ruf = randomUniformForest(Xtrain, Ytrain, mtry = 4)

## (predict and) Post-process :
# bostonNew.predAll.ruf = predict(bostonNew.ruf, Xtest, type = "all")

# bostonNew.postProcessPred.ruf = postProcessingVotes(bostonNew.ruf, 
# predObject = bostonNew.predAll.ruf)

## Assess new predictions :
# statsModel = model.stats(bostonNew.postProcessPred.ruf, Ytest, regression = TRUE)

## Convergence : grow more trees
# bostonNew.moreTrees.ruf = rUniformForest.grow(bostonNew.ruf, Xtrain, ntree = 100)
}

\concept{ensemble methods}
\concept{supervised learning}
\concept{unsupervised learning}
\concept{classification}
\concept{regression}
\concept{Random Forests}
\concept{Random Uniform Forests}
\concept{clustering}
\concept{dimension reduction}
\concept{variable importance}