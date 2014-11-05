\name{randomUniformForest}
\alias{randomUniformForest}
\alias{randomUniformForest.default}
\alias{randomUniformForest.formula}
\alias{plot.randomUniformForest}
\alias{print.randomUniformForest}
\alias{summary.randomUniformForest}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Random Uniform Forests for Classification, Regression and Unsupervised Learning
}
\description{
  Ensemble model for classification, regression and unsupervised learning, based on a forest of unpruned and randomized binary decision trees. Unlike \emph{Breiman's Random Forests}, each tree is grown by sampling, \emph{with replacement}, a set of variables at each node. Each cut-point is generated randomly, according to the \emph{continuous Uniform distribution between two random points of each candidate variable}. Optimal random node is, then, selected among many full random nodes by maximizing information gain (classification) or minimizing 'L2' (or 'L1') distance (regression). Unlike \emph{Extremely Randomized Trees}, data are either bootstrapped or sub-sampled for each tree. Random Uniform Forests are aimed to lower correlation between trees, to offer a deep analysis of variable importance and to allow native distributed and incremental learning. The unsupervised mode introduces clustering and dimension reduction, using a three-layer engine (dissimilarity matrix, Multidimensional Scaling and k-means or hierarchical clustering). 
}
\usage{
\method{randomUniformForest}{formula}(formula, data = NULL, \dots)
\method{randomUniformForest}{default}(X, Y = NULL, xtest = NULL, ytest = NULL, 
	ntree = 100,
	mtry = ifelse(bagging, ncol(X), floor(4/3*ncol(X))),
	nodesize = 1,
	maxnodes = Inf,
	depth = Inf,
	depthcontrol = NULL,
	regression = ifelse(is.factor(Y), FALSE, TRUE),
	replace = ifelse(regression, FALSE, TRUE),
	OOB = TRUE,
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
	subsamplerate = ifelse(regression, 0.7, 1),
	importance = TRUE,
	bagging = FALSE,
	unsupervised = FALSE,
	unsupervisedMethod = c("uniform univariate sampling", 
	"uniform multivariate sampling", "with bootstrap"),
	classwt = NULL,
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
	rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L2", "L1"),
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	threads = "auto",
	parallelpackage = "doParallel",
	\dots)	
	\method{print}{randomUniformForest}(x, \dots)
	\method{summary}{randomUniformForest}(object, maxVar = 30, border = NA, \dots)
	\method{plot}{randomUniformForest}(x, threads = "auto", \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{maxVar}{
	    maximum number of variables to plot and print when summarizing  a randomUniformForest object.
	}
	\item{border}{
	    positive integer value or NA. Change color of the borders when plotting variable importance. By default, NA, which disables border.
	}
	\item{x, object}{
	    an object of class randomUniformForest.
	}
	\item{data}{
	    in case of formula, a data frame, or matrix, containing the variables (including response) and their values.
	}
	\item{X, formula}{
	    a data frame, or matrix, of predictors, or a formula describing the model to be fitted. Note that, it is strongly recommended to avoid formula when using options or with large samples.
	}
	\item{Y}{
	   a response vector. If it is a factor, classification is assumed, otherwise regression is computed.
	}
	\item{xtest}{
	     a data frame or matrix (like X) containing predictors for the test set.
	}
	\item{ytest}{
	     responses for the test set, if provided.
	}
	\item{ntree}{
	     number of trees to grow. Default value is 100. Do not set it too small.
	}
	\item{mtry}{
	     number of variables randomly sampled with replacement as candidates at each split. Default value is floor(4/3*ncol(X)) unless 'bagging' or 'randomfeature' options are specified. One can also set mtry = "random". For regression, increasing 'mtry' value usually leads to better accuracy. Note that mtry = 1 will lead to a \emph{purely uniformly random forest}. 'mtry' has also an effect when assessing variable importance.
	}
	\item{nodesize}{
	     minimal size of terminal nodes. Default value is 1 (for both classification and regression) and usually produce best results, as it reduces bias when trees are fully grown. Variance is increased, but that is exactly what Random Uniform Forests need. 
	}
	\item{maxnodes}{
	 	  maximal number of nodes for each tree. Default value is 'Inf', growing trees to maximum size. Random number of nodes is allowed, setting option to "random". 
	}
	\item{depth}{
	     depth of each tree. By default, Trees are fully grown. Maximum depth, for a balanced tree, is floor(log(n)/log(2)) where n = nrow(X). Stumps are not allowed, hence smallest depth is 3. Note that 'depth' has an effect when assessing variable importance. Enabling, in conjunction, 'depthcontrol' activates a deeper competition between nodes in order to reduce the loss of accuracy induced by the tree's depth.
	}
	\item{depthcontrol}{
	    an integer, beginning at 1. Let algorithm controls the growth of each tree, letting the optimization criterion depend on the number of nodes as the tree is growing. More precisely, the option activates an internal measure against which algorithm is competing. 'depthcontrol' usually works well for large samples and regression and is strongly recommended each time one wants to control both speed (lower values) and accuracy (higher values). 
	}
	\item{regression}{
	     only needed if either classification or regression has to be set explicitly. Otherwise, model checks if 'Y' is a factor (classification), or not (regression) before computing task. If Y is not a factor and one wants to do classification, must be set to FALSE. 
	}
	\item{replace}{
	     if TRUE, sampling of cases is done with replacement. By default, TRUE for classification, FALSE for regression.
	}
	\item{OOB}{
	     if replace is TRUE, then if OOB is TRUE, "Out-of-bag" evaluation is done, resulting in an estimate of generalization (and mean squared) error and bounds. OOB option has overhead on computing time, but it is one of the most useful option.
	}
	\item{BreimanBounds}{
	    if TRUE, computes all theoretical properties provided by Breiman (2001), since Random Uniform Forests inherit of Random Forests properties. For classification, it gives the two bounds of prediction error, average correlation between trees, strength and standard deviation of strength. In classsification, Breiman's bound should act as an upper bound; if not then overfitting is likely to happen. 
		A special case where the first Breiman's bound does not work is when classes are imbalanced. In such case, the second bound can be used since it overrides imbalanced classes and is the \emph{upper bound of prediction error}. It could be loose but, if there are enough trees and data, will be strongly reduced using options. For regression, model returns an estimate of the forest theoretical prediction error, its upper bound, mean prediction error of a tree, average correlation between trees residuals and expected squared bias. The estimate of theoretical prediction error of the forest (and the upper bound) is not a bound of test (mean squared) error. Note that for multi-class problems or large files, 'BreimanBounds' option requires a lot of computing time.
	}
	\item{subsamplerate}{
	     value is the rate of sub-sampling (Bulhmann et Yu, 2002) for training sample. By default, 0.7 for regression (1, e.g. no sub-sampling in classification). If 'replace' is TRUE, 'subsamplerate' can be set to values greater than 1. For regression, if only accuracy is needed, setting 'replace' to FALSE and 'subsamplerate' to 1, may improve results in some cases but OOB evaluation (and Breiman's bounds) will be lost.
	}
	\item{importance}{
	     should importance of predictors be assessed? By default, TRUE. 
	}
	\item{bagging}{
	     if TRUE, \emph{Bagging} (Breiman, 1996) of random uniform decision trees is done. Useful to compare "Bagging of random uniform decision trees" and usual "Bagging of trees". For regression, can give sometimes better results than sampling, with replacement, variables.
	}
	\item{unsupervised}{
	    unsupervised learning mode, following Breiman's ideas. Note that one has to call the second stage of unsupervised learning, see \code{\link{unsupervised.randomUniformForest}} to obtain a full object that will allow  clustering.
	}
	\item{unsupervisedMethod}{
	    method that has to be used to turn unsupervised problem into a supervised one. Note that 'unsupervisedMethod' uses either one (then bootstrap will not happen) or two arguments, the second one always being "with bootstrap" allowing, then, to use bootstrap. 
	}
	\item{classwt}{
	     for classification only. Prior of the classes. Need not add up to one. Useful for imbalanced classes. Note that if one wants to compute many forests and combine them, with 'classwt' enabled for only few of them, all others forests must have 'classwt' enabled, leading to classwt = rep(1, 'number of classes') for forests that one do not need to compute weights.
	}
	\item{oversampling}{
	     for classification, a scalar between 1 and -1 for over or under-sampling of minority or majority class, setting by the value of  'targetclass'. For example, if set to -0.3, and 'targetclass' set to 1, then first class (assumed to be the majority class) will be undersampled by a proportion of 30 percent. if set to 0.5, and 'targetclass' set to 2, then second class (assumed to be the minority class) will be  oversampled by a proportion of 50 percent. In all cases, size of original matrix will not be modified, hence 'oversampling' option is implicitly associated with undersampling. If one wants to force true oversampling, 'rebalanced sampling' is the alternative. Note that for both classification and regression, 'oversampling' is the value that can be used to perturb responses values (see 'outputperturbationsampling') and then must be set between 0.05 and 1. Let it to its default value if one wants to run 'outputperturbationsampling' automatically, especially if classification is done with 'oversampling' option.
	}
	\item{targetclass}{
	     for classification only. Which class (given by its subscript, e.g. 1 for first class) should be targeted by 'oversampling' or 'outputperturbationsampling' option ? 
	}
	\item{outputperturbationsampling}{
	     if TRUE, let model applies a random perturbation over responses vector. For classification, 'targetclass' must be set to the class  that will be perturbed. By default 5 percent of the values will be perturbed, but more is allowed (up to 100 percent) using 'oversampling' option. 'outputperturbationsampling' is aimed to reduce correlation between trees (residuals) and can be monitored using Breiman's bounds. It also works better with many trees and may help to protect from overfitting, since model does no longer use response values but a random gaussian variable with the same mean than response but not the same variance.
	}
	\item{rebalancedsampling}{
	     for classification only. Can be set to TRUE or to a vector containing the desired sample size for each class. If TRUE, model builds samples where all classes are equally distributed, leading to exactly balanced classes, by either oversampling or undersampling. If vector, size of train data will change, according to the sum of the values in vector.
	}
	\item{featureselectionrule}{
	    which optimization criterion should be chosen for growing trees ? By default, model uses "entropy" (in classification) to compute information gain function. If set to "random", model chooses randomly between Gini criterion and entropy for each node of each tree. For regression, sum of squared residuals ("L2") or sum of absolute residuals ("L1") is allowed, or "random".
	}
	\item{randomcombination}{
	     vector containing features index and, eventually, weight(s) for (random) combination of features. For example, if a combination of feature 1 and feature 2 is desired with a weight of 0.2 for the first, then randomcombination = c(1, 2, 0.2). If a combination of feature 3 and feature 4 is needed in the same time with a weight of 0.5, then randomcombination = c(1, 2, 3, 4, 0.2, 0.5). For full random combinations, one just put features in a vector, and then, in our example, randomcombination = c(1, 2, 3, 4). Useful sometimes to both reduce correlation between trees and increase their average individual strength (e.g. accuracy).
	}
	\item{randomfeature}{
	     if TRUE, a forest of totally randomized trees (e.g. purely random forest) will be grown. In this case, there is no optimization. Useful as a base result for forest of randomized trees.
	}
	\item{categoricalvariablesidx}{
	     which variables should be considered as categorical? By default, value is NULL, and then categorical variables remain in the same approach as continuous ones. If 'X' is a data frame, value can be set to "all", in which case model will automatically identify categorical variables and will use a different scheme for growing trees using these variables. If 'X' is a matrix, one has to set a vector containing subscripts of categorical variables. For example, if feature 1 and feature 2 are categorical, then categoricalvariablesidx = c(1,2). Note that using formula will automatically build dummies and throw them to the model. Current engine for categorical features in random uniform forests does not use dummies and is also designed to work with variable that take discrete values. 
	}
	\item{na.action}{
	    how to deal with NA data? By default, na.action = "fastImpute", using rough replacement with median or most frequent values. If speed is not required, na.action = "accurateImpute", can lead to better results using model to impute NA values. But, one can not not expect a fast imputation. na.action = "omit", just omit NA values and it's the only option when using formula. Note that many, accurate and fast, models (softImpute, missMDA, missForest, Amelia, ...) are available on CRAN to impute missing data. 
	}
	\item{logX}{
	     applies logarithm transformation to all predictors whose values are strictly positive, and ignores the others.
	}
	\item{classcutoff}{
	     for classification only. Change proportion of votes needed to get majority. First value of vector is the name of the class (between quotes) that has to be assessed. Second value is a kind of weight needed to get majority. For example, in a problem which classes "bad" and "good", and classcutoff = c("bad", 0.4), it means that class "bad" needs 'Cte/0.4' times more votes than class "good" to be the majority class when predicting a new observation, where Cte = 0.5 and all votes sum to the number of trees. 	
	}
	\item{threads}{
	     compute model in parallel for computers with many cores. Default value is 'auto', letting model run on all logical cores minus 1. User can set 'threads' to any value greater than 1. Note that, in Windows, logical cores consume same memory than physical ones, but will not speed up computation linearly with the number of logical cores. Note that it is strongly recommended to use only one thread for small samples.
	}
	\item{parallelpackage}{
	    which parallel back-end to use for computing parallel tasks ? By default and for ease of use, 'doParallel' is the package retained for now. Should not be modified. It has the great advantage to allow killing task, e.g. pushing the 'Stop' button without freezing R. Note that in this case and, at least, in Windows, one may need to manually terminate processes, using the 'Task manager' in order to avoid them occupying uselessly cores.  
	}
	\item{\dots}{
	  not currently used.
	}
}

\details{
 Random Uniform Forests are inspired by Bagging and Breiman's Random Forests (tm) but have many differences at theoretical and algorithmic levels. Random Uniform Forests build many randomized and unpruned trees and the four main differences with Random Forests are :\cr
 - sampling, with replacement, features\cr
 - sub-sampling data, in the case of regression,\cr
 - generating random cut-points according to the Uniform distribution, i.e. cut-points usually not belong to data but are virtual points drawn between the minimum and the maximum, or between two random points, of each candidate variable at each node, using the continuous Uniform distribution, since all points are (or will be always converted to) numeric values.\cr
 - optimization criterion. Maximizing Information gain is preferably used for classification. For regression, sum of squared (or absolute) residuals is computed for each candidate node (region) then, for each sampled feature, the metrics are summed for each pair of complementary nodes. The chosen pair is the one that reaches the minimum. More precisely, in regression only sums are involved and only in the candidate nodes (not in the current one).\cr
The enumeration above leads to a tree that is grown using global optimization, for the current partition, to select each node. Sampling features with replacement increases the competition between nodes, in order to limit variance, especially in the regression case where prediction error depends more on the model than in the classification case.

Others differences also appear at the node level. Like Random Forests, classification is done by majority vote, and regression by averaging trees outputs but : \cr
- trees are explicitly designed to have an average low bias, while trying to tame the happening of an increasing variance, and are thus optimized to reach a high level of randomness. The forest maintains the bias and reduces variance, since variance of the forest is approximatively (in regression) the product of average correlation between trees residuals and average variance of trees. This leads to the same scheme for prediction error. Note that the decrease in correlation can not be obtained in the same time than a decreasing variance. The main work is to decrease correlation faster than the growth of variance. Note also that low correlation is mandatory to reach convergence, especially in regression where average correlation tends to be high. Type vignette("randomUniformForestsOverview", package = "randomUniformForest") in the R prompt to get a summary of technical details.
 
 Others main features, thanks to Breiman's ideas, to the ensemble structure and to the Bayesian framework, are :\cr
 - some other paradigms of ensemble learning (like Bagging) using options,\cr
 - all Breiman's bounds,\cr
 - post-processing votes in order to lower MSE by lowering bias, see \code{\link{postProcessingVotes}},\cr
 - deep analysis of variable importance and selection, see \code{\link{importance.randomUniformForest}} and \code{\link{partialImportance}},\cr
 - partial dependencies, opening the way to extrapolation, see \code{\link{partialDependenceOverResponses}} and \code{\link{partialDependenceBetweenPredictors}},\cr
 - visualization tools and tables to help for a better interpretation,\cr
 - missing values imputation, see \code{\link{fillNA2.randomUniformForest}}, \cr
 - treatment of imbalanced classes,\cr
 - cost sensitive learning,\cr
 - quantile regression,\cr
 - prediction and confidence intervals, see \code{\link{bCI}}, \cr 
 - internal MapReduce paradigm for large datasets that can fit in memory, see \code{\link{rUniformForest.big}},\cr
 - incremental learning for large datasets that can not fit in memory, see \code{\link{rUniformForest.combine}}, \cr
 - distributed learning, allowing to run many different models on different data (sharing, at least, some features) on many computers and combine them in a single one, in different manners, for predictions. Note that one has to carefully manage the \emph{i.i.d.} assumption in order to see convergence happen.\cr
 - unsupervised learning and dimension reduction, see \code{\link{unsupervised.randomUniformForest}}.
  
 In particular, incremental learning is native, since the model uses random cut-points, and one can remove (see \code{\link{rm.trees}}), duplicate or add trees (but not modify) at each step of the incremental process.\cr
 Note that the model is not allowing results to be reproducible using the set.seed() function. One reason is that many (including essential) options run at the tree (or node) level in order to decrease correlation and use many random seeds internally. Since convergence is the primal property of Random Forests, for the same training sample, even if results will slightly vary, one has to consider the OOB estimate and Breiman's upper bound (in classification) as the main guarantees. They are effective only under the \emph{i.i.d.} assumption. If enough data is available, one can derive an OOB bound, giving conditions where test error would remain below OOB estimate and, consequently, below Breiman's bounds (see vignette).
 
 Note that speed is currently not at the state of art for small datasets, due to the majority R code and some constant overhead that seems coming from the parallelism. However some of the critical parts of the algorithm are written in C++, thanks to the Rcpp package. For large datasets the gap is greatly reduced, due to shortcuts added to the R code and increased randomness. That is the case when the dimension is getting high, or for regression. A great speed-up can also be achieved with the 'depth' option (for values close to 10), or the 'maxnode' one, at the cost of some loss in accuracy.
}

\value{
  An object of class randomUniformForest, which is a list with the following components:
  
  \item{forest }{list of tree objects, OOB objects (if OOB = TRUE), variable importance objects (if importance = TRUE). }
  \item{predictionObject }{if 'xtest' is not NULL, prediction objects.}
  \item{errorObject }{statistics about errors of model.}
  \item{forestParams }{almost all parameters of the model.}
  \item{classes }{original labels of response vector in cas of classification.}
  \item{logX }{TRUE, if logarithm transformation has been called.}
  \item{y }{training responses.}
  \item{variablesNames }{vector of variables names.}
  \item{call }{the original call to \code{randomUniformForest}.}
}

\references{
Biau, G., Devroye, L., Lugosi, G., 2008. \emph{Consistency of random forests and other averaging classifiers}. The Journal of Machine Learning Research 9, 2015-2033.

Breiman, L, 1996. \emph{Heuristics of instability and stabilization in model selection}. The Annals of Statistics 24, no. 6, 2350-2383.

Breiman, L., 1996. \emph{Bagging predictors}. Machine learning 24, 123-140.

Breiman, L., 2001. \emph{Random Forests}, Machine Learning 45(1), 5-32.

Breiman, L., Friedman, J.H., Olshen, R.A., Stone, C., 1984. \emph{Classification and Regression Trees}. New York: Chapman and Hall.

Ciss, S., 2014. PhD thesis: \emph{Forets uniformement aleatoires et detection des irregularites aux cotisations sociales}. Universite Paris Ouest Nanterre, France. In french.\cr
English title : \emph{Random Uniform Forests and irregularity detection in social security contributions}.\cr
Link : https://www.dropbox.com/s/q7hbgeafrdd8qtc/Saip_Ciss_These.pdf?dl=0  

Ho, T.K., 1998. \emph{The random subspace method for constructing decision forests}. IEEE Transactions on Pattern Analysis and Machine Intelligence 20, 832-844.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{predict.randomUniformForest}}, \code{\link{rUniformForest.big}}, \code{\link{rUniformForest.combine}}, \code{\link{rUniformForest.grow}}, \code{\link{importance.randomUniformForest}}, \code{\link{rm.trees}},
 \code{\link{roc.curve}}, \code{\link{fillNA2.randomUniformForest}}, \code{\link{getTree.randomUniformForest}}, \code{\link{unsupervised.randomUniformForest}} 
}

\examples{
# NOTE : use option 'threads = 1' (disabling parallel processing) to speed up computing 
# for small samples, since parallel processing is useful only for computationally intensive tasks

###### Part One : quick guide

## not run
#### Classification 
# data(iris)
# iris.ruf <- randomUniformForest(Species ~ ., data = iris, threads = 1)
# iris.ruf ## or print(iris.ruf)

## plot OOB error: 
# plot(iris.ruf, threads = 1)

## print and plot (global) variable importance and some statistics about trees:
# summary(iris.ruf)

#### Regression

## Note that when formula is used, missing values are automatically deleted and dummies
## are built for categorical features
# data(airquality)
# ozone.ruf <- randomUniformForest(Ozone ~ ., data = airquality, threads = 1)
# ozone.ruf

## plot OOB error: 
# plot(ozone.ruf, threads = 1)

## Bagging
# ozone.bagging.ruf <- randomUniformForest(Ozone ~ ., data = airquality,
# bagging = TRUE, threads = 1)

## Ensemble of totally randomized trees, e.g. purely random forest
# ozone.prf <- randomUniformForest(Ozone ~ ., data = airquality, randomfeature = TRUE, threads = 1)


#### Common case: use X, as a matrix or data frame, and Y, as a response vector, 

#### Classification : iris data, training and testing
data(iris)

## define random train and test sample. "Species" is the response vector
# iris.train_test <- init_values(iris[,-which(colnames(iris) == "Species")], iris$Species,
# sample.size = 1/2)

## iris train and test samples
# iris.train = iris.train_test$xtrain
# species.train = iris.train_test$ytrain
# iris.test = iris.train_test$xtest
# species.test = iris.train_test$ytest

## iris train and test modelling
# iris.train_and_test.ruf <- randomUniformForest(iris.train, species.train,
# xtest = iris.test, ytest = species.test, threads = 1)

# iris.train_and_test.ruf

## Balanced sampling : equal sample size for all labels
# iris.train_and_test.balancedsampling.ruf <- randomUniformForest(iris.train, species.train,
# xtest = iris.test, ytest = species.test, rebalancedsampling = TRUE, threads = 1)
								
###### Part Two : Summarized case studies (remove comments to run)

#### Classification : Wine Quality data
## http://archive.ics.uci.edu/ml/datasets/Wine+Quality
## We use red wine quality file : data have 1599 observations, 12 variables and 6 classes.
 
# data(wineQualityRed)
# wineQualityRed.data = wineQualityRed

## class and observations
# Y = wineQualityRed.data[,"quality"]
# X = wineQualityRed.data[, -which(colnames(wineQualityRed.data) == "quality")]

## First look : train model with default parameters (and retrieve estimates)
## call it standard model.
# wineQualityRed.std.ruf <- randomUniformForest(X, as.factor(Y), threads = 2)

## see OOB evaluation and parameters
# wineQualityRed.std.ruf 

## see statistics about the forest and global variable importance
# summary(wineQualityRed.std.ruf)

## But some labels do not have enough observations to assess variable importance
## merging class 3 and 4. Merging class 7 and 8 to get enough data.
# Y[Y == 3] = 4
# Y[Y == 8] = 7

## make Y as a factor, change names and get a summary
# Y = as.factor(Y)
# levels(Y) = c("3 or 4", "5", "6", "7 or 8")
# table(Y)

## learn a new model to get a better view on variable importance
## note : Y is now a factor, the model will catch the learning task as a classification
# wineQualityRed.new.ruf <- randomUniformForest(X, Y)
# wineQualityRed.new.ruf 

## global variable importance is more consistent
# summary(wineQualityRed.new.ruf)

## plot OOB error (needs some computing)
# plot(wineQualityRed.new.ruf, threads = 2)

## go deeper in assessing variable importance, using a high level of interaction
# importance.wineQualityRed <- importance(wineQualityRed.new.ruf, Xtest = X, maxInteractions = 6)
									
## visualize : global importance, importance based on interactions, 
## importance based on labels, partial dependencies for all influential variables 
## (loop over the prompt to get others partial dependencies)
## get more points, using option whichOrder = "all", default option.

# plot(importance.wineQualityRed, Xtest = X, whichOrder = "first")

## look at some specific labels from a (very) local view point
## what features for a very good wine (class 7 or 8) ?
# pImportance.wineQualityRed.class7or8 <- partialImportance(X, importance.wineQualityRed, 
# whichClass = "7 or 8", nLocalFeatures = 6)
											
## but how do they act ?
## get it feature after feature, recalling partial dependence
## and considering feature at the first order 
## assuming the feature is the most important, at least for the class one need to assess.

# pDependence.wineQualityRed.totalSulfurDioxide <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "total.sulfur.dioxide", 
# whichOrder = "first", outliersFilter = TRUE)
											
## see what happen then for "alcohol"
# pDependence.wineQualityRed.alcohol <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "alcohol", 
# whichOrder = "first", outliersFilter = TRUE)

#### Regression : Auto MPG 
## http://archive.ics.uci.edu/ml/datasets/Auto+MPG
## 398 observations, 8 variables, missing values
## Variable to predict : "mpg", miles per gallon 

## data(autoMPG)
# autoMPG.data = autoMPG

# Y = autoMPG.data[,"mpg"]
# X = autoMPG.data[,-which(colnames(autoMPG.data) == "mpg")]

## remove "car name" which is a variable with unique ID (car models)
# X = X[, -which(colnames(X) == "car name")]

## train the model and get OOB evaluation
# autoMPG.ruf <- randomUniformForest(X, Y)

## assess variable importance
# importance.autoMPG <- importance(autoMPG.ruf, Xtest = X)
# plot(importance.autoMPG, Xtest = X)

## what are the features that lead to a lower consumption (and high mpg)?
# pImportance.autoMPG.high <- partialImportance(X, importance.autoMPG, 
# threshold = mean(Y), thresholdDirection = "high", nLocalFeatures = 6)
											
## Look at "weight" and "acceleration" dependence
## note that option perspective = TRUE allows a 3D representation
# pDependence.autoMPG.weightAndAcceleration <- 
# partialDependenceBetweenPredictors(X, importance.autoMPG, c("weight", "acceleration"),
# whichOrder = "all", perspective = FALSE, outliersFilter =  TRUE)

##dtFW
}

\keyword{classification}
\keyword{regression}
\keyword{tree}
\keyword{ensemble}
\keyword{machine}
\keyword{incremental}
\keyword{learning}