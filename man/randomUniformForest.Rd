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
  Ensemble model for classification, regression and unsupervised learning, based on a forest of unpruned and randomized binary decision trees. Unlike \emph{Breiman's Random Forests}, each tree is grown by sampling, \emph{with replacement}, a set of variables before splitting each node. Each cut-point is generated randomly, according to the \emph{continuous Uniform distribution between two random points of each candidate variable or using its whole current support}. Optimal random node is, then, selected among many full random ones by maximizing Information Gain (classification) or minimizing a distance (regression), 'L2' (or 'L1'). Unlike \emph{Extremely Randomized Trees}, data are either \emph{bootstrapped or sub-sampled for each tree}. From the theoretical side, Random Uniform Forests are aimed to lower correlation between trees and to offer a deep analysis of variable importance. The unsupervised mode introduces clustering and dimension reduction, using a three-layer engine: dissimilarity matrix, Multidimensional Scaling (or Spectral decomposition) and k-means (or hierarchical clustering). From the practical side, Random Uniform Forests are designed to provide a complete analysis of (un)supervised problems and to allow native distributed and incremental learning.
}
\usage{
\method{randomUniformForest}{formula}(formula, data = NULL, subset = NULL, \dots)
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
	na.action = c("veryFastImpute", "fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	subset = NULL,
	usesubtrees = FALSE,
	threads = "auto",
	parallelpackage = "doParallel",
	\dots)	
	\method{print}{randomUniformForest}(x, biasCorrection = FALSE, \dots)
	\method{summary}{randomUniformForest}(object, maxVar = 30, border = NA, \dots)
	\method{plot}{randomUniformForest}(x, threads = "auto", \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{maxVar}{
	    maximum number of variables to plot and print when summarizing a randomUniformForest object.
	}
	\item{border}{
	    positive integer value or NA. Change colour of the borders when plotting variable importance. By default, NA, which disables border.
	}
	\item{biasCorrection}{
		if TRUE, bias correction that will be applied to model errors.
	}
	\item{x, object}{
	    an object of class randomUniformForest.
	}
	\item{data}{
	    in case of formula, a data frame, or matrix, containing the variables (including response) and their values.
	}
	\item{subset}{	
		an index vector indicating which rows should be used.
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
	     number of variables randomly sampled with replacement as candidates at each split. Default value is floor(4/3*ncol(X)) unless 'bagging' or 'randomfeature' options are specified. For regression, increasing 'mtry' value usually leads to better accuracy. Note that mtry = 1 will lead to a \emph{purely uniformly random forest}. 'mtry' has also an effect when assessing variable importance. Random 'mtry' is allowed setting mtry = "random". 
	}
	\item{nodesize}{
	     minimal size of terminal nodes. Default value is 1 (for both classification and regression) and usually produce best results, as it reduces bias when trees are fully grown. Variance is increased, but that is exactly what Random Uniform Forests need. Random 'nodesize' is allowed, setting option to "random". For each tree 'nodesize' will take a random value between 1 and 50, using increments of 5.
	}
	\item{maxnodes}{
	 	  maximal number of nodes for each tree. Default value is 'Inf', growing trees to maximum size. Random number of nodes is allowed, setting option to "random". 
	}
	\item{depth}{
	    depth of each tree. By default, Trees are fully grown. Maximum depth, for a balanced tree, is log(n)/log(2).
		In regression it will usually be the case. Stumps are not allowed, hence smallest depth is 3. Note that 'depth' has an effect when assessing variable importance. Enabling, in conjunction, 'depthcontrol' activates a deeper competition between nodes in order to reduce the loss of accuracy induced by the tree's depth. Random depth is allowed, setting option to "random". 
	}
	\item{depthcontrol}{
	    an integer, beginning at 1. Let algorithm controls the growth of each tree, letting the optimization criterion depend on the number of nodes as the tree is growing. More precisely, the option activates an internal measure against which algorithm is competing. 'depthcontrol' usually works well in conjunction with 'depth' option for large samples and regression, and is strongly recommended each time one wants to control both speed (lower values) and accuracy (higher values). 
	}
	\item{regression}{
	     only needed if either classification or regression has to be set explicitly. Otherwise, model checks if 'Y' is a factor (classification), or not (regression) before computing task. If 'Y' is not a factor and one wants to do classification, must be set to FALSE. 
	}
	\item{replace}{
	     if TRUE, sampling of cases is done with replacement. By default, TRUE for classification, FALSE for regression.
	}
	\item{OOB}{
	     if replace is TRUE, then if OOB is TRUE, "Out-of-bag" evaluation is done, resulting in an estimate of generalization (and mean squared) error and bounds. OOB option has an overhead on computing time especially for a large number of trees and large datasets, but it is one of the most useful option.
	}
	\item{BreimanBounds}{
	    if TRUE, computes all theoretical properties provided by Breiman (2001), since Random Uniform Forests inherit of Random Forests properties. For classification, it gives the two bounds of prediction error, average correlation between trees, strength and standard deviation of strength. For regression, model returns an estimate of the forest theoretical prediction error, its upper bound, mean prediction error of a tree, average correlation between trees residuals and expected squared bias. Note that for multi-class problems or large files, 'BreimanBounds' option requires a lot of computing time. For more details see 'Notes' below.
	}
	\item{subsamplerate}{
	     value is the rate of sub-sampling (Bulhmann et Yu, 2002) for training sample. By default, 0.7 for regression (1, e.g. no sub-sampling, in classification). If 'replace' is TRUE, 'subsamplerate' can be set to values greater than 1. For regression, if only accuracy is needed, setting 'replace' to FALSE and 'subsamplerate' to 1, may improve results in some cases but OOB evaluation (and Breiman's bounds) will be lost.
	}
	\item{importance}{
	     should importance of predictors be assessed? By default, TRUE. Note that it is strongly recommended to set it to FALSE for large datasets (then enable it for a subset of the dataset) since it is a intensive task.
	}
	\item{bagging}{
	     if TRUE, \emph{Bagging} (Breiman, 1996) of random uniform decision trees is done. Useful to compare "Bagging of random uniform decision trees" and usual "Bagging of trees". For regression, can give sometimes better results than sampling, with replacement, variables at each node.
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
	     for classification, a scalar between 1 and -1 for over or undersampling of minority or majority class, given by the value of 'targetclass'. For example, if set to 'oversampling = -0.3', and 'targetclass = 1', then the first class (assumed to be the majority class) will be undersampled by a proportion of 30 percent. if set to 0.5, and 'targetclass' set to 2, then second class (assumed to be the minority class) will be oversampled by a proportion of 50 percent. In all cases, size of original matrix will not be modified, hence 'oversampling' option is implicitly associated with undersampling. If one wants to force true oversampling or undersampling, 'rebalanced sampling' is the alternative. Note that for both classification and regression, 'oversampling' is also the value that can be used to perturb response values (see 'outputperturbationsampling'). Hence, in classification one should avoid to use 'outputperturbationsampling' and 'oversampling' together.
	}
	\item{targetclass}{
	     for classification only. Which class (given by its subscript, e.g. 1 for first class) should be targeted by 'oversampling' or 'outputperturbationsampling' option ? 
	}
	\item{outputperturbationsampling}{
	     if TRUE, let model applies a random perturbation over responses vector. For classification, 'targetclass' must be set to the class  that will be perturbed. By default 5 percent of the values will be perturbed, but more is allowed (more than 100 percent, with bootstrap) using the 'oversampling' option. 'outputperturbationsampling' is aimed to reduce correlation between trees (residuals) and can be monitored using Breiman's bounds. It also works better with many trees and may help to protect from overfitting, since model does no longer use training labels (or values) but, in regression, a random Gaussian variable with similar (but possibly different) mean than response but different variance. Note that each tree generates its own Gaussian variable whose parameters vary.
	}
	\item{rebalancedsampling}{
	     for classification only. Can be set to TRUE or to a vector containing the desired sample size for each class. If TRUE, model builds samples where all classes are equally distributed, leading to exactly balanced classes, by either oversampling or undersampling. If vector, size of train data will change, according to the sum of the values in the vector. If the number of minority class cases asked in the vector is greater than the original one, sampling with replacement is done for these cases.
	}
	\item{featureselectionrule}{
	    which optimization criterion should be chosen for growing trees ? By default, model uses "entropy" (in classification) to compute the Information Gain function. If set to "random", model chooses randomly between Gini criterion and entropy for each node of each tree. For regression, sum of squared residuals ("L2") or sum of absolute residuals ("L1") is allowed, or "random".
	}
	\item{randomcombination}{
	     vector containing features index and, eventually, weight(s) for (random) combination of features. For example, if a combination of feature 1 and feature 2 is desired with a weight of 0.2 for the first, then randomcombination = c(1, 2, 0.2). If a combination of feature 3 and feature 4 is needed in the same time with a weight of 0.5, then randomcombination = c(1, 2, 3, 4, 0.2, 0.5). For full random combinations, one just put features in a vector, and then, in our example, randomcombination = c(1, 2, 3, 4). Useful sometimes to both reduce correlation between trees and increase their average individual strength.
	}
	\item{randomfeature}{
	     if TRUE, a forest of totally randomized trees (e.g. purely random forest) will be grown. In this case, there is no optimization. Useful as a base result for forest of randomized trees.
	}
	\item{categoricalvariablesidx}{
	     which variables should be considered as categorical? By default, value is NULL, and then categorical variables remain in the same approach as continuous ones. If 'X' is a data frame, option can be set to "all", in which case the model will automatically identify categorical variables and will use a different scheme for growing trees using these variables (see the Details section). If 'X' is a matrix, one has to set a vector containing subscripts of categorical variables. For example, if feature 1 and feature 2 are categorical, then categoricalvariablesidx = c(1,2). Note that using formula will automatically build dummies and throw them to the model. Current engine for categorical features in random uniform forests does not use dummies and is also designed to work with variable that take discrete values. Note that assessing categorical variable increases the computation time.
	}
	\item{na.action}{
	    how to deal with NA data? By default, na.action = "veryFastImpute", using rough replacement with median, mean or most frequent values. If speed is not required, na.action = "accurateImpute", can lead to better results using model to impute NA values. na.action = "omit", just omit NA values and it's the only option when using formula. Note that many, accurate and fast, models (softImpute, missMDA, missForest, Amelia, ...) are available on CRAN to impute missing data. Note that "accurateImpute" calls the \code{\link{fillNA2.randomUniformForest}} function, hence it is recommended to impute data outside of the model, especially  for large datasets.
	}
	\item{logX}{
	     applies logarithm transformation to all predictors whose values are strictly positive, and ignores the others.
	}
	\item{classcutoff}{
	    for classification only. Change proportion of votes needed to get majority. First value of vector is the name of the class (between quotes) that has to be assessed. Second value is a kind of weight needed to get majority. For example, in a problem with classes "bad" and "good", and 'classcutoff = c("bad", 0.4)', it means that class "bad" needs 'Cte/0.4' times more votes than class "good" to be the majority class when predicting the label of a new observation, where Cte = 0.5 and sum of all votes equals to the number of trees. Note that the option can be monitored with the OOB evaluation.	
	}
	\item{usesubtrees}{
		allows to grow trees with an arbitrary depth (default one being half of the maximum depth), that are, in a second step, updated by extending all non pure nodes. The option is mainly designed to reduce memory footprint at the cost of an increasing computing time, for, at least, the same expected accuracy than in the standard case. Note that it seems also more robust when training real life and large datasets, especially when the average depth of trees is far from the theoretical one. 
	}
	\item{threads}{
	     compute model in parallel for computers with many cores. Default value is 'auto', letting model run on all logical cores minus 1. User can set 'threads' to any value greater than 1. Note that, in Windows, logical cores consume same memory than physical ones, but will not speed up computation linearly with the number of logical cores. Note also that it is strongly recommended to use only one thread for small samples.
	}
	\item{parallelpackage}{
	    which parallel back-end to use for computing parallel tasks ? By default and for ease of use, 'doParallel' is the package retained for now. Should not be modified. It has the great advantage to allow killing task, e.g. pushing the 'Stop' button without freezing R. Note that in this case and, at least, in Windows, one may need to manually terminate processes, using the 'Task manager' in order to avoid them occupying uselessly cores.  
	}
	\item{\dots}{
	  not currently used.
	}
}

\details{
Random Uniform Forests have been designed to provide a complete analysis of supervised and unsupervised problems.
From pre-processing to predictions, interpretation and summary of results. Hence, many functions and options have been built around the algorithm in order to let the user simply put the data as they come (matrix or data frame), getting basis objects that come with usual (and statistical) measures of assessment for all critical parts. 

Random Uniform Forests are inspired by \emph{Bagging} and Breiman's \emph{Random Forests} (tm) but have many differences at theoretical and algorithmic levels. They build many randomized and unpruned binary decision trees and the four main differences with Random Forests are:\cr
 - sampling \emph{with replacement} a set of features to generate each candidate node,\cr
 - subsampling data, in the case of regression,\cr
 - generating \emph{random cut-points according to the Uniform distribution}, i.e. cut-points usually not belong to data but are virtual points drawn between the minimum and the maximum, or between two random points, of each candidate variable at each node, using the continuous Uniform distribution, since all points are (or will always be converted to) numeric values.\cr
 - optimization criterion. Maximizing the Information Gain is preferably used for classification. For regression, sum of squared (or absolute) residuals is computed for each candidate node (region), using for each a random sampled pair of feature and cut-point. Then the metrics are summed for each pair of complementary nodes. The chosen pair is the one that reaches the minimum over all current candidates pairs, and define the optimal (and random) child node. More precisely, \emph{in regression only sums are involved and only in the candidate nodes (not in the current one)}. Note that it could (while not clear) also be the case for Random Forests.
 
The enumeration above leads to a large and deep tree that is grown using global optimization, for the current partition, to select each node. Sampling features, with replacement, increases the competition between nodes, in order to limit variance, especially in the regression case where prediction error depends more on the model than in the classification case.

Others differences also appear at the node level. Like Random Forests, classification is done by majority vote, and regression by averaging trees outputs but: \cr
- trees can be updated with streaming data (currently disabled for further tests),\cr
- forests with different parameters and data can be combined to form one forest,\cr
- trees are explicitly designed to have an average low bias, while trying to tame the happening of an increasing variance, and are thus optimized to reach a high level of randomness. The forest maintains the bias and reduces variance, since variance of the forest is approximatively (in regression) the product of average correlation between trees and average variance of trees. This leads to a similar scheme for the prediction error bounded by 'average correlation between trees residuals' x 'average variance of trees residuals'. Note that the trend in decreasing correlation can not be obtained in the same time than a decreasing variance. The main work is to \emph{decrease correlation faster than the growth of variance}. Low correlation is mandatory to reach convergence and prevent overfitting, especially in regression where average correlation tends to be high (Ciss, 2015a, 2015c). One may experiment it since the model produces, by default, all Breiman's bounds and their details.
 
 Others main features, thanks to Breiman's ideas, to the ensemble structure and to the Bayesian framework, are:\cr
 - some other paradigms of ensemble learning (like Bagging of random uniform decision trees or ensemble of totally randomized trees) using options,\cr
 - functions to manipulate and plot trees, see \code{\link{getTree.randomUniformForest}} and friends\cr
 - all Breiman's bounds,\cr
 - (internal) pre-processing in order to handle (almost) any matrix or data frame,\cr
 - post-processing votes in order to lower MSE by reducing bias, see \code{\link{postProcessingVotes}},\cr
 - changing majority vote, using options ('classcutoff')\cr
 - \emph{output perturbation sampling}, lowering more the correlation, replacing completely (for regression and for each tree) the training vector of responses by an independent random Gaussian one with similar mean but different variance,\cr 
 - deep analysis of variable importance and selection, see \code{\link{importance.randomUniformForest}} and \code{\link{partialImportance}},\cr
 - partial dependencies, opening the way to extrapolation, see \code{\link{partialDependenceOverResponses}} and \code{\link{partialDependenceBetweenPredictors}},\cr
 - visualization tools and tables for almost any essential function, \cr
 - generic function to assess results, see \code{\link{model.stats}},\cr
 - generic cross-validation function, see \code{\link{generic.cv}}, \cr
 - missing values imputation, see \code{\link{rufImpute}},\cr
 - many methods for imbalanced classes, using options ('oversampling', 'rebalancedsampling', 'classwt', 'classcutoff', 'usesubtrees') or \code{\link{reSMOTE}}\cr
 - cost-sensitive learning, using options 'classwt' (which is dual) and friends\cr
 - native handling of categorical variables using a randomization mechanism at the node level. More precisely, the algorithm selects for each candidate node, randomly and before the splitting process, two values. The first one keeps its position along the variable while the second replaces, temporarily, all others values of the variable. This leads to a binary variable that can be treated like a numerical one. After the splitting, the variable recovers its original values. Since cut-points are almost virtual and random (a cut-point is not a point of the training sample), one has just to take care that the random splitting would not weaker the variable.\cr
 - quantile regression, see \code{\link{predict.randomUniformForest}}\cr
 - (new methods for) prediction and confidence intervals, see \code{\link{bCI}}, \cr
 - unsupervised learning, see \code{\link{unsupervised.randomUniformForest}} and friends:\cr
 - dimension reduction, using MDS or Spectral decomposition, see \code{\link{unsupervised.randomUniformForest}},\cr
 - dynamic clustering, allowing to split/merge/modify on-the-fly cluster(s), see \code{\link{modifyClusters}} and friends,\cr
 - visualization, allowing to display, manipulate and assess clusters in an easy way, see \code{\link{unsupervised.randomUniformForest}} plot and print methods,\cr
 - variable importance for clusters, inherited from the supervised case,\cr
 - cluster analysis, closing unsupervised learning in a compact and granular view, see \code{\link{clusterAnalysis}},\cr
 - native parallelism, thanks to the \emph{parallel}, \emph{doParallel} and \emph{foreach} packages,\cr
 - internal MapReduce paradigm for large datasets that can fit in memory, see \code{\link{rUniformForest.big}},\cr
 - incremental learning for large datasets that can not fit in memory, see \code{\link{rUniformForest.combine}}, \cr
 - distributed learning, allowing to run many different models on different data (sharing, at least, some features) on many computers and combine them in a single one, in different manners, for predictions. See \code{\link{rUniformForest.combine}} examples. Note that one has to carefully manage the \emph{i.i.d.} assumption in order to see convergence happen.
   
 In particular, incremental learning is native, since the model uses random cut-points, and one can remove, duplicate, add or modify/update trees at each step of the incremental process, see \code{\link{rm.trees}}.\cr
 The model is not allowing results to be exactly reproducible using the set.seed() function. One reason is that many (including essential) options run at the tree (or node) level in order to decrease correlation and use many random seeds internally. Since \emph{convergence is the primal property of Random Forests}, for the same enough large training sample, even if results will slightly vary, one has to consider the OOB estimate and Breiman's upper bound (in classification) as the main guarantees. They are effective only under the \emph{i.i.d.} assumption. If enough data is available, one can derive OOB bounds (Ciss, 2015c), leading the test error to be bounded by the OOB error, the latter
 itself bounded by the Breiman's bounds.
 
 Note that speed is currently not at the state-of-the-art for small datasets, due to the majority of R code and some constant overhead that seems coming from the parallelism. However some of the critical parts of the algorithm are written in C++, thanks to the Rcpp package. For large datasets the gap is greatly reduced, thanks to shortcuts added to the R code and increased randomness. That is the case when the dimension is getting high, or for regression. A great speed-up can also be achieved with 'depth', 'maxnodes', 'subsamplerate', 'mtry' (large values or with mtry = 1), 'randomfeature' (in combination with mtry), 'rebalancedsampling' (in classification) options or by combining (see \code{\link{rUniformForest.combine}}) many forests built upon chunks of data. All these tools will usually have a cost, loss in accuracy, depending on the dataset and the task. 
}

\value{
  An object of class randomUniformForest, which is a list with the following components:
  
  \item{forest }{list of tree objects, OOB objects (if OOB = TRUE), variable importance objects (if importance = TRUE). }
  \item{predictionObject }{if 'xtest' is not NULL, prediction objects.}
  \item{errorObject }{statistics about errors of the model.}
  \item{forestParams }{almost all parameters of the model.}
  \item{classes }{original labels of response vector in case of classification.}
  \item{logX }{TRUE, if logarithm transformation has been called.}
  \item{y }{training responses.}
  \item{variablesNames}{vector of variables names.}
  \item{call }{the original call to \code{\link{randomUniformForest}}.}
}

\note{
- About Breiman's bounds: In classification, Breiman's bound should act as an upper bound; if not then overfitting is likely to happen.  A special case where the first Breiman's bound does not work is when classes are imbalanced. In such case, the second bound can be used since it overrides imbalanced classes and is the \emph{upper bound of prediction error}. It could be loose but, if there are enough trees and data, will be strongly reduced using options. The estimate of theoretical prediction error of the forest (and the upper bound) is not a bound of test (mean squared) error.
}

\references{
 Amit, Y., Geman, D., 1997. Shape Quantization and Recognition with Randomized Trees. \emph{Neural Computation} 9, 1545-1588.

Biau, G., Devroye, L., Lugosi, G., 2008. Consistency of random forests and other averaging classifiers. \emph{The Journal of Machine Learning Research} 9, 2015-2033.

Breiman, L, 1996. Heuristics of instability and stabilization in model selection. \emph{The Annals of Statistics} 24, no. 6, 2350-2383.

Breiman, L., 1996. Bagging predictors. \emph{Machine learning} 24, 123-140.

Breiman, L., 2001. Random Forests, \emph{Machine Learning} 45(1), 5-32.

Breiman, L., 2001. Statistical Modeling: The Two Cultures (with comments and a rejoinder by the author). \emph{Statistical Science} 16, no. 3, 199-231.

Breiman, L., Friedman, J.H., Olshen, R.A., Stone, C., 1984. \emph{Classification and Regression Trees}. New York: Chapman and Hall.

Ciss, S., 2014. PhD thesis: \emph{Forets uniformement aleatoires et detection des irregularites aux cotisations sociales}. Universite Paris Ouest Nanterre, France. In french.\cr
English title : \emph{Random Uniform Forests and irregularity detection in social security contributions}.\cr
Link : https://www.dropbox.com/s/q7hbgeafrdd8qtc/Saip_Ciss_These.pdf?dl=0

Ciss, S., 2015a. \emph{Random Uniform Forests}. Preprint. hal-01104340.
 
Ciss, S., 2015b. \emph{Variable Importance in Random Uniform Forests}. Preprint. hal-01104751. 

Ciss, S., 2015c. \emph{Generalization Error and Out-of-bag Bounds in Random (Uniform) Forests}. Preprint. hal-01110524.

Dietterich, T.G., 2000. \emph{Ensemble Methods in Machine Learning}, in : Multiple Classifier Systems, Lecture Notes in Computer Science. Springer Berlin Heidelberg, pp. 1-15.

Efron, B., 1979. Bootstrap Methods: Another Look at the Jackknife. \emph{The Annals of Statistics} 7, 1-26.   

Ho, T.K., 1998. The random subspace method for constructing decision forests.  \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence} 20, 832-844.

Vapnik, V.N., 1995. \emph{The nature of statistical learning theory}. Springer-Verlag New York, Inc., New York, NY, USA.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{predict.randomUniformForest}}, \code{\link{rUniformForest.big}}, \code{\link{rUniformForest.combine}}, \code{\link{rUniformForest.grow}}, \code{\link{importance.randomUniformForest}}, \code{\link{rm.trees}},
 \code{\link{roc.curve}}, \code{\link{rufImpute}}, \code{\link{getTree.randomUniformForest}}, \code{\link{unsupervised.randomUniformForest}}, \code{\link{clusteringObservations}}  
}

\examples{
## not run
## NOTE : use option 'threads = 1' (disabling parallel processing) to speed up computing 
## for small samples, since parallel processing is useful only for computationally 
## intensive tasks

###### PART ONE : QUICK GUIDE

#### Classification 

# data(iris)
# iris.ruf <- randomUniformForest(Species ~ ., data = iris, threads = 1)

## Regular companions (from 1 to 18):
## 1 -  model, parameters, statistics:
# iris.ruf ## or print(iris.ruf)

## 2 - OOB error: 
# plot(iris.ruf, threads = 1)

## 3 - (global) variable importance, some statistics about trees:
# summary(iris.ruf)

#### Regression

## NOTE: when formula is used, missing values are automatically deleted and dummies
## are built for categorical features

# data(airquality)
# ozone.ruf <- randomUniformForest(Ozone ~ ., data = airquality, threads = 1)
# ozone.ruf

## plot OOB error: 
# plot(ozone.ruf, threads = 1)

## 4 - Alternative modelling:
## 4.1 bagging like:
# ozone.bagging <- randomUniformForest(Ozone ~ ., data = airquality,
# bagging = TRUE, threads = 1)

## 4.2 Ensemble of totally randomized trees, e.g. purely random forest:
# ozone.prf <- randomUniformForest(Ozone ~ ., data = airquality, 
# randomfeature = TRUE, threads = 1)

## 4.3 Extremely randomized trees like:
# ozone.ETlike <- randomUniformForest(Ozone ~ ., data = airquality, 
# subsamplerate = 1, replace = FALSE, bagging = TRUE, mtry = floor((ncol(airquality)-1)/3),
# threads = 1)

## Common case: use X, as a matrix or data frame, and Y, as a response vector, 

#### Classification : iris data, training and testing

# data(iris)

## define random training and test sample :
## "Species" is the response vector

# set.seed(2015)
# iris.train_test <- init_values(iris[,-which(colnames(iris) == "Species")], iris$Species,
# sample.size = 1/2)

## training and test samples:
# iris.train = iris.train_test$xtrain
# species.train = iris.train_test$ytrain
# iris.test = iris.train_test$xtest
# species.test = iris.train_test$ytest

## 5 - training and test (or validation) modelling:
# iris.train_test.ruf <- randomUniformForest(iris.train, species.train,
# xtest = iris.test, ytest = species.test, threads = 1)

## 6 - all-in-one results:
# iris.train_test.ruf

## 7 - Alternative modelling: imbalanced classes
## balanced sampling (for example): equal sample size for all labels

# iris.train_test.balancedsampling.ruf <- randomUniformForest(iris.train, species.train,
# xtest = iris.test, ytest = species.test, rebalancedsampling = TRUE, threads = 1)
								
###### PART TWO : SUMMARIZED CASE STUDY

#### Classification : Wine Quality data
## http://archive.ics.uci.edu/ml/datasets/Wine+Quality
## We use 'red wine quality' file : data have 1599 observations, 12 variables and 6 classes.
 
# data(wineQualityRed)
# wineQualityRed.data = wineQualityRed

## class and observations
# Y = wineQualityRed.data[, "quality"]
# X = wineQualityRed.data[, -which(colnames(wineQualityRed.data) == "quality")]

## First look : train model with default parameters (and retrieve estimates)
# wineQualityRed.std.ruf <- randomUniformForest(X, as.factor(Y))
# wineQualityRed.std.ruf 

## (global) Variable Importance:
# summary(wineQualityRed.std.ruf)

## But some labels do not have enough observations to assess variable importance
## merging class 3 and 4. Merging class 7 and 8 to get enough labels.
# Y[Y == 3] = 4
# Y[Y == 8] = 7

## make Y as a factor, change names and get a summary
# Y = as.factor(Y)
# levels(Y) = c("3 or 4", "5", "6", "7 or 8")
# table(Y)

## learn a new model to get a better view on variable importance
## NOTE: Y is now a factor, the model will catch the learning task as a classification one
# wineQualityRed.new.ruf <- randomUniformForest(X, Y)
# wineQualityRed.new.ruf 

## global variable importance is more consistent
# summary(wineQualityRed.new.ruf)

## plot OOB error (needs some computing)
# plot(wineQualityRed.new.ruf)

## 8 - alternative Modelling: use subtrees (small trees, extended then reassembled)
## may change something, depending on data
# wineQualityRed.new.ruf <- randomUniformForest(X, Y, usesubtrees = TRUE)

## 9 - deep variable importance:
## 9.1 - interactions are granular: use more for consistency, or less to see primary information
## 9.2 - a table is printed with details
# importance.wineQualityRed <- importance(wineQualityRed.new.ruf, Xtest = X, maxInteractions = 6)
									
## 10 - visualization: 
## 10.1 - global importance, interactions, importance based on interactions, 
## importance based on labels, partial dependencies for all influential variables 
## (loop over the prompt to get others partial dependencies)
## 10.2 - get more points, using option whichOrder = "all", default option.

# plot(importance.wineQualityRed, Xtest = X, whichOrder = "first")

## 11 - Cluster analysis: (if quick answers are needed)
## Note: called 'cluster' since it was first designed for clustering
## 11.1 - choose the granularity : components, maximum features, (as) categorical ones
## 11.2 - get a compact view
## 11.3 - see how importance is explaining the data
# analysis.wineQualityRed = clusterAnalysis(importance.wineQualityRed, X, components = 3, 
# maxFeatures = 3, clusteredObject = wineQualityRed.new.ruf, categorical = NULL, OOB = TRUE)

## 11.4 - interpretation: 
## Numerical features average: a good wine has much less volatile acidity,
## much more citric acid, ... than a wine of low quality.
## Most influential features:  while volatile.acidity seems to be important,...
## (Component frequencies:) ..., all variables must be taken into account, since information
## provided by the most important ones does not enough cover the whole available information.

## 11.5 - Complementarity:
## One does not forget to look plot of importance function. clusterAnalysis( )
## is a summarized view of the former and should not have contradictory terms
## but, eventually, complementary ones.

## 12 - Partial importance: (local) variable importance per class
## which features for a very good wine (class 7 or 8) ?
## Note: in classification, partial importance is almost the same than "variable importance over
## labels", being more local but they have different interpretations. The former is exclusive.
# pImportance.wineQualityRed.class7or8 <- partialImportance(X, importance.wineQualityRed, 
# whichClass = "7 or 8", nLocalFeatures = 6)
											
## 13 - Partial dependencies: how response relies to each variable or a pair of ones?
## 13.1 - admit options.
## get it feature after feature, recalling partial dependence and considering feature 
## at the first order assuming the feature is the most important, 
## at least for the class one need to assess.

# pDependence.wineQualityRed.totalSulfurDioxide <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "total.sulfur.dioxide", 
# whichOrder = "first", outliersFilter = TRUE)

## 13.2 - Look for the second order (assuming the feature is the second most important, 
## at least for the class one need to assess).
# pDependence.wineQualityRed.totalSulfurDioxide <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "total.sulfur.dioxide", 
# whichOrder = "second", outliersFilter = TRUE)

## 13.3 - Look at all orders: no assumptions, simply look the average effect 
# pDependence.wineQualityRed.totalSulfurDioxide <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "total.sulfur.dioxide", 
# whichOrder = "all", outliersFilter = TRUE)
											
## see what happens then for "alcohol" (more points using option 'whichOrder = "all"')
# pDependence.wineQualityRed.alcohol <- partialDependenceOverResponses(X, 
# importance.wineQualityRed, whichFeature = "alcohol",  
# whichOrder = "first", outliersFilter = TRUE)

## 13.4 - Translate interactions into dependence : pair of features
## is interaction leading to the same class (underlying structure)?
## is dependence linear ? 
## for which values of the pair is the dependence the most effective ?
 
# pDependence.wineQualityRed.sulfatesAndVolatileAcidity <- partialDependenceBetweenPredictors(X, 
# importance.wineQualityRed, c("sulphates", "volatile.acidity"), 
# whichOrder = "all", outliersFilter = TRUE)

#### Regression : Auto MPG 
## http://archive.ics.uci.edu/ml/datasets/Auto+MPG
## 398 observations, 8 variables,
## Variable to predict : "mpg", miles per gallon 

# data(autoMPG)
# autoMPG.data = autoMPG

# Y = autoMPG.data[,"mpg"]
# X = autoMPG.data[,-which(colnames(autoMPG.data) == "mpg")]

## remove "car name" which is a variable with unique ID (car models)
# X = X[, -which(colnames(X) == "car name")]

## train the default model and get OOB evaluation
# autoMPG.ruf <- randomUniformForest(X, Y)

## assess variable importance (ask more points with 'maxInteractions' option)
## NOTE: importance strongly depends on 'ntree' and 'mtry' parameters
# importance.autoMPG <- importance(autoMPG.ruf, Xtest = X)

## 14 - Dependence on most important predictors: marginal distribution of the response
## over each variable
# plot(importance.autoMPG, Xtest = X)

## 15 - Extrapolation:
## recalling partial dependencies and getting points
## NOTE : points are the result of the forest classifier and not of the training responses
# pDependence.autoMPG.weight <- partialDependenceOverResponses(X, importance.autoMPG,
# whichFeature = "weight", whichOrder = "all", outliersFilter = TRUE)

## 16 - Visualization again: view as discrete values
# visualize again 'model year' as a discrete variable and not as a continuous one 
# pDependence.autoMPG.modelyear <- partialDependenceOverResponses(X, importance.autoMPG,
# whichFeature = "model year", whichOrder = "all", maxClasses = 30)

## 16 - Partial importance for regression: see important variables only for a part 
## of response values
## what are the features that lead to a lower consumption (and high mpg)?
# pImportance.autoMPG.high <- partialImportance(X, importance.autoMPG, 
# threshold = mean(Y), thresholdDirection = "high", nLocalFeatures = 6)
											
## 17 - Partial dependencies between covariates: 
## look at "weight" and "acceleration" dependence 
# pDependence.autoMPG.weightAndAcceleration <- 
# partialDependenceBetweenPredictors(X, importance.autoMPG, c("weight", "acceleration"),
# whichOrder = "all", perspective = FALSE, outliersFilter = TRUE)

## 18 - More visualization: 3D (looking to the prompt to start animation)
## Note: requires some computation
# pDependence.autoMPG.weightAndAcceleration <- 
# partialDependenceBetweenPredictors(X, importance.autoMPG, c("weight", "acceleration"),
# whichOrder = "all", perspective = TRUE, outliersFilter = FALSE)

##dtFW
}

\concept{machine learning}
\concept{supervised learning}
\concept{classification}
\concept{regression}
\concept{decision tree}
\concept{ensemble methods}