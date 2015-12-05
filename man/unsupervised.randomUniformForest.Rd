\name{unsupervised.randomUniformForest}
\alias{unsupervised.randomUniformForest}
\alias{unsupervised}
\alias{plot.unsupervised}
\alias{print.unsupervised}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Unsupervised Learning with Random Uniform Forests
}
\description{
   Unsupervised mode of Random Uniform Forests, allowing clustering, dimension reduction, visualization and variable importance, using a three-layer engine: dissimilarity matrix, Multidimensional Scaling (MDS) and k-means or hierarchical clustering. The unsupervised mode does not require the number of clusters to be known, thanks to the gap statistic, and inherit of main algorithmic properties of the supervised mode, hence allowing (almost) any type of variable.   
}
\usage{
\method{unsupervised}{randomUniformForest}(object,
	baseModel = c("proximity", "proximityThenDistance",  "importanceThenDistance"),
	endModel = c("MDSkMeans", "MDShClust", "MDS"),
	endModelMetric = NULL,
	samplingMethod = c("uniform univariate sampling", 
	"uniform multivariate sampling", "with bootstrap"),
	MDSmetric = c("metricMDS", "nonMetricMDS"),
	proximityMatrix = NULL,
	outliersFilter = FALSE, 
	Xtest = NULL, 
	predObject = NULL, 
	metricDimension = 2, 
	coordinates = c(1,2),
	bootstrapReplicates = 100,
	clusters = NULL,
	maxIters = NULL,
	importanceObject = NULL,
	maxInteractions = 2,
	reduceClusters = FALSE, 
	maxClusters = 5,
	mapAndReduce = FALSE,
	OOB = FALSE,
	subset = NULL, 
	seed = 2014,
	uthreads = "auto",
	\dots)	
	\method{print}{unsupervised}(x, \dots)
	\method{plot}{unsupervised}(x, importanceObject = NULL, xlim = NULL, ylim = NULL, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x, object}{
	    an object of class randomUniformForest or a matrix (or data frame) from which the unsupervised mode will begin. If matrix or data frame, a full unsupervised learning object will be modelled, beginning by the computation of a randomUniformForest object (following Breiman's ideas). Note that in this latter case, data must be provided using 'Xtest' option. 
	}
	\item{xlim}{
		vector of 2 points or NULL. Plotting option for zooming in the graph.
	}
	\item{ylim}{
		vector of 2 points or NULL. Plotting option for zooming in the graph.
	}
	\item{baseModel}{
	    'baseModel' defines the way that algorithm will compute dissimilarity matrix. If 'proximity', a proximities matrix of observations will be computed. It can be square, if there is enough memory, or compressed up to a low dimension matrix, using trees and biggest proximities. The matrix is, then, sent to the Multidimensional scaling (MDS) function. If 'proximitythenDistance', the MDS function computes a distance first before applying scaling. If 'importanceThenDistance', then an instance of the importance.randomUniformForest() function will be used. This latter is an alternative of proximities matrix and is useful for difficult problems or when dimension of the problem is high. Note that for all cases, the matrix outputted can  be compressed up to a 2 dimensional matrix without loosing much information. For large data sets, compression automatically happens.
	}
	\item{endModel}{
	   in all cases, dimension reduction is always done. The MDS process is one of the two engines to achieve the task and cannot be overrided. It will send its output to K-means algorithm, if 'MDSkmeans', or to hierarchical clustering, if 'MDShClust'.
	}
	\item{endModelMetric}{
	   the metric or algorithm that has to be used when computing 'endModel'. See 'algorithm' in kmeans() function.
	}
	\item{samplingMethod}{
	    method that has to be used to turn unsupervised problem into a supervised one. 'samplingMethod' uses either one (then bootstrap will not happen) or two arguments, the second one always being "with bootstrap" allowing, then, to use bootstrap. Note that 'samplingMethod' is identical to 'unsupervisedMethod' in \code{\link{randomUniformForest}}. See Details section for more information.
	}
	\item{MDSmetric}{
	    one of the metric or non-metric to achieve MDS. See cmdscale() function and isoMDS() one of the MASS R package.
	}
	\item{proximityMatrix}{
	    proximities matrix, of size 'n x n', can be provided, coming from a previous unsupervised.randomUniformForest object or can be external.
	}
	\item{outliersFilter}{
	    if TRUE, outliers in the MDS output will be filtered to achieve clustering. Currently experimental.
	}
	\item{Xtest}{
	 	if one provides a randomUniformForest object for first argument, then 'Xtest' must be filled with a test sample in order to achieve unsupervised learning. In this case the test sample will be clustered with the same classes as in the randomUniformForest object. For example, one can consider 'Xtest' option as a way to add more training cases independently to the randomUniformForest classifier. The main paradigm of this situation is when only a few cases are labelled while the others are not and can (or do) not be used as a validation set. If first argument of unsupervised.randomUniformForest() function is already a test sample, then 'Xtest' is not mandatory but one can, however, put here another test sample in order to know how the model will react. More precisely, first argument will act, as a training sample, in all the layers of the unsupervised learning model while the sample in 'Xtest' will only have effects on the high-level layers, namely the MDS function and the clustering algorithm.
	}
	\item{predObject}{
	    one may provide here a full prediction object, coming from the randomUniformForest predict function (with option type='all') and any data set. Two cases are possible. The first one involves a supervised problem (coming with a randomUniformForest at first argument of the unsupervised.randomUniformForest() function) that one wants to cluster the predictions. Hence, one will need to put the dataset using 'Xtest' and its predictions in 'predObject'. The second case arises when one wants to speed up computations or exploit the space with a single supervised randomUniformForest model. Any new dataset (to put in 'Xtest') and its prediction object will be assessed by the unsupervised mode partially (using 'MDS' option) or totally. The output will be able to update the model in a supervised manner.
	}
	\item{metricDimension}{
	    for the MDS process, the dimension of the final data. Actually, a 'n x p' matrix (or data frame) used for unsupervised learning will end with a 'n x 2' matrix to be clustered. Hence, visualization will be made easy  and one will not have to worry about features, since the whole process will match full details which can be sent back to the randomUniformForest algorithm which can be reused to analyse both data and features. Another essential argument is the wish to avoid dealing constantly with large data sets, especially when combining two unsupervised objects. Usually a 'metricDimension' higher than 2 is useless.
	}
	\item{coordinates}{
	    MDS coordinates to plot. First and second are usually enough since they concentrate almost all the information and are not correlated.
	}
	\item{bootstrapReplicates}{
	    bootstrap replicates to be sent to the gap statistic. This latter is used to automatically find the number of clusters (if 'MDSkmeans' is called). See the clusGap() function in the cluster R package.
	}
	\item{clusters}{
	    number of clusters to find, if one knows it.
	}
	\item{maxIters}{
	    number of iterations of the K-means algorithm. See kmeans() function.
	}
	\item{importanceObject}{
	    if a object of class importance is provided and 'baseModel' option is 'importanceThenDistance', then the dissimilarity matrix will
		be computed using it. Useful if current data are not enough clean and one has a former trustfully randomUniformForest object. Note that dissimilarity matrix based on importance objects are recommended only when 'proximity' or 'proximityThenDistance' have some problems.
	}
	\item{maxInteractions}{
	    maximum number of interactions when computing an object of class importance. More means lots of details (and possibly noise). Less means
		faster computation and possibly less reliable results. 
	}
	\item{reduceClusters}{
	    experimental approach to reduce the number of clusters. One should use the \code{\link{modifyClusters}} function, a user-friendly and reversible algorithm to modify the number of clusters.
	}
	\item{maxClusters}{
	    maximum number of clusters to find. Used in conjunction with the clusGap() function from the R package "cluster".
	}
	\item{mapAndReduce}{
	    for large data sets, this option maps (internally) the whole unsupervised learning process in order to decrease memory footprint and to reduce computation time. To avoid computing a 'n x n' proximity matrix, only a subsample of the data is first used. After then the sampled proximity matrix is sent to the MDS process that generates MDS points. These ones are learned and all others MDS points will be predictions. Note that if the data have more than 2000 rows then 'mapAndReduce' will be partially enabled in order to avoid memory leaks while trying to preserve as much as possible the loss of accuracy that comes from the prediction task.
	}
	\item{OOB}{
	    see \code{\link{clusteringObservations}}
	}
	\item{subset}{			
		an index vector indicating which rows should be used.
	}
	\item{seed}{
		seed value to allow convergence (but not reproducibility since random Uniform forests are stochastic but converge). See Details section for more information.
	}
	\item{uthreads}{
	    allows multi-core computation using, by default and for the current computer, all logical cores minus 1 for each task that can be done using many cores in parallel. 
	}
	\item{\dots}{
	  allow all options of \code{\link{randomUniformForest}} to be passed to the unsupervised.randomUniformForest() function, e.g. ntree, mtry,  nodesize, ...
	}
}

\details{
 The unsupervised mode of Random Uniform Forests is designed to provide dimension reduction, clustering and a full analysis of features and observations. The process uses a tree-layer engine built around a randomUniformForest object.

 It can be summarized using the following chain :\cr
 RandomUniformForest object --> dissimilarity matrix --> multidimensional scaling --> clustering algorithm --> clusters --> that can be computed into an object of class importance.randomUniformForest. This latter is, optionally, used to analyse features, their links with clusters and the links between observations, features and clusters.
 
 First step involves Breiman's ideas. Since Random Uniform Forests inherit of all properties of Random Forests, they are able to implement the key concepts provided by Breiman for the unsupervised case : \cr
 - create a synthetic data set by scrambling independently (for example uniformly) each column of the original dataset, \cr
 - merge both synthetic and original dataset and give each one a label, \cr
 - run Random (Uniform) Forests to this new dataset and retrieve OOB errors, \cr
 - the lower the errors, the more clustering will be easy. If error is too high, say close to 50 percent, then there is (almost) no hope \cr
 - once the forest classifier is built, one can now move to the second step.
 
 The second step used proximities matrix. If data are a 'n x p' matrix (or data frame) then proximities matrix will have a 'n x n' size, meaning that for each observation, one will search, for each tree, which other observation falls in the same terminal node, then increase the proximity of the two observations by one. At the end, all observations and all trees will have been processed, leading to a proximities matrix that will be normalized by the number of trees. Since it can be very large, it is possible (but currently disabled due to the compromise that is needed with large datasets between accuracy and computing time) to use a 'n x B' proximities matrix, where 'B' is the number of trees. A further step is to use a 'n x q' (biggest) proximities matrix, where 'q' can be as small as 2 in order to compress the data to their maximum.
 
 Once proximities matrix (or dissimilarity matrix using for example '1 - proximities') has been computed, the third step is to enter in the MDS process. MDS is needed for dimension reduction (if not already happened) and mainly to generate decorrelated components in which the points will reside. The two first components are usually enough to get good visualization. Unlike PCA, they are used only to achieve the best possible separation between points. Note that we allow proximities to be transformed into distances since we found that it produces sometimes better results.
 
 The next step concludes the three-layer engine by calling a clustering algorithm, preferably K-means, to partition the MDS points since we use the features only in the early phase of the algorithm. Hence, we manipulate coordinates and points in a new space where MDS is the rule that matters. K-means algorithm is simply a way to provide a measurable way of the clusters which already exist, since clustering is almost done earlier. Note that the number of clusters is automatically adjusted by the gap statistic. If not enough, clusters structure can be instantaneously modified (see \code{\link{modifyClusters}}), letting the silhouette coefficient (Rousseeuw) tell the last word.
 
 The unsupervised learning is then partially achieved. An object with all tools is generated and can be returned in order to be assessed by the randomUniformForest algorithm that will provide a deeper analysis and visualization tools : \cr
 - what are the important features, e.g. the most discriminant ? \cr
 - what are their interactions, e.g. how features relate to the whole clustering scheme ? \cr
 - Partial dependencies. \cr
 - Links between observations and features. \cr
 - ... \cr
 
The unsupervised model is turned into a supervised one using the \code{\link{as.supervised}} function, with all options one need to call for the supervised case, then doing the analysis by the \code{\link{importance.randomUniformForest}} function. This last step is essential if one wants to access to the full details. Moreover, data are now turned into a training sample for the algorithm. If new data become available, one may use the supervised case or the unsupervised one. The former has the great advantage to have a complexity in O(B*n*log(n)). Indeed, due to the lot of computation involved, the unsupervised case requires much more time than a simple call to K-means algorithm but also provides more details.\cr 
In comparison with a faster algorithm and for a small dataset, the unsupervised mode of Random Uniform Forests
always provides :\cr
- dimension reduction that leads to a visualization in the MDS space (in two dimensions),\cr
- a resulted clustering scheme that processes MDS points, thus integrating all features in the same representation,\cr
- an object that can be used to evaluate and analyse most important features, their interactions and their links with observations, according to the clustering scheme.

When going toward large datasets, the unsupervised mode becomes hybrid with a fully unsupervised mode for a subsample of data and a supervised mode that learns MDS points and predict them for the remaining rows of the sample. This allows to strongly reduce computation time and to make the resulted object recyclable.\cr
- Hence, new data can be learned incrementally (combining or updating objects) and/or separately to the former object.
- The main argument resides in the fact that the combination of Random Uniform Forests and MDS space allows data to be dynamically clustered as the time is going (i.e. new data are coming) or in the space, by changing, on the fly, the clustering representation.

Note that the whole engine is stochastic, with almost no possibility of reproducibility using the set.seed() function. However since random Uniform Forests converge, seed has been added to the second primary part of the unsupervised mode, the creation of the synthetic dataset. It means that most of the work to achieve a good clustering representation is devoted to the randomUniformForest part, allowing one to assess parameters independently for each layer and look for the ones that have the main effect on the clustering. 
}

\value{
  An object of class unsupervised, which is a list with the following components:
  
  \item{proximityMatrix}{the resulted dissimilarity matrix.}
  \item{MDSModel}{the resulted Multidimensional scaling model.}
  \item{unsupervisedModel}{the resulted unsupervised model with clustered observations in unsupervisedModel$cluster.}
  \item{largeDataLearningModel}{ if the dataset is large, the resulted model that learned a sample of the MDS points, and predicted others points.}
  \item{gapStatistics}{if K-means algorithm has been called, the results of the gap statistic. Otherwise NULL.}
  \item{rUFObject}{Random Uniform Forests object.}
  \item{nbClusters}{Number of clusters found.}
  \item{params}{options of the model.}
}

\references{
 Abreu, N., 2011. Analise do perfil do cliente Recheio e desenvolvimento de um sistema promocional. Mestrado em Marketing, ISCTE-IUL, Lisbon
 
 Breiman and Cutler web site : http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
   
 Cox, T. F., Cox, M. A. A., 2001. Multidimensional Scaling. Second edition. \emph{Chapman and Hall}. 
 
 Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika} 53, 325-328.
 
 Kaufman, L., Rousseeuw, P.J., 1990. Finding Groups in Data: An Introduction to Cluster Analysis (1 ed.). New York: John Wiley. 
 
 Lloyd, S. P., 1957, 1982. Least squares quantization in PCM. Technical Note, Bell Laboratories. Published in 1982 in \emph{IEEE Transactions on Information Theory} 28, 128-137.
 
 Murtagh, F., Legendre, P., 2013. Ward's hierarchical agglomerative clustering method: which algorithms implement Ward's criterion? \emph{Journal of Classification} (in press).
 
 Tibshirani, R., Walther, G., Hastie, T., 2001. Estimating the number of data clusters via the Gap statistic. \emph{Journal of the Royal Statistical Society} B, 63, 411-423.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{modifyClusters}}, \code{\link{mergeClusters}},
 \code{\link{clusteringObservations}}, \code{\link{as.supervised}}
}

\examples{
## not run
## 1 - the famous iris dataset

## load data
# data(iris)

## run unsupervised modelling, removing labels :
## Default options, letting the 'gap statistic' find the number of clusters, 
## except for 'baseModel' which is slightly more efficient 
## with the "proximityThenDistance" argument

# iris.rufUnsupervised = unsupervised.randomUniformForest(iris[,-5], 
# baseModel = "proximityThenDistance", seed = 2014, threads = 1)

## view a summary
# iris.rufUnsupervised

## one may assess the gap statistic by calling the 'modifyClusters( )' function,  
## increasing or decreasing the number of clusters and looking the variations 
## of the silhouette coefficient.
## For example, if 4 clusters are found (since we know there are 3) :

# iris.rufUnsupervised2 = modifyClusters(iris.rufUnsupervised, decreaseBy = 1)

## plot clusters 
# plot(iris.rufUnsupervised)

## 2 - Full example with details
## Wholesale customers data (UCI machine learning repository)

# URL = "http://archive.ics.uci.edu/ml/machine-learning-databases/00292/"
# datasetName = "Wholesale%20customers%20data.csv"
# wholesaleCustomers = read.csv(paste(URL, datasetName, sep =""))

## modelling, letting the algorithm deal with all problems :
## categorical features, number of clusters, dimension reduction, visualization,
## variable importance, links between features and observations,...

# wholesaleCustomers.rufUnsupervised = unsupervised.randomUniformForest(wholesaleCustomers, 
# nodesize = 10, bagging = TRUE, ntree = 200, categoricalvariablesidx = "all")

## assess quality of the clustering :
## (and change eventually model parameters, e.g. 'baseModel' or 'endModel', 
## running again the model to get a better clustering one, looking the average silhouette 
## or the distance between clusters)

# wholesaleCustomers.rufUnsupervised

## visualization : at first, only clusters 
# plot(wholesaleCustomers.rufUnsupervised)

## but, we may need more :
## get details, turning first the model in a supervised one

# wholesaleCustomers.rufSupervised = as.supervised(wholesaleCustomers.rufUnsupervised, 
# wholesaleCustomers, bagging = TRUE, ntree = 200, 
# nodesize = 10, categoricalvariablesidx = "all")

## Is the learning efficient (using OOB evaluation) ?
# wholesaleCustomers.rufSupervised

## get variable importance, leading to a full analysis and visualization
## while limiting interactions of variables to 3 orders

# wholesaleCustomers.importance = importance(wholesaleCustomers.rufSupervised, 
# Xtest = wholesaleCustomers, maxInteractions = 3)

## a - visualize : features, interactions, partial dependencies, features in clusters
## NOTE : tile window in the R menu to see all plots. Loop over the prompt to see
## all matched partial dependencies

# plot(wholesaleCustomers.importance, Xtest = wholesaleCustomers)

## we get global variable importance (information gain), interactions, partial dependencies,
## and variable importance over labels. See vignette for more details.

## b - more visualization : (another look on 'variable importance over labels')
# featuresCluster1 = partialImportance(wholesaleCustomers, wholesaleCustomers.importance, 
# whichClass = 1)

## c - visualization : clusters and most important features
# plot(wholesaleCustomers.rufUnsupervised, importanceObject = wholesaleCustomers.importance)

## d - table : see individual links between observations and features
## the table show each observation with its associated features 
## and their frequencies of occurrence

# featuresAndObs = as.data.frame(wholesaleCustomers.importance$localVariableImportance$obs)
# frequencyFeaturesIdx = grep("Frequency", colnames(featuresAndObs))
# featuresNames = apply(featuresAndObs[,-c(1,frequencyFeaturesIdx)], 2, 
# function(Z) colnames(wholesaleCustomers)[Z])
# featuresAndObs[,-c(1,frequencyFeaturesIdx)] = featuresNames

# head(featuresAndObs)

## NOTE : since features are almost in monetary units, one may assess clusters 
## by looking the sum of all features per cluster and turn the problem 
## into a 'revenues per cluster and feature' one that can be linked 
## with the clustering process and visualization tools.

## first, merge outliers and retrieve clusters
# Class = mergeOutliers(wholesaleCustomers.rufUnsupervised)

## then add classes
# wholesaleCustomersClusterized = cbind(wholesaleCustomers, Class)

## finally compute revenues per cluster and feature.
## Note that this view may give more insights on how the algorithm clusters data.

# revenuePerClusterAndFeature = 
# aggregate(wholesaleCustomersClusterized[,-c(1,2,9)], list(Class), sum)

## see results
# revenuePerClusterAndFeature

## revenuePerCluster : leading to know where and how more work might happen...
# rowSums(revenuePerClusterAndFeature[,-1])
}

\keyword{unsupervised}
\keyword{learning}
\keyword{clustering}
\keyword{dimension}
\keyword{reduction}