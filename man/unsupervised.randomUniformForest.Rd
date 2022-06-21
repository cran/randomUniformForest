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
   Unsupervised mode of Random Uniform Forests is designed to provide, in all cases, clustering, dimension reduction, easy visualization, deep variable importance, relations between observations, variables and clusters. It also comes with two specific points: easy assessment (cluster analysis) and dynamic clustering, allowing to change on-the-fly any clustering shape. A three-layer engine is used: dissimilarity matrix, Multidimensional Scaling (MDS) or Spectral decomposition, and k-means or hierarchical clustering. The unsupervised mode does not require the number of clusters to be known, thanks to the gap statistic, and inherits of main algorithmic properties of the supervised mode, allowing (almost) any type of variable.   
}
\usage{
\method{unsupervised}{randomUniformForest}(object,
	baseModel = c("proximity", "proximityThenDistance", "importanceThenDistance"),
	endModel = c("MDSkMeans", "MDShClust", "MDS", "SpectralkMeans"),
	endModelMetric = NULL,
	samplingMethod = c("uniform univariate sampling", 
	"uniform multivariate sampling", "with bootstrap"),
	MDSmetric = c("metricMDS", "nonMetricMDS"),
	proximityMatrix = NULL,
	fullProximityMatrix = TRUE,
	sparseProximities = FALSE,
	outliersFilter = FALSE, 
	Xtest = NULL, 
	predObject = NULL, 
	metricDimension = 2,
	distance = c("euclidean", "maximum", "manhattan", "canberra", "binary",
	"minkowski"),
	coordinates = c(1,2),
	bootstrapReplicates = 100,
	clusters = NULL,
	maxIters = NULL,
	nstart = 1,
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
	\method{plot}{unsupervised}(x, importanceObject = NULL, xlim = NULL, ylim = NULL, coordinates = NULL, \dots)
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
	    'baseModel' defines the way that algorithm will compute dissimilarity matrix. If 'proximity', a proximities matrix of observations will be computed. It can be square, if there is enough memory, or compressed up to a low dimension matrix, using trees and biggest proximities (currently disabled). By default, the matrix is sent to the Multidimensional scaling (MDS) function. If 'proximitythenDistance', the MDS function computes a distance first before applying scaling. If 'importanceThenDistance', then an instance of the importance.randomUniformForest() function will be used. This latter is an alternative of proximities matrix and is useful for difficult problems or when dimension of the problem is high. Note that for all cases, the matrix outputted can  be compressed up to a 2 dimensional matrix without loosing much information. For large data sets, compression automatically happens.
	}
	\item{endModel}{
	   in all cases, dimension reduction is always done. The MDS process is one of the two engines to achieve the task and cannot be overrided (except if one specifies 'SpectralkMeans' in the option). It will send its output to K-means algorithm, if 'MDSkmeans', or to hierarchical clustering, if 'MDShClust'.
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
	\item{fullProximityMatrix}{
	    if TRUE, compute the 'n x n' proximities, otherwise compute the 'n x B' ones, where 'B' is the number of trees.
	}
	\item{sparseProximities}{
	    if TRUE, proximities are filled with 1 (which happens when an observation falls in the same terminal node than another, whatever the number of times is) and 0 (an observation does not fall in any terminal node than another). If FALSE, Breiman's formulation is used. Depending on the data, both can be useful.		
	}
	\item{outliersFilter}{
	    if TRUE, outliers in the MDS output will be filtered to achieve clustering. Currently experimental.
	}
	\item{Xtest}{
	 	if one provides a randomUniformForest object for first argument, then 'Xtest' must be filled with a test sample in order to achieve unsupervised learning. In this case the test sample will be clustered with the same classes as in the randomUniformForest object. For example, one can consider 'Xtest' option as a way to add more training cases independently to the randomUniformForest classifier. The main paradigm of this situation is when only a few cases are labelled while the others are not and can (or do) not be used as a validation set. If first argument of unsupervised.randomUniformForest() function is already a test sample, then 'Xtest' is not mandatory but one can, however, put here another test sample in order to know how the model will react. More precisely, first argument will act, as a training sample, in all the layers of the unsupervised learning model while the sample in 'Xtest' will only have effects on the high-level layers, namely the MDS (or spectral) function and the clustering algorithm.
	}
	\item{predObject}{
	    one may provide here a full prediction object, coming from the \code{\link{predict.randomUniformForest}} function (with option type='all') and any data set. Two cases are possible. The first one involves a supervised problem (coming with a randomUniformForest at first argument of the unsupervised.randomUniformForest() function) that one wants to cluster the predictions. Hence, one will need to put the dataset using 'Xtest' and its predictions in 'predObject'. The second case arises when one wants to speed up computations or exploit the space with a single supervised randomUniformForest model. Any new dataset (to put in 'Xtest') and its prediction object will be assessed by the unsupervised mode partially (using 'MDS' option) or totally. The output will be able to update the model in a supervised manner.
	}
	\item{metricDimension}{
	    for the MDS or spectral process, the dimension of the final representation. By default, a 'n x p' matrix (or data frame) used for unsupervised learning will end with a 'n x 2' matrix, with decorrelated components, that will be clustered. Using a low dimension avoids issues with large data sets, especially when combining two unsupervised objects. Usually a 'metricDimension' higher than 2 is useless except for spectral decomposition,
		where its value should be assessed.
	}
	\item{distance}{
		the chosen distance from which proximities will be transformed if MDS... is used as 'endModel'.
	}
	\item{coordinates}{
	    MDS (or eigenvectors of the spectral decomposition) coordinates to plot. First and second are usually enough since they concentrate almost all the information and are not correlated. In the Spectral case, coordinates should be chosen with care (using the plot method to assess them) since they change a lot the representation.
		Note that number of chosen coordinates should be small in order to avoid long computation for large datasets.
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
	\item{nstart}{
	    see kmeans() function.
	}
	\item{importanceObject}{
	    if a object of class importance is provided and 'baseModel' option is 'importanceThenDistance', then the dissimilarity matrix will be computed using it. Useful if current data are not enough clean and one has a former trustfully randomUniformForest object. Note that dissimilarity matrix based on importance objects are recommended only when 'proximity' or 'proximityThenDistance' have some issues.
	}
	\item{maxInteractions}{
	    maximum number of interactions when computing an object of class importance. More means lots of details (and possibly noise). Less means faster computation and possibly less reliable results. 
	}
	\item{reduceClusters}{
	    experimental approach to reduce the number of clusters. One should use the \code{\link{modifyClusters}} function, a user-friendly and reversible algorithm to modify the number of clusters.
	}
	\item{maxClusters}{
	    maximum number of clusters to find. Used in conjunction with the clusGap() function from the R package "cluster".
	}
	\item{mapAndReduce}{
	    for large data sets, this option maps (internally) the whole unsupervised learning process in order to decrease memory footprint and to reduce computation time. To avoid computing a 'n x n' proximity matrix, only a subsample of the data is first used. After then the sampled proximity matrix is sent to the MDS process that generates MDS points. These ones are learned and all others MDS points will be predictions. Note that if the data have more than 10000 rows then 'mapAndReduce' will be partially enabled and will try to preserve as much as possible the loss of accuracy that comes from the prediction task.
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
 The unsupervised mode of Random Uniform Forests is designed to provide dimension reduction, clustering, visualization  and a full analysis of features and observations. The process uses a tree-layer engine built around a randomUniformForest object. One may view the whole task as a kind of PCA. The main purpose of unsupervised Random Uniform Forests is to be highly versatile, allowing to assess datasets.

 It can be summarized using the following chain :\cr
 Unsupervised randomUniformForest object --> dissimilarity matrix --> multidimensional scaling or spectral decomposition --> clustering algorithm --> clusters --> that can be turn into a supervised object, see \code{\link{as.supervised}} and analysed in depth with \code{\link{importance.randomUniformForest}}, with many details, and \code{\link{clusterAnalysis}}, which is both compact and granular. 
 
 First step involves Breiman's ideas. Since Random Uniform Forests inherit of all properties of Random Forests, they are able to implement the key concepts provided by Breiman for the unsupervised case : \cr
 - create a synthetic data set by scrambling independently (for example uniformly) each column of the original dataset, \cr
 - merge both synthetic and original dataset and give each one a label, \cr
 - run Random (Uniform) Forests to this new dataset and retrieve OOB errors, \cr
 - the lower the OOB errors, the more clustering will be easy. If error is too high, say more than 50 percent, then there is (almost) no hope, \cr
 - once the forest classifier is built, one can now move to the second step.
 
 The second step used proximities matrix. If data are a 'n x p' matrix (or data frame) then proximities matrix will have a 'n x n' size, meaning that for each observation, one will search, for each tree, which other observation falls in the same terminal node, then increase the proximity of the two observations by one. At the end, all observations and all trees will have been processed, leading to a proximities matrix that will be normalized by the number of trees. Since it can be very large, it is possible (but currently disabled due to the compromise that is needed with large datasets between accuracy and computing time) to use a 'n x B' proximities matrix, where 'B' is the number of trees. A further step is to use a 'n x q' (biggest) proximities matrix, where 'q' can be as small as 2 in order to compress the data to their maximum. Note that proximities matrix has a second mode, for difficult cases, which produces a binary matrix.
 In our implementation, we considered all cases where n > 10,000 as problems that must be predicted, since the proximities matrix will, otherwise, have size greater than 800 Mo. In this case, learning points imply to learn MDS (or spectral) ones as responses and the data with n <= 10,000 as predictors. Then predict the remaining data.  
 
 Once proximities matrix (or dissimilarity matrix using for example '1 - proximities') has been computed, the third step is to enter in the MDS (or spectral) process. MDS is needed for dimension reduction (if not already happened) and mainly to generate decorrelated components in which the points will reside. The two first components are usually enough to get good visualization. Unlike PCA, they are used only to achieve the best possible separation between points. Note that we allow proximities to be transformed into distances since we found that it produces sometimes better results.
 In the case of spectral decomposition, two eigenvectors are usually enough to get a good clustering scheme, and the algorithm allows to arbitrary choose them, while the default option limits the number of eigenvectors.
 
 The next step concludes the three-layer engine by calling a clustering algorithm, preferably K-means, to partition the MDS (or spectral) points since we use the features only in the early phase of the algorithm. Hence, we manipulate coordinates and points in a new space where MDS (or spectral) is the rule that matters. K-means algorithm is simply a way to provide a measurable way of the clusters which already exist, since clustering is almost done earlier. Note that the number of clusters can be automatically adjusted by the gap statistic (Tibshirani, Walther, Hastie, 2001). If not enough, clusters structure can be instantaneously modified, merged or split (see \code{\link{modifyClusters}}, \code{\link{mergeClusters}}, \code{\link{splitClusters}}) letting the silhouette coefficient (Rousseeuw, 1987) have the last word.
 
 The unsupervised learning is then partially achieved. An object with all tools is generated and can be returned in order to be assessed by the \code{\link{randomUniformForest}} algorithm that will provide a deeper analysis and visualization tools : \cr
 - most influential features,\cr
 - interactions,\cr
 - partial dependencies,\cr
 - dimension reduction,\cr
 - visualization,\cr
 - cluster analysis,\cr
 - ... \cr
 
The unsupervised model is turned into a supervised one using the \code{\link{as.supervised}} function, with all options one need to call for the supervised case, then doing the analysis by the \code{\link{importance.randomUniformForest}} function to get full details. Moreover, data are now turned into a training sample for the algorithm. If new data become available, one may use the supervised case or the unsupervised one using \code{\link{update.unsupervised}}. The former has the great advantage to have a complexity, when predicting, in O(B*n*log(n)), where 'B' is the number of trees. Indeed, due to the lot of computation involved, the unsupervised case requires much more time than a simple call to K-means algorithm but also provides much more details and more abilities to cluster.\cr 

When going toward large datasets, the unsupervised mode becomes hybrid with a fully unsupervised mode for a random subsample of the data and a supervised mode that learns MDS (or spectral) points and predict them for the remaining rows of the sample. This allows to strongly reduce computation time and to make the resulted object recyclable. New data can be learned incrementally (combining or updating objects) and/or separately to the former object. 

Note that the whole engine is stochastic, with almost no possibility of exact reproducibility using the set.seed() function. However, since Random Uniform Forests converge, seed has been added to the second part of the first layer, the creation of the synthetic dataset. It means that most of the work to achieve a good clustering representation is devoted to the \code{\link{randomUniformForest}}  part, allowing one to assess parameters independently for each layer and look for the ones that have the main effect on the clustering. Moreover, the stochastic nature is the main argument of unsupervised learning in Random Uniform Forests, since it is first expected that random models produce random effects. As a consequence, the model is first self-consistent. Efficiency of results can be assessed independently, using \code{\link{clusterAnalysis}} which makes the bridge between clustering and original data.
}

\value{
  An object of class unsupervised, which is a list with the following components:
  
  \item{proximityMatrix}{the resulted dissimilarity matrix.}
  \item{MDSModel}{the resulted Multidimensional scaling or spectral decomposition model.}
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
 
 Ng, A. Y., Jordan, M. I., Weiss, Y., 2002. On spectral clustering: Analysis and an algorithm. \emph{Advances in neural information processing systems}, 2, 849-856.
 
 Rousseeuw, P. J., 1987. Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis. \emph{Computational and Applied Mathematics} 20, 53-65.
 
 Tibshirani, R., Walther, G., Hastie, T., 2001. Estimating the number of data clusters via the Gap statistic. \emph{Journal of the Royal Statistical Society} B, 63, 411-423.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{modifyClusters}}, \code{\link{mergeClusters}}, \code{\link{splitClusters}}
 \code{\link{clusteringObservations}}, \code{\link{as.supervised}}, 
 \code{\link{update.unsupervised}}, \code{\link{combineUnsupervised}},
 \code{\link{clusterAnalysis}}
}

\examples{
## not run
#### A - Overview : the famous iris dataset

## load data
# data(iris)

## run unsupervised modelling, removing labels :
## Default options, letting the 'gap statistic' find the number of clusters.

## Note: the stochastic nature of the algorithm acts as an hidden parameter 
## whose one has to take care of.

## default (MDS) clustering
# iris.uruf = unsupervised.randomUniformForest(iris[,-5])

## for iris data, default option is not stable.
 
## Note: correlation between covariates play a key role on stability and clustering.
## Neutralize correlation requires to choose splitting variable at random.

## Increase 'nodesize' to get more stability
# iris.uruf2 = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2)

## would be better and stable (if repeated)

## Regular companions (from point 1 to point 14):
## 1 - Summary with raw assessment
# iris.uruf2

## 2 - Plot in the new space in 2D
# plot(iris.uruf2)

## 3 - Modify clusters (increasing or decreasing the number of clusters 
## and looking the variations of the silhouette coefficient).
## For example, if 4 clusters are found (since we know there are 3) :

# iris.uruf4Clusters = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2, 
# clusters = 4)
# iris.uruf3Clusters = modifyClusters(iris.uruf4Clusters, decreaseBy = 1)

## 4 - (or) Merge clusters if there are too many 
## merge the second and the third

## Note: one can play with modify/merge/split and 'plot' many times in order to improve
## the clustering
# iris.urufmerge = mergeClusters(iris.urufnsupervised4Clusters, c(1,2))  

## 5 - Alternatives : try Spectral clustering
## 'coordinates' specify which eigenvectors have to be sent to the clustering algorithm
## 'metricDimension' states for the maximal number of eigenvectors to use.
## Try 'coordinates = c(2,3)' and assess it by visualization

# iris.urufSpectral = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2, 
# endModel = "SpectralkMeans", metricDimension = 3, coordinates = c(2,3))

# plot(iris.urufSpectral)

## or try them all and choose/visualize the best representation
# iris.urufSpectral = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2, 
# endModel = "SpectralkMeans", metricDimension = 5, coordinates = c(1:5))

# plot(iris.urufSpectral, coordinates = c(1,2))
# plot(iris.urufSpectral, coordinates = c(2,3))
## ...

#### B - Full example with details
## Wholesale customers data (UCI machine learning repository)

# URL = "http://archive.ics.uci.edu/ml/machine-learning-databases/00292/"
# datasetName = "Wholesale%20customers%20data.csv"
# wholesaleCustomers = read.csv(paste(URL, datasetName, sep =""))

## modelling : three ways  
## - let the algorithm deal with all layers: categorical features, number of clusters,...
## - control the first layer to get stability (proximities matrix) and let the
## algorithm control the others.
## - control all manually.

## first way
# wholesaleCustomers.uruf.1 = unsupervised.randomUniformForest(wholesaleCustomers)

## assess quality of the clustering :
# wholesaleCustomers.uruf.1

## second way (we keep it) : use 'bagging' option to let each feature define a candidate node,
## use 'sparseProximities' that provides a better separation

## third way (seems the best) : spectral clustering, specifying clusters and changing seed
## to stabilize the clustering scheme when modelling is repeated.

# wholesaleCustomers.uruf.3 = unsupervised.randomUniformForest(wholesaleCustomers, 
# ntree = 500, sparseProximities = TRUE, endModel = "SpectralkMeans", bagging = TRUE, 
# categorical = c("Channel", "Region"), clusters = 3, seed = 2016)

## Speed up computation: less trees + increasing (minimal) node size + do not use logical cores
# wholesaleCustomers.uruf.3 = unsupervised.randomUniformForest(wholesaleCustomers, 
# ntree = 100, nodesize = 10, sparseProximities = TRUE, endModel = "SpectralkMeans", 
# bagging = TRUE, categorical = c("Channel", "Region"), clusters = 3, seed = 2016)

## Note: 'Channel' and 'Region' are categorical and differ from other products (consumer goods)
## At first, we keep all, seeking links between the context of consumption and the products.

# wholesaleCustomers.uruf.2 = unsupervised.randomUniformForest(wholesaleCustomers, 
# ntree = 500, sparseProximities = TRUE, bagging = TRUE, categorical = c("Channel", "Region"))

## Regular companions :
## 6 - Assess the randomUniformForest object (low OOB error is usually better but not a rule)
# wholesaleCustomers.uruf.2$rUF

## 7 - Look for influential variables (before clustering)
# summary(wholesaleCustomers.uruf.2$rUF)

## assess quality of the clustering and remove/add/merge clusters to see if things are better
# wholesaleCustomers.uruf.2
# plot(wholesaleCustomers.uruf.2)

## 8 - Get details : turning first the model into a supervised one
# wholesaleCustomers.ruf.2 = as.supervised(wholesaleCustomers.uruf.2, wholesaleCustomers, 
# ntree = 500, categorical = c("Channel", "Region"), BreimanBounds = FALSE)

## 9 - Assess the 'Learning clusters' process
# wholesaleCustomers.ruf.2

## 10 - Get variable importance : leading to a full analysis and visualization
# wholesaleCustomers.importance = importance(wholesaleCustomers.ruf.2, 
# Xtest = wholesaleCustomers, maxInteractions = 8)

## Visualize all: features, interactions, partial dependencies,...
## Note: tile window in the R menu to see all plots. Loop over the prompt to see
## all matched partial dependencies

# plot(wholesaleCustomers.importance, Xtest = wholesaleCustomers, whichOrder = "all")
## For details, type vignette("VariableImportanceInRandomUniformForests", 
## package = "randomUniformForest")
## Note : Variable importance strongly depends on the model and the chosen clusters

## 11 - more visualization: (another look on 'variable importance over labels')
## for each cluster, we can see which variables matter and their order.
# featuresCluster1 = partialImportance(wholesaleCustomers, 
# wholesaleCustomers.importance, whichClass = 1)
# ...

## 12 - Refining analysis : two possible levels, with and without Region (and Channel)

## 12.a - visualizing clusters and features: keeping Region and Channel
# plot(wholesaleCustomers.uruf.2, importanceObject = wholesaleCustomers.importance)

## 12.b - Removing Region and Channel to better assess others products:
## one may want to lead an analysis independent to the context

## If Region and/or Channel matter:
## - Solution 1: reach stability by first trying and repeating full random models, 
## leading to a new model (with variables of interest) and the former one. 
## - Solution 2: clustering with all features (as above), 
## assessment on variables of interest handled by the importance() function.
## - Solution 3: clustering with all features (as above), clusterAnalysis() function 
## is able to provide a granular view. We choose this solution.

## 13 - Cluster analysis: aggregated links between observations, features and clusters
## clusterAnalysis mainly shows what is happening when looking deeper in details

## 13.a - first, look inside clusters 
# wholesaleCustomersFinalSummary.ruf = clusterAnalysis(wholesaleCustomers.importance, 
# wholesaleCustomers, components = 3, maxFeatures = 4)

## 13.b - Numerical and categorical features aggregation
## same function, more options
## Note : here features are valuable, revenue for each. Hence one can exploit that.

# wholesaleCustomersFinalSummary.ruf = clusterAnalysis(wholesaleCustomers.importance, 
# wholesaleCustomers, components = 3, maxFeatures = 4, 
# clusteredObject = wholesaleCustomers.uruf.2, categorical = c("Channel", "Region"))

## 14 - Conciliate analysis
## clusterAnalysis( ) provides both influential features (leading to a better clustering)
## and valuable ones. Due to the purpose of the task, both types can be different
## but should give insights about the next task. 
## Suppose one wants to maximize revenues. Where could one put the efforts ?

## First, look to the global point of view (clusters) : 
# revenuePerCluster = 
# rowSums(wholesaleCustomersFinalSummary.ruf$numericalFeaturesAnalysis[,-c(1,2)])

## Then, look to the local point of view :
## for example, influential features in one cluster might be the ones that one need
## to develop in  another one... avoiding brute force techniques

## Notes:
## 1 - cluster analysis must be consistent with results of importance( ) function.
## If not, then something is going wrong in the modelling.
## 2 - To complete an analysis, one should take care of the whole process.
## For example, changing the clustering may change a lot the analysis.
}

\concept{unsupervised learning}
\concept{clustering}
\concept{dimension reduction}
\concept{variable importance}