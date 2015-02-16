\name{clusterAnalysis}
\alias{clusterAnalysis}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cluster (or classes) analysis of importance objects. 
}
\description{
   Provides a full analysis of clustered objects in a compact and granular representation. More precisely, observations, features and clusters are analysed in a same scheme, leading to unify and interpret all results of the unsupervised mode in one way. The function prints at most 5 tables and is designed to be an extension of importance object and also
   works in the supervised mode.
}
\usage{
clusterAnalysis(object, X, 
components = 2, 
maxFeatures = 2, 
clusteredObject = NULL, 
categorical = NULL,
OOB = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
	    an object of class 'importance'.
	}
	\item{X}{
	   the original data frame (or matrix) used to compute 'importance'.
	}
	\item{components}{
	  number of components to use for assessing cluster (or class). It has the same purpose than eigenvectors in PCA.
	  Note that the number must not excess the value of 'maxInteractions' in \code{\link{importance.randomUniformForest}}.
	}
	\item{maxFeatures}{
	  number of components to print for the visibility of analysis.
	}
	\item{clusteredObject}{
	  in the case of unsupervised learning, the clustered object coming from \code{\link{unsupervised.randomUniformForest}} or it could be class predictions or labels of the training set (whose, for both, size should be equal to the number of rows of 'X') or an object of class randomUniformForest
	  from which the function will search either OOB predictions or test set predictions, the former having priority.
	}
	\item{categorical}{
	  which variables (by names or subscripts) should be considered as categorical? It is currently recommended to let
	  the algorithm find it.
	}
	\item{OOB}{
	  if TRUE and clusteredObject contains an OOB object, compute the analysis using the OOB classifier. Useful when there are no test labels. 
	}
}

\details{
	'clusterAnalysis( )' function is designed to provide a straightforward interpretation of a clustering problem or, in classification, the relation found by the model between features, observations and labels.
	The main ingredient of the function comes from an object of class importance, which gives the local importance and its derivatives (see Ciss, 2015b). For each observation, one records the features, over all trees, that comes in
	first, second, ..., and so on position in the terminal nodes where falls the observation. One also records the frequency (number of times divided by number of trees) at which each of the recorded features appears. This gives
	the first outputted table ('featuresAndObs' see below), with all observations (in rows) and in columns their class or cluster, their position (from first one to last one) and frequency. Note that a position is a local component
	, named 'localVariable', e.g. 'localVariable1' means that for the i-th observation, the recorded variable was the one 
	that had the highest frequency in the terminal node where the observation has been dropped.
	
	The resulted table is then aggregated, providing in rows, the cluster (or the class) predicted by the model and
	in columns a component (with the same purpose than an eigenvector), that aggregates variables recorded for all observations of the dataset. We call the table 'clusterAnalysis' since it connects clusters (or class) with dependent variables on each component. The table is granular since one may choose (using options) both number of components and number of variables that will be displayed allowing an easy viewing on variables that matter. In the table, variables in each component are displayed in an decreased order of influence and are usually co-dependent. 
	
	The importance (between 0 and 1) of each component is displayed in the third table 'componentAnalysis' allowing the user to define how many components are needed to explain the clusters. More means more variables and the assistance of object of class importance, friends and methods (plot, print, summary) will be needed. Less means easy, but possibly partial, interpretation. Note that the 'clusterAnalysis( )' function is, partially, the summarized view of \code{\link{importance.randomUniformForest}}. 
	
	The three resulted tables comes directly from the modelling and, in classification, do not use training labels
	but predictions (like clusters are predictions computed by the model). One may want to to know how the modelling is relying on the true nature of data. That is the main interest of the 'clusterAnalysis( )' function. The latter unifies modelling and observed points. The last tables, named 'numericalFeaturesAnalysis' and/or 'categoricalFeaturesAnalysis', aggregates (by a sum, mean, standard deviation or most frequent state) the values of each feature for each predicted cluster or class.\cr
	1 - If the clustering scheme is pertinent, then the table must show differences between clusters, by looking to the aggregated values of features.\cr 
	2 - It must provide insights, especially in a cost-benefit analysis, about the features one needs to take care in order to have more knowledge of each cluster (or class).\cr 
	3 - In a general case and in Random Uniform Forests, the modelling may not match, even is the clusters are well separated, either the true distribution (number of cases per cluster or class) or the true effects (relations between covariates, relations in each cluster,...) present in the data. The 'clusterAnalysis( )' function, unifying all provided tables, must ensure that:\cr
		3.a - true, or at least main or possible, effects are found,\cr
		3.b - the approximated distribution generated is enough close to the true one.\cr
	
	From the points above, the main argument is that Random Uniform Forests are almost stochastic. Hence random effects are first expected. If repeated modelling lead to similar results (clusters and within relations), then results are, at least, self-consistent and clusterAnalysis( ) provides a way to look how the modelling enlighten the data.
}


\value{
  A list of 3 to 5 data frames with the following components:
   \item{featuresAndObs}{ a data frame with on rows the observations and on columns a local 'component' (localVariable) representing a meta-variable that accumulates dependent covariates by their names. Components are sorted by level of interactions (from highest to the lowest) then after for each local component is quantified by the influence (named LocalVariableFrequency, from 0 to 1) it gets on each observation over all trees. Note that the data frame is called from 'object' and the dimension depends to the value of 'maxInteractions' when calling \code{\link{importance.randomUniformForest}} .}
   \item{clusterAnalysis}{ from the data frame above, an aggregated data frame giving a more or less, but compact, granular view of components elements (i.e. variables). Note that for each component, within variables have a co-dependence. Each row is a cluster. }
   \item{componentAnalysis}{ from the data frame on top, an aggregated data frame giving a more or less, but compact, granular view of component frequencies similar to a "percentage of variance explained" of an eigen value.
   Each row is a cluster. }
   \item{numericalFeaturesAnalysis}{ from the data frame on top and 'X', an aggregated data frame giving sum of values of each numerical variable within each cluster. }
   \item{categoricalFeaturesAnalysis}{ from the data frame on top and the data, an aggregated data frame giving the most frequent value of each categorical variable within each cluster.}
}

\references{
 Ciss, S., 2015b. Variable Importance in Random Uniform Forests. hal-1104751
 https://hal.archives-ouvertes.fr/hal-01104751
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\seealso{
 \code{\link{unsupervised.randomUniformForest}},  \code{\link{importance.randomUniformForest}}  
}

\examples{
## not run
## load iris data
# data(iris)

## A - Clustering
## run unsupervised modelling, removing labels and committing 4 clusters
# iris.rufUnsupervised = unsupervised.randomUniformForest(iris[,-5], mtry = 1, nodesize = 2)

## as supervised
# iris.rufUnsup2sup = as.supervised(iris.rufUnsupervised, iris[,-5], mtry = 1, nodesize = 2)

## importance
# iris.ruf.importance = importance(iris.rufUnsup2sup, Xtest = iris[,-5], maxInteractions = 4)

## cluster analysis
# iris.ruf.clusterAnalysis = clusterAnalysis(iris.ruf.importance, iris[,-5], components = 3)

## or
# iris.ruf.clusterAnalysis = clusterAnalysis(iris.ruf.importance, iris[,-5], components = 3, 
# clusteredObject = iris.rufUnsupervised)

## B - Classification (using OOB evaluation)
# iris.ruf = randomUniformForest(Species ~., data = iris, mtry = 1, nodesize = 2, threads = 1)

## importance
# iris.ruf.importance = importance(iris.ruf, Xtest = iris[,-5], maxInteractions = 4)

## analysis
# iris.ruf.clusterAnalysis = clusterAnalysis(iris.ruf.importance, iris[,-5], components = 3, 
# OOB = TRUE, clusteredObject = iris.ruf)
}

\concept{Variable Importance analysis}