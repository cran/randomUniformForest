\name{reSMOTE}
\alias{reSMOTE}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  REplication of a Synthetic Minority Oversampling TEchnique for highly imbalanced datasets
}
\description{
  Produce new samples of the minority class by creating synthetic data of the original ones by randomization. After that, the new dataset contains all original data + synthetic data labelled with the minority class. Main argument is derived from Breiman's ideas on the construction of synthetic data for the unsupervised mode of random (Uniform) Forests. The new dataset preserves distribution of the covariates of the original dataset (at least with default 'samplingMethod' option).
}
\usage{
reSMOTE(X, Y, 
	majorityClass = 0, 
	increaseMinorityClassTo = NULL, 
	samplingMethod = c("uniform univariate sampling", 
		"uniform multivariate sampling", "with bootstrap"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{X}{
	    a matrix (training sample) of numeric values.
	}
	\item{Y}{
	   labels (vector or factor) of the training sample.
	}
	\item{majorityClass}{
	  value(s) of the majority class(es).
	}
	\item{increaseMinorityClassTo}{
	  the percentage of minority class one would want to reach in the original sample for a less difficult classification task. Should be less than 50 percent, but enough more than the original percentage of the minority class. 
	}
	\item{samplingMethod}{
	  which method to use for sampling from data ? Note that one may use any of the two first method with the 'with bootstrap' option, for example 'samplingMethod = c("uniform univariate sampling",  "with bootstrap")' or 'samplingMethod = "uniform univariate sampling"'.
	}
}

\details{
	The main purpose of reSMOTE( ) function is to increase AUC. The function is specifically designed for
	random Uniform Forests with possible conjunction of options 'classwt', 'rebalancedsampling' and/or 'oversampling', also designed to treat imbalanced classes. Since ensemble models tend to choose majority class too quickly, reSMOTE introduces a perturbation that forces the algorithm to learn data more slowly. Note that the wanted proportion of minority class ('increaseMinorityClassTo' option) plays a key role and is probably algorithm and data dependent.
	Note also that reSMOTE() function is its in early phase and is mainly designed for very imbalanced datasets and when others options find their limits.	
}

\value{
  A matrix or a data frame of the new training sample embedding original training sample and synthetic training sample.
  Labels are put in the last column, with the column name noted "Class" (if Y is a factor) or "Response" (if not).
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}


\examples{
## generate a skewed and not easy dataset to learn
# set.seed(2014)
# n = 1000
# p = 100
# X = simulationData(n, p)
# X = fillVariablesNames(X)
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]
# Y = as.factor(ifelse(rule > (sample(5,1)*(mean(rule) - 3*epsilon1*X[,7])), "+","-"))

## distribution of labels
# distributionOfY = table(Y)

## ratio of majority/minority sample
# max(distributionOfY)/min(distributionOfY)

## generate train and test sets, experiencing concept drift

# set.seed(2014)
# train_test = init_values(X, Y, sample.size = 5/10)
# X1 = train_test$xtrain
# Y1 = train_test$ytrain
# X2 = train_test$xtest
# Y2 = train_test$ytest

# table(Y1)
# table(Y2)

## learn and predict with base model
# baseModel.ruf = randomUniformForest(X1, Y1, xtest = X2, ytest = Y2)
# baseModel.ruf

## reSMOTE train sample : we choose 'increaseMinorityClassTo = 0.85'.
## Changing its value has a strong effect on accuracy (hence calibration is needed)
# X1Y1.resmoted = reSMOTE(X1, Y1, majorityClass = "-", increaseMinorityClassTo = 0.85)

## learn and predict with reSMOTEd model
# reSMOTEdModel.ruf = randomUniformForest(X1Y1.resmoted[,1:p], X1Y1.resmoted[,"Class"], 
# xtest = as.data.frame(X2), ytest = Y2)

## compare both models, looking the AUC (and AUPR, if misclassified positives cases 
## have a much higher cost than misclassified negative ones)

# baseModel.ruf
# reSMOTEdModel.ruf

## EXTENSIONS : increasing the number of trees and/or changing the majority 
## for the voting mechanism may be a costless alternative for improving AUC.

## MORE TREES + CLASS REWEIGHTING ('classwt')
## setting 'classwt' to 'c(1,7)' means that the minority class "+" needs
## 7.14% (0.5/7) of the votes to get majority : more efficient here

# baseModelWeights.ruf = randomUniformForest(X1, Y1,
# xtest = X2, ytest = Y2, ntree = 500, classwt = c(1,7), OOB = FALSE)

## update reSMOTE with class reweighting also works... but less efficiently here
## except if cost of false negative cases is very high

# reSMOTEdModelWeights.ruf = randomUniformForest(X1Y1.resmoted[,1:p], 
# X1Y1.resmoted[,"Class"], xtest = as.data.frame(X2), ytest = Y2, 
# ntree = 500, classwt = c(1, 1.05), OOB = FALSE)

# baseModelWeights.ruf
# reSMOTEdModelWeights.ruf

## UNDERSAMPLING OF THE MAJORITY CLASS ('oversampling')
## option 'oversampling = -0.85' removes randomly and for each tree 85 percent
## of the majority class then applies an oversampling of the minority class 
## to preserve sample size before the bootstrap sampling.

# baseModelUndersampling.ruf = randomUniformForest(X1, Y1, xtest = X2, ytest = Y2, 
# ntree = 500, oversampling = -0.85, targetclass = 1, OOB = FALSE)

# baseModelUndersampling.ruf
}

\keyword{imbalanced}
\keyword{classes}
\keyword{skew}
\keyword{oversampling}