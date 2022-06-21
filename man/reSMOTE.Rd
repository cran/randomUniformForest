\name{reSMOTE}
\alias{reSMOTE}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  REplication of a Synthetic Minority Oversampling TEchnique for highly imbalanced datasets
}
\description{
  Produce new samples of the minority class by creating synthetic data of the original ones by randomization. After that, the new dataset contains all original data + synthetic data labelled with the minority class. Main argument is derived from Breiman's ideas on the construction of synthetic data for the unsupervised mode of random (Uniform) Forests. The new dataset preserves the distribution of the original dataset covariates(at least with default 'samplingMethod' option).
}
\usage{
reSMOTE(X, Y, 
	majorityClass = 0, 
	increaseMinorityClassTo = NULL,
	conditionalTo = NULL,
	samplingFromGaussian = FALSE,  
	samplingMethod = c("uniform univariate sampling", 
		"uniform multivariate sampling", "with bootstrap"),
	maxClasses = max(10, which.is.factor(X[,, drop = FALSE], 
	count = TRUE)), 
	seed = 2014)
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
	  label(s) of the majority class(es).
	}
	\item{increaseMinorityClassTo}{
	  proportion of minority class observations that one wants to reach.
	}
	\item{conditionalTo}{
	  the targeted variable (usually the response) that one needs the new (artificial) observations to be sampled.
	  More precisely, synthetic dataset will be generated conditionally to each target value.
	}
	\item{samplingFromGaussian}{
	  not used for now.
	}
	\item{samplingMethod}{
	  which method to use for sampling from data? Note that one may use any of the two first methods with the 'with bootstrap' option, for example 'samplingMethod = c("uniform univariate sampling",  "with bootstrap")' or 'samplingMethod = "uniform univariate sampling"'.
	}
	\item{maxClasses}{
	  To be documented.
	}
	\item{seed}{
	  seed to use when sampling distribution.
	}
}

\details{
	The main purpose of the reSMOTE( ) function is to increase AUC. The function is specifically designed for
	random Uniform Forests with possible conjunction of options 'classwt', 'rebalancedsampling' and/or 'oversampling', also designed to treat imbalanced classes. Since ensemble models tend to choose majority class too quickly, reSMOTE introduces a perturbation that forces the algorithm to learn data more slowly. Note that the wanted proportion of minority class ('increaseMinorityClassTo' option) plays a key role and is probably model and data dependent.
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
## not run
## A - generate a skewed dataset to learn
# set.seed(2014)
# n = 1000
# p = 100
# X = simulationData(n, p)
# X = fillVariablesNames(X)
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]
# Y = as.factor(ifelse(rule > (3*(mean(rule) - 4/3*epsilon1*X[,7])), "-","+"))

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

## reSMOTE train sample : we choose 'increaseMinorityClassTo = 0.3'.
# X1Y1.resmoted = reSMOTE(X1, Y1, majorityClass = "-",
# conditionalTo = Y1, increaseMinorityClassTo = 0.3)

## learn and predict with reSMOTEd model
# reSMOTEdModel.ruf = randomUniformForest(X1Y1.resmoted[,1:p], X1Y1.resmoted[,"Class"], 
# xtest = X2, ytest = Y2)

## compare both models, looking the AUC (and AUPR, if misclassified positives cases 
## have a much higher cost than misclassified negative ones)

# baseModel.ruf
# reSMOTEdModel.ruf

## B - German credit dataset 
# URL = "http://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/"
# dataset = "german.data-numeric"
# germanCreditData = read.table(paste(URL, dataset, sep =""))

# Y = germanCreditData[,25]
# X = germanCreditData[,-25]

# set.seed(2014)
# train_test = init_values(X, Y, sample.size = 5/10)
# X1 = train_test$xtrain
# Y1 = train_test$ytrain
# X2 = train_test$xtest
# Y2 = train_test$ytest

# 1 - default modelling: AUC around 0.67; AUPR around 0.56

# germanCredit.baseModel.ruf = randomUniformForest(X1, as.factor(Y1), 
# xtest = X2, ytest = as.factor(Y2), ntree = 200, importance = FALSE)
# germanCredit.baseModel.ruf

## 2 - reSMOTEd modelling: AUC around 0.70; AUPR around 0.60; 
## Note : reSMOTE is not consistent with OOB evaluation.

# p = ncol(X1)
# X1Y1.resmoted = reSMOTE(as.true.matrix(X1), as.factor(Y1), majorityClass = "1", 
# conditionalTo = Y1, increaseMinorityClassTo = 0.6)
# X1.new = X1Y1.resmoted[,1:p]
# Y1.new = X1Y1.resmoted[,"Class"]
# germanCredit.reSMOTEd.ruf2 = randomUniformForest(X1.new, as.factor(Y1.new), 
# xtest = X2, ytest = as.factor(Y2), ntree = 200, importance = FALSE)

## 2.1 - try reSMOTE + Random Forests: AUC around 0.70; AUPR around 0.58;

# library(randomForest) 
# germanCredit.reSMOTEd.rf = randomForest(X1.new, as.factor(Y1.new), 
# xtest = X2, ytest = as.factor(Y2))
# statsRF = model.stats(as.factor(as.numeric(germanCredit.SMOTEd.rf$test$pred)), as.factor(Y2))

## 3 - use benchmark: SMOTE, implemented in 'unbalanced' R package

# require(unbalanced)
## SMOTE : generate 'perc.over/100' new instance for each rare instance.
## 'perc.under/100' is the number of instances of majority class for each smoted observation.
## Generate smoted observations to have a similar distribution (but more examples) 
## of the labels than the one generated by 'reSMOTE()'.

# XX = ubBalance(X1, as.factor(Y1-1), type = "ubSMOTE", perc.over = 400, perc.under = 100)

## 3.1 apply to randomUniformForest: AUC around 0.68; AUPR around 0.55;

# germanCredit.SMOTEd.ruf = randomUniformForest(XX$X, as.factor(XX$Y), 
# xtest = X2, ytest = as.factor(Y2-1), ntree = 200, importance = FALSE, BreimanBounds = FALSE)
 
## 3.2 apply to randomForest: AUC around 0.68; AUPR around 0.57;
 
# germanCredit.SMOTEd.rf = randomForest(XX$X, as.factor(XX$Y),
# xtest = X2, ytest = as.factor(Y2-1))
# statsModel = model.stats(as.factor(as.numeric(germanCredit.SMOTEd.rf$test$pred)), 
# as.factor(Y2))

## 4 - Stabilize : AUC around 0.71; AUPR around 0.62;
## a - generate a large number of minority class examples
## b - use rebalanced sampling in the model : equal size for each class 
## so that size of the new training sample equals size of the original training sample
## c - increase minimum number of observations in terminal nodes (nodesize)
## d - increase number of trees.

# X1Y1.resmoted = reSMOTE(as.true.matrix(X1), as.factor(Y1), majorityClass = "1", 
# conditionalTo = Y1, increaseMinorityClassTo = 1)
# X1.new = X1Y1.resmoted[,1:p]
# Y1.new = X1Y1.resmoted[,"Class"]

# germanCredit.reSMOTEd.ruf2 = randomUniformForest(X1.new, as.factor(Y1.new), 
# xtest = X2, ytest = as.factor(Y2), rebalancedsampling = c(250,250), nodesize = 10,
# ntree = 1000, importance = FALSE, OOB = FALSE)

# 4.1 - Apply to Random Forests : AUC around 0.72; AUPR around 0.62;
# germanCredit.reSMOTEd.rf = randomForest(X1.new, as.factor(Y1.new), 
# xtest = X2, ytest = as.factor(Y2), sampsize = c(250,250), nodesize = 10, ntree = 1000)
# statsModel = model.stats(as.factor(as.numeric(germanCredit.reSMOTEd.rf$test$pred)), 
# as.factor(Y2))
}

\keyword{imbalanced}
\keyword{classes}
\keyword{skew}
\keyword{oversampling}