\name{model.stats}
\alias{model.stats}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Common statistics for a vector (or factor) of predictions and a vector (or factor) of responses
}
\description{
   Given a vector of predictions and a vector of responses, provide some statistics and plots like AUC, AUPR, confusion matrix, F1-score, geometric mean, residuals, mean squared and mean absolute error.
}
\usage{
model.stats(predictions, responses, regression = FALSE, OOB = FALSE, plotting = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{predictions}{
	    a vector (or factor, if classification) of predictions.
	}
	\item{responses}{
	    a vector (or factor, if classification) of responses of the same length than 'predictions'.
	}
	\item{regression}{
		if FALSE, considered arguments are treated as a classification task.
	}
	\item{OOB}{
		if TRUE, expects 'prediction' to be an object of class randomUniformForest (with option 'OOB' enabled)
		in order to assess OOB predictions.
	}
	\item{plotting}{
		if TRUE, displays graphics. Set it to FALSE in the case of a regression with large datasets.
	}
}


\value{
  print and plot metrics.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}

\examples{
## not run
## Classification : synthetic data

# set.seed(2014)
# n = 1000
# p = 100 
# X = simulationData(n, p)
# X = fillVariablesNames(X)
# epsilon1 = runif(n,-1,1)
# epsilon2 = runif(n,-1,1)
# rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]
# Y = as.factor(ifelse(rule > mean(rule), "+","-")) 

# training and test sets

# train_test = init_values(X, Y, sample.size = 1/2)
# X1 = train_test$xtrain
# Y1 = train_test$ytrain
# X2 = train_test$xtest
# Y2 = train_test$ytest

# train model
# synth.ruf = randomUniformForest(X1, as.factor(Y1))

# evaluates OOB predictions
# statsOOB.pred.synth.ruf = model.stats(synth.ruf, as.factor(Y1), OOB = TRUE)

# predict 
# pred.synth.ruf = predict(synth.ruf, X2)

# statistics : produces also two plots
# stats.pred.synth.ruf = model.stats(pred.synth.ruf, as.factor(Y2))

# or, trick, do all in two lines
# synth.ruf = randomUniformForest(X1, as.factor(Y1), xtest = X2, ytest = as.factor(Y2))
# stats.pred.synth.ruf = model.stats(synth.ruf, as.factor(Y2))

## regression : synthetic data
# Y = rule
# Y1 = Y[train_test$train_idx]
# Y2 = Y[train_test$test_idx]

# synth.ruf = randomUniformForest(X1, Y1)
# statsOOB.pred.synth.ruf = model.stats(synth.ruf, Y1, OOB = TRUE, regression = TRUE)
# pred.synth.ruf = predict(synth.ruf, X2)
# stats.pred.synth.ruf = model.stats(pred.synth.ruf, Y2, regression = TRUE)
}