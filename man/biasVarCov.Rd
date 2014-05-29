\name{biasVarCov} 
\alias{biasVarCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bias-Variance-Covariance Decomposition
}

\description{
 Bias-Variance-Covariance decomposition for Mean Squared Error (MSE) or test error in binary classification, between a response vector and its estimate, over the test sample. For every estimate, based on training examples, MSE on the test sample, between values of the response and values of the estimate, can be decomposed in noise (variance of the response), squared bias (between response and estimate), variance (of the estimate) and covariance (of the response and the estimate). Same decomposition arrives for binary classification with responses in \{0, 1\}.
}

\usage{
biasVarCov(predictions, target, regression = FALSE, idx = 1:length(target))
}

\arguments{
   \item{predictions}{
    vector of predictions for the test sample.
}
 \item{target}{
    response vector for the test sample.
}
 \item{regression}{
   if TRUE, decomposition is done for regression. If FALSE, for classification.
}
 \item{idx}{
    not currently used.
}
}

\value{
a list with the following components :
   \item{MSE, predError}{mean squared error or test error for the test sample.}
  \item{squaredBias}{squared bias.}
  \item{predictionsVar}{variance of the estimate.}
  \item{predictionsTargetCov}{covariance between estimate and response.}
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
n = 100;  p = 10
# simulate 'p' gaussian vectors with random parameters between -10 and 10.
X <- simulationData(n,p)

# make a rule to create response vector
epsilon1 = runif(n,-1,1)
epsilon2 = runif(n,-1,1)
rule = 2*(X[,1]*X[,2] + X[,3]*X[,4]) + epsilon1*X[,5] + epsilon2*X[,6]

# Classification
Y <- as.factor(ifelse(rule > mean(rule), 1, 0))

# define random train and test sample of (X,Y)
trainTestIdx <- cut(sample(1:nrow(X), nrow(X)), 2, labels= FALSE)

# training sample
Xtrain = X[trainTestIdx == 1,]
Ytrain = Y[trainTestIdx == 1]

# test sample
Xtest = X[trainTestIdx == 2,]
Ytest = Y[trainTestIdx == 2]

# learning model
synthDataset.ruf <- randomUniformForest(Xtrain, Ytrain, threads = 1, 
ntree = 20, BreimanBounds = FALSE)
synthDataset.pred.ruf <- predict(synthDataset.ruf, Xtest)

# test error
testError = 1 - sum(synthDataset.pred.ruf == Ytest)/length(Ytest)
testError

# test error decomposition
testError.dec <- biasVarCov(synthDataset.pred.ruf, Ytest, regression = FALSE)

# Regression
Y = rule

# training sample
Xtrain = X[trainTestIdx == 1,]
Ytrain = Y[trainTestIdx == 1]

# test sample
Xtest = X[trainTestIdx == 2,]
Ytest = Y[trainTestIdx == 2]

# learning model
synthDataset.ruf <- randomUniformForest(Xtrain, Ytrain, threads = 1,
ntree = 20, BreimanBounds = FALSE)
synthDataset.pred.ruf <- predict(synthDataset.ruf, Xtest)

# MSE
MSE <- sum((synthDataset.pred.ruf - Ytest)^2)/length(Ytest)

# test error decomposition
MSE.dec <- biasVarCov(synthDataset.pred.ruf, Ytest, regression = TRUE)
}

\concept{MSE decomposition}
\concept{generalization error decomposition}