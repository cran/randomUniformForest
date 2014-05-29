\name{init_values} 
\alias{init_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Training and validation samples from data
}

\description{
 Draw training and test samples from data. Samples can be accessed by subsctioting original data or by their own references.
}

\usage{
init_values(X, Y = NULL, sample.size = 0.5, 
data.splitting = "ALL", 
unit.scaling = FALSE, 
scaling = FALSE, 
regression = FALSE)
}

\arguments{
   \item{X}{
    a matrix or dataframe to be splitted in training and validation sample
}
 \item{Y}{
    a response vector for the observed data.
}
 \item{sample.size}{
   size of the needed training sample in proportion of the nulber of observations in original data.
}
 \item{data.splitting}{
    not currently used.
}
\item{unit.scaling}{
    if TRUE, scale all data in X between 0 and 1, if they are all positive, or between -1 and 1.
}
\item{scaling}{
	if TRUE, centers and scales data, so each variable willhave mean 0 abd variance 1.
}
 \item{regression}{
     if TRUE and scaling  = TRUE, Y will also be scaled.
}
}

\value{
a list with the following components :
   \item{xtrain}{a matrix or data frame representing the training sample.}
  \item{ytrain}{a response vector representing the training responses according to the training sample.}
  \item{xtest }{a matrix or data frame representing the validation sample.}
  \item{ytest }{a response vector representing the validation responses according to the validation sample.}
  \item{train_idx}{subscripts of the training sample.}
  \item{test_idx}{subscripts of the validation sample.}
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr}
}

\examples{
data(iris)
Y <- iris$Species
X <- iris[,-which(colnames(iris) == "Species")]
trainingAndValidationsamples <- init_values(X, Y, sample.size = 0.5)

Xtrain = trainingAndValidationsamples$xtrain
Ytrain = trainingAndValidationsamples$ytrain
Xvalid = trainingAndValidationsamples$xtest
Yvalid = trainingAndValidationsamples$ytest
}