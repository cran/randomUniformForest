\name{model.stats}
\alias{model.stats}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Common statistics for a vector of predictions and a vector of responses
}
\description{
   Given a vector of predictions and a vector of responses, provide some statistics and plots like AUC, AUPR, confusion matrix, F1-score, geometric mean, residuals, mean squared and mean absolute error.
}
\usage{
model.stats(predictions, responses, regression = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{predictions}{
	    a vector (or factor) of predictions.
	}
	\item{responses}{
	    a vector (or factor) of responses of the same length than 'predictions'.
	}
	\item{ regression}{
		if FALSE, considered arguments are treated as a classification task.
	}
}


\value{
  print and plot metrics and put them in a list.
}

\author{
  Saip Ciss \email{saip.ciss@wanadoo.fr} 
}