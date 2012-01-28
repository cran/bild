\name{bild-class}
\Rdversion{1.1}
\docType{class}
\alias{bild-class}
\title{Class "bild" for results of a maximum likelihood estimation}
\description{This class encapsulates results of a maximum likelihood procedure.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("bild", ...)}, 
but most often as the result of a call to \code{\link{bild}}.}
\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{Object of class \code{"matrix"}. Estimated parameters.}
    \item{\code{se}:}{Object of class \code{"matrix"}. Standard errors of estimated parameters.}
    \item{\code{covariance}:}{Object of class \code{"matrix"}. Covariance of estimated parameters.}
    \item{\code{correlation}:}{Object of class \code{"matrix"}. Correlation of estimated parameters.}
    \item{\code{log.likelihood}:}{Object of class \code{"numeric"}. The value of the log likelihood.}
    \item{\code{message}:}{Object of class \code{"integer"}. A character string giving any additional information returned by the optimizer, or NULL. 
    See \code{\link[stats]{optim}} for details.}
    \item{\code{n.cases}:}{Object of class \code{"numeric"}. Number of individual profiles used in the optimization procedure.}
    \item{\code{ni.cases}:}{Object of class \code{"numeric"}. Number of individual profiles in the dataset.}
    \item{\code{aic}:}{Object of class \code{"numeric"}. The Akaike information criterion for a fitted model object.}
    \item{\code{residuals}:}{Object of class \code{"numeric"}. The residuals of estimated parameters.}
    \item{\code{s.residuals}:}{Object of class \code{"numeric"}. The residuals of estimated parameters summed over the individual profile.}
    \item{\code{ind.probability}:}{Object of class \code{"numeric"}. The transitions probabilities.}
    \item{\code{prob.matrix}:}{Object of class \code{"matrix"}. The matrix of transitions probabilities.}
    \item{\code{Fitted}:}{Object of class \code{"numeric"}. The fitted values for the estimated parameters.}
    \item{\code{Fitted.av}:}{Object of class \code{"numeric"}. }
    \item{\code{Time}:}{Object of class \code{"numeric"}. Vector of time points.}
    \item{\code{model.matrix}:}{Object of class \code{"matrix"}. The model matrix.}
    \item{\code{y.matrix}:}{Object of class \code{"matrix"}. The matrix of response values.}
    \item{\code{subset.data}:}{Object of class \code{"data.frame"}. The data subset if considered.}
    \item{\code{y.av}:}{Object of class \code{"numeric"}. The average of the response value over an individual profile.}
    \item{\code{f.value}:}{Object of class \code{"factor"}. Indicates the \code{aggregation} factor if present.}
    \item{\code{call}:}{Object of class \code{"language"}. The call to \code{"bild"}.} 
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x="bild", y="missing")}: Plots six type of plots.}
    \item{show}{\code{signature(object="bild")}: Display object briefly.}
    \item{summary}{\code{signature(object="bild")}: Generate object summary.}
    \item{getAIC}{\code{signature(object="bild")}: Returns a numeric value corresponding to the AIC of the fitted model.}
    \item{getLogLik}{\code{signature(object="bild")}: Returns a numeric value corresponding to the log-Likelihood of the fitted model.}
	 }
}
\keyword{class}

