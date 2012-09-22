\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,bild-method}
\title{Methods for Function summary in package "bild"}
\description{ Summarize objects}
\usage{
\S4method{summary}{bild}(object, cov=FALSE, cor=FALSE)}
\arguments{ 
  \item{object}{an object of class \code{\link[=bild-class]{bild}}.}
  \item{cov}{if set to TRUE prints the matrix of covariances between parameters estimates. The default is FALSE.}
  \item{cor}{if set to TRUE prints the matrix of correlations between parameters estimates. The default is FALSE.}  
   }  
\details{ Computes and returns a list of summary statistics of the fitted linear model given a \code{bild} object, 
using the components (list elements) "call" and "terms" from its argument, plus 
depending on the structure of the dependence model chosen, the table for the estimates of coefficients will appear 
\code{log.psi1} if the dependence structure of the process corresponds to a first-order 
Markov chain, or both \code{log.psi1} and \code{log.psi2} if the dependence structure of the process corresponds 
to a second-order Markov chain. \code{log.psi1} is the log-odds ratio between adjacent observations and 
\code{log.psi2} is the logarithm of conditional odds ratio for observations separated by one time point. 
If the structure of the dependence model chosen includes the random intercept (models "\code{indR}", "\code{MC1R}" and "\code{MC2R}") 
the estimate of the random effect (\code{omega}) will also appear where \eqn{\omega=\log\sigma^2}{\omega=log(\omega^2)}.}
\section{Methods}{
\describe{
\item{\code{signature(object = "ANY")}:}{Generic function.}
\item{\code{signature(object = "bild")}:}{Generate a summary as an object of class 
\code{\link[=summary.bild-class]{summary.bild}}, 
containing information about the matched call to \code{bild}, the number of profiles in the data, 
the number of profiles used in the fit, the log-likelihood, the AIC,
a table with estimates, asymptotic SE, t-values and p-values,  
the estimated correlation and variance-covariance matrix for the estimated parameters if the user wishes, 
and a message giving additional information returned by the optimizer.}
}}
\keyword{methods}

