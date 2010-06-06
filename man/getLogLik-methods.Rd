\name{getLogLik-methods}
\docType{methods}
\alias{getLogLik-methods}
\alias{getLogLik,bild-method}
\title{Extract Log-Likelihood}
\description{Methods for function \code{getLogLik} extracting the Log-Likelihood 
for one fitted model object from class \code{\link[=bild-class]{bild}}.}
\usage{
\S4method{getLogLik}{bild}(object)
}
\arguments{
\item{object}{an object of class \code{\link[=bild-class]{bild}}.}}
\value{Returns a numeric value corresponding to the log-Likelihood of the fitted model.}
\section{Methods}{
\describe{
    \item{\code{signature(object="bild")}:}{Returns a numeric value corresponding to the log-Likelihood of the fitted model.}
}}
\keyword{methods}
