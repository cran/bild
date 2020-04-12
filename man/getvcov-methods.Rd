\name{getvcov-methods}
\docType{methods}
\alias{getvcov-methods}
\alias{getvcov,bild-method}
\title{Methods for function \code{getvcov}}
\description{Extract the variance-covariance matrix of a fitted model object from class \code{\link[=bild-class]{bild}}. }

\usage{
\S4method{getvcov}{bild}(object)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=bild-class]{bild}}.}
    }  

\section{Methods}{
\describe{
\item{\code{signature(object="bild")}:}{Returns a numeric value corresponding to the variance-covariance matrix of the fixed effect estimates of the fitted model 
\code{\link{bild}} object.} 
}}

\keyword{methods}

