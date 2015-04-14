\name{summary.bild-class}
\Rdversion{1.1}
\docType{class}
\alias{summary.bild-class}
\title{Class "summary.bild", Summary of "bild" Objects}
\description{Extract of \code{\link{bild}} object.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("summary.bild", ...)}, 
but most often by invoking \code{summary} on an \code{bild} object. They contain values meant for printing by \code{show}.}
\section{Slots}{
\describe{
    \item{\code{coefficients}:}{Object of class \code{"matrix"}. Estimated parameters.}
    \item{\code{se}:}{Object of class \code{"matrix"}. Standard errors of estimated parameters.}
    \item{\code{covariance}:}{Object of class \code{"matrix"}. Covariance of estimated parameters.}
    \item{\code{correlation}:}{Object of class \code{"matrix"}. Correlation of estimated parameters.}
    \item{\code{log.likelihood}:}{Object of class \code{"numeric"}. The value of the log likelihood.}
    \item{\code{message}:}{Object of class \code{"integer"}. A character string giving any additional information returned 
    by the optimizer, or NULL. See \code{\link[stats]{optim}} for details.}
    \item{\code{n.cases}:}{Object of class \code{"numeric"}. Number of individual profiles used in the optimization procedure.}
    \item{\code{ni.cases}:}{Object of class \code{"numeric"}. Number of individual profiles in the dataset.}
    \item{\code{aic}:}{Object of class \code{"numeric"}. The Akaike information criterion for a fitted model object.}
    \item{\code{call}:}{Object of class \code{"language"}. The \code{call} that generated \code{bild} object.} 
}
}
\section{Extends}{Class \code{"\link[=bild-class]{bild}"}, directly.}
\section{Methods}{
\describe{
 \item{show}{\code{signature(object = "summary.bild")}: Pretty-prints object.}
 }
}

\seealso{ \code{\link{bild}}, \code{\link{bild-class}} }
\keyword{class}
