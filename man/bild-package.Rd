\name{bild-package}
\alias{bild-package}
\docType{package}
\title{bild: a package for BInary Longitudinal Data}
\description{
\code{bild} performs logistic regression for binary longitudinal
  data, allowing for serial dependence among observations from a given
  individual and a random intercept term. Estimation is via maximization
  of the exact likelihood of a suitably defined model. Missing values and 
  unbalanced data are allowed, with some restrictions. 
}
\details{This package contains functions to perform the fit of parametric models via likelihood method for binary
longitudinal data using "S4" classes and methods as implemented in the \code{methods} package.}
\author{M. Helena \enc{Gonçalves}{Goncalves}, M. \enc{Salomé}{Salome} Cabral and Adelchi Azzalini}
\references{Azzalini, A. (1994). Logistic regression for autocorrelated data with application to repeated measures. 
\emph{Biometrika}, 81, 767-775. Amendment: (1997) vol. 84, 989.

\enc{Gonçalves}{Goncalves}, M. Helena (2002) \emph{ Likelihood methods for discrete longitudinal data}. PhD thesis, 
Faculty of Sciences, University of Lisbon.
  
\enc{Gonçalves}{Goncalves}, M. Helena and Azzalini, A. (2008). Using Markov chains for marginal modelling of binary longitudinal data 
in an exact likelihood approach. \emph{Metron}, vol LXVI, 2, 157-181.
  
M. Helena \enc{Gonçalves}{Goncalves},  M. \enc{Salomé}{Salome} Cabral and Adelchi Azzalini (2012). 
The R Package \code{bild} for the Analysis of Binary Longitudinal Data. \emph{Journal of Statistical Software}, 46(4), 1-17.}

\section{Acknowledgments}{We would like to thank the CRAN team for help in the fine tuning of the Fortran code.}

\keyword{package}

\seealso{\code{\link{bild-class}}, \code{\link{bild}}, \code{\link[methods]{Methods}}}      

