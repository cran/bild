\name{bildControl}
\alias{bildControl}
\title{Auxiliary for Controlling "bild" Fitting}
\description{Auxiliary function as user interface for \code{bild} fitting }
\usage{bildControl(maxit = 100, abstol = 1e-006, reltol = 1e-006)}
\arguments{ 
 \item{maxit}{maximum number of iterations.}      
 \item{abstol}{absolute convergence tolerance.}     
 \item{reltol}{relative convergence tolerance.}
 }      
 \details{          
   See R documentation of \code{\link[stats]{optim}} for details of standard default values 
   for the remaining options not cosidered in \code{bildControl}.              
  }    
   
\value{A list with the arguments as components.}    

\seealso{\code{\link{bild-class}}, \code{\link[stats]{optim}}}

\keyword{function}
