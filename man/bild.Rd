\name{bild}
\alias{bild}
\title{ fit of parametric models for binary longitudinal data via likelihood method}
\description{\code{bild} performs the fit of parametric models via likelihood method. 
Serial dependence and random intercept are allowed according to the stochastic model chosen. 
Missing values and unbalanced data are automatically accounted for computing the 
likelihood function.}
\usage{bild(formula = formula(data), data, time, id, subSET, aggregate = FALSE, 
       start = NULL, trace = FALSE, dependence="ind", method = "BFGS", 
       control = bildControl(), integrate = bildIntegrate())}
\arguments{ 
  \item{formula}{a description of the model to be fitted of the form response~predictors}
  \item{data}{a \code{data} frame containing the variables in the formula. NA values are allowed. 
  If data is missing, an error message is produced. See "Details".}
  \item{time}{a string that matches the name of the \code{time} variable in data. By default, the program expects a variable named \code{time} 
   to be present in the \code{data.frame}, otherwise the name of the variable playing the role of time must be declared by assigning \code{time} here. }
  \item{id}{a string that matches the name of the \code{id} variable in \code{data}. By default, the program expects a variable named \code{id} 
   to be present in the \code{data.frame}, otherwise the name of the variable playing the role of \code{id} must be declared by assigning \code{id} here. }
  \item{subSET}{an optional expression indicating the subset of the rows of \code{data} that should be 
    used in the fit. All observations are included by default.}
  \item{aggregate}{a string that permits the user identify the factor to be used in the \code{\link{plot-methods}}.}
  \item{start}{a vector of initial values for the nuisance parameters of the likelihood. The dimension of the vector is according 
  to the structure of the dependence model.}
  \item{trace}{logical flag: if TRUE, details of the nonlinear optimization are printed. By default the flag is set to FALSE.}
  \item{dependence}{expression stating which \code{dependence} structure should be used in the fit. The default is "ind".
  According to the stochastic model chosen serial dependence and random effects are allowed. 
  There are five options: "\code{ind}" (independence), "\code{MC1}" (first order Markov Chain), "\code{MC2}" (second order Markov Chain),  
  "\code{MC1R}" (first order Markov Chain with random intercept) or "\code{MC2R}" (second order Markov Chain with random intercept).}
  \item{method}{The \code{method} to be used in the optimization process: "\code{BFGS}","\code{CG}", "\code{L-BFGS-B}" and "\code{SANN}". 
  The default is "\code{BFGS}". 
  See \code{\link[stats]{optim}} for details.}
  \item{control}{a list of algorithmic constants for the optimizer \code{optim}. 
  See R documentation of \code{optim.control} for details and possible control options. By default, \code{bild} sets the maximum number 
  of iterations (\code{maxit}) equal to 100, the absolute convergence tolerance (\code{abstol}) and the relative 
  convergence tolerance (\code{rel.tol}) equal to 1e-6 and uses the \code{\link[stats]{optim}} standard default values for the remaining options.}
  \item{integrate}{a list of algorithmic constants for the computation of a definite integral using a Fortran-77 subroutine. See "Details".}      
}      
   \section{Background}{ 
   Assume that each subject of a given set has been observed at number of 
   successive time points. For each subject and for each time point, a binary 
   response variable, taking value 0 and 1, and a set of covariates are 
   recorded.  
   The underlying methodology builds a logistic regression model for the
   probability that the response variable  takes value 1
   as a function of the  covariates, taking into account that successive 
   observations from the same individual cannot be assumed to be independent. 

   The basic model for serial dependence is of Markovian type of the first order
   (denoted \code{MC1} here), suitably constructed so that the logistic regression
   parameters maintain the same meaning as in ordinary logistic regression for
   independent observations. The serial dependence parameter is the logarithm of 
   the odds-ratio between probabilities of adjacent observations, which is 
   assumed to be constant for all adjacent pairs, and it is denoted here 
   \code{log.psi1}.

   An extension of this formulation allows a Markovian dependence of the second 
   order, denoted \code{MC2} here.  In this case there are two parameters which 
   regulate serial dependence: \code{log.psi1} as before and \code{log.psi2} 
   which is the analogous quantity for observations which are two time units apart, 
   conditionally on the intermediate value.

   Individual random effects can be incorporated in the form of a random
   intercept term of the linear predictor of the logistic regression,
   assuming a normal distribution of mean 0 and variance \eqn{\sigma^2}, 
   parameterized as \eqn{\omega=\log(\sigma^2)}{\omega=log(\sigma^2)}. 
   The combination of serial Markov dependence with a random intercept corresponds here  
   to the dependence structures \code{MC1R} and \code{MC2R}."
    
   Original sources of the above formulation are given by Azzalini (1994), as
   for the first order Markov dependence, and by \enc{Gonçalves}{Goncalves} (2002) and 
   \enc{Gonçalves}{Goncalves} and Azzalini (2008) for the its extensions.
   }
\details{
\code{data} are contained in a \code{data.frame}. Each element of the \code{data} argument must be identifiable by a name. 
The simplest situation occurs when all subjects are observed at the same time points. 
The response variable represent the individual profiles of each subject, it is expected
a variable in the \code{data.frame} that identifies the correspondence of each component of the response variable to the subject that it belongs, 
by default is named \code{id} variable. It is expected a variable named \code{time} to be present in the \code{data.frame}. 
If the \code{time} component has been given a different name, this should be declared. 
The \code{time} variable should identify the time points that each individual profile has been observed.
 
 When it is expected that all subjects in one experiment to be observed at the same time points, but in practice some of the subjects were
 not observed in some of the scheduled occasions, NA values can then be inserted in the response variable. 
 If a response profile is replicated several times, a variable called \code{counts} must be created accordingly. 
 This vector is used for weighting the response profile indicating for each individual profile the number of times that is replicated. 
 The vector \code{counts} must repeat the number of the observed replications of each individual profile as many times as the number of observed time 
 points for the correspondent profile. The program expect such vector to be named \code{counts}. 
 If each profile has been observed only once, the construction of the vector \code{counts} is not required. 

\code{subSET} is an optional expression indicating the subset of \code{data} that should be 
    used in the fit. This is a logical statement of the type 
    \code{variable 1} == "a" & \code{variable 2} > x  
    which identifies the observations to be selected. All observations are included by default.    



For the models with random intercept \code{MC1R} and \code{MC2R}, 
\code{bild} compute integrals based on a Fortran-77 subroutine package 
\code{QUADPACK}. For some data sets, when the dependence structure has 
a random intercept term, the user could have the need to do a specification 
of the \code{integrate} argument list changing  
the integration limits in the \code{bildIntegrate} function. 
The \code{\link{bildIntegrate}} is an auxiliary function for controlling \code{bild} 
fitting. See the example of \code{\link{locust}} data.   


}

\value{An object of class \code{\link[=bild-class]{bild}}.}
\references{Azzalini, A. (1994). Logistic regression for autocorrelated data with application to repeated measures.  
\emph{Biometrika}, 81, 767-775. Amendment: (1997) vol. 84, 989.

\enc{Gonçalves}{Goncalves}, M. Helena (2002). \emph{ Likelihood methods for discrete longitudinal data}. PhD thesis, 
Faculty of Sciences, University of Lisbon.
  
\enc{Gonçalves}{Goncalves}, M. Helena and Azzalini, A. (2008). Using Markov chains for marginal modelling of binary longitudinal data 
in an exact likelihood approach. \emph{Metron}, vol LXVI, 2, 157-181.
  
\enc{Gonçalves}{Goncalves} MH, Cabral MS and Azzalini A (2012). 
The R Package \code{bild} for the Analysis of Binary Longitudinal Data. \emph{Journal of Statistical Software}, 46(9), 1-17.}
  
\author{M. Helena \enc{Gonçalves}{Goncalves}, M. \enc{Salomé}{Salome} Cabral and Adelchi Azzalini}

\seealso{\code{\link{bild-class}}, \code{\link{bildControl}}, \code{\link{bildIntegrate}}, \code{\link[stats]{optim}}}

\examples{\donttest{ 
#####  data= airpollution, dependence="MC2R"
str(airpollution)

air2r <- bild(wheeze~age+smoking, data=airpollution, trace=TRUE, time="age",
    aggregate=smoking, dependence="MC2R")

summary(air2r)
getAIC(air2r)
getLogLik(air2r)
plot(air2r) 

####  data=muscatine, dependence="MC2R"
str(muscatine)
                                            
# we decompose the time effect in orthogonal components
muscatine$time1 <- c(-1, 0, 1)
muscatine$time2 <- c(1, -2, 1)

musc2 <- bild(obese~(time1+time2)*sex, data=muscatine, time="time1", 
        aggregate=sex, trace=TRUE, dependence="MC2")
        
summary(musc2)
getAIC(musc2)
getLogLik(musc2)

}}

\keyword{function}
