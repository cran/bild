\name{muscatine}
\alias{muscatine}
\docType{data}
\title{Muscatine}
\description{This example is a subset of data from the Muscatine Coronary Risk Factor Study, a longitudinal 
study of coronary risk factors in school children from Muscatine (Iowa, USA). }
\usage{data(muscatine)}
\format{
  A data frame with 156 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{identifies de number of the individual profile. This vector contains observation of 52 individuals.}
    \item{\code{obese}}{a numeric vector that identify the obesity status (1="yes", 0="no") of a child at each occasion. }
    \item{\code{sex}}{a factor with levels \code{1} for "female" and \code{0} for "male". }
    \item{\code{time}}{a numeric vector (1,2,3) indicating the observed time points.}
    \item{\code{counts}}{a numeric vector indicating the number of times that each profile is replicated.}
  }
}
\details{The data set presented by Fitzmaurice, Laird and Lipsitz (1994) contains records 
on 1014 children who were 7-9 years old in 1977 and were examined in 1977, 1979 and 1981. 
Height and weight were measured in each survey year, and those with relative weight greater 
than 110% of the median weight in their respective stratum were classified as obese. 
The binary response of interest is whether the child is obese (1) or not (0). However, many 
data records are incomplete, since not all children participate in all the surveys. This data set 
was also analyzed by Azzalini (1994).}
\source{ Fitzmaurice, G. M.,  Laird, N. M. and Lipsitz, S. R. (1994). Analyzing incomplete 
longitudinal binary responses: a likelihood based approach. \emph{Biometrics}, 38, 602-612.}
\references{Azzalini, A. (1994). Logistic regression for autocorrelated data with 
application to repeated measures. \emph{Biometrika}, 81, 767-775.}
\examples{
str(muscatine)

# we decompose the time effect in orthogonal components
muscatine$time1 <- c(-1, 0, 1)
muscatine$time2 <- c(1, -2, 1)

# second order Markov Chain without random effects
musc2 <- bild(obese~(time1+time2)*sex, data=muscatine, time="time1", 
        aggregate=sex, trace=TRUE, dependence="MC2")

summary(musc2)
getAIC(musc2)
getLogLik(musc2)
}
\keyword{datasets}

