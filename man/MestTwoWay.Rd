\name{MestTwoWay}
\alias{MestTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Way ANOVA for M-estimators
}
\description{
\code{MestTwoWay} computes a two-way ANOVA for M-estimators.
}

\usage{MestTwoWay(formula, data, estimator = c("mom_est", "onestep_est", "median"), 
  nboot = 500, distance = c("mahalanobis", "projected"), seed = 123, alpha = 0.05, 
  na.rm = TRUE, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the two factors.}
  \item{data}{a data frame containing the variables in \code{formula}.} 
  \item{estimator}{a character string to select one of the estimators. "mom_est": the modified one-step estimator of location based on Huber's Psi, "onestep_est": one-step estimator of location using Huber's Psi, "median": median.}
  \item{nboot}{the number of bootstrap samples.}
  \item{distance}{the specified distance. Default is set to "mahalanobis".}
  \item{seed}{a seed number for the reproducibility of results. Default is set to 123.}
  \item{alpha}{the level of significance to assess the statistical difference. Default is set to alpha = 0.05.} 
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}
}

%\note{
%The tests available from this R library are based on two Generalized P-value approch, for Two-Way ANOVA under unequal variances and cell frequencies. The first test, the gPQ, is an extension Li et al. (2011), and the second test, gPB, is a generalized test that is numerically equivalent to the Parametric Bootstrap test derived by Xu et al. (2013). The gPQ test tends to assure the intended size of the test, but somewhat conservative, especially when the sample sizes are small.  
%
%The gPB test tends to exceed the intended size of the test. Hence, the gPB is recommended for situations of small sample sizes, and gPQ otherwise. 
%}

\value{

A list containing the following components:
\item{output}{a data frame of output.}
\item{alpha}{the level of significance to assess the statistical difference.} 
\item{method}{the method name.}
\item{data}{a data frame containing the variables in which NA values (if exist) are removed.}
\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the two factors.}

}

\references{

Wilcox, R.R. (2011). Introduction to robust estimation and hypothesis testing. Academic press.

}

\author{
Muhammed Ali Yilmaz, Osman Dag
}

\examples{

library(twowaytests)
data(alveolar)

MestTwoWay(cell ~ ovalbumin*treatment, data = alveolar)

out <- MestTwoWay(cell ~ ovalbumin*treatment, data = alveolar, verbose = FALSE)
paircompTwoWay(out)

out <- MestTwoWay(cell ~ treatment*ovalbumin, data = alveolar, verbose = FALSE)
paircompTwoWay(out)

}

\seealso{
\code{\link[WRS2]{t2way}} 
 
}

\keyword{functions}
