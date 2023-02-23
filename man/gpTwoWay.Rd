\name{gpTwoWay}
\alias{gpTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Way ANOVA under Heteroscedasticity
}
\description{
\code{gpTwoWay} computes a two-way ANOVA for main effects and interaction effect under heteroscedasticity.
}

\usage{gpTwoWay(formula, data, method = c("gPB","gPQ"), seed = 123, alpha = 0.05, 
  na.rm = TRUE, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the two factors.}
  \item{data}{a data frame containing the variables in \code{formula}.} 
  \item{method}{a character string to select the method. "gPB": Parametric Bootstrap based Generalized Test, "gPQ": Generalized Pivotal Quantity based Generalized Test.} 
  \item{seed}{a seed number for the reproducibility of results. Default is set to 123.}
  \item{alpha}{the level of significance to assess the statistical difference. Default is
set to alpha = 0.05.} 
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}


}

\note{
These tests available from this R library are based on two Generalized P-value approaches, for Two-Way ANOVA under unequal variances and cell frequencies. The first test, the gPQ, is an extension Li et al. (2011), and the second test, gPB, is a generalized test that is numerically equivalent to the Parametric Bootstrap test derived by Xu et al. (2013). The gPQ test tends to assure the intended size of the test, but somewhat conservative, especially when the sample sizes are small.  

The gPB test tends to exceed the intended size of the test. Hence, the gPB is recommended for situations of small sample sizes, and gPQ otherwise. 
}


\value{

A list containing the following components:
\item{output}{a data frame of output.}
\item{alpha}{the level of significance to assess the statistical difference.} 
\item{method}{the selected method used in generalized test.}
\item{data}{a data frame containing the variables in which NA values (if exist) are removed.}
\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the two factors.}

}

\references{

Ananda, M.M., Dag, O., Weerahandi, S. (2022). Heteroscedastic two-way ANOVA under constraints. \emph{Communications in Statistics-Theory and Methods}, 1-16.

}


\author{
Sam Weerahandi, Osman Dag, Malwane Ananda
}


\examples{


###Example 1

library(twowaytests)
data(alveolar)

\donttest{

# to use Parametric Bootstrap based Generalized Test
gpTwoWay(cell ~ ovalbumin*treatment, data = alveolar, method = "gPB")

# to use Generalized Pivotal Quantity based Generalized Test
gpTwoWay(cell ~ ovalbumin*treatment, data = alveolar, method = "gPQ")

out <- gpTwoWay(cell ~ ovalbumin*treatment, data = alveolar, method = "gPB")
paircompTwoWay(out)

out <- gpTwoWay(cell ~ treatment*ovalbumin, data = alveolar, method = "gPB")
paircompTwoWay(out)
}
}
\keyword{functions}
