\name{aovTwoWay}
\alias{aovTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Way ANOVA 
}
\description{
\code{aovTwoWay} computes a two-way ANOVA for main effects and interaction effect.
}

\usage{aovTwoWay(formula, data, alpha = 0.05, na.rm = TRUE, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the two factors.}
  \item{data}{a data frame containing the variables in \code{formula}.} 
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
Toutenburg, H., Shalabh, S., Shalabh, H. (2002). Statistical analysis of designed experiments. Springer.
}

\author{
Muhammed Ali Yilmaz, Osman Dag
}

\examples{

library(twowaytests)
data(alveolar)

aovTwoWay(cell ~ ovalbumin*treatment, data = alveolar)

out <- aovTwoWay(cell ~ ovalbumin*treatment, data = alveolar)
paircompTwoWay(out)

out <- aovTwoWay(cell ~ ovalbumin*treatment, data = alveolar, alpha = 0.10)
paircompTwoWay(out)

out <- aovTwoWay(cell ~ treatment*ovalbumin, data = alveolar, alpha = 0.10)
paircompTwoWay(out)
}

\keyword{functions}
