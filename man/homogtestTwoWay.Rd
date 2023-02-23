\name{homogtestTwoWay}
\alias{homogtestTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance Homogeneity Tests for Two-Way Layout
}
\description{
\code{homogtestTwoWay} performs variance homogeneity tests including Levene, Bartlett, Fligner-Killeen tests for two-way layout.
}

\usage{homogtestTwoWay(formula, data, method = c("Levene", "Bartlett", "Fligner"), 
  alpha = 0.05, na.rm = TRUE, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} the corresponding two factors.}
  \item{data}{a data frame containing the variables in \code{formula}} 
  \item{method}{a character string to select one of the variance homogeneity tests. "Levene": Levene's test, "Bartlett": Bartlett's test, "Fligner": Fligner-Killeen test.} 
  \item{alpha}{the level of significance to assess variance homogenity. Default is set to alpha = 0.05.}
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list containing the following components:
\item{statistic}{the corresponding test statistic.}
\item{parameter}{the parameter(s) of the approximate corresponding distribution of the test statistic. The corresponding distribution is F distribution for Levene's test, Chi-square disribution for Bartlett's test and Fligner-Killeen test.}
\item{p.value}{the p-value of the test.}

}


\author{
Osman Dag, Malwane Ananda, Sam Weerahandi
}

\seealso{
\code{\link[car]{leveneTest}} \code{\link[stats]{bartlett.test}} \code{\link[stats]{fligner.test}} 
 
}



\examples{

library(twowaytests)
data(alveolar)

homogtestTwoWay(cell ~ ovalbumin*treatment, data = alveolar)
homogtestTwoWay(cell ~ ovalbumin*treatment, data = alveolar, method = "Bartlett")

}

\keyword{functions}
