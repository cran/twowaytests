\name{nortestTwoWay}
\alias{nortestTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normality Tests for Two-Way Layout
}
\description{
\code{nortestTwoWay} performs normality tests including Shapiro-Wilk, Shapiro-Francia, Kolmogorov-Smirnov, Anderson-Darling, Cramer-von Mises, Pearson Chi-square tests, and also assess the normality of each group through plots for two-way layout.
}

\usage{nortestTwoWay(formula, data, method = c("SW", "SF", "LT", "AD", "CVM", "PT"), 
  alpha = 0.05, plot = c("qqplot", "histogram"), na.rm = TRUE, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} the corresponding two factors.}
  \item{data}{a data frame containing the variables in the formula \code{formula}} 
  \item{method}{a character string to select one of the normality tests. "SW": Shapiro-Wilk test, "SF": Shapiro-Francia test, "LT": Lilliefors (Kolmogorov-Smirnov) test, "AD": Anderson-Darling test, "CVM": Cramer-von Mises test, "PT": Pearson Chi-square test.} 
  \item{alpha}{the level of significance to assess normality. Default is set to alpha = 0.05.}
  \item{plot}{a character string to select one of the plots including qqplot and histogram. The red line is the density line of normal distribution. If plot is set to NULL, the plot is not returned.}
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A data frame gives the test results for the normality of each subset via corresponding normality.

}


\author{
Osman Dag, Sam Weerahandi, Malwane Ananda
}

\seealso{
\code{\link[twowaytests]{gplotTwoWay}} \code{\link[stats]{shapiro.test}} 
 
}



\examples{
library(twowaytests)
data(alveolar)
nortestTwoWay(count ~ ovalbumin*treatment, data = alveolar)
}

\keyword{functions}
