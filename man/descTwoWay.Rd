\name{descTwoWay}
\alias{descTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Descriptive Statistics for Two-Way Layout
}
\description{
\code{descTwoWay} produces basic descritive statistics including sample size, mean, standard deviation, median, minimum value, maximum value, 25th quantile, 75th quantile, skewness, kurtosis, the number of missing value for two-way layout.
}

\usage{descTwoWay(formula, data)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives two factors.}
  \item{data}{a data frame containing the variables in the formula \code{formula}} 
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a data.frame of output.

}



\author{
Osman Dag, Sam Weerahandi, Malwane Ananda 
}



\examples{

library(twowaytests)
data(alveolar)
descTwoWay(count ~ ovalbumin*treatment, data = alveolar)
}

\keyword{functions}
