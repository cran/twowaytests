\name{gplotTwoWay}
\alias{gplotTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-and-Whisker Plots and Error Bars for Two-Way Layout
}
\description{
\code{gplotTwoWay} produce box-and-whisker plots and error bars for two-way layout.
}

\usage{gplotTwoWay(formula, data, type = c("errorbar", "boxplot"), xlab = NULL, 
  ylab = NULL, title = NULL, legend.title = NULL, width = NULL, 
  option = c("se", "sd"), na.rm = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the corresponding two factors.}
  \item{data}{a data frame containing the variables in the formula \code{formula}} 
  \item{type}{a character string to select one of the plots. "errorbar": error bar, "boxplot": box-and-whisker plot.} 
  \item{xlab}{a label for the x axis, defaults to a description of x.} 
  \item{ylab}{a label for the y axis, defaults to a description of y.} 
  \item{title}{a main title for the plot.} 
  \item{legend.title}{a legend title for the plot.}
  \item{width}{a numeric giving the width of the little lines at the tops and bottoms of the error bars (defaults to 0.20) and the width of the boxes for box-and-whisker plots (defaults to 0.75).} 
  \item{option}{a character string to select one of the options to draw error bars with standard error or standard deviation. "se": standard error, "sd": standard deviation. Defaults to "se".} 
  \item{na.rm}{a logical indicating whether NA values should be stripped before the computation proceeds.}

}

\details{

The upper whisker of box-and-whisker plots extends from the hinge to the highest value that is within 1.5 * IQR of the hinge, where IQR is the inter-quartile range. The lower whisker extends from the hinge to the lowest value within 1.5 * IQR of the hinge. Data out of the ends of the whiskers are outliers and plotted as points.

}


\author{
Osman Dag, Malwane Ananda, Sam Weerahandi
}

\examples{

library(twowaytests)
data(alveolar)

# to obtain error bar (mean +- standard deviation) 
gplotTwoWay(count ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "sd")

# to obtain error bar (mean +- standard error) 
gplotTwoWay(count ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "se")

# to obtain box-and-whisker plot 
gplotTwoWay(count ~ ovalbumin*treatment, data = alveolar, type = "boxplot")
}

\keyword{functions}
