\name{gplotTwoWay}
\alias{gplotTwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-and-Whisker, Violin Plots and Error Bars for Two-Way Layout
}
\description{
\code{gplotTwoWay} produces box-and-whisker, violin plots, and error bars for two-way layout.
}

\usage{gplotTwoWay(formula, data, type = c("errorbar", "violin", "boxplot"), 
  color_manual = NULL, back_color = FALSE, xlab = NULL, ylab = NULL, 
  title = NULL, legend.title = NULL, width = NULL, option = c("sd", "se"), 
  na.rm = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} gives the corresponding two factors.}
  \item{data}{a data frame containing the variables in \code{formula}.} 
  \item{type}{a character string to select one of the plots. "errorbar": error bar, "violin": violin plot, and "boxplot": box-and-whisker plot.} 
  \item{color_manual}{a vector of colors. A palette can also be defined with \code{\link[wesanderson]{wes_palette}}. Default is set to "FantasticFox1" available in \code{\link[wesanderson]{wes_palette}}.} 
  \item{back_color}{a logical to make background grey or white.}   
  \item{xlab}{a label for the x axis, defaults to a description of x.} 
  \item{ylab}{a label for the y axis, defaults to a description of y.} 
  \item{title}{a main title for the plot.} 
  \item{legend.title}{a legend title for the plot.}
  \item{width}{a numeric giving the width of the little lines at the tops and bottoms of the error bars (defaults to 0.35) and the width of the boxes for box-and-whisker plots (defaults to 0.75).} 
  \item{option}{a character string to select one of the options to draw error bars with standard error or standard deviation. "se": standard error, "sd": standard deviation. Defaults to "se".} 
  \item{na.rm}{a logical indicating whether NA values should be stripped before the computation proceeds.}

}

\details{

The upper whisker of box-and-whisker plots extends from the hinge to the highest value that is within 1.5 * IQR of the hinge, where IQR is the inter-quartile range. The lower whisker extends from the hinge to the lowest value within 1.5 * IQR of the hinge. Data out of the ends of the whiskers are outliers and plotted as points.

}


\author{
Osman Dag, Muhammed Ali Yilmaz
}

\examples{

library(twowaytests)
data(alveolar)

# to obtain error bar (mean +- standard deviation) 
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "sd")

# to change factor places in graph
gplotTwoWay(cell ~ treatment*ovalbumin, data = alveolar, type = "errorbar", option = "sd")

# to obtain error bar (mean +- standard error) 
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "se")

# to obtain box-and-whisker plot 
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "boxplot")

# to obtain violin plot 
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "violin")

# to define specific colors
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "sd", 
color_manual = c("#999999", "#E69F00", "#56B4E9"))

# to change color palette
library(wesanderson)
gplotTwoWay(cell ~ ovalbumin*treatment, data = alveolar, type = "errorbar", option = "sd", 
color_manual = wes_palette(name = "GrandBudapest1", n = 3))



}

\keyword{functions}
