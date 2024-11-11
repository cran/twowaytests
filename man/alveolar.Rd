\name{alveolar}
\alias{alveolar}
\docType{data}
\title{
Alveolar Cell Count Data
}
\description{
Alveolar cell count data were collected by Lacroix et al. (2002). 
}
\usage{data(alveolar)}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{ovalbumin}}{a factor with levels: no (ovalbumin-nonsensitized); yes (ovalbumin-sensitized)}
    \item{\code{treatment}}{a factor with levels: acetaldehyde; air; benzaldehyde}
    \item{\code{cell}}{Alveolar count x 10^6}
  }
}
%%  \details{
%%  ~~ If necessary, more details than the __description__ above %%  ~~
%%  }


\references{

Lacroix, G., Tissot, S., Rogerieux, F., Beaulieu, R., Cornu, L., Gillet, C., ..., Bois, F. Y. (2002). Decrease in ovalbumin-induced pulmonary allergic response by benzaldehyde but not acetaldehyde exposure in a guinea pig model. \emph{Journal of Toxicology and Environmental Health, Part A}, \bold{65:14}, 995-1012.

}


\examples{

library(twowaytests)
data(alveolar)

}
\keyword{datasets}
