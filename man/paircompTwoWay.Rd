\name{paircompTwoWay}
\alias{paircompTwoWay}
\alias{paircompTwoWay.default}
\alias{paircompTwoWay.twt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Comparisons for Two Way Layout
}
\description{
\code{paircompTwoWay} is a generic function for pairwise comparisons by adjusting p-values.
}

\usage{
\method{paircompTwoWay}{twt}(x, adjust.method = c("bonferroni", "holm", "hochberg", "hommel", "BH", 
  "BY", "fdr", "none"), verbose = TRUE, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{twt} object.}
  \item{adjust.method}{Method for adjusting p values (see \code{\link[stats]{p.adjust}}). Default is set to "bonferroni".} 
  \item{verbose}{a logical for printing output to R console.}
  \item{...}{Additional arguments affecting multiple comparisons of groups in two-way independent designs.}

}

\details{
The paircompTwoWay function makes pairwise comparisons depending on significance of interaction term. It calculates p-values with corresponding same test if the interaction term is not statistically significant. Otherwise, it calculates p-values with \code{\link[onewaytests]{st.test}} after two-way ANOVA, \code{\link[onewaytests]{wt.test}} after two-way ANOVA under heteroscedasticity and \code{\link[onewaytests]{mw.test}} after other tests. Then, it adjusts p-values with \code{\link[stats]{p.adjust}}.
}


\value{
Returns a data.frame of output.
}



\author{
Muhammed Ali Yilmaz, Osman Dag}



\examples{

library(twowaytests)
data(alveolar)

out <- aovTwoWay(cell ~ ovalbumin*treatment, data = alveolar)
paircompTwoWay(out)

out <- aovTwoWay(cell ~ treatment*ovalbumin, data = alveolar)
paircompTwoWay(out)

out <- aovTwoWay(cell ~ ovalbumin*treatment, data = alveolar, alpha = 0.10)
paircompTwoWay(out)

out <- aovTwoWay(cell ~ treatment*ovalbumin, data = alveolar, alpha = 0.10)
paircompTwoWay(out)
paircompTwoWay(out, adjust.method = "hochberg")

}

\keyword{functions}
