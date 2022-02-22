\name{pMAP}
\alias{pMAP}
\title{Modified Alpha Power Distribution}
\description{Cumulative distribution function of MAP distribution for any base distribution}
\value{An object of the same length as q, giving the cdf values computed at q}
\usage{pMAP(q,par,distr,lower.tail = TRUE, log.p = FALSE)}
\arguments{

  \item{q}{vector of quantiles}
  \item{par}{named vector or list of the distribution parameters}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{lower.tail}{logical; if TRUE, probability is P(X <= x), otherwise, P(X > x).}
  \item{log.p}{logical; if TRUE, probabilitiy p is given as log(p).}
  }
\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed Hussein}
\examples{pMAP(q=2,par=list(alpha=1,beta=1.5,rate=3),distr="exp")
pMAP(c(2,3,4),par=c(scale=2,shape=2,alpha=3,beta=2),distr="weibull")}

\keyword{modified qlpha power distribution, family of distributions}
