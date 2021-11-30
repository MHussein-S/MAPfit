\name{dMAP}
\alias{dMAP}
\title{Modified Alpha Power Distribution}
\description{Probability density function of MAP distribution for any base distribution}
\value{An object of the same length as x, giving the pdf values computed at x }
\usage{dMAP(x,par,distr, log = FALSE)}
\arguments{

  \item{x}{vector of quantiles}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{log}{logical; if TRUE, the values of log(pdf), at x, are given.}
  }
\references{}
\author{Mohamed Hussein}
\examples{dMAP(x=2,par=list(alpha=1,beta=1.5,rate=3),distr="exp")
dMAP(c(2,3,4),par=c(scale=2,shape=2,alpha=3,beta=2),distr="weibull")}

\keyword{modified qlpha power distribution, family of distributions}


