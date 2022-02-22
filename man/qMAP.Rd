\name{qMAP}
\alias{qMAP}
\title{Modified Alpha Power Exponential Distribution}
\description{Quantile function of MAP distribution for any base distribution}
\value{An object of the same length as p, giving the quantile function values computed at p}
\usage{qMAP(p,par,distr,lower=0,upper,lower.tail=TRUE,log.p=FALSE)}
\arguments{
  \item{p}{scalar or vector of probabilities.}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{lower.tail}{logical; if TRUE, probability are P(X <= x), otherwise, P(X > x).}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed Hussein}
\examples{qMAP(p=0.5,c(rate=0.2,alpha=2,beta=2),"exp",upper=100)
  qMAP(0.25,c(rate=0.2,alpha=2,beta=2),"exp",upper=100)
  qMAP(0.75,c(alpha=2,beta=2),"norm",upper=100)
  qMAP(log(0.5),c(alpha=5,beta=1,rate=0.5),"exp",lower=0,upper=100,log.p=TRUE)
  qMAP(0.10,c(alpha=1.5,beta=1,rate=0.5),"exp",lower=0,upper=100,lower.tail=FALSE)
  qMAP(p=0.5,c(rate=0.2,alpha=2,beta=2),"exp",upper=100)
}
\keyword{modified alpha power distribution, family of distributions}
