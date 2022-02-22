\name{rMAP}
\alias{rMAP}
\title{Modified Alpha Power Distribution}
\description{Generates random deviates from MAP distribution for any base distribution using acceptance rejection criteria}
\value{Random deviates from the modified alpha power distribution, length of the result is determined by n}
\usage{rMAP(n,par,distr)}
\arguments{

  \item{n}{number of random variates to be generated}
  \item{par}{named vector or list of the distribution parameters}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  }

\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed Hussein}
\examples{x<-rMAP(n=500,par=c(rate=2,alpha=3,beta=2),distr="exp")
rMAP(n=100,par=c(scale=2,shape=2,alpha=3,beta=2),distr="weibull")}
\keyword{modified qlpha power distribution, family of distributions}
