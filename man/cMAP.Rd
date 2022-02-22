\name{cMAP}
\alias{cMAP}
\title{Modified Alpha Power Distribution}
\description{The rth moment about the mean of the modefied alpha power distribution for any base distribution}
\value{An object of the same length as x, giving the rth centeral moment}
\usage{cMAP(r,par,distr)}
\arguments{

  \item{r}{moment of order r}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  }

\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}

\author{Mohamed Hussein}
\examples{cMAP(r=2,par=c(rate=2,alpha=3,beta=2),distr="exp")
cMAP(r=3,par=c(rate=0.2,alpha=1,beta=2),distr="exp")
cMAP(r=3,par=c(scale=2,shape=2,alpha=3,beta=2),distr="weibull")
}
\keyword{modified qlpha power distribution, family of distributions}
