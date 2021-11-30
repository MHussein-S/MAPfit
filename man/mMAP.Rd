\name{mMAP}
\alias{mMAP}
\title{Modified Alpha Power Distribution}
\description{The rth moment about the origin of the modefied alpha power distribution for any base distribution}
\value{An object of the same length as x, giving the rth moment}
\usage{mMAP(r,par,distr)}
\arguments{

  \item{r}{moment of order r}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  }

\references{}
\author{Mohamed Hussein}
\examples{mMAP(r=2,par=c(rate=2,alpha=3,beta=2),distr="exp")
mMAP(r=3,par=c(rate=0.2,alpha=1,beta=2),distr="exp")
mMAP(r=1,par=c(scale=2,shape=2,alpha=3,beta=2),distr="weibull")
}
\keyword{modified qlpha power distribution, family of distributions}
