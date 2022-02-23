\name{modeMAP}
\alias{modeMAP}
\title{Modified Alpha Power Distribution}
\description{Mode of MAP distribution by direct maximization of log(g(x))}
\value{Approximation of the Mode of MAP distribution}
\usage{modeMAP(start=0,parm,distr,lower=0,upper)}
\arguments{
  \item{start}{Initial values for the parameters to be maximized over.}
  \item{parm}{named vector or list of the distribution parameters.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{upper}{the upper end point of the interval to be searched.}
  }
\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed
Hussein}
\examples{modeMAP(parm=c(alpha=2,beta=2,rate=0.2),distr="exp",upper=1000)
modeMAP(parm=c(alpha=3,beta=2,rate=1),distr="exp",upper=1000)
}

\keyword{modified alpha power distribution, family of distributions}


