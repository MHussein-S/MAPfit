\name{modeMAPE}
\alias{modeMAPE}
\title{Modified Alpha Power Exponential Distribution}
\description{Mode of MAPE distribution by direct maximization of log(g(x))}
\value{Approximation of the Mode of MAPE distribution}
\usage{modeMAPE(start=0,parm,lower=0,upper)}
\arguments{
  \item{start}{Initial values for the parameters to be maximized over.}
  \item{parm}{named vector or list of the distribution parameters.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  }
\references{}
\author{Mohamed
Hussein}
\examples{modeMAPE(parm=c(alpha=2,beta=2,rate=0.2),upper=1000)
modeMAPE(parm=c(alpha=3,beta=2,rate=1),upper=1000)
}

\keyword{modified alpha power distribution, family of distributions}


