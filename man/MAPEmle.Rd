\name{MAPEmle}
\alias{MAPEmle}
\title{Modified Alpha Power Distribution}
\description{MLE of parameters, KS test statistic and p-value, AIC,BIC, HQIC, -log-likelihood function and convergence status for MAPE distribution}
\value{Maximum likelihood estimates,Kolmogorov Smirnov test statistic and p-value, Akaike Information Criterion,Bayesian Information Criterion, Hannan-Quinn information criterion,  the negative log-likelihood function and convergence status for Modified Alpha-Powrer distribution.}
\usage{MAPEmle(data,start)}
\arguments{
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{start}{named vector or list of the starting values of the distribution parameters}
}

\references{}
\author{Mohamed Hussein}
\examples{
#{
data<-rMAP(n=500,par=c(rate=2,alpha=3,beta=2),distr="exp")
MAPEmle(data,start=c(1.1,1.1,1.1))
#}
}
\keyword{modified alpha power distribution, family of distributions}
