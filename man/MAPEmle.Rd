\name{MAPEmle}
\alias{MAPEmle}
\title{Modified Alpha Power Distribution}
\description{MLE of parameters, KS test statistic and p-value, AIC,BIC, HQIC, -log-likelihood function and convergence status for MAPE distribution}
\value{Maximum likelihood estimates, standard error, 95\% CI, Kolmogorov Smirnov test statistic and p-value, Akaike Information Criterion,Bayesian Information Criterion, Hannan-Quinn information criterion,  the negative log-likelihood function and convergence status for Modified Alpha-Powrer distribution.}
\usage{MAPEmle(data,start)}
\arguments{
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{start}{named vector or list of the starting values of the distribution parameters}
}

\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed Hussein}
\examples{
#{
data<-rMAP(n=500,par=c(rate=2,alpha=3,beta=2),distr="exp")
MAPEmle(data,start=c(1.1,1.1,1.1))
#}
}
\keyword{modified alpha power distribution, family of distributions}
