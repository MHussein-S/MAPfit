\name{plotMAP}
\alias{plotMAP}
\title{Modified Alpha Power Distribution}
\description{Plots the empirical and theoretical density, CDF and P-P plot}
\value{Graph of empirical versus theoretical cdf, pdf}
\usage{plotMAP(data, distr, para, histo = TRUE, breaks = "default", demp = TRUE)}
\arguments{

  \item{data}{a vector of data values for which the distribution was fitted}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{para}{named vector or list of the distribution parameters}
  \item{histo}{A logical to plot the histogram using the hist function}
  \item{breaks}{If "default" the histogram is plotted with the function hist with its default
breaks definition. Else breaks is passed to the function hist.}
  \item{demp}{A logical to plot the empirical density on the first plot (alone or superimposed on
the histogram depending of the value of the argument histo) using the density
function.}
  }
\references{Hussein, M.; Elsayed, H.; Cordeiro, G.M. A New Family of Continuous Distributions: Properties and Estimation. Symmetry 2022, 14, 276. https://doi.org/10.3390/sym14020276}
\author{Mohamed Hussein}
\examples{
#{
data<-rMAP(n=500,par=c(rate=2,alpha=3,beta=2),distr="exp")
plotMAP(data,distr = "exp",para=c(rate=2,alpha=3,beta=2))
plotMAP(data,distr = "exp",para=c(rate=2,alpha=3,beta=2),histo=FALSE)
#}
}
\keyword{modified qlpha power distribution, family of distributions}
