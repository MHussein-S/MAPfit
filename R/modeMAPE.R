# Mode of MAPE by 1- Dmax: direct maximization of log(g(x))
#2- rootF: Root finding of d/dx(log(g(x)))
modeMAPE<-function(start=0,parm,lower=0,upper,method="Dmax")
{
  alpha<-parm["alpha"]
  beta<-parm["beta"]
  rate<-parm["rate"]
  if (method=="Dmax")
    {
    objfn <- function(param,x)
      {
      Logg <- log(alpha)+log(rate)-log(alpha*beta-1)+log(beta)*(1-exp(-rate*x))^2-log(alpha)*exp(-rate*x)-rate*x+log(log(alpha)+2*(1-exp(-rate*x))*log(beta))
      return(-Logg)
      }
    opt <- optim(par = start,fn=objfn, param = parm,method = "Brent",lower = 0,upper = upper)
    Mode<-opt$par
    }
  else if (method=="rootF")
    {
    fn<-function (x,param)
      {
      CDF<-pexp(x,rate = rate)
      PDF<-dexp(x,rate)
      df<--rate^2*exp(-rate*x)
      dLogg <- 2*CDF*PDF*log(beta)+PDF*log(alpha)+df/PDF
      dLogg<-dLogg+(2*PDF*log(beta))/(log(alpha)+2*CDF*log(beta))
      }
    uni <- uniroot(f=fn, param=parm,lower = 0,upper = upper)
    Mode<-uni$root
  }
  }
