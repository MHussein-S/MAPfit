# Mode of MAPE by direct maximization of log(g(x))
modeMAPE<-function(start=0,parm,lower=0,upper)
  {
  alpha<-parm["alpha"]
  beta<-parm["beta"]
  rate<-parm["rate"]
  objfn <- function(param,x)
    {
    Logg <- log(alpha)+log(rate)-log(alpha*beta-1)+log(beta)*(1-exp(-rate*x))^2-log(alpha)*exp(-rate*x)-rate*x+log(log(alpha)+2*(1-exp(-rate*x))*log(beta))
    return(-Logg)
    }
  opt <- optim(par = start,fn=objfn, param = parm,method = "Brent",lower = lower,upper = upper)
  opt$par
  }

