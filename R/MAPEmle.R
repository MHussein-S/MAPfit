MAPEmle<-function (data,start)
{
  #total gradient of log-likelihood function of MAPE
  grlnlMAP <- function(par, obs, ...)
  {
    #gradient of log-likelihood function of MAPE "individual contribution"
    grdMAP <- function(par,x)
    {
      alpha<-par[1]
      beta<-par[2]
      rate<-par[3]
      dalpha<-1/alpha-beta/(alpha*beta-1)-1/alpha*exp(-rate*x)+1/alpha*(1/(log(alpha)+2*(1-exp(-rate*x))*log(beta)))
      dbeta<--alpha/(alpha*beta-1)+1/beta*(1-exp(-rate*x))^2+2*(1-exp(-rate*x))/beta*(1/(log(alpha)+2*(1-exp(-rate*x))*log(beta)))
      drate<-1/rate+2*log(beta)*x*exp(-rate*x)*(1-exp(-rate*x))+log(alpha)*x*exp(-rate*x)-x+2*log(beta)*x*exp(-rate*x)/(log(alpha)+2*(1-exp(-rate*x))*log(beta))
      res<-c(dalpha, dbeta,drate)
      return(res)
    }
    -rowSums(sapply(obs, function(x) grdMAP(par=par,x)))
  }
  # -log-likelihood function of MAPE distribution
  objfn <- function(parm,obs,...)
  {
    n<-length(obs)
    alpha<-parm[1]
    beta<-parm[2]
    rate<-parm[3]
    if ((alpha < 1) | (beta < 1))
      stop("MAP distribution not defined for alpha and/or beta <1")
    else if (alpha*beta==1)
      stop("MAP distribution not defined for alpha*beta=1")
    LL <- n*log(alpha)+n*log(rate)-n*log(alpha*beta-1)+log(beta)*sum((1-exp(-rate*obs))^2)-log(alpha)*sum(exp(-rate*obs))-rate*sum(obs)+sum(log(log(alpha)+2*(1-exp(-rate*obs))*log(beta)))
    return(-LL)
  }
  #variance covarince matrix
  vCmatrix <- function(par, obs)
  {
    n<-length(obs)
    alpha<-par[1]
    beta<-par[2]
    lambda<-par[3]
    domn<-log(alpha)+2*(1-exp(-lambda*obs))*log(beta)
    val1<-sum(1/domn)
    val2<-sum(1/domn^2)
    
    d2alpha<- -n/alpha^2 + n*beta^2/(alpha*beta-1)^2 + 1/alpha^2*sum(exp(-lambda*obs)) - 1/alpha^2*val2 - 1/alpha*val1
    
    d2alpha_beta<- n/(alpha*beta-1)^2 - 2/(alpha*beta)*sum((1-exp(-lambda*obs))/domn^2)
    
    d2alpha_lambda<-1/alpha*sum(obs*exp(-lambda*obs)) - 2*log(beta)/alpha*sum(obs*exp(-lambda*obs)/domn^2)
    
    d2beta<-n*alpha^2/(alpha*beta-1)^2 - 1/beta^2*sum((1-exp(-lambda*obs))^2) -4/beta^2*sum((1-exp(-lambda*obs))^2/domn^2) - 2/beta^2*sum((1-exp(-lambda*obs))/domn)
    
    d2beta_lambda<-2/beta*sum(obs*exp(-lambda*obs)*(1-exp(-lambda*obs))) - 4*log(beta)/beta*sum(obs*exp(-lambda*obs)*(1-exp(-lambda*obs))/domn^2)+2/beta*sum(obs*exp(-lambda*obs)/domn)
    
    d2lambda<- -n/lambda^2 + 2*log(beta)*sum(obs^2*exp(-2*lambda*obs)) - 2*log(beta)*sum(obs^2*exp(-lambda*obs)*(1-exp(-lambda*obs))) - 2*log(beta)*sum(obs^2*exp(-lambda*obs)) - 4*(log(beta))^2*sum(obs^2*exp(-2*lambda*obs)/domn^2) - 2*log(beta)*sum(obs^2*exp(-lambda*obs)/domn)
    
    I<-matrix(c(-d2alpha,-d2alpha_beta,-d2alpha_lambda,-d2alpha_beta,-d2beta,-d2beta_lambda,-d2alpha_lambda,-d2beta_lambda,-d2lambda),nrow = 3,ncol = 3,byrow = T)
    solve(I)
    
  }
  npar <- 3
  Mat <- diag(3)
  colnames(Mat) <- c("alpha","beta","rate")
  rownames(Mat) <- paste0("constr", 1:3)
  initconstr <- Mat %*% start - c(1,1,0)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=objfn,grad = grlnlMAP, obs = data , ui=Mat,ci = c(1,1,0),hessian = T)
  if (is.null(names(opt$par)))
    names(opt$par) <- c("alpha","beta","rate")
  #---
  A<-vCmatrix(opt$par,data)
  prop_sigma<-sqrt(diag(A))
  upper<-opt$par+qnorm(0.975)*prop_sigma
  lower<-opt$par+qnorm(0.025)*prop_sigma
  #---
  loglik <- -opt$value
  k<-npar
  n<-length(data)
  AIC<--2*loglik+2*k
  BIC<- -2*loglik+ k*log(n)
  HQIC<- -2*loglik+2*k*log(log(n))
  res=cbind(opt$par,prop_sigma,lower,upper)
  colnames(res)=c("MLE","Std. Err", "Inf. 95% CI","Sup. 95% CI")
  res1=cbind(AIC,BIC, HQIC, opt$value)
  colnames(res1)=c("AIC","BIC","HQIC", "-log(Likelihood)")
  rownames(res1)=c("")
  yemp<-ecdf(data)
  ytheo<-pMAP(sort(data),par=opt$par,'exp' )
  KS<-ks.test(yemp(data),ytheo)
  res2=cbind(KS$statistic,KS$p.value)
  colnames(res2)=c("KS Statistic","KS p-value")
  rownames(res2)=c("")
  res3=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res3)=c("")
  rownames(res3)=c("")
  list("Estimates"=res,"Measures"=res1,"Kolmogorov-Smirnov Test"=res2,"Convergence Status"=res3)
}