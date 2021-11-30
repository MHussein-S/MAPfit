# Generate random deviates from MAP distribution for any base distribution

rMAP<-function(n,par,distr)
{
  if (is.null(names(par)))
    stop("'par' must be a named list or vector")
  ddistname <- paste("d", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  rdistname <- paste("r", distr, sep = "")
  if (!exists(ddistname, mode = "function"))
    stop("The ", ddistname, " function must be defined")
  if (!is.list(par))
    par<-as.list(par)
  #exclude alpha from par
  alphapar<-match("alpha",names(par))
  if(!is.na(alphapar))
  {
    alpha<-par$alpha
    par<-par[-alphapar] #exclude alpha parameter from par
  }
  else
    stop(" 'alpha' parameter not defined")
  #exclude beta from par
  betapar<-match("beta",names(par))
  if(!is.na(betapar))
  {
    beta<-par$beta
    par<-par[-betapar] #exclude beta parameter from par
  }
  else
    stop(" 'beta' parameter not defined")
  argdistname <- names(formals(ddistname))
  m <- match(names(par), argdistname)
  if (any(is.na(m)))
  {
    stop(" 'par' must specify names which are arguments to ", distr)
  }
  else if ((alpha < 1) | (beta < 1))
  {
    stop("MAP distribution not defined for alpha and/or beta <1")
  }
  else if (alpha*beta==1)
    stop("MAP distribution not defined for alpha*beta=1")
  x<-c()
  while (length(x)<n)
  {
    #1.	Generate a rv Y distributed as F
    Y<-do.call(rdistname, c(n=1, as.list(par)))
    # 2.	Generate U (independent of Y)
    U<-runif(1)
    C<-alpha*beta*(log(alpha)+2*log(beta))/(alpha*beta-1)
    #--
    CDF<-do.call(pdistname, c(list(Y), as.list(par)))
    PDF<-do.call(ddistname, c(list(Y), as.list(par)))
    MAPPDF<-1/(alpha*beta-1)*beta^(CDF^2)*alpha^CDF*PDF*(log(alpha)+2*CDF*log(beta))
    #--
    criteria<-MAPPDF/(C*PDF)
    if (U<=criteria)
      x<-c(x,Y)
  }
  return(x)
}
