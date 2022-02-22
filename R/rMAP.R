# Generate random deviates from MAP distribution for any base distribution

rMAP<-function(n,par,distr)
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  pdistname<- paste("p", distr, sep = "")
  rdistname<- paste("r", distr, sep = "")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " distribution is not defined"))
  alphapar<-match("alpha",names(par))
  betapar<-match("beta",names(par))
  if(is.na(alphapar))
    stop(" 'alpha' parameter is not defined")
  if(is.na(betapar))
    stop(" 'beta' parameter is not defined")
  args <- names(formals(ddistname))
  alpha<-par$alpha
  beta<-par$beta
  distparn<-setdiff(names(par),c("alpha","beta"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specifies names of parameters which are not valid for ",ddistname)
  if ((alpha < 1) | (beta < 1))
  {
    stop("MAP distribution not defined for alpha and/or beta <=1")
  }
  else if (alpha*beta==1)
    stop("MAP distribution not defined for alpha*beta=1")
  
  ##########################
  x<-c()
  while (length(x)<n)
  {
    #1.	Generate a rv Y distributed as F
    Y<-do.call(rdistname, c(n=1, as.list(distpar)))
    # 2.	Generate U (independent of Y)
    U<-runif(1)
    C<-alpha*beta*(log(alpha)+2*log(beta))/(alpha*beta-1)
    #--
    CDF<-do.call(pdistname, c(list(Y), as.list(distpar)))
    PDF<-do.call(ddistname, c(list(Y), as.list(distpar)))
    MAPPDF<-1/(alpha*beta-1)*beta^(CDF^2)*alpha^CDF*PDF*(log(alpha)+2*CDF*log(beta))
  ##########################
    criteria<-MAPPDF/(C*PDF)
    if (U<=criteria)
      x<-c(x,Y)
  }
  return(x)
}
