mMAP<-function(r,par,distr)
{
    if(!is.list(par))
      par<- as.list(par)
    if (is.null(names(par)))
      stop("'par' must be a named list")
    ddistname <- paste("d", distr, sep = "")
    qdistname <- paste("q", distr, sep = "")
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
  
  ########################## integrand ##########################
  integrand <- function(y,alpha,beta,qdistname,par,r)
  {
    L1<-do.call(qdistname, c(list(y), as.list(par)))
    L1<-beta^(y^2)*alpha^y*(log(alpha)+2*y*log(beta))*L1^r
    L1<-L1/(alpha*beta-1)
  }
  #-------------------------------------------
  a<-integrate(integrand, lower = 0, upper =1,alpha=alpha,beta=beta,qdistname=qdistname ,par=distpar,r=r)
  return(a$value)
}


