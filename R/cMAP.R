# r-th central moment for any base distribution
cMAP<-function(r,par,distr)
{
  if (!is.list(par))
    par<-as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
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
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  if (!exists(ddistname, mode = "function"))
    paste("The " , distr, " disribution", "is not defined")
  parnm <- names(par)
  args <- names(formals(ddistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
  {
    stop(" 'par' must specify names which are arguments to ", distr)
  }
  else if ((alpha < 1) | (beta < 1))
  {
    stop("MAP distribution not defined for alpha and/or beta <=1")
  }
  else if (alpha*beta==1)
    stop("MAP distribution not defined for alpha*beta=1")
  #------------------ integrand ---------------
  integrand <- function(y,alpha,beta,distr,par,r)
  {#integrand
    qdistname <- paste("q", distr, sep = "")
    L1<-do.call(qdistname, c(list(y), as.list(par)))
    L1<-beta^(y^2)*alpha^y*(log(alpha)+2*y*log(beta))*(L1-mMAP(1,c(alpha=alpha,beta=beta,par),distr))^r
    L1<-L1/(alpha*beta-1)
  }#integrand
  #-------------------------------------------
  a<-integrate(integrand, lower = 0, upper =1,alpha=alpha,beta=beta,distr=distr ,par=par,r=r)
  return(a$value)
}
