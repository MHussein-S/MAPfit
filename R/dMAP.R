#pdf of MAP for any base distribution
dMAP<-function (x,par,distr, log = FALSE)
{
  if(is.null(names(par)))
    stop(" 'par' must be a named list.")
  ddistname <- paste("d", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  argdistname <- names(formals(ddistname))
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
  m <- match(names(par), argdistname)
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
  fbase<-do.call(ddistname, c(list(x), as.list(par)))
  Fbase<-do.call(pdistname, c(list(x), as.list(par)))
  d<-1/(alpha*beta-1)*beta^(Fbase^2)*alpha^Fbase*fbase*(log(alpha)+2*Fbase*log(beta))
  if (log)
    d <- log(d)
  return(d)
}
#------------------------------------------------------
