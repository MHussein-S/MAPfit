dMAP<-function (x,par,distr, log = FALSE)
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  
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
  fbase<-do.call(ddistname, c(list(x), as.list(distpar)))
  Fbase<-do.call(pdistname, c(list(x), as.list(distpar)))
  d<-1/(alpha*beta-1)*beta^(Fbase^2)*alpha^Fbase*fbase*(log(alpha)+2*Fbase*log(beta))
  if (log)
    d <- log(d)
  return(d)
}
#------------------------------------------------------
