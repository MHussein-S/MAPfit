#Quantile function of MAP distribution
qMAP<-function(p,par,distr,lower=0,upper,lower.tail=TRUE,log.p=FALSE)
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
  pdistname <- paste("p", distr, sep = "")
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
  if (log.p==TRUE)
    p<-exp(p)
  if (lower.tail==FALSE)
    p<-1-p
  qfun<-function (x,p,alpha,beta,distpar,pdistname)
  {
    CDF<-do.call(pdistname, c(list(x), as.list(distpar)))
    MAPCDF<-((beta^(CDF^2)*alpha^CDF-1)/(alpha*beta-1))
    MAPCDF-p
    }
  #--
  res<-c()
  for (pi in p)
    {
    uni <- uniroot(f=qfun, p=pi, alpha=alpha,beta=beta, distpar=par, pdistname=pdistname,f.lower = -pi,f.upper = pi,lower = lower,upper = upper)$root
    res<-c(res,uni)
    }
  return(res)
  }


