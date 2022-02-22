pMAP<-function (q,par,distr,lower.tail = TRUE, log.p = FALSE )
{
  if(!is.list(par))
    par<- as.list(par)
  if (is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  qdistname <- paste("q", distr, sep = "")
  pdistname<- paste("p", distr, sep = "")
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
  ######################################
  Fbase<-do.call(pdistname, c(list(q), as.list(distpar)))
  p<-(beta^(Fbase^2)*alpha^Fbase-1)/(alpha*beta-1)
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  return(p)
}
