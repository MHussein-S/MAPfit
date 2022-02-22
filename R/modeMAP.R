# Mode of MAPE by direct maximization of log(g(x))
modeMAP<-function(start=0,parm,distr,lower=0,upper)
{
  if(is.null(names(parm))|any(names(parm)==""))
    stop("'parm' must be a named numeric vector of the form 'c(name=val,name=val,...)'")
  if (!is.character(distr))
    stop("distr must be a character string naming the baseline distribution")
  setpar<-function(ddistname,allpar)
  {
    if(!is.list(allpar))
      allpar<- as.list(allpar)
    if (!exists(ddistname, mode="function"))
      stop(paste("The ", ddistname, " distribution is not defined"))
    args <- names(formals(ddistname))
    alpha<-allpar$alpha
    beta<-allpar$beta
    
    extrpar<-list(alpha=alpha,beta=beta)
    distparn<-setdiff(names(allpar),c("alpha","beta"))
    distpar<-allpar[distparn]
    m <- match(distparn,args)
    if (any(is.na(m)))
      stop("you specifies names of parameters which are not valid for ",ddistname)
    return(list(extrpar=extrpar,distpar=distpar))
  }
  objfn <- function(param,x,distr)
    {
    ddistname <- paste("d",distr,sep="")
    pdistname <- paste("p",distr,sep="")
    parset<-setpar(ddistname,allpar=as.list(param))
    alpha<-parset$extrpar$alpha
    beta<-parset$extrpar$beta
    distpar<-parset$distpar
    ddistname <- paste("d",distr,sep="")
    pdistname <- paste("p",distr,sep="")
    F<-do.call(pdistname,c(list(x),as.list(distpar)))
    f<-do.call(ddistname,c(list(x),as.list(distpar)))
    Logg<--log(alpha*beta-1)+F^2*log(beta)+F*log(alpha)+log(f)+log(log(alpha)+2*F*log(beta))
    return(-Logg)
  }
  opt <- optim(par = start,fn=objfn, param = parm,distr=distr,method = "Brent",lower = lower,upper = upper)
  
  opt$par
}







