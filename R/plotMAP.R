#Empirical and theoretical: dens, CDF
#and P-P plot
plotMAP<-function (data, distr, para, histo = TRUE, breaks = "default", demp = TRUE)
{
  if (missing(data) || !is.vector(data, mode = "numeric"))
    stop("data must be a numeric vector")
  if ((missing(distr) & !missing(para)) || (!missing(distr) &  missing(para)))
    stop("distr and para -must defined")
  if (!histo & !demp)
    stop("one the arguments histo and demp must be put to TRUE")
  s <- sort(data)
  obsp<-ppoints(s) #observed probabilities
  if (length(s) != length(obsp))
    stop("problem when computing probabilities.")
  n <- length(data)
  if (missing(distr))
  {
    par(mfrow = c(1, 2))
    if (histo)
    {
      if (demp)
      {
        if (breaks == "default")
          h <- hist(data, freq = FALSE, xlab = "Data", main = "Empirical density")
        else
          h <- hist(data, freq = FALSE, xlab = "Data", main = "Empirical density", breaks = breaks)
        lines(density(data)$x, density(data)$y, lty = 1, col	= "red")
      }
      else
      {
        if (breaks == "default")
          h <- hist(data, freq = FALSE, xlab = "Data", main ="Histogram")
        else
          h <-hist(data, freq = FALSE, xlab = "Data", main ="Histogram", 	breaks	= breaks)
      }
    }
    else
    {
      h <- hist(data, freq = FALSE, plot = FALSE)

      plot(density(data)$x, density(data)$y, lty = 1, col = "red", type = "l", xlab = "Data", main = paste("Empirical density"), ylab = "Density")

    }
    plot(s, obsp, main = paste("Cumulative distribution"), xlab = "Data", xlim = c(min(h$breaks), max(h$breaks)), ylab = "CDF")
  }
  else #distr not missing
  {
    if (is.null(names(para)))
      stop("'para' must be a named list or vector")
    ddistname <- paste("d", distr, sep = "")
    pdistname <- paste("p", distr, sep = "")
    if (!exists(ddistname, mode = "function"))
      stop("The ", ddistname, " function must be defined")
    if (!is.list(para))
      para<-as.list(para)
    #exclude alpha from para
    alphapar<-match("alpha",names(para))
    if(!is.na(alphapar))
    {
      alpha<-para$alpha
      para<-para[-alphapar] #exclude alpha parameter from para
    }
    else
      stop(" 'alpha' parameter not defined")
    #exclude beta from para
    betapar<-match("beta",names(para))
    if(!is.na(betapar))
    {
      beta<-para$beta
      para<-para[-betapar] #exclude beta parameter from para
    }
    else
      stop(" 'beta' parameter not defined")
    argdistname <- names(formals(ddistname))
    m <- match(names(para), argdistname)
    if (any(is.na(m)))
    {
      stop(" 'para' must specify names which are arguments to ", distr)
    }
    else if ((alpha < 1) | (beta < 1))
    {
      stop("MAP distribution not defined for alpha and/or beta <1")
    }
    else if (alpha*beta==1)
      stop("MAP distribution not defined for alpha*beta=1")
    #-------------
    par(mfrow = c(2, 2))
    if (breaks == "default")
      h <- hist(data, plot = FALSE)
    else
    h <- hist(data, breaks = breaks, plot = FALSE)
    xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)

    CDF<-do.call(pdistname, c(list(xhist), as.list(para)))
    PDF<-do.call(ddistname, c(list(xhist), as.list(para)))
    MAPCDF<-(beta^(CDF^2)*alpha^CDF-1)/(alpha*beta-1)
    MAPPDF<-1/(alpha*beta-1)*beta^(CDF^2)*alpha^CDF*PDF*(log(alpha)+2*CDF*log(beta))
    yhist <- MAPPDF
    if (length(yhist) != length(xhist))
      stop("problem when computing densities.")
    ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(max(h$density),max(density(data)$y)))
    if (histo)
    {
      hist(data, freq = FALSE, xlab = "Data", ylim = c(0, ymax), breaks = h$breaks, main = paste("Empirical and theoretical dens."),ylab= "Density", xlim = c(min(h$breaks), max(h$breaks)))
      lines(xhist,yhist , lty = 2, col = "black", type = "l", xlab = "Data", main = paste("Empirical and theoretical dens."), ylab= "Density", 	xlim = c(min(h$breaks), max(h$breaks)))
      if (demp)
      {
       lines(stats::density(data)$x, stats::density(data)$y, lty = 1, col = "red")
      legend("topright", bty = "n", lty = c(2,1), col = c("black", "red"), legend = c("theoretical","empirical"), bg = "white", cex = 0.7)
      }
    }
    else
    {
      plot(xhist,yhist , lty = 2, col = "black", type = "l", xlab = "Data", main = paste("Empirical and theoretical dens."), ylab= "Density", 	xlim = c(min(h$breaks), max(h$breaks)))
      if (demp)
      {
        lines(stats::density(data)$x, stats::density(data)$y, lty = 1, col = "red",xlim = c(min(h$breaks), max(h$breaks)))
        legend("topright", bty = "n", lty = c(2,1), col = c("black", "red"), legend = c("theoretical","empirical"), bg = "white", cex = 0.7)
      }
    }
    xmin <- min(h$breaks)
    xmax <- max(h$breaks)
    plot(s, obsp, main = paste("Empirical and theoretical CDFs"),  xlab = "Data", 	ylab = "CDF", xlim = c(xmin, xmax)) #Empirical CDF

     sfin <- seq(xmin, xmax, by = (xmax - xmin)/1000)
    CDF_fin<-do.call(pdistname, c(list(sfin), as.list(para))) #Theoretical CDF
    PDF_fin<-do.call(ddistname, c(list(sfin), as.list(para)))
    MAPCDF_fin<-(beta^(CDF_fin^2)*alpha^CDF_fin-1)/(alpha*beta-1)
    MAPPDF_fin<-1/(alpha*beta-1)*beta^(CDF_fin^2)*alpha^CDF_fin*PDF_fin*(log(alpha)+2*CDF_fin*log(beta))
    graphics::lines(sfin, MAPCDF_fin, lty = 1, col = "red")
    CDF_s<-do.call(pdistname, c(list(s), as.list(para)))
    PDF_s<-do.call(ddistname, c(list(s), as.list(para)))
    MAPCDF_s<-(beta^(CDF_s^2)*alpha^CDF_s-1)/(alpha*beta-1)
    MAPPDF_s<-1/(alpha*beta-1)*beta^(CDF_s^2)*alpha^CDF_s*PDF_s*(log(alpha)+2*CDF_s*log(beta))
    if (length(MAPCDF_s) != length(obsp))
      stop("problem when computing probabilities.")
    plot(MAPCDF_s, obsp, main = "P-P plot", xlab = "Theoretical probabilities", ylab = "Empirical probabilities")

    abline(0, 1, col="red")
  }
}


