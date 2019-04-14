#MMRR models for landscape genetic analysis of IBE and IBR in Emel et al. submitted
#Sarah Emel

library(rdist)
library(vegan)
library(ade4)
library(psych)

#FUNCTIONS
#function to convert list of pairwise distances to matrix
listtomatrix <- function(x){
  an <- as.character(c(1:20))
  M <- array(0, c(length(an), length(an)), list(an, an))
  i <- match(data1$pop1, an)
  j <- match(data1$pop2, an)
  M[cbind(i,j)] <- M[cbind(j,i)] <- x
  y<-data.matrix(M)
  return(y)
}

#MMRR function
MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

#unfold converts the lower diagonal elements of a matrix into a vector
#unfold is called by MMRR
unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  return(x)
}


#GET IBD AND IBR
data1<-read.csv(file="Sang_IBD_IBR.csv")

#IBD: euclidean distance
IBD<-listtomatrix(data1$eucdist)

#IBR: top models across all populations
IBRlcpall<-listtomatrix(data1$topLCPall114)
IBRcsall<-listtomatrix(data1$topCSall106)

#GET GENETIC DISTANCE
#genetic distance
fst<-read.csv(file="FSToutliersremoved.csv",header=T)
Y<-data.matrix(fst)

#IBE AND MMRR
#read in soil data
data<-read.csv(file="envfinalall.csv") #site means
#Make individual names characters (not factors)
data$Pop <- as.character(data$Pop)
row.names(data)<-data$Pop #changing population names column to row names
data<-data[,-1] #remove population column
names(data) <- c("Cluster", "UTM_N", "UTM_E", "K", "P", "pH", "Ca_Mg", "Mg", "Cr", "Cu", "Ni","Pb", "Zn")
str(data) #look at the structure of the dataframe

#Confirm that genotypes and environmental data are in the same order:
identical(colnames(fst), rownames(data))

pairs.panels(data[,4:13], scale = F, cex.cor = 1)
#K is highly correlated with Cu, so remove it. We also want to remove Pb

pred <- data[,c(5:11,13)] #These are the variables to test in 'all'


#####################################

#All variables
#PCA
pc<-prcomp(~., data=pred,scale=T) #PCA 
summary(pc) #prints SD, prop. of var., and cumulative proportion for each PC axis
pc
plot(pc) #plots histogram of variance explained by each axis
biplot(pc) #plots location of variables on first two PCs

#PLOT THE PCA
#adapted from http://r.789695.n4.nabble.com/Color-points-in-princomp-plot-td4689621.html
myBiplot <- function (x, choices = 1L:2L, scale = 1,
                      pc.biplot = FALSE, var.axes = TRUE,
                      type = "t",
                      col,
                      col.arrows = "#FF0000",
                      col.text = "#000000",
                      cex = rep(par("cex"), 2),
                      expand = 1, 
                      xlabs = NULL, ylabs = NULL,
                      xlim = NULL, ylim = NULL, 
                      main = NULL, sub = NULL,
                      xlab = NULL, ylab = NULL, 
                      arrow.len = 0.1,
                      ...
)

{
  if (length(choices) != 2L) 
    stop("length of choices must be 2")
  if (!length(scores <- x$x)) 
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
         domain = NA)
  if (is.complex(scores)) 
    stop("biplots are not defined for complex PCA")
  
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  
  if (scale < 0 || scale > 1) 
    warning("'scale' is outside [0, 1]")
  if (scale != 0) 
    lam <- lam^scale
  else lam <- 1
  if (pc.biplot) 
    lam <- lam/sqrt(n)
  
  y <- t(t(x$rotation[, choices]) * lam)
  x <- t(t(scores[, choices])/lam)  #note that from here on
  #x is no longer the PC object
  #originally pased into the function
  n <- nrow(x)
  p <- nrow(y)
  
  if (missing(xlabs)) {
    xlabs <- dimnames(x)[[1L]]
    if (is.null(xlabs)) 
      xlabs <- 1L:n
  }
  xlabs <- as.character(xlabs)
  dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
  
  if (missing(ylabs)) {
    ylabs <- dimnames(y)[[1L]]
    if (is.null(ylabs)) 
      ylabs <- paste("Var", 1L:p)
  }
  ylabs <- as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
  
  if (length(cex) == 1L) 
    cex <- c(cex, cex)
  
  unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
                                  abs(max(x, na.rm = TRUE)))
  rangx1 <- unsigned.range(x[, 1L])
  rangx2 <- unsigned.range(x[, 2L])
  rangy1 <- unsigned.range(y[, 1L])
  rangy2 <- unsigned.range(y[, 2L])
  
  if (missing(xlim) && missing(ylim)) 
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
  else if (missing(xlim)) 
    xlim <- rangx1
  else if (missing(ylim)) 
    ylim <- rangx2
  
  ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  on.exit(par(op))
  op <- par(pty = "s")
  if (!is.null(main)) 
    op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
  
  #first, plot scores - either normally, or as row labels
  if (type == "p") {
    plot(x, type = type, xlim = xlim, ylim = ylim, col = col, 
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
  }
  else if (type == "t") {
    plot(x, type = "n", xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    text(x, xlabs, cex = cex[1L], col = col.text, ...)
    
  }
  else if (type == "n") {  #plot an empty frame
    plot(x, type = type, xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
  }
  
  par(new = TRUE)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * 
         ratio, xlab = "", ylab = "", col = col.arrows, ...)
  axis(3, col = col.arrows, ...)
  axis(4, col = col.arrows, ...)
  box(col = "#000000")
  text(y, labels = ylabs, cex = cex[2L], col = col.arrows, ...)
  if (var.axes) 
    arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = col.arrows, 
           length = arrow.len)
  #now replot into xlim, ylim scaled by lam, to reset par("usr")
  #to the  correct values needed for subsequent application of points(),
  #text() etc.
  par(new = TRUE)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  plot(0, type = "n", xlim = xlim * lam[1], ylim = ylim * lam[2], 
       xlab = '', ylab = '', sub = '', main = '', xaxt='n', yaxt='n', axes=FALSE)
  
  invisible()
  
}

myBiplot(pc, choices=1:2, cex = 0.7, col.text="#000000") #slightly nicer than default biplot

#calculate IBE from first 2 PCA axes
first2axes<-cbind(pc$x[,1],pc$x[,2])
first2axes<-as.matrix(first2axes)
IBE.all<-rdist(first2axes,"euclidean",p=2L)
IBE.all<-data.matrix(IBE.all)

#MMRR

#IBE.all
X<-list(IBE.all)
MMRR1<-MMRR(Y=Y,X=X,nperm=10000)
MMRR1

#top lcp all model, IBE.all
X<-list(IBRlcpall)
MMRR1<-MMRR(Y=Y,X=X,nperm=10000)
MMRR1

#top CS all model, IBE.all
X<-list(IBRcsall,IBE.all)
MMRR1<-MMRR(Y=Y,X=X,nperm=10000)
MMRR1