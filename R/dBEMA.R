#' Estimate the number of latent factors in B.
#'
#' For this to work properly, the input coefficient matrix 'Stand.B' should be properly standardized.
#'
#' @param Stand.B Coefficient matrix.
#' @param N Sample size.
#' @param alpha Quantile of bulk eigenvalues.
#' @param n.avg Number of simulations to estimate the distribution of bulk eigenvalues.
#' @return The number of latent factors.
#' @export
dBEMA <- function(Stand.B, N, alpha=0.2, n.avg=10){  #l contains the eigenvalues
  if (nrow(Stand.B)<ncol(Stand.B)) {
    p <- nrow(Stand.B)
    S <- ncol(Stand.B)
    l <- svd(1/S*Stand.B%*%t(Stand.B),nu=0,nv=0)$d
  } else {
    S <- nrow(Stand.B)
    p <- ncol(Stand.B)
    l <- svd(1/S*t(Stand.B)%*%Stand.B,nu=0,nv=0)$d
  }
  n <- N
  rm(Stand.B,N)
  if (p>n) {l <- l[1:n]}
  o=optimize(loss,interval=c(0.1, 5),l=l,alpha=alpha,p=p,n=n,S=S,n.avg=n.avg)
  a=o$minimum
  L=matrix(data=0,nrow=500,ncol=p)
  for (i in 1:500){
    x1 <- matrix(rnorm(n*p),nrow=p,ncol=n)*sqrt(rgamma(p,a,1))
    if (!is.null(S)) {x1 <- sweep(x = x1, MARGIN = 2, STATS = sqrt(RMTstat::rmp(n = n, ndf = S, pdim = n)), FUN = "*")}
    if (p<=n){
      L[i,]=svd(x1%*%t(x1),nu=0,nv=0)$d/n
    } else {
      L[i,1:n]=svd(t(x1)%*%x1,nu=0,nv=0)$d/n
    }
  }

  L.quant <- t(apply(L,2,function(x){c(mean(x),quantile(x,c(0.1,0.5,0.75,0.9)))}))
  colnames(L.quant) <- c("Mean","10%","50%","75%","90%")
  #for (i in 1:p){
  #  l1[i]=mean(L[,i])
  #  l2[i]=quantile(L[,i],0.9)
  #  l3[i]=quantile(L[,i],0.1)
  #  l4[i]=quantile(L[,i],0.5)
  #}
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~L.quant[k,1]-1)$coef[[1]]
  plot(l)
  lines(L.quant[,1]*s1,col='red')
  lines(L.quant[,4]*s1,col='orange')
  lines(L.quant[,2]*s1,col='blue')
  lines(L.quant[,3]*s1,col='violet')
  abline(max(L.quant[,4]*s1),0,col='orange')
  K.hat <- sum(l>max(L.quant[,4]*s1))
  return(list(nFactors=K.hat,L=L*s1,L.quant=L.quant*s1,alpha=alpha,a=a,b=s1,eigs=l))
}

loss <- function(proposal,l,alpha,p,n,S=NULL,n.avg=10){
  L=matrix(data=0,nrow=n.avg,ncol=p)
  for (i in 1:n.avg){
    #wi=rWishart(n, p, diag(rgamma(p,proposal,1)))/n
    x1 <- matrix(rnorm(n*p),nrow=p,ncol=n)*sqrt(rgamma(p,proposal,1))
    if (!is.null(S)) {x1 <- sweep(x = x1, MARGIN = 2, STATS = sqrt(RMTstat::rmp(n = n, ndf = S, pdim = n)))}
    #x1=mvrnorm(n,rep(0,p),diag(rgamma(p,proposal,1)))
    if (p<=n){
      L[i,]=svd(x1%*%t(x1),nu=0,nv=0)$d/n
    }
    if (p>n){
      L[i,1:n]=svd(t(x1)%*%x1,nu=0,nv=0)$d/n
    }

  }
  l1=colMeans(L)
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  l1=s1*l1
  return(sum((l1-l)[k]^2))
}

Plot.Eigs <- function(out) {
  L.quant <- out$L.quant
  l <- out$eigs
  plot(l,xlab="Eigenvalue index",ylab="Eingevalue")
  lines(L.quant[,1],col='red')
  lines(L.quant[,4],col='orange')
  lines(L.quant[,2],col='blue')
  lines(L.quant[,3],col='violet')
  abline(max(L.quant[,4]),0,col='orange')
  return(sum(l>max(L.quant[,4])))
}

Plot.Diff <- function(Eigs, a, alpha, n.avg=10) {
  L=matrix(data=0,nrow=n.avg,ncol=p)
  for (i in 1:n.avg){
    #wi=rWishart(n, p, diag(rgamma(p,proposal,1)))/n
    x1 <- matrix(rnorm(n*p),nrow=p,ncol=n)*sqrt(rgamma(p,a,1))
    #x1=mvrnorm(n,rep(0,p),diag(rgamma(p,proposal,1)))
    if (p<=n){
      L[i,]=svd(x1%*%t(x1),nu=0,nv=0)$d/n
    }
    if (p>n){
      L[i,1:n]=svd(t(x1)%*%x1,nu=0,nv=0)$d/n
    }

  }
  l1=colMeans(L)
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  l1=s1*l1
  plot(l1[k],l[k],xlab="Est",ylab="True")
  abline(a=0,b=1,col="red")
}
