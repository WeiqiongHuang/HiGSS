#' Simulate metabolite GWAS data from the factor model
#'
#' @param IDs A data.frame object specifying the pathway structure. It will be randomly generated if unspecified.
#' @param J Number of superpathways used to generate IDs. Default is 5.
#' @param K Number of latent factors. Default is 10.
#' @param N Sample size. Default is 1000.
#' @param S Number of genetic variants. Default is 1e4.
#' @param p.l.spike Probability of generating spike pathways.
#' @param p.out Probability of generating outliers.
#' @param common.spikes A vector of subpathway names. All entries of L belonging to those subpathway are set as 0. Default is NULL.
#' @param regenerate.paras A list containing parameters to regenerate L. Default is NULL.
#' @param fix.L.scaling Default is FALSE. The function will automatically scale L to make the signal strength compariable to noise level. If it's set as TRUE, 'L.scaling' must also be specified.
#' @param L.scaling User specificed scalar to scale L.
#' @param G.model Distribution of non-zero genetic effects. Can take values of either 'gaussian' or 'laplace'.
#' @param p.G.nonzero Probability of generating non-zero entries in G. Default is 0.005.
#' @param p.Delta.0 Probability of non-genetically regulated metabolites in Delta. Default is 0.2.
#' @return A list containg matrices in mtGWAS factor analysis model.
#' @export
simuFactorModel <- function(IDs=NULL,J=5,K=10,N=1000,S=1e4,G.model=c("gaussian","laplace"),p.l.spike=0.05,p.G.nonzero=0.005,p.out = 0.01,p.Delta.0=0.2,common.spikes = NULL,regenerate.L=FALSE,regenerate.paras = NULL,fix.L.scaling = FALSE, L.scaling = NULL){
  if(is.null(IDs)){
    #simulate J superpathways
    IDs <- simuPathway(J)
  }
  m <- nrow(IDs)

  #simulate L
  if(regenerate.L==TRUE){
    L <- matrix(data = NA,nrow = m,ncol = K)
    for (k in 1:K) {
      L[,k] <- rnorm(n = m,mean = regenerate.paras$L.parameters[[k]]$mu,regenerate.paras$L.parameters[[k]]$sd)
    }
    L <- L[,regenerate.paras$column.orders]
  }else{
    p.0 <- rbeta(n = K,shape1 = 2,shape2 = 1)
    L <- matrix(data = NA,nrow = m,ncol = K)
    L.parameters <- list()
    for (k in 1:K){
      simu.L <- rHDP(IDs = IDs,gamma = 4,alpha0 = 2,p.00 = 0.05,b.l = -6,b.u = 6,p.0 = p.0[k],p.out = p.out,tau = 10,sigma.alpha = 3,sigma.beta = 2,common.spikes = common.spikes)
      L.parameters[[k]] <- simu.L
      L[,k] <- simu.L$L
    }
    column.orders <- order(colSums(L^2),decreasing = T)
    L <- L[,column.orders]
    regenerate.paras <- list(L.parameters,column.orders);names(regenerate.paras) <- c("L.parameters","column.orders")
  }
  #simulate G
  sigSNPs <- LaplacesDemon::rbern(n = S*K,prob = p.G.nonzero)
  sigma.alt <- runif(n = K,min = 6,max = 9)
  if(G.model=="gaussian"){
    G <- matrix(data=ifelse(sigSNPs==1,rnorm(n = S*K,mean = 0,sd = sigma.alt[k]/sqrt(N)),0), nrow = S, ncol = K,byrow = T)
  }else if(G.model=="laplace"){
    G <- matrix(data=ifelse(sigSNPs==1,rlaplace(n = S*K,mu = 0,sigma = sigma.alt[k]/sqrt(2*N)),0), nrow = S, ncol = K,byrow = T)
  }
  X <- simuGenotype(N = N,S = S)$X
  X <- X-rowMeans(X)%*%t(rep(1,N))
  Dss = rowSums(X^2)
  Xi <- matrix(data = rnorm(n = N*K,mean = 0,sd = 1),nrow = N,ncol = K)
  C <- t(X)%*%G+Xi

  #scale L
  SigmaSquare <- rgamma(n = m,shape = 1,rate = 1)# mean = 1
  E <- matrix(data = rnorm(n = N*m,mean = 0,sd = rep(sqrt(SigmaSquare),N)),nrow = m,ncol = N)
  if(fix.L.scaling==FALSE){
    L <- L/sqrt(eigen(L%*%t(L))[[1]][K]*N/eigen(E%*%t(E))[[1]][1]/1.1)
  }else{
    L <- L*L.scaling
  }


  #simulate Delta
  Delta <- simuDelta(S = S,IDs = IDs,Dss = Dss,p0 =p.Delta.0,SigmaSquare = SigmaSquare,gamma = 5,alpha = 5)
  #max(rowSums(Delta!=0))

  #Y=LC+DeltaX+E
  Y <- as.matrix(L%*%t(C))+Delta%*%X+E
  #row-centering
  Y <- Y-rowMeans(Y)%*%t(rep(1,N))
  B.hat <- Y%*%t(X);B.hat <- sapply(1:S, function(s) B.hat[,s]/sqrt(Dss[s]))

  result <- list(IDs,L,G,X,Xi,C,Delta,E,Y,B.hat,regenerate.paras);names(result) <- c("metabolites","L","G","X","Xi","C","Delta","E","Y","B.hat","regenerate.paras")
  return(result)
}



#' Simulate a super-pathway and sub-pathway file
#'
#' @details
#' Only number of super-pathways is required to be specified, the number of subpathways and total number of metabolites will be randomly generated.
#'
#' @param J Number of super-pathways to simulate
#' @return A data frame of three coloumns, namely the superpathway, subpathway and metabolite indices.
#' @examples IDs <- simuPathway(5)
#' @export
simuPathway <- function(J){
  SUPER_PATHWAY <- vector()
  SUB_PATHWAY <- vector()
  for (j in 1:J) {
    G <- rpois(n = 1,lambda = 8)
    for (g in 1:G) {
      mg <- ifelse(rpois(n = 1,lambda = 10)>13,rpois(n = 10,lambda = 30),rpois(n = 5,lambda = 10))
      SUB_PATHWAY <- c(SUB_PATHWAY, rep(paste0(j,":",g),mg))
      SUPER_PATHWAY <- c(SUPER_PATHWAY, rep(j,mg))
    }
  }
  IDs <- (cbind(SUPER_PATHWAY,SUB_PATHWAY,c(1:length(SUPER_PATHWAY))))
  colnames(IDs) <- c("SUPER_PATHWAY","SUB_PATHWAY","metabolite")
  IDs <- as.data.frame(IDs)
  return(IDs)
}

#simulate L from HDP
#IDs: contains pathway information
#gamma,alpha0: concentration parameters
#p.00: prob of spike
#b.l,b.u: upper and lower bounds for outliers
#p.0: prob of non-spike but 0 mean
#p.out: prob of outliers
#tau: parameter of hourse shoe prior
#sigma.alpha,sigma.beta: parameters of gamma prior for variances
rHDP <- function(IDs,gamma,alpha0,p.00,b.l,b.u,p.0,p.out,tau,sigma.alpha,sigma.beta,common.spikes = NULL){
  J <- length(unique(IDs$SUPER_PATHWAY))
  L <- vector();mu <- vector();sd <- vector()
  M.k <- vector()#number of tables serving dish k
  K.g <- vector()#dish assignments of subpathways
  theta <- matrix(data = NA,nrow = 1,ncol = 2);theta <- theta[-1,]
  sup <- unique(IDs$SUPER_PATHWAY)
  is.spike <- vector()
  is.out <- vector()
  for (j in 1:J) {
    T.j <- matrix(data = NA,nrow = 1,ncol = 2);T.j <- T.j[-1,]#dish assignments of tables
    sub <- unique(IDs$SUB_PATHWAY[IDs$SUPER_PATHWAY== sup[j]])
    for (g in 1:length(sub)) {
      mg <- sum( IDs$SUB_PATHWAY== sub[g] & IDs$SUPER_PATHWAY== sup[j])
      if(sub[g]%in%common.spikes){
        spike <- 1
      }else{
        spike <- LaplacesDemon::rbern(n = 1,prob = p.00)
      }
      if(spike==1){
        L <- c(L,rep(0,mg));mu <- c(mu,rep(0,mg));sd <- c(sd,rep(0,mg))
        K.g <- c(K.g,0)
        #is.spike <- c(is.spike,rep(1,mg))
        is.spike <- c(is.spike,1)
        is.out <- c(is.out,rep(0,mg))
      }else{
        #is.spike <- c(is.spike,rep(0,mg))
        is.spike <- c(is.spike,0)
        t.g <- sample(x = 1:(length(T.j[,1])+1),size = 1,prob = c(T.j[,1],alpha0))
        if(t.g == length(T.j[,1])+1){#new table
          T.j <- rbind(T.j,c(0,0))
          T.j[t.g,1] <- 1
          #assign dish to new table
          k.t <- sample(x = 1:(length(M.k)+1),size = 1,prob = c(M.k,gamma))
          if(k.t==length(M.k)+1){#new dish
            M.k[k.t] <- 1
            new.mu <- ifelse(LaplacesDemon::rbern(n = 1,prob = p.00)==1,0, rnorm(n = 1,mean = 0,sd = abs(rcauchy(n = 1,location = 0,scale = 1))*tau) )
            while (abs(new.mu)>max(abs(b.l),abs(b.u))) {
              new.mu <- ifelse(LaplacesDemon::rbern(n = 1,prob = p.00)==1,0, rnorm(n = 1,mean = 0,sd = abs(rcauchy(n = 1,location = 0,scale = 1))*tau) )
            }
            new.sigma <- sqrt(rgamma(n = 1,shape = sigma.alpha,rate = sigma.beta))
            theta <- rbind(theta,c(new.mu,new.sigma))
          }else{
            M.k[k.t] <-M.k[k.t] + 1
          }
          T.j[t.g,2] <- k.t
        }else{
          T.j[t.g,1] <- T.j[t.g,1]+1
        }
        k.g <- T.j[t.g,2]
        K.g <- c(K.g,k.g)
        mu.g <- theta[k.g,1];sigma.g <- theta[k.g,2]
        outliers <- rbinom(n = mg,size = 1,prob = p.out)
        L.out <- ifelse(LaplacesDemon::rbern(n = mg,prob = 0.5)==1,runif(n = mg,min = b.l,max = b.l+2),runif(n = mg,min = b.u-2,max = b.u))
        L.in <- rnorm(n = mg,mean = mu.g,sd = sigma.g)
        L.g <- ifelse(outliers==1,L.out,L.in)
        L.g <- ifelse((L.g-mu.g+3*sigma.g)*(mu.g+3*sigma.g-L.g)>0,L.in,L.out)

        is.out <- c(is.out,(L.g-mu.g+3*sigma.g)*(mu.g+3*sigma.g-L.g)<0)
        L <- c(L,L.g);mu <- c(mu,rep(mu.g,mg));sd <- c(sd,rep(sigma.g,mg))

      }
    }
  }
  result <- list(L,mu,sd,theta,is.spike,is.out,K.g);names(result) <- c("L","mu","sd","Theta","spikes","outliers","clusters")
  return(result)
}
rGM <- function(S,K){
  C <- matrix(nrow = S,ncol = K)
  for (k in 1:K) {
    C[,k] <- ifelse(LaplacesDemon::rbern(n = S,prob = sample(x = c(15:30),size = 1)/S)==1, rnorm(n = S,mean = 0,sd = 2), rnorm(n = S,mean = 0,sd = 1) )
  }
  return(C)
}

#simulate independent genotype matrix
simuGenotype <- function(N,S){
  X <- matrix(data = NA,nrow = S,ncol = N)
  MAF <- runif(n = S,min = 0.05,max = 0.95)
  for (s in 1:S) {
    X[s,] <- rbinom(n = N,size = 2,prob = MAF[s])
  }
  result <- list(X,MAF);names(result) <- c("X","MAF")
  return(result)
}

#simulate Delta
#S: number of SNPs
#IDs: contains pathway information
#gamma,alpha: concentration parameters
#Dss=diag(XTX)
#p.0: prob of pi_m=0
#SigmaSquare: diagnoal elements of noise covariance matrix Sigma
simuDelta <- function(S,IDs,gamma=7,alpha=5,Dss,p0=0.2,SigmaSquare){
  phi2 <- 25
  # Pi <- c(0,99)/10000
  # H.Pi <- c(p0,(1-p0)*(exp(-20*Pi[-1])/sum(exp(-20*Pi[-1]))))
  Pi <- exp(seq(log(1/10000),log(99/10000),length.out =99))
  H.Pi <- (1-p0)*(exp(-300*Pi)/sum(exp(-300*Pi)))
  Pi <- c(0,Pi);H.Pi <- c(p0,H.Pi)
  #plot(x = Pi,y = H.Pi)

  G <- length(unique(paste0(IDs$SUPER_PATHWAY,":",IDs$SUB_PATHWAY)))#number
  M.k <- matrix(data = 0,nrow = G,ncol = length(Pi));colnames(M.k) <- Pi
  Delta <- matrix(data = 0,nrow = nrow(IDs),ncol = S)
  for (j in 1:G) {
    #j=2
    T.j <- matrix(data = NA,nrow = 1,ncol = 3);T.j <- T.j[-1,]
    Metabolites.j <- which(paste0(IDs$SUPER_PATHWAY,":",IDs$SUB_PATHWAY)==unique(paste0(IDs$SUPER_PATHWAY,":",IDs$SUB_PATHWAY))[j])
    for (m in Metabolites.j) {
      #sample T
      #m <- 1
      if(nrow(T.j)==0){
        #sample theta for T.new
        theta.t <- sample(x = Pi,size = 1,prob = LaplacesDemon::rdirichlet(n = 1,alpha = colSums(M.k)+gamma*H.Pi))
        T.j <- rbind(T.j,c(nrow(T.j)+1,theta.t,1))
        M.k[j,which(Pi==theta.t)] <- M.k[j,which(Pi==theta.t)]+1
        theta.m <- theta.t
      }else{
        T.new <- LaplacesDemon::rbern(n = 1,prob = sum(T.j[,3])/(sum(T.j[,3])+alpha))
        if(T.new){
          theta.t <- sample(x = Pi,size = 1,prob = LaplacesDemon::rdirichlet(n = 1,alpha = colSums(M.k)+gamma*H.Pi))
          T.j <- rbind(T.j,c(nrow(T.j)+1,theta.t,1))
          M.k[j,which(Pi==theta.t)] <- M.k[j,which(Pi==theta.t)]+1
          theta.m <- theta.t
        }else{
          T.m <- sample(x = T.j[,1],size = 1,prob = T.j[,3]/sum(T.j[,3]))
          T.j[which(T.j[,1]==T.m),3] <- T.j[which(T.j[,1]==T.m),3]+1
          theta.m <- T.j[which(T.j[,1]==T.m),2]
        }
      }
      if(theta.m!=0){
        #theta.m <- 1-0.999
        Delta[m,] <- ifelse( LaplacesDemon::rbern(n = S,prob = theta.m)==1,rnorm(n = S,mean = 0,sd = sqrt(phi2*SigmaSquare[m]/Dss  )) ,0 )
      }
    }
  }
  return(Delta)
}


