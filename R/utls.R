estimate_mode <- function(x) {
  if(length(x)==1){
    return(x)
  }else{
    d <- density(x)
    return(d$x[which.max(d$y)])
  }
}

addSpace <- function(x,p,w){
  y <- x[1:(p[1]-1)]
  gap.centers <-length(y)+w/2
  y <- c(y,rep(NA,w))
  for (i in 1:length(p)) {
    if(i<length(p)){upper = p[i+1]-1}else{upper = length(x)}
    y <- c(y,x[p[i]:upper])
    if(i<length(p)){
      gap.centers <- c(gap.centers,length(y)+w/2)
      y <- c(y,rep(NA,w))
    }
  }
  return(list(y,gap.centers))
}

matrix.up <- function(x){
  if(is.matrix(x)){
    x[col(x) > row(x)]
  }else{
    x
  }
}


#' Combining the subpathways sharing parameters with large probability to form clusters
#' @param w The similarity matrix. It's rows and columns should be named.
#' @param method The method used to measure the similarity between one element and a group of elements. Can take values of "average", "min" and "max". Default is "max".
#' @param include.singleton Should singleton subpathways be included in estimated clusters. Default is False.
#' @param p.singleton If 'include.singleton' is True and 'p.singleton' is provided, only subpathways that do not share parameter with other pathways with large probability are considered as true singletons.
#' @return A list of vectors, each vector contains the elements in the estimated cluster.
#' @export
graphicalCluster <- function(w,method="max",threshold,include.singleton=F,p.singleton=NULL){
  clusters <- list()
  x <- w;r=nrow(x)
  while (r>0) {
    if(nrow(x)==1){
      r <- 0
    }else if(nrow(x)==2){
      if(x[1,2]>=threshold){
        clusters <- append(clusters,list(rownames(x)))
      }else{
        if(include.singleton==T){
          if(is.null(p.singleton)){
            clusters <- append(clusters,list(rownames(x)[1]));clusters <- append(clusters,list(rownames(x)[2]))
          }else{
            if(p.singleton[1]>=threshold){
              clusters <- append(clusters,list(rownames(x)[1]))
            }
            if(p.singleton[2]>=threshold){
              clusters <- append(clusters,list(rownames(x)[2]))
            }
          }
        }

      }
      r <- 0
    }else{
      cluster <- vector()
      continue <- T
      while(continue==T){
        if(length(cluster)==0){
          if(max( matrix.up(x) )>=threshold){
            pairs <- which(x == max( matrix.up(x) ) , arr.ind = TRUE)
            pairs <- pairs[pairs[,1]!=pairs[,2],]
            init.pair <- pairs[1,]
            cluster <- c(cluster ,rownames(x)[init.pair])
          }else{
            continue <- F
          }
        }else if(is.vector(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])){
          if(method=="average"){
            if(mean(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])>=threshold){cluster <- rownames(x)}
          }else if(method=="min"){
            if(min(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])>=threshold){cluster <- rownames(x)}
          }else if(method=="max"){
            if(max(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])>=threshold){cluster <- rownames(x)}
          }
          continue <- F
        }else{
          if(method=="average"){
            indix <- which.max(matrixStats::rowMeans(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)]))
            best.pathway <- names(indix)
            distance <- matrixStats::rowMeans(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])[indix]
          }else if(method=="min"){
            indix <- which.max(apply(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)], 1, FUN = min))
            best.pathway <- names(indix)
            distance <- matrixStats::rowMins(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])[indix]
          }else if(method=="max"){
            indix <- which.max(apply(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)], 1, FUN = max))
            best.pathway <- names(indix)
            distance <- matrixStats::rowMaxs(x[-which(rownames(x)%in%cluster),which(rownames(x)%in%cluster)])[indix]
          }

          if(distance>=threshold){
            cluster <- c(cluster,best.pathway)
          }else{
            continue <- F
          }
        }
      }
      pathway.to.remove <- which(rownames(x)%in% cluster)
      x <- x[-pathway.to.remove,-pathway.to.remove]
      r <- ifelse(is.matrix(x) ,nrow(x),0 )
      if(length(cluster)>0){clusters <- append(clusters,list(cluster))}

    }
  }
  if(include.singleton==T){
    if(is.null(p.singleton)){
      if(length(clusters)==0){
        singletons <- rownames(w)
      }else{
        singletons <- rownames(w)[rowSums(sapply(clusters, function(cluster){rownames(w)%in%cluster}))==0]
      }
    }else{
      if(length(clusters)==0){
        singletons <- rownames(w)[p.singleton >= threshold]
      }else{
        singletons <- rownames(w)[rowSums(sapply(clusters, function(cluster){rownames(w)%in%cluster}))==0]
        singletons <- singletons[singletons%in% rownames(w)[p.singleton >= threshold]]
      }
    }
    if(length(singletons)>0){
      for (singleton in singletons) {
        clusters <- append(clusters,singleton)
      }
    }
  }
  return(clusters)
}

#' Calculates the coverage rate of estimated clusters
#' @param clusters.true A list of vectors, each vector contains the elements in the true cluster.
#' @param clusters.result A list of vectors, each vector contains the elements in the estimated cluster.
#' @return The coverage rate of estimated clusters.
#' @export
coverageRate <- function(clusters.true,clusters.result){
  return( sum(sapply(clusters.true, function(G.j){ max( sapply(clusters.result, function(R.i){length(intersect(R.i,G.j))/length(G.j)}) )*length(G.j) }))/sum( sapply(clusters.true, function(G.j){ length(G.j) }) ) )
}

#' Calculates the inclusion rate of estimated clusters
#' @param clusters.true A list of vectors, each vector contains the elements in the true cluster.
#' @param clusters.result A list of vectors,  each vector contains the elements in the estimated cluster.
#' @return The inclusion rate of estimated clusters.
#' @export
inclusionRate <- function(clusters.true,clusters.result){
  return(sum(sapply(clusters.result, function(R.i){ max( sapply(clusters.true, function(G.j){length(intersect(R.i,G.j))/length(R.i)}) )*length(R.i) }))/sum( sapply(clusters.result, function(R.i){ length(R.i) }) ))
}
