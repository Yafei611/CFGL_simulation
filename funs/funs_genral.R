

cal_sample_cov <- function(dat){
  if (!is.list(dat)) stop("input should be a list")
  nt <- length(dat)
  S <- list()
  for (i in 1:nt){
    S[[i]] <- cov(dat[[i]])
  }
  return(S)
}


cal_sparsity <- function(theta){
  if (!is.list(theta)) stop("input should be a list")
  nt <- length(theta)
  sps_s <- rep(0,nt)
  for (i in 1:nt){
    temp <- theta[[i]]
    diag(temp) <- 0
    sps_s[i] <- sum(temp!=0)
  }
  sps_to <- sum(sps_s)
  return(list(sp_s=sps_s, sp_to=sps_to))
}


cal_loglik <- function(theta,S,n,mutipiler=1e250){
  l <- length(theta)
  loglik<-0
  for (k in 1:l){
    loglik <- loglik +n[k]/2 * ( (log(det(rbind(theta[[k]][1,]*mutipiler,theta[[k]][-1,]))) -log(mutipiler)) 
                                 - sum(diag(S[[k]] %*% theta[[k]])) )
  }
  return(loglik)
}

cal_AIC <- function(theta,loglik){
  l <- length(theta)
  pn <- 0
  for (k in 1:l){
    pn <- pn + sum(theta[[k]][upper.tri(theta[[k]],diag = F)]!=0)
  }
  aic <- 2*pn - 2*loglik
  return(aic)
}

cal_BIC <- function(theta,loglik,n){
  l <- length(theta)
  pn <- 0
  for (k in 1:l){
    pn <- pn + log(n[k])*sum(theta[[k]][upper.tri(theta[[k]],diag = F)]!=0)
  }
  bic <- pn - 2*loglik
  return(bic)
}


cal_BIC2 <- function(theta,loglik,n,scrmat,mindiff=1e-5){
  pn <- 0
  for (k in 1:3){
    theta[[k]] <- theta[[k]][upper.tri(theta[[k]],diag = F)]
    scrmat[[k]] <- scrmat[[k]][upper.tri(scrmat[[k]],diag = F)]
    pn <- pn + log(n[k])*sum(theta[[k]]!=0)
  }
  
  n.mg12 <- sum( (abs(theta[[1]]-theta[[2]])<mindiff) & scrmat[[1]] & (abs(theta[[1]])>0) )
  n.mg13 <- sum( (abs(theta[[1]]-theta[[3]])<mindiff) & scrmat[[2]] & (abs(theta[[1]])>0) )
  n.mg23 <- sum( (abs(theta[[2]]-theta[[3]])<mindiff) & scrmat[[3]] & (abs(theta[[2]])>0) )
  n.mg123 <- sum( (abs(theta[[2]]-theta[[3]])<mindiff) & (abs(theta[[1]]-theta[[3]])<mindiff) & scrmat[[3]] & (abs(theta[[1]])>0) )
  pn <- pn - n.mg12 - n.mg13 - n.mg23 + n.mg123
  
  bic <- pn - 2*loglik
  return(bic)
}





get_SSE  <- function(dm,tm,include.diag=F){
  sse <- list()
  for (i in 1:length(tm)){
    if (include.diag==F) {
      diag(dm[[i]]) <- 0
      diag(tm[[i]]) <- 0
    }
    sse[[i]] <- sum((dm[[i]]-tm[[i]])^2)
    names(sse)[i] <- paste("t",i,sep="")
  }
  sse$all <- sum(sse$t1+sse$t2+sse$t3)
  return(sse)
}

get_FPR_TPR  <- function(dm,tm){
  ln <- length(dm)
  tab <- matrix(0,nr<-ln+1,nc=6)
  colnames(tab) <- c("tpr","fpr","tpn","fpn","tp","tn")
  rownames(tab) <- c(paste("t",1:ln,sep=""),"all")
  tab <- as.data.frame(tab)
  for (i in 1: ln){
    diag(tm[[i]]) = 0
    diag(dm[[i]]) = 0
    p = tm[[i]]!=0; sum.p = sum(p)
    n = tm[[i]]==0; sum.n = sum(n)
    dp = dm[[i]]!=0; 
    tpn = sum(dp&p); tpr = tpn/sum.p
    fpn = sum(dp&n); fpr = fpn/sum.n
    tab[i,] <- c(tpr,fpr,tpn,fpn,sum.p,sum.n)
  }
  tab[ln+1,3:6] <- colSums(tab[1:ln,])[3:6]
  tab[ln+1,1:2] <- c(tab$tpn[ln+1]/tab$tp[ln+1],tab$fpn[ln+1]/tab$tn[ln+1]) 
  return(tab)
}

get_ROC <- function(dm,tm){
  #diag not included !!!!!!!
  tp.absv = NULL;
  fp.absv = NULL;
  ntp = 0;
  nfp = 0;
  sum.p = 0;
  sum.n = 0;
  
  for (i in 1: length(dm)){
    p = tm[[i]]!=0; diag(p) = 0; sum.p =sum.p+sum(p)
    n = tm[[i]]==0; sum.n =sum.n+sum(n)
    dp = dm[[i]]!=0; diag(dp) = 0
    tp = dp&p; tp.absv = c(tp.absv,abs(dm[[i]][tp]))
    fp = dp&n; fp.absv = c(fp.absv,abs(dm[[i]][fp]))
  }
  
  max.ntp = length(tp.absv)
  
  temp = rep(0,sum.p)
  #temp[1:max.ntp] = sort(tp.absv,decreasing = T)
  if (max.ntp>0) temp[1:max.ntp] = sort(tp.absv,decreasing = T)
  nfp.vec = apply(matrix(temp,nc=1),1,function(x) sum(fp.absv>=x))
  ntp.vec = c(0:max.ntp,rep(max.ntp,sum.p-max.ntp))[-1]
  
  tpr.vec = ntp.vec/sum.p
  fpr.vec = nfp.vec/sum.n
  return(list(tp=ntp.vec,fp=nfp.vec,tpr=tpr.vec,fpr=fpr.vec))
}


#####

get_n_olp_2t <- function(mat1,mat2,diag.include=F){
  if (!diag.include){
    diag(mat1) <- 0
    diag(mat2) <- 0
  }
  n <- sum(mat1!=0&mat2!=0)
  return(n)
}

get_n_olp2 <- function(mat1,mat2,mat3,diag.include=F){
  if (!diag.include){
    diag(mat1) <- 0
    diag(mat2) <- 0
    diag(mat3) <- 0
  }
  n <- sum(mat1!=0&mat2!=0&mat3==0)
  return(n)
}

get_n_olp3 <- function(mat1,mat2,mat3,diag.include=F){
  if (!diag.include){
    diag(mat1) <- 0
    diag(mat2) <- 0
    diag(mat3) <- 0
  }
  n <- sum(mat1!=0&mat2!=0&mat3!=0)
  return(n)
}

theta2rmat <- function(theta,top_edge=NULL,min_edge=0,keep.diag=F){
  temp0 <- NULL
  temp1 <- NULL
  temp2 <- NULL
  
  rmat <- list()
  rmat2 <- list()
  
  for (i in 1:length(theta)) {
    temp <- diag(theta[[i]])
    rmat[[i]] <- -theta[[i]]/sqrt(temp%*%t(temp))
  }
  
  if (!is.null(top_edge)){
    for (i in 1:length(rmat)) {
      temp0 <- rmat[[i]]
      diag(temp0) <- 0
      temp1 <- c(temp1,abs(as.vector(temp0)))
    }
    temp1 <- temp1[which(temp1!=0)]
    min_edge0 <- temp1[which(rank(-temp1,ties.method = "random")==top_edge)]
    if (min_edge0>min_edge) print("Min_edge was overrided")
    min_edge <- max(min_edge0,min_edge)
  }
  
  print(paste("min_edge is",min_edge))
  
  for (i in 1:length(rmat)) {
    temp2 <- rmat[[i]]
    diag(temp2) <- 0
    temp2[abs(temp2)<min_edge] <- 0
    rmat2[[i]] <- temp2
    if (keep.diag) diag(rmat2[[i]]) <- diag(rmat[[i]])
  }
  return(rmat2)
}



















































