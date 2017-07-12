
library("matrixcalc")
library("mvtnorm")
library("igraph")

crt_net <- function(nodes.n,net.a,net.m){
  # generate scale free net
  nt <- sample_pa(n = nodes.n,power = net.a,m = net.m,directed = FALSE)
  output <- get.adjacency(nt,names = FALSE,edges = FALSE, sparse = FALSE)
  return(output)
}

crt_mat_A0 <- function(netm,edge.rng){
  # generate A mat, net edge value from uniform
  p <- dim(netm)[1]
  non0.loc <- (netm==1)*upper.tri(netm)==1
  edge.n <- sum(sum(non0.loc))
  edge.v <- runif(n = edge.n,min = edge.rng[1],max =  edge.rng[2])
  edge.v <- edge.v*sign(runif(edge.n,-1,1))
  A0 <- matrix(0,nc=p,nr=p)
  A0[non0.loc] <- edge.v
  A0 <- (A0 + t(A0))
  diag(A0) = 1
  return(A0)
}

crt_mat_U <- function(A0,mag){
  p <- dim(A0)[1]
  non0.loc <- upper.tri(A0)*(abs(A0)>0)==1
  edge.n <- sum(sum(non0.loc))
  temp2 <- runif(n = edge.n,min = -1,max = 1)
  edge.v <- sign(temp2)*(abs(temp2)*mag+mag)
  U <- matrix(0,nc=p,nr=p)
  U[non0.loc] <- edge.v
  U <- (U + t(U))
  return(U)
}

crt_mat_A1 <- function(A0,U=0){
  # make A pd
  egv1 <- min(eigen(A0)$value)
  egv2 <- NULL
  if (is.matrix(U)) {egv2 <- min(eigen(A0+U)$value)}
  delta <- abs(min(egv1,egv2))+0.05
  A1 <- A0+U+delta*diag(dim(A0)[1])
  return(A1)
}


crt_mat_corm <- function(A1){
  # create cor matrix since diag=1
  p = dim(A1)[1]
  inv.A1 = solve(A1)
  diag.inv.A1 = diag(inv.A1)
  corm = inv.A1/sqrt(matrix(diag.inv.A1,nc=p,nr=p,byrow = T)*matrix(diag.inv.A1,nc=p,nr=p,byrow = F))
  return(corm)
}


simu_net <- function(mod.size,mod.type,net.a,net.m,edge.rng,mag){
  mi <- 1
  lb <- rep(0,sum(mod.size))
  for (i in 1:length(mod.size)){
    lb[mi:(mi+mod.size[i]-1)] <- i
    mi <- mi+mod.size[i]
  }
  
  p <- length(lb)
  t <- dim(mod.type)[2]+1
  
  matA.all <- list()
  corm.all <- list()
  pcsm.all <- list()
  
  for (i in 1:t){
    matA.all[[i]] <- matrix(0,nc=p,nr=p)
    corm.all[[i]] <- matrix(0,nc=p,nr=p)
    pcsm.all[[i]] <- matrix(0,nc=p,nr=p)
  }
  
  for (modi in 1:max(lb)){
    
    mod.area <- which(lb==modi)
    is.pd <- FALSE
    
    if (sum(mod.type[modi,]=="nm")==2) {
      temp <- diag(mod.size[modi])
      
      matA.all[[1]][mod.area,mod.area] <- temp
      corm.all[[1]][mod.area,mod.area] <- temp
      pcsm.all[[1]][mod.area,mod.area] <- temp
      
      matA.all[[2]][mod.area,mod.area] <- temp
      corm.all[[2]][mod.area,mod.area] <- temp
      pcsm.all[[2]][mod.area,mod.area] <- temp
      
      matA.all[[3]][mod.area,mod.area] <- temp
      corm.all[[3]][mod.area,mod.area] <- temp
      pcsm.all[[3]][mod.area,mod.area] <- temp
    }
    
    else{
      while (!is.pd){
        netm <- crt_net(mod.size[modi],net.a[modi],net.m[modi])
        A0 <- crt_mat_A0(netm,edge.rng)
        A1 <- crt_mat_A1(A0)
        is.pd <- is.positive.definite(A1)
        if (!is.pd) print(paste("non-pd occurs: modi =",modi))
      }
      
      corm <- crt_mat_corm(A1)
      pcsm <- solve(corm)
      pcsm[abs(pcsm)<1e-8]=0
      
      matA.all[[1]][mod.area,mod.area] <- A1
      corm.all[[1]][mod.area,mod.area] <- corm
      pcsm.all[[1]][mod.area,mod.area] <- pcsm
      
      for (ti in 2:t){
        
        is.pd <- FALSE
        
        if (mod.type[modi,ti-1]=="ss"){
          matA.all[[ti]][mod.area,mod.area] <- A1
          corm.all[[ti]][mod.area,mod.area] <- corm
          pcsm.all[[ti]][mod.area,mod.area] <- pcsm
        }
        
        if (mod.type[modi,ti-1]=="sd"){
          while (!is.pd){
            A1.2 <- crt_mat_A1(A0,U = crt_mat_U(A0,mag))
            is.pd <- is.positive.definite(A1.2)
            if (!is.pd) print(paste("non-pd occurs: modi =",modi))
          }
          corm2 <- crt_mat_corm(A1.2)
          pcsm2 <- solve(corm2)
          pcsm2[abs(pcsm2)<1e-8]=0
          
          matA.all[[ti]][mod.area,mod.area] <- A1.2
          corm.all[[ti]][mod.area,mod.area] <- corm2
          pcsm.all[[ti]][mod.area,mod.area] <- pcsm2
          rm("A1.2","corm2","pcsm2")
        }
        
        if (mod.type[modi,ti-1]=="dd"){
          while (!is.pd){
            netm2 <- crt_net(mod.size[modi],net.a[modi],net.m[modi])
            A0.2 <- crt_mat_A0(netm2,edge.rng)
            A1.2 <- crt_mat_A1(A0.2)
            is.pd <-is.positive.definite(A1.2)
            if (!is.pd) print(paste("non-pd occurs: modi =",modi))
          }
          
          corm2 <- crt_mat_corm(A1.2)
          pcsm2 <- solve(corm2)
          pcsm2[abs(pcsm2)<1e-8]=0
          
          matA.all[[ti]][mod.area,mod.area] <- A1.2
          corm.all[[ti]][mod.area,mod.area] <- corm2
          pcsm.all[[ti]][mod.area,mod.area] <- pcsm2
          rm("A0.2","A1.2","corm2","pcsm2")
        }
        
        if (mod.type[modi,ti-1]=="10"){
          matA.all[[ti]][mod.area,mod.area] <- diag(dim(netm)[1])
          corm.all[[ti]][mod.area,mod.area] <- diag(dim(netm)[1])
          pcsm.all[[ti]][mod.area,mod.area] <- diag(dim(netm)[1])
        }
      }
    }
  }
  return(list(mod.lb = lb, pcsm=pcsm.all,corm=corm.all,matA=matA.all))
}


simu_expr <- function(covm,sampn){
  expr=list()
  for (i in 1:length(covm)){
    expr[[i]] <- rmvnorm(n = sampn[i],sigma=covm[[i]])
    names(expr)[i] <- paste("t",i,sep="")
  }
  return(expr)
}





