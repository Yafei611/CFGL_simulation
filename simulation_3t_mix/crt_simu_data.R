
setwd("d:/R_work/jgl/simulation_3t_mix/");
rm(list=ls(all=TRUE));
source("../funs/funs_simu_pcsm.R")
source("../funs/funs_genral.R")
source("../funs/funs_diff_pcsm.R")
source("../funs/funs_JGLS.R")

load("data/mar17_8s/smconfig.Rdata")

rdseed <- 5
scr.mat.diffcut <- 0.01

l <- length(smconfig)

for (i in 1:l){
  settings$smconfig.id <- i
  sm <- smconfig[[settings$smconfig.id]]
  tnet <- list()
  simu_data <- list()
  scrmats.true <- list()
  
  for (k in 1:settings$simu.times){
    set.seed(k*5)
    t.net <- simu_net(sm$mod.size,sm$mod.type,sm$net.a,sm$net.m,sm$edge.rng,sm$mag)
    
    scrmats.true[[k]] <- list(abs(t.net$pcsm[[1]]-t.net$pcsm[[2]])<scr.mat.diffcut,
                              abs(t.net$pcsm[[1]]-t.net$pcsm[[3]])<scr.mat.diffcut,
                              abs(t.net$pcsm[[2]]-t.net$pcsm[[3]])<scr.mat.diffcut)
                              
    expr <- simu_expr(covm = t.net$corm,sampn = sm$sampn)
    expr <- centering_data(expr)
    simu_data[[k]] <- expr
    tnet[[k]] <- t.net
    print(paste("S",i,k))
  }
  save(file=paste("data/mar17_8s/simu_data/simu_data_",i,".Rdata",sep=""),sm, tnpar, settings, tnet, simu_data, scrmats.true)
}

