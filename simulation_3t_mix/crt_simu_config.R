
setwd("d:/R_work/jgl/simulation_3t_mix/");
rm(list=ls(all=TRUE));

settings <- list()
settings$verbose <- F
settings$simu.times <- 20
settings$tag <- "mar17_8s"

########

tnpar <- list()
tnpar$lam1          <- NULL
tnpar$lam2          <- NULL
tnpar$penalize.diag <- T
tnpar$w.alpha       <- 0.40
tnpar$s.selected    <- NULL


########

smconfig <- list()

##
## sm1
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(100,100,100)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[1]] <- sm

##
## sm2
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(100,100,100)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[2]] <- sm

##
## sm3
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(100,100,100)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[3]] <- sm

##
## sm4
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "ss","ss","ss","ss" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(100,100,100)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[4]] <- sm



##
## sm5
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(50,50,50)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[5]] <- sm

##
## sm6
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(50,50,50)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[6]] <- sm

##
## sm7
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "dd","dd","dd","dd" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(50,50,50)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[7]] <- sm

##
## sm8
sm <- list()
sm$mod.size <- rep(50,12)
sm$mod.type <- matrix(c("nm","nm",  "ss","ss",   "ss","ss","ss","ss",   "sd","sd","sd","sd",
                        "nm","nm",  "ss","ss",   "sd","sd","sd","sd",   "ss","ss","ss","ss" ),nc=2)

sm$net.a    <- rep(1,12)
sm$net.m    <- rep(1,12)
sm$edge.rng <- c(0.2,1) 
sm$sampn    <- c(50,50,50)
sm$mag      <- 2.5*sqrt(log(sum(sm$mod.size))/(sum(sm$sampn)/2))*2

smconfig[[8]] <- sm


dir.create(paste("data/",settings$tag,sep=""))
save(file=paste("data/",settings$tag,"/smconfig.Rdata",sep=""),settings,smconfig,tnpar)


