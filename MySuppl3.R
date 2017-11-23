## beginning state A and ending state A:
source("HKY.R")
frq <- c(0.2,0.3,0.3,0.2)
RateM <- matrix(c(-1,0.6,0.2,0.2,
                  0.6,-1,0.2,0.2,
                  6,6,-20,8,
                  0.3,0.3,0.4,-1),4,byrow=T)
BgnSt <- 4 ; EndSt <- 3
## Time interval and number of intermediate points:
BgnTm <- 0.01 ; EndTm <- 3 ; nSim <- 30
TmV <- seq(BgnTm,EndTm,len=nSim)
## Number of independent samples:
nSmpl <- 100
##----------------------------------------------------------
#pdf("SupplFig.pdf",width=12.0,height=12.0)
par(mfrow=c(1,2))
##----------------------------------------------------------
## Firstly consider REJECTION sampling  
## 
cat("Rejection sampling","\n")
source("Rejection.R")
## Acceptance probability
p.acc <- AccPrbV(BgnSt,EndSt,RateM,TmV)
## Number of expected iterations per sample
n.sbst <- NsbstV(BgnSt,EndSt,RateM,TmV)
## Total number of expected iterations 
n.oper <- n.sbst/p.acc
## Times for rejection sampling
cpu.r.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitR","SmplR")
for (sim in 1:nSim){
  cat(sim)
  cpu.r <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.r <- cpu.r+RejectSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.r.mat[sim,] <- cpu.r
}
cat("\n")
## Linear regression coefficient
x1 <- 1/p.acc
x2 <- n.oper
alpha1 <- 1/(t(x1)%*%x1)*x1%*%cpu.r.mat[,1]
alpha2 <- 1/(t(x2)%*%x2)*x2%*%cpu.r.mat[,2]
yR <- rowSums(cpu.r.mat)
## Plot
plot(TmV,yR,ylim=c(0,0.7),
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=1)
title("Begin state T and end state C")
points(TmV,alpha1*x1+alpha2*x2,type="l",col="red")

##----------------------------------------------------------
## Secondly consider DIRECT sampling  
## (corresponding to last two plots in Fig6 and Fig1, right plot)

cat("Direct sampling","\n")
source("Direct.R")

## Times for direct sampling
cpu.d.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitD","SmplD")
for (sim in 1:nSim){
  cat(sim)
  cpu.d <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.d <- cpu.d+DirectSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.d.mat[sim,] <- cpu.d
}
cat("\n")
## Number of expected iterations
n.cnd.sbst <- cndNsbstV(BgnSt,EndSt,RateM,TmV)
## Linear regression coefficient
mnInit <- mean(cpu.d.mat[,1])
x <- n.cnd.sbst
alpha <- 1/(t(x)%*%x)*x%*%cpu.d.mat[,2]

## Plot
yI <- rowSums(cpu.d.mat)
points(TmV,yI,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=2)
points(TmV,mnInit+alpha*x,type="l",col="red")

##----------------------------------------------------------
## Thirdly consider UNIFORMIZATION  
## (corresponding to last two plots in Fig7 and Fig1, right plot)
cat("Uniformization","\n")
source("Uniform.R")

## Times for uniform sampling
cpu.u.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitU","SmplU")
for (sim in 1:nSim){
  cat(sim)
  cpu.u <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.u <- cpu.u+UniformSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.u.mat[sim,] <- cpu.u
}
cat("\n")
## Number of expected iterations
n.virt.sbst <- vNsbstV(BgnSt,EndSt,RateM,TmV)
## Linear regression coefficient
mnInit <- mean(cpu.u.mat[,1])
x <- n.virt.sbst
alpha <- 1/(t(x)%*%x)*x%*%cpu.u.mat[,2]
yU <- rowSums(cpu.u.mat)
## Plots
points(TmV,yU,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=4)
points(TmV,mnInit+alpha*x,col="red",type="l")

legend("topleft",legend=c("Rejection sampling","Direct sampling","Uniformization"),pch=c(1,2,4))
##-----------------------------------------------
#dev.off()

##-----------------------------------------------
##-----------------------------------------------
##-----------------------------------------------

## beginning state A and ending state G:
source("HKY.R")
frq <- c(0.2,0.3,0.3,0.2)
RateM <- HKYRate(2,c(0.2,0.3,0.3,0.2))
BgnSt <- 3 ; EndSt <- 4
## Time interval and number of intermediate points:
BgnTm <- 0.01 ; EndTm <- 3 ; nSim <- 30
TmV <- seq(BgnTm,EndTm,len=nSim)
## Number of independent samples:
nSmpl <- 100
##----------------------------------------------------------
#pdf("SupplFig.pdf",width=12.0,height=12.0)
##----------------------------------------------------------
## Firstly consider REJECTION sampling  
## 
cat("Rejection sampling","\n")
source("Rejection.R")
## Acceptance probability
p.acc <- AccPrbV(BgnSt,EndSt,RateM,TmV)
## Number of expected iterations per sample
n.sbst <- NsbstV(BgnSt,EndSt,RateM,TmV)
## Total number of expected iterations 
n.oper <- n.sbst/p.acc
## Times for rejection sampling
cpu.r.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitR","SmplR")
for (sim in 1:nSim){
  cat(sim)
  cpu.r <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.r <- cpu.r+RejectSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.r.mat[sim,] <- cpu.r
}
cat("\n")
## Linear regression coefficient
x1 <- 1/p.acc
x2 <- n.oper
alpha1 <- 1/(t(x1)%*%x1)*x1%*%cpu.r.mat[,1]
alpha2 <- 1/(t(x2)%*%x2)*x2%*%cpu.r.mat[,2]
yR <- rowSums(cpu.r.mat)
## Plot
plot(TmV,yR,ylim=c(0,0.1),
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=1)
title("Begin state C and end state T")
points(TmV,alpha1*x1+alpha2*x2,type="l",col="red")

##----------------------------------------------------------
## Secondly consider DIRECT sampling  
## (corresponding to last two plots in Fig6 and Fig1, right plot)

cat("Direct sampling","\n")
source("Direct.R")

## Times for direct sampling
cpu.d.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitD","SmplD")
for (sim in 1:nSim){
  cat(sim)
  cpu.d <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.d <- cpu.d+DirectSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.d.mat[sim,] <- cpu.d
}
cat("\n")
## Number of expected iterations
n.cnd.sbst <- cndNsbstV(BgnSt,EndSt,RateM,TmV)
## Linear regression coefficient
mnInit <- mean(cpu.d.mat[,1])
x <- n.cnd.sbst
alpha <- 1/(t(x)%*%x)*x%*%cpu.d.mat[,2]

## Plot
yI <- rowSums(cpu.d.mat)
points(TmV,yI,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=2)
points(TmV,mnInit+alpha*x,type="l",col="red")

##----------------------------------------------------------
## Thirdly consider UNIFORMIZATION  
## (corresponding to last two plots in Fig7 and Fig1, right plot)
cat("Uniformization","\n")
source("Uniform.R")

## Times for uniform sampling
cpu.u.mat <- matrix(0,nrow=nSim,ncol=2)
colnames(cpu.r.mat) <- c("InitU","SmplU")
for (sim in 1:nSim){
  cat(sim)
  cpu.u <- c(0,0)
  for (smpl in 1:nSmpl){
    cpu.u <- cpu.u+UniformSampl(BgnSt,EndSt,RateM,TmV[sim])$ptm
  }
  cpu.u.mat[sim,] <- cpu.u
}
cat("\n")
## Number of expected iterations
n.virt.sbst <- vNsbstV(BgnSt,EndSt,RateM,TmV)
## Linear regression coefficient
mnInit <- mean(cpu.u.mat[,1])
x <- n.virt.sbst
alpha <- 1/(t(x)%*%x)*x%*%cpu.u.mat[,2]
yU <- rowSums(cpu.u.mat)
## Plots
points(TmV,yU,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=4)
points(TmV,mnInit+alpha*x,col="red",type="l")

legend("topleft",legend=c("Rejection sampling","Direct sampling","Uniformization"),pch=c(1,2,4))
##-----------------------------------------------
#dev.off()
