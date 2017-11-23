## beginning state A and ending state A:
source("stationary.R")
#mat=matrix(c(-2,1,1,1,-2,1,1,1,-2),3)
#mat=matrix(c(-11,10,1,1,-2,1,1,10,-11),3,byrow=T)
#mat=matrix(c(-20,10,10,1,-2,1,1,1,-2),3,byrow=T)
#mat=matrix(c(-11,10,1,10,-11,1,1,1,-2),3,byrow=T)
mat=matrix(c(-20,10,10,10,-11,1,1,10,-11),3,byrow=T)
eigen=solve1(mat)[1]
RateM=solve1(mat)[[2]]
#RateM <- HKYRate(2,c(0.2,0.3,0.3,0.2))
BgnSt <- 1 ; EndSt <- 1
## Time interval and number of intermediate points:
#BgnTm <- 0.01 ; EndTm <- 30 ; nSim <- 100
BgnTm <- 0.01 ; EndTm <- 3 ; nSim <- 30
TmV <- seq(BgnTm,EndTm,len=nSim)
## Number of independent samples:
nSmpl <- 100
##----------------------------------------------------------
#pdf("SupplFig.pdf",width=12.0,height=12.0)
par(mfrow=c(1,3))
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
x1R <- 1/p.acc
x2R <- n.oper
alpha1R <- 1/(t(x1R)%*%x1R)*x1R%*%cpu.r.mat[,1]
alpha2R <- 1/(t(x2R)%*%x2R)*x2R%*%cpu.r.mat[,2]
yR <- rowSums(cpu.r.mat)
## Plot


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
mnInitI <- mean(cpu.d.mat[,1])
xI <- n.cnd.sbst
alphaI <- 1/(t(xI)%*%xI)*xI%*%cpu.d.mat[,2]

## Plot
yI <- rowSums(cpu.d.mat)


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
mnInitU <- mean(cpu.u.mat[,1])
xU <- n.virt.sbst
alphaU <- 1/(t(xU)%*%xU)*xU%*%cpu.u.mat[,2]
yU <- rowSums(cpu.u.mat)
## Plots
ymax=max(c(yR,yI,yU))

plot(TmV,yR,ylim=c(0,ymax+0.05),
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=1)
title("Begin state X and end state X")
points(TmV,alpha1R*x1R+alpha2R*x2R,type="l",col="red")

points(TmV,yI,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=2)
points(TmV,mnInitI+alphaI*xI,type="l",col="red")

points(TmV,yU,
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=4)
points(TmV,mnInitU+alphaU*xU,col="red",type="l")

legend("topleft",legend=c("Rejection sampling","Direct sampling","Uniformization"),pch=c(1,2,4),
       cex=0.8)
##-----------------------------------------------
#dev.off()

##-----------------------------------------------
##-----------------------------------------------
##-----------------------------------------------

## beginning state A and ending state G:
source("HKY.R")
frq <- c(0.2,0.3,0.3,0.2)
#RateM <- HKYRate(2,c(0.2,0.3,0.3,0.2))
BgnSt <- 1 ; EndSt <- 2
## Time interval and number of intermediate points:
#BgnTm <- 0.01 ; EndTm <- 3 ; nSim <- 30
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
x1R <- 1/p.acc
x2R <- n.oper
alpha1R <- 1/(t(x1R)%*%x1R)*x1R%*%cpu.r.mat[,1]
alpha2R <- 1/(t(x2R)%*%x2R)*x2R%*%cpu.r.mat[,2]
yR <- rowSums(cpu.r.mat)
## Plot


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
mnInitI <- mean(cpu.d.mat[,1])
xI <- n.cnd.sbst
alphaI <- 1/(t(xI)%*%xI)*xI%*%cpu.d.mat[,2]

## Plot
yI <- rowSums(cpu.d.mat)


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
mnInitU <- mean(cpu.u.mat[,1])
xU <- n.virt.sbst
alphaU <- 1/(t(xU)%*%xU)*xU%*%cpu.u.mat[,2]
yU <- rowSums(cpu.u.mat)
## Plots
ymax=max(c(yR,yI,yU))

plot(TmV,yR,ylim=c(0,0.05+ymax),
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=1)
title("Begin state X and end state Y")
points(TmV,alpha1R*x1R+alpha2R*x2R,type="l",col="red")

points(TmV,yI,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=2)
points(TmV,mnInitI+alphaI*xI,type="l",col="red")

points(TmV,yU,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=4)
points(TmV,mnInitU+alphaU*xU,col="red",type="l")

legend("topleft",legend=c("Rejection sampling","Direct sampling","Uniformization"),pch=c(1,2,4),
       cex=0.8)
##-----------------------------------------------
#dev.off()



BgnSt <- 1 ; EndSt <- 3
## Time interval and number of intermediate points:
#BgnTm <- 0.01 ; EndTm <- 30 ; nSim <- 100

TmV <- seq(BgnTm,EndTm,len=nSim)
## Number of independent samples:
nSmpl <- 100
##----------------------------------------------------------
#pdf("SupplFig.pdf",width=12.0,height=12.0)
par(mfrow=c(1,2,3))
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
x1R <- 1/p.acc
x2R <- n.oper
alpha1R <- 1/(t(x1R)%*%x1R)*x1R%*%cpu.r.mat[,1]
alpha2R <- 1/(t(x2R)%*%x2R)*x2R%*%cpu.r.mat[,2]
yR <- rowSums(cpu.r.mat)
## Plot


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
mnInitI <- mean(cpu.d.mat[,1])
xI <- n.cnd.sbst
alphaI <- 1/(t(xI)%*%xI)*xI%*%cpu.d.mat[,2]

## Plot
yI <- rowSums(cpu.d.mat)


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
mnInitU <- mean(cpu.u.mat[,1])
xU <- n.virt.sbst
alphaU <- 1/(t(xU)%*%xU)*xU%*%cpu.u.mat[,2]
yU <- rowSums(cpu.u.mat)
## Plots
ymax=max(c(yR,yI,yU))

plot(TmV,yR,ylim=c(0,ymax+0.05),
     xlab="Time",ylab="Total CPU",
     cex.axis=2,cex.lab=1.4,pch=1)
title("Begin state X and end state Z")
points(TmV,alpha1R*x1R+alpha2R*x2R,type="l",col="red")

points(TmV,yI,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=2)
points(TmV,mnInitI+alphaI*xI,type="l",col="red")

points(TmV,yU,
       xlab="Time",ylab="Total CPU",
       cex.axis=2,cex.lab=1.4,pch=4)
points(TmV,mnInitU+alphaU*xU,col="red",type="l")

legend("topleft",legend=c("Rejection sampling","Direct sampling","Uniformization"),pch=c(1,2,4),
       cex=0.8)
##-----------------------------------------------
#dev.off()
