## This R script reproduces Figure 5 in
## Hobolth and Stone:
## Efficient simulation from finite-state, continuous time
## Markov chains with incomplete observations
##----------------------------------------------------------
## GY acceptance probabilities
##----------------------------------------------------------
## Example 2: GY rate matrix
source("Rejection.R")
source("GY.R")
#pdf("GY-Accp.pdf",width=10.0,height=4.0)
par(mfrow=c(1,2))
##----------------------------------------------------------
## Case a: EndSt=AAG
##----------------------------------------------------------
#RateM <- GYRate(2,0.01,codon.freq)
source("HKY.R")
RateM <- HKYRate(2,c(0.2,0.3,0.3,0.2))
BgnSt <- 1 ; EndSt <- 1  ## 2=AAC, 3=AAG, 61=TTT
BgnTm <- 0.01 ; EndTm <- 4 ; nSim <- 100
TmV <- seq(BgnTm,EndTm,len=nSim)
nSmpl <- 250
p.acc <- AccPrbV(BgnSt,EndSt,RateM,TmV)
## Plot
plot(TmV,p.acc,type="l",
     xlab="Time",ylab="Acceptance probability",
     cex.axis=2,cex.lab=1.4,ylim=c(0,1)) #,log="y")
lines(TmV,AccPrbAppxV(BgnSt,EndSt,RateM,TmV),lty="dashed")
lines(TmV,rep(0.2,nSim),lty="dashed")
title("Begin state A and end state A")
##----------------------------------------------------------
## Case b: EndSt=AAC
##----------------------------------------------------------
RateM <- HKYRate(2,c(0.2,0.3,0.3,0.2))
BgnSt <- 1 ; EndSt <- 2  ## 2=AAC, 3=AAG, 61=TTT
p.acc <- AccPrbV(BgnSt,EndSt,RateM,TmV)
## Plot
plot(TmV,p.acc,type="l",
     xlab="Time",ylab="Acceptance probability",
     cex.axis=2,cex.lab=1.4,ylim=c(0,1)) #,log="y")
lines(TmV,AccPrbAppxV(BgnSt,EndSt,RateM,TmV),lty="dashed")
lines(TmV,rep(0.3,nSim),lty="dashed")
title("Begin state A and end state G")

