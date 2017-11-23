##--------------------------------------------------------------
##--------------------------------------------------------------
## i) REJECTION SAMPLING
## ii) ACCEPTANCE PROBABILITY (and a first-order approximation)
## iii) EXPECTED NUMBER OF JUMPS
##--------------------------------------------------------------
## i) REJECTION SAMPLING
##--------------------------------------------------------------
## Name: RejectSampl
## Date: July 13, 2007
## Author: Asger Hobolth
## Input: BgnSt: Begin state
##        EndSt: End state
##        RateM: Rate matrix
##        Tm: Time
## Output: List consisting of
##         St: States
##         Tm: Times
##         (including states BgnSt and EndSt at times 0 and Tm)
##         ptm: Process time
##         (divided into time for initialization and time in while loop)  
##--------------------------------------------------------------
## SimpleFrwdSampl: Simple forward sampling
SimpleFrwdSampl <- function(BgnSt,RateM,Tm){
  ptm1 <- proc.time()[1]
  nSt <- nrow(RateM)  ## Size of state space
  StSp <- 1:nSt       ## State space
  ## Choose state at time=0 from initial distribution
  X <- vector(mode="numeric")  ## X holds states of path
  X[1]<- BgnSt                 ## beginning state is BgnSt
  T <- vector(mode="numeric")  ## T holds times of path
  T[1] <- 0                    ## beginning time is 0
  ## Simulate waiting time before jump from exponential distribution
  T[2] <- rexp(1,-RateM[X[1],X[1]])
  X[2] <- sample(StSp[-X[1]],1,prob=RateM[X[1],-X[1]])
  ## Simulate states and waiting times similarly as above
  ## until time is larger than Tm
  cnt <- 2 ## counter
  ptm2 <- proc.time()[1]
  while (T[cnt] < Tm) {
    T[cnt+1] <- T[cnt]+rexp(1,-RateM[X[cnt],X[cnt]])
    X[cnt+1] <- sample(StSp[-X[cnt]],1,prob=RateM[X[cnt],-X[cnt]])
    cnt <- cnt+1
  }
  ptm3 <- proc.time()[1]
  ## Output state changes and corresponding times
  Path <- list()
  Path$Tm <- c(T[1:(cnt-1)],Tm)
  Path$St <- c(X[1:(cnt-1)],X[cnt-1])
  Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
  return(Path)
}
##--------------------------------------------------------------
## ConstrFrwdSampl: Constrained forward sampling
## (Constraint: At least one jump)
ConstrFrwdSampl <- function(BgnSt,RateM,Tm){
  ptm1 <- proc.time()[1]
  nSt <- nrow(RateM)  ## Size of state space
  StSp <- 1:nSt       ## State space
  X <- vector(mode="numeric")  ## X holds states of path
  X[1] <- BgnSt                ## beginning state is BgnSt
  T <- vector(mode="numeric")  ## T holds times of path
  T[1] <- 0                    ## beginning time is 0
  ## Simulate waiting time before first jump from
  ## truncated exponential distribution
  RateAwayBgn <- -RateM[BgnSt,BgnSt]
  T[2] <- -log(1-runif(1)*(1-exp(-Tm*RateAwayBgn)))/RateAwayBgn
  X[2] <- sample(StSp[-BgnSt],1,prob=RateM[BgnSt,-BgnSt])
  ## Simulate states and waiting times
  ## until time is larger than Tm
  cnt <- 2  ## counter
  ptm2 <- proc.time()[1]
  while (T[cnt] < Tm) {
    T[cnt+1] <- T[cnt]+rexp(1,-RateM[X[cnt],X[cnt]])
    X[cnt+1] <- sample(StSp[-X[cnt]],1,prob=RateM[X[cnt],-X[cnt]])
    cnt <- cnt+1
  }
  ptm3 <- proc.time()[1]
  ## Output state changes and corresponding times
  Path <- list()
  Path$Tm <- c(T[1:(cnt-1)],Tm)
  Path$St <- c(X[1:(cnt-1)],X[cnt-1])
  Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
  return(Path)
}
##--------------------------------------------------------------
## Output: ptm:
##        Total time spent on initialization and in while loop
RejectSampl <- function(BgnSt,EndSt,RateM,Tm){
  if (BgnSt==EndSt){
    p.acc <- AccPrbSame(BgnSt,RateM,Tm)
    if (p.acc < 1e-3) {
      SimPath <- list()
      SimPath$Tm <- c(0,NA,Tm)
      SimPath$St <- c(BgnSt,NA,EndSt)
      SimPath$ptm <- c(NA,NA)
      return(SimPath)
    }
    SimPath <- SimpleFrwdSampl(BgnSt,RateM,Tm)
    ptm <- SimPath$ptm
    while (SimPath$St[length(SimPath$St)]!=EndSt){
      SimPath <- SimpleFrwdSampl(BgnSt,RateM,Tm)
      ptm <- ptm+SimPath$ptm
    }
    SimPath$ptm <- ptm
    return(SimPath)
  }
  if (BgnSt!=EndSt){
    p.acc <- AccPrbDiff(BgnSt,EndSt,RateM,Tm)
    if (p.acc < 1e-3) {
      SimPath <- list()
      SimPath$Tm <- c(0,NA,Tm)
      SimPath$St <- c(BgnSt,NA,EndSt)
      SimPath$ptm <- c(NA,NA)
      return(SimPath)
    }
    SimPath <- ConstrFrwdSampl(BgnSt,RateM,Tm)
    ptm <- SimPath$ptm
    while (SimPath$St[length(SimPath$St)]!=EndSt){
      SimPath <- ConstrFrwdSampl(BgnSt,RateM,Tm)
      ptm <- ptm+SimPath$ptm
    }
    SimPath$ptm <- ptm
    return(SimPath)
  }
}
##-----------------------------------------------------------
## ii) ACCEPTANCE PROBABILITY
##-----------------------------------------------------------
## Acceptance probability when BgnSt and EndSt are the same
##-----------------------------------------------------------
AccPrbSame <- function(BgnSt,RateM,Tm){
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  PrbM <- U%*%diag(exp(Tm*Lam))%*%InvU
  return(PrbM[BgnSt,BgnSt])
}
AccPrbSameV <- function(BgnSt,RateM,TmV){
  len <- length(TmV)
  AccPrbV <- rep(0,len)
  for (i in 1:len) AccPrbV[i] <- AccPrbSame(BgnSt,RateM,TmV[i])
  return(AccPrbV)
}
##------------------------------------------------------------
## Acceptance probability when BgnSt and EndSt are different
##------------------------------------------------------------
AccPrbDiff <- function(BgnSt,EndSt,RateM,Tm){
  nSt <- nrow(RateM)
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  ## Numbers needed to determine if a first jump is made
  PrbM <- U%*%diag(exp(Tm*Lam))%*%InvU
  BgnRate <- -RateM[BgnSt,BgnSt]
  PrbBgnBgn <- PrbM[BgnSt,BgnSt]
  ## The q_a+lam_j=0 case only occurs under artificial circumstances
  ## and is neglected
  if (BgnRate %in% Lam) cat("A j exists such that q_a+lam_j=0","\n")
  JVec <- (exp(Tm*Lam)-exp(-Tm*BgnRate))/(Lam+BgnRate)
  KVec <-
  rowSums(U*matrix(rep(InvU[,EndSt]*JVec,nSt),nrow=nSt,byrow=TRUE))
  AccPrb <- (1/(1-exp(-Tm*BgnRate)))*sum(RateM[BgnSt,-BgnSt]*KVec[-BgnSt])  
  return(AccPrb)
}
AccPrbDiffV <- function(BgnSt,EndSt,RateM,TmV){
  len <- length(TmV)
  AccPrbV <- rep(0,len)
  for (i in 1:len) AccPrbV[i] <- AccPrbDiff(BgnSt,EndSt,RateM,TmV[i])
  return(AccPrbV)
}
##------------------------------------------------------------
## Acceptance probability 
##------------------------------------------------------------
AccPrbV <- function(BgnSt,EndSt,RateM,TmV){
  if (BgnSt==EndSt) return(AccPrbSameV(BgnSt,RateM,TmV))
  if (BgnSt!=EndSt) return(AccPrbDiffV(BgnSt,EndSt,RateM,TmV))
}
##------------------------------------------------------------
## Acceptance probability approximation
##------------------------------------------------------------
AccPrbAppxV <- function(BgnSt,EndSt,RateM,TmV){
  if (BgnSt!=EndSt){
    nSt <- nrow(RateM)
    StSp <- 1:nSt
    Qa <- -RateM[BgnSt,BgnSt]
    Qab <- RateM[BgnSt,EndSt]
    Qb <- -RateM[EndSt,EndSt]
    nonab <- StSp[-c(BgnSt,EndSt)]
    intcpt <- Qab/Qa
    slope1 <- sum(RateM[BgnSt,nonab]*RateM[nonab,EndSt])/2/Qa
    slope2 <- -Qab*Qb/2/Qa
    return(intcpt+(slope1+slope2)*TmV)
  }
  if (BgnSt==EndSt){
    Qa <- -RateM[BgnSt,BgnSt]
    return(1-Qa*TmV)
  }
}
##------------------------------------------------------------
## iii) EXPECTED NUMBER OF JUMPS
##------------------------------------------------------------
## Expected number of jumps 
## conditional on beginning state (note: NOT!! ending state)
## and used when BgnSt and EndSt are the same
##------------------------------------------------------------
NsbstSame <- function(BgnSt,RateM,Tm){
  nSt <- nrow(RateM)
  StSp <- 1:nSt
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  ## Integrals
  IVec <- rep(0,nSt)
  for (j in 1:nSt){
    IVec[j] <- ifelse(identical(all.equal(Lam[j],0),TRUE),
                      Tm,
                      (exp(Lam[j]*Tm)-1)/Lam[j])
  }
  ## Time spent in a state
  TmVec <- rep(0,nSt)
  for (i in 1:nSt){
    TmVec[i] <- sum(U[BgnSt,]*InvU[,i]*IVec)
  }
  ## Expected number of jumps
  N <- 0
  for (i in 1:nSt){
    for (j in StSp[-i]){
      N <- N+RateM[i,j]*TmVec[i]
    }
  }
  return(N)
}
NsbstSameV <- function(BgnSt,RateM,TmV){
  len <- length(TmV)
  NV <- rep(0,len)
  for (i in 1:len) NV[i] <- NsbstSame(BgnSt,RateM,TmV[i])
  return(NV)
}
##------------------------------------------------------------
## Expected number of jumps 
## conditional on beginning state (note: NOT!! ending state)
## and used when BgnSt and EndSt different
##------------------------------------------------------------
## Time spent in states after first jump [entries 1:nSt] and in
## BgnSt before jump [entry (nSt+1)]
RejTime <- function(BgnSt,RateM,Tm,NewSt){
  nSt <- nrow(RateM)
  StSp <- 1:nSt
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  ## Rate matrix entry 
  QBgn <- -RateM[BgnSt,BgnSt]
  ##----------------------------------
  ## Time spent in BgnSt before jump
  ##----------------------------------
  TmBgn <- 1/(1-exp(-Tm*QBgn))*
    (-exp(-Tm*QBgn)/QBgn-Tm*exp(-Tm*QBgn)+1/QBgn)
  ##----------------------------------
  ## Time spent in states after jump
  ##----------------------------------
  ## Integrals
  JVec <- rep(0,nSt)
  for (j in 1:nSt){
    if (all.equal(Lam[j],0)==TRUE){
      #cat("Lam=0 for j=",j,"\n")
      JVec[j] <-
        (exp(-Tm*QBgn)/QBgn)*(Tm*exp(Tm*QBgn)-exp(Tm*QBgn)/QBgn+1/QBgn)
      }
    if (all.equal(Lam[j],QBgn)==TRUE){
      #cat("Lam=Qj for j=",j,"\n")
      JVec[j] <- (1/Lam[j])*(
        (exp(Tm*Lam[j])/QBgn)*(1-exp(-Tm*QBgn))-Tm )
    }
    if (all.equal(Lam[j],0)!=TRUE & all.equal(Lam[j],QBgn)!=TRUE){
      #cat("Lam!=Qj for j=",j,"\n")
      JVec[j] <- (1/Lam[j])*
        (exp(Tm*Lam[j])/QBgn*(1-exp(-Tm*QBgn))-
        (1/(Lam[j]-QBgn))*(exp((Lam[j]-QBgn)*Tm)-1))
    }
  }
  TmInSt <- QBgn/(1-exp(-Tm*QBgn))*colSums(U[NewSt,]*InvU*JVec)
  return(c(TmInSt,TmBgn))
}
## Total number of substitutions
NsbstDiff <- function(BgnSt,RateM,Tm){ 
  nSt <- nrow(RateM)
  StSp <- 1:nSt
  nonBgnSt <- StSp[-BgnSt]
  QBgn <- -RateM[BgnSt,BgnSt]
  TmInSt <- rep(0,nSt+1)
  for (NewSt in nonBgnSt){
    TmInSt <- TmInSt+
      RejTime(BgnSt,RateM,Tm,NewSt)*RateM[BgnSt,NewSt]/QBgn
  }
  nExpSbst <- 1-sum(TmInSt[StSp]*diag(RateM))
  return(nExpSbst)
}
NsbstDiffV <- function(BgnSt,RateM,TmV){
  len <- length(TmV)
  NV <- rep(0,len)
  for (i in 1:len) {
    cat("evaluating no",i,"out of",len,"\n")
    NV[i] <- NsbstDiff(BgnSt,RateM,TmV[i])
  }
  return(NV)
}
##------------------------------------------------------------------
## Expected number of sbst
##------------------------------------------------------------------
NsbstV <- function(BgnSt,EndSt,RateM,TmV){
  if (BgnSt==EndSt) return(NsbstSameV(BgnSt,RateM,TmV))
  if (BgnSt!=EndSt) return(NsbstDiffV(BgnSt,RateM,TmV))
}
