##---------------------------------------------------------------
##---------------------------------------------------------------
## i) UNIFORM SAMPLING
## ii) EXPECTED NUMBER OF JUMPS -- INCLUDING VIRTUAL
##---------------------------------------------------------------
## i) UNIFORM SAMPLING
##---------------------------------------------------------------
## Name: UniformSampl
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
##         (divided into time for initialization and time for sampling)
##---------------------------------------------------------------
UniformSampl <- function(BgnSt,EndSt,RateM,Tm){
  ptm1 <- proc.time()[1]
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  PrbBgnEnd <- (U%*%diag(exp(Tm*Lam))%*%InvU)[BgnSt,EndSt]
  ## Determine max diagonal entry and construct transition matrix
  nSt <- nrow(RateM)
  Mx <- max(-diag(RateM))
  TransM <- diag(rep(1,times=nSt))+RateM/Mx
  TransMn <- TransM
  ## TransMn is n'th power of TransM
  ##------------------------------------------------------------
  ## Simulate number of jumps
  rU <- runif(1)
  ifelse(BgnSt==EndSt,cum <- dpois(0,Mx*Tm)/PrbBgnEnd,cum <- 0)
  notExceed <- TRUE
  if (cum>rU) { notExceed <- FALSE }
  nJmp <- 0
  ptm2 <- proc.time()[1]
  while (notExceed){
    nJmp <- nJmp+1 
    prb <- dpois(nJmp,Mx*Tm)*TransMn[BgnSt,EndSt]/PrbBgnEnd
    cum <- cum+prb
    if (cum>rU) notExceed <- FALSE
    ## Update transition matrices
    ## TransArr holds n'th power of TransM
    if (nJmp==1) TransArr <- array(TransM,c(nSt,nSt,1))
    if (nJmp!=1) TransArr <- array(c(TransArr,TransMn),c(nSt,nSt,nJmp))
    TransMn <- TransMn%*%TransM
  }
  #cat("nJmp:",nJmp,"\n")
  ptm3 <- proc.time()[1]
  ##--------------------------------------------------------
  ## if (nJmp==0): done
  if (nJmp==0){
    Path <- list()
    Path$St <- c(BgnSt,EndSt)
    Path$Tm <- c(0,Tm)
    Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
    return(Path)
  }
  ## if (nJmp==1)
  if (nJmp==1){
    ## if virtual jump: done
    if (BgnSt==EndSt){
      Path <- list()
      Path$St <- c(BgnSt,EndSt)
      Path$Tm <- c(0,Tm)
      Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
      return(Path)
    }
    ## if true jump: done
    if (BgnSt!=EndSt){
      Path <- list()
      Path$St <- c(BgnSt,EndSt,EndSt)
      Path$Tm <- c(0,Tm*runif(1),Tm)
      Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
      return(Path)
    }
  }
  ## Case (nJmp >= 2):
  ## Simulate jumping times
  JmpTmV <- Tm*sort(runif(nJmp))
  ## Simulate states (last state always EndSt)
  JmpStV <- rep(0,(nJmp-1))
  Prb1 <- TransM[BgnSt,]
  for (i in 1:(nJmp-1)){
    Prb2Mat <- TransArr[,,(nJmp-i)]
    Prb2 <- Prb2Mat[,EndSt]
    JmpStV[i] <-
      sample(1:nSt,size=1,replace=TRUE,round(Prb1*Prb2,digits=10))
    Prb1 <- TransM[JmpStV[i],]
  }
  ptm3 <- proc.time()[1]
  JmpStV <- c(JmpStV,EndSt)
  ## Remove virtual substitutions
  TrueSub <- c(BgnSt,JmpStV[1:nJmp])!=c(JmpStV[1:nJmp],EndSt)
  Path <- list()
  Path$St <- c(BgnSt,JmpStV[TrueSub],EndSt)
  Path$Tm <- c(    0,JmpTmV[TrueSub],Tm)
  Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
  return(Path)
}
##---------------------------------------------------------------
## ii) EXPECTED NUMBER OF JUMPS -- INCLUDING VIRTUAL
##---------------------------------------------------------------
## Expected number of jumps -- including virtual
##---------------------------------------------------------------
vNsbst <- function(BgnSt,EndSt,RateM,Tm){
  nSt <- nrow(RateM)
  StSp <- 1:nSt
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  ## Probability matrix
  PrbM <- U%*%diag(exp(Tm*Lam))%*%InvU
  PrbBE <- PrbM[BgnSt,EndSt]
  Mx <- max(-diag(RateM))
  TransM <- diag(rep(1,times=nSt))+RateM/Mx
  vN <- Mx*Tm*((TransM%*%PrbM)[BgnSt,EndSt])/PrbBE
  return(vN)
}
vNsbstV <- function(BgnSt,EndSt,RateM,TmV){
  len <- length(TmV)
  NV <- rep(0,len)
  for (i in 1:len) NV[i] <- vNsbst(BgnSt,EndSt,RateM,TmV[i])
  return(NV)
}
