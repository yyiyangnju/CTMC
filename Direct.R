##---------------------------------------------------------------
##---------------------------------------------------------------
## i) DIRECT SAMPLING
## ii) EXPECTED NUMBER OF JUMPS
##     CONDITIONAL ON BEGINNING AND ENDING STATE
##---------------------------------------------------------------
## i) DIRECT SAMPLING
##---------------------------------------------------------------
## Name: DirectSampl
## Date: July 6, 2007
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
##---------------------------------------------------------------
DirectSampl <- function(BgnSt,EndSt,RateM,Tm){
  ptm1 <- proc.time()[1]
  nSt <- nrow(RateM) ## Size of state space
  StSp <- 1:nSt      ## State space
  X <- vector(mode="numeric") ## X holds states of path
  X[1] <- BgnSt               ## beginning state is BgnSt
  T <- vector(mode="numeric") ## T holds times of path
  T[1] <- 0                   ## beginning time is 0
  EndTm <- Tm
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  ## Numbers needed to determine if a first jump is made
  PrbM <- U%*%diag(exp(Tm*Lam))%*%InvU
  BgnRate <- -RateM[BgnSt,BgnSt]
  PrbBgnBgn <- PrbM[BgnSt,BgnSt]
  ## Determine if a jump is made and update accordingly.
  cnt <- 1  ## count number of jumps
  ptm2 <- proc.time()[1]
  while (BgnSt!=EndSt | runif(1)>exp(-BgnRate*Tm)/PrbBgnBgn) {
    ## When calculating integrals below I assume q_a+lam_j is non-zero.
    ## The q_a+lam_j=0 case only occurs under artificial circumstances
    ## and is neglected
    if (BgnRate %in% Lam) cat("A j exists such that q_a+lam_j=0","\n")  
    ##--------------------------------------
    ## Simulate new state
    ##--------------------------------------
    PrbBgnEnd <- PrbM[BgnSt,EndSt]
    JVec <- (exp(Tm*Lam)-exp(-Tm*BgnRate))/(Lam+BgnRate)
    PrbV <- (RateM[BgnSt,]/PrbBgnEnd)*
      rowSums(U*matrix(rep(InvU[,EndSt]*JVec,nSt),nrow=nSt,byrow=TRUE))
    NewSt <- sample(StSp[-BgnSt],size=1,prob=pmax(0,PrbV[-BgnSt]))
    ##--------------------------------------
    ## Simulate (conditional) waiting time 
    ##--------------------------------------
    ## Define cumulative density function
    CDF.Fct <- function(tm){
      KVec.tm <- exp(Tm*Lam)/(Lam+BgnRate)*(1-exp(-(Lam+BgnRate)*tm))
      CDF.tm <- (RateM[BgnSt,NewSt]/PrbBgnEnd/PrbV[NewSt])*
        sum(U[NewSt,]*InvU[,EndSt]*KVec.tm)
      return(CDF.tm)
    }
    ## Waiting time
    rU <- runif(1)
    RootFct <- function(tm,rU) CDF.Fct(tm)-rU
    uroot <- uniroot(RootFct,interval=c(0,Tm),rU=rU)
    ## cat("number of iterations",uroot$iter,"\n")
    Wt <- uroot$root
    ##----------------------------------------------------------
    ## Update states and times of sample path
    X[cnt+1] <- NewSt
    T[cnt+1] <- T[cnt]+Wt
    ##-----------------------------------------------------------
    ## Repeat procedure with new BgnSt state and new time interval
    cnt <- cnt+1
    Tm <- Tm-Wt
    BgnSt <- NewSt
    ## Numbers needed to determine if a new jump is made
    PrbM <- U%*%diag(exp(Tm*Lam))%*%InvU
    BgnRate <- -RateM[BgnSt,BgnSt]
    PrbBgnBgn <- PrbM[BgnSt,BgnSt]
  }
  ptm3 <- proc.time()[1]
  ## Output state changes and corresponding times
  Path <- list()
  Path$St <- c(X,EndSt)
  Path$Tm <- c(T,EndTm)
  Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
  return(Path)
}
##---------------------------------------------------------------------
## ii) EXPECTED NUMBER OF JUMPS
##     CONDITIONAL ON BEGINNING AND ENDING STATE
##---------------------------------------------------------------------
JFct <- function(Lam,Tm){
  nSt <- length(Lam)
  JMat <- matrix(0,nSt,nSt)
  expTmLam <- exp(Tm*Lam)
  for (i in 1:nSt){
    for (j in 1:nSt){
      JMat[i,j] <- ifelse(identical(all.equal(Lam[i],Lam[j]),TRUE),
                          Tm*expTmLam[i],
                          (expTmLam[i]-expTmLam[j])/(Lam[i]-Lam[j]))
    }
  }
  return(JMat)
}
cndNsbst <- function(BgnSt,EndSt,RateM,Tm){
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
  ## Integrals
  JMat <- JFct(Lam,Tm)
  ## Expected number of jumps
  N <- 0
  for (i in 1:nSt){
    for (j in StSp[-i]){
      N <- N+RateM[i,j]*
        sum( U[BgnSt,]*InvU[,i]*colSums(U[j,]*InvU[,EndSt]*JMat) )/PrbBE
    }
  }
  return(N)
}
cndNsbstV <- function(BgnSt,EndSt,RateM,TmV){
  len <- length(TmV)
  NV <- rep(0,len)
  for (i in 1:len){
    cat("evaluating no",i,"out of",len,"\n")
    NV[i] <- cndNsbst(BgnSt,EndSt,RateM,TmV[i])
  }
  return(NV)
}
