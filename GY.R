##------------------------------------------
## Goldman and Yang rate matrix function
##------------------------------------------
#codon.mat <- read.table(file="~/raleigh/Projects/ExactSmpl/programs/codon.mat",header=TRUE)
codon.mat <- read.table(file="codon.mat",header=TRUE)
aa <- codon.mat[,"aa"]
codon.freq <- codon.mat[,"frq"]
##---------------------------------------------------
## Order: ACGT
purine <- c(TRUE,FALSE,TRUE,FALSE)
##--------------------------------------------------
## Determine stop codons
stop.codon <- rep(FALSE,64) ## stop codon
stop.codon[which(aa=="Stp")] <- TRUE
## Define types of substitutions
syn.ts <- matrix(rep(0,64^2),nrow=64,ncol=64)    ## syn transi
syn.tv <- matrix(rep(0,64^2),nrow=64,ncol=64)    ## syn transv
nonsyn.ts <- matrix(rep(0,64^2),nrow=64,ncol=64) ## nonsyn transi
nonsyn.tv <- matrix(rep(0,64^2),nrow=64,ncol=64) ## nonsyn transv
## Row index determined by (i1,i2,i3)
## Col index determined by (j1,j2,j3)
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      rw <- (i1-1)*16+(i2-1)*4+i3
      ## Change in first position
      for (j1 in 1:4){
        cl <- (j1-1)*16+(i2-1)*4+i3
        if (j1!=i1){
          ts <- ifelse(purine[i1]==purine[j1],1,0)
          syn <- ifelse(aa[rw]==aa[cl],1,0)
          syn.ts[rw,cl] <- syn*ts
          syn.tv[rw,cl] <- syn*(1-ts)
          nonsyn.ts[rw,cl] <- (1-syn)*ts
          nonsyn.tv[rw,cl] <- (1-syn)*(1-ts)
        }
      }
      ## Change in second position
      for (j2 in 1:4) {
        cl <- (i1-1)*16+(j2-1)*4+i3
        if (j2!=i2){
          ts <- ifelse(purine[i2]==purine[j2],1,0)
          syn <- ifelse(aa[rw]==aa[cl],1,0)
          syn.ts[rw,cl] <- syn*ts
          syn.tv[rw,cl] <- syn*(1-ts)
          nonsyn.ts[rw,cl] <- (1-syn)*ts
          nonsyn.tv[rw,cl] <- (1-syn)*(1-ts)
        }
      }
      ## Change in third position
      for (j3 in 1:4) {
        cl <- (i1-1)*16+(i2-1)*4+j3
        if (j3!=i3){
          ts <- ifelse(purine[i3]==purine[j3],1,0)
          syn <- ifelse(aa[rw]==aa[cl],1,0)
          syn.ts[rw,cl] <- syn*ts
          syn.tv[rw,cl] <- syn*(1-ts)
          nonsyn.ts[rw,cl] <- (1-syn)*ts
          nonsyn.tv[rw,cl] <- (1-syn)*(1-ts)
        }
      }
    }
  }
}
##-------------------------------------------------------------
## Define GY matrix
## Input: kappa: ts/tv rate ratio
##        omega: dn/ds ratio
##        codon.freq: vector of length 64
##-------------------------------------------------------------
GYRate <- function(kappa,omega,codon.freq){
  GY.r <- kappa*syn.ts+syn.tv+omega*kappa*nonsyn.ts+omega*nonsyn.tv
  GY.r <- GY.r%*%diag(codon.freq)
  ## Remove stopcodons from rate matrix
  GY.r <- GY.r[!stop.codon,!stop.codon]
  if (identical(all.equal(sum(codon.freq[!stop.codon]),1),FALSE)) {
    cat("codon.freq does not add to one","\n")
  }
  ## Scaling constraint
  scl <- sum(codon.freq[!stop.codon]*GY.r)
  GY.r <- GY.r/scl
  ## Get the diagonal right
  for (i in 1:61) GY.r[i,i] <- -sum(GY.r[i,-i])
  return(GY.r)
}
