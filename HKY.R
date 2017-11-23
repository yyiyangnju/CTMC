##----------------------------------------------------------
## HKY rate matrix function
##----------------------------------------------------------
HKYRate <- function(kappa,frq){
  SymM <- matrix(c(0,kappa,1,1,
                   kappa,0,1,1,
                   1,1,0,kappa,
                   1,1,kappa,0),
                 nrow=4,ncol=4,byrow=TRUE)
  FrqM <- matrix(rep(frq,4),byrow=TRUE,nrow=4,ncol=4)
  RateM <- SymM*FrqM
  ## Get the diagonal right
  for (i in 1:4) RateM[i,i] <- -sum(RateM[i,-i])
  ## Scale such that one substitution is expected in one time unit 
  RateM <- RateM/sum(-diag(RateM)*frq)
  return(RateM)
}
##------------------------------------------------------------
