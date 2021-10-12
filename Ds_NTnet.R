Ds_NTnet <- function(Prime, N, bound){
   # gp set methods in NT_net.
   # Prime is the row vector.
   # produce N dots dispersed well-proportionedly all the whole area of a given closed region.
   # bound is the boundary of a given closed region.
  
  sPrime <- matrix(sqrt(Prime),nrow=1,ncol=length(Prime))
  line <- matrix(seq(1,N,1),nrow=N,ncol=1)
  tVec <- kronecker(line,sPrime)
  cVec <- tVec-floor(tVec)
  temp <- matrix(rep(1,N),nrow=N,ncol=1)
  aAdd <- kronecker(temp, matrix(bound[1,],nrow=1,ncol=length(Prime)))
  bAdd <- kronecker(temp,matrix((bound[2,]-bound[1,]),nrow=1,ncol=length(Prime)))
  ntData <- aAdd + bAdd*cVec
  
  return(ntData)
  
}

