## MEST
## Computing the p-value of the MEST for indentifying the association of gene-gene
## or gene-environment interaction in case-control studies.
## outcome: The binary response variable coded as 0 (controls) and 1 (cases). No default;
## snp1: SNP variable coded as 0-0.5-1.  No default;
## covar: All the covariates of interest included in the model as main effects;
## snp2: SNP variable that interacts with snp1. The default is NULL. If not, it will test the gene-gene interaction;
## envir: Environment variable that interacts with snp1. The default is NULL. If not, it will test the gene-environment interaction; 
## Note that snp2 and envir could not be both NULL.
## ifint: Test the joint effects or interactive effects. The default is T, that represents interactive effects; 
## B: iterations.

library(MASS)
library(survival)
library(mvtnorm)


MEST.Pval <- function(outcome, snp1, covar, snp2=NULL, envir=NULL, ifint=T, B){
  
  dex <- which(!is.na(snp1))
  snp1 <- snp1[dex]
  outcome <- outcome[dex]
  covar <- covar[dex,]
  len <- ncol(covar)
  
  delta <- 10^(-2)
  Prime <- c(2,3)
  num <- length(outcome)
  
  # ## single-marker analysis
  # 
  # if(is.null(snp2) & is.null(envir)){
  #   if(sum(is.na(snp1))!=0){
  #     warning("Missing values exist in SNPs.")
  #   }else{
  #       y <- as.matrix(outcome)
  #       geno <- as.matrix(snp1)
  #       if(!is.null(covar)){
  #         cova2 <- as.matrix(covar)
  #       }
  #       p.max3 <- MAX3(y, geno, covariates=cova2, Score.test=T, Wald.test=F, rhombus.formula=T)$p.value
  #   }
  # }
  
  ## gene-gene or gene-environment
  
  dataper <- data.frame(outcome, covar)
  # glm.mc <- glm(outcome~covar,data=dataper,family="binomial")
  glm.mc <- glm(outcome~.,data=dataper,family="binomial")
  betamc <- matrix(glm.mc$coefficients)
  # m.design <- matrix(cbind(rep(1,num),covar),num,len+1)
  m.design <- as.matrix(cbind(rep(1,num),covar))
  prmc <- exp(m.design%*%betamc)/(1+exp(m.design%*%betamc))
  
  if(!is.null(snp2) & !is.null(envir)){
    warning("Invalid input")
  }
  
  ## gene-gene
  
  else if(!is.null(snp2) & is.null(envir)){
    
    snp2 <- snp2[dex]
    
    if(length(which(snp1==1))<=1|length(which(snp2==1))<=1){
      warning("rare risk allele!")
    }
    
    ## test interaction
    
    else if(ifint){
      
      ## permutation 
      
      mestmc <- c()
      
      for (itermc in 1:B){

        cat('########### the',itermc,'round...###########\n')

        outcomemc <- rep(NA,num)
        for(i in 1:num){
          outcomemc[i] <- rbinom(1,1,prmc[i])
        }
        
        # parameter space
        tmc <- 1
        shrinkratemc <- 0.5
        Nmc <- 300
        boundmc <- matrix(c(0,1,0,1),2,2)
        ethetempmc <- matrix(NA,nrow=100,ncol=2)
        while (max((boundmc[2,]-boundmc[1,])/2)>=delta) {
          try(
            {
              if(tmc>1){
                ntnetmc <- rbind(Ds_NTnet(Prime,Nmc,boundmc),ethetempmc[tmc-1,])
              }else{
                ntnetmc <- Ds_NTnet(Prime,Nmc,boundmc)
              }
              if(length(which(ntnetmc[,1]<10^(-2)|ntnetmc[,2]<10^(-2)))!=0){
                ntnetmc <- ntnetmc[-which(ntnetmc[,1]<10^(-2)|ntnetmc[,2]<10^(-2)),]
              }
              tempmc <- apply(ntnetmc,1,funcmc1,dataymc=outcomemc, snp1=snp1, snp2=snp2, covar=covar)
              mesttempmc <- max(tempmc)
              ethetempmc[tmc,] <- ntnetmc[which(tempmc==mesttempmc),]
              aboundmc <- apply(rbind((ethetempmc[tmc,]-shrinkratemc*((boundmc[2,]-boundmc[1,])/2)),boundmc[1,]),2,max)
              bboundmc <- apply(rbind((ethetempmc[tmc,]+shrinkratemc*((boundmc[2,]-boundmc[1,])/2)),boundmc[2,]),2,min)
              boundmc <- rbind(aboundmc,bboundmc)
              
              Nmc <- 100
              tmc <- tmc+1
              shrinkratemc <- shrinkratemc^tmc
            }, silent = TRUE
            )
        }
        mestmc[itermc]=mesttempmc
      }
      
      ## compute p-value

      # parameter space
      t <- 1
      shrinkrate <- 0.5
      N <- 300
      bound <- matrix(c(0,1,0,1),2,2)
      ethetemp <- matrix(NA,nrow=100,ncol=2)
      while (max((bound[2,]-bound[1,])/2)>=delta) {
        try(
          {
            if(t>1){
              ntnet <- rbind(Ds_NTnet(Prime,N,bound),ethetemp[t-1,])
            }else{
              ntnet <- Ds_NTnet(Prime,N,bound)
            }
            if(length(which(ntnet[,1]<10^(-2)|ntnet[,2]<10^(-2)))!=0){
              ntnet <- ntnet[-which(ntnet[,1]<10^(-2)|ntnet[,2]<10^(-2)),]
            }
            temp <- apply(ntnet,1,func1, outcome=outcome, snp1=snp1, snp2=snp2, covar=covar)
            mesttemp <- max(temp)
            ethetemp[t,] <- ntnet[which(temp==mesttemp),]
            abound <- apply(rbind((ethetemp[t,]-shrinkrate*((bound[2,]-bound[1,])/2)),bound[1,]),2,max)
            bbound <- apply(rbind((ethetemp[t,]+shrinkrate*((bound[2,]-bound[1,])/2)),bound[2,]),2,min)
            bound <- rbind(abound,bbound)
            
            N <- 100
            t <- t+1
            shrinkrate <- shrinkrate^t
          }, silent = TRUE
          )
      }
      meststat <- mesttemp
  
      p.value <- sum(abs(mestmc)>=meststat)/B
      return(p.value)
    }
    
    ## test joint effects
    
    else{
      
      ## permutation 
      
      mestmc <- c()
      for (itermc in 1:B){

        cat('########### the',itermc,'round...###########\n')

        outcomemc <- rep(NA,num)
        for(i in 1:num){
          outcomemc[i] <- rbinom(1,1,prmc[i])
        }
        
        # parameter space
        tmc <- 1
        shrinkratemc <- 0.5
        Nmc <- 1000
        boundmc <- matrix(c(0,1,0,1),2,2)
        ethetempmc <- matrix(NA,nrow=100,ncol=2)
        while (max((boundmc[2,]-boundmc[1,])/2)>=delta) {
          try(
          {
            if(tmc>1){
              ntnetmc <- rbind(Ds_NTnet(Prime,Nmc,boundmc),ethetempmc[tmc-1,])
            }else{
              ntnetmc <- Ds_NTnet(Prime,Nmc,boundmc)
            }
            if(length(which(ntnetmc[,1]<10^(-2)|ntnetmc[,2]<10^(-2)))!=0){
              ntnetmc <- ntnetmc[-which(ntnetmc[,1]<10^(-2)|ntnetmc[,2]<10^(-2)),]
            }
            tempmc <- apply(ntnetmc,1,funcmc2,dataymc=outcomemc, snp1=snp1, snp2=snp2, covar=covar)
            mesttempmc <- max(tempmc)
            ethetempmc[tmc,] <- ntnetmc[which(tempmc==mesttempmc),]
            aboundmc <- apply(rbind((ethetempmc[tmc,]-shrinkratemc*((boundmc[2,]-boundmc[1,])/2)),boundmc[1,]),2,max)
            bboundmc <- apply(rbind((ethetempmc[tmc,]+shrinkratemc*((boundmc[2,]-boundmc[1,])/2)),boundmc[2,]),2,min)
            boundmc <- rbind(aboundmc,bboundmc)
            
            Nmc <- 120
            tmc <- tmc+1
            shrinkratemc <- shrinkratemc^tmc
          }, silent = TRUE
          )
        }
        mestmc[itermc] <- mesttempmc
      }
      
      ## compute p-value
      
      # parameter space
      t <- 1
      shrinkrate <- 0.5
      N <- 1000
      bound <- matrix(c(0,1,0,1),2,2)
      ethetemp <- matrix(NA,nrow=100,ncol=2)
      while (max((bound[2,]-bound[1,])/2)>=delta) {
        try(
          {
            if(t>1){
              ntnet <- rbind(Ds_NTnet(Prime,N,bound),ethetemp[t-1,])
            }else{
              ntnet <- Ds_NTnet(Prime,N,bound)
            }
            if(length(which(ntnet[,1]<10^(-2)|ntnet[,2]<10^(-2)))!=0){
              ntnet <- ntnet[-which(ntnet[,1]<10^(-2)|ntnet[,2]<10^(-2)),]
            }
            temp <- apply(ntnet,1,func2, outcome=outcome, snp1=snp1, snp2=snp2, covar=covar)
            mesttemp <- max(temp)
            ethetemp[t,] <- ntnet[which(temp==mesttemp),]
            abound <- apply(rbind((ethetemp[t,]-shrinkrate*((bound[2,]-bound[1,])/2)),bound[1,]),2,max)
            bbound <- apply(rbind((ethetemp[t,]+shrinkrate*((bound[2,]-bound[1,])/2)),bound[2,]),2,min)
            bound <- rbind(abound,bbound)
            
            N <- 120
            t <- t+1
            shrinkrate <- shrinkrate^t
          }, silent = TRUE
          )
      }
      meststat <- mesttemp
      
      p.value <- sum(abs(mestmc)>=meststat)/B
      return(p.value)
    }
  }
  
  ## gene-environment
  
  else if(is.null(snp2) & !is.null(envir)){
    
    envir <- envir[dex]
    theta <- seq(0.1,1,0.01)
    
    # dummy variable
    z <- matrix(0,nrow=num,ncol=2)
    z[which(snp1==0.5),2] <- 1
    z[which(snp1==1),1] <- 1
    
    eex <- envir
    dum1 <- z[,1]
    dum2 <- z[,2]
    dum1e <- z[,1]*envir
    dum2e <- z[,2]*envir
    
    if(length(which(snp1==1))<=1){
      warning("rare risk allele!")
    }
    
    ## test interaction
    
    else if(ifint){
      
      ## permutation
      
      mestmc <- c()
      for (itermc in 1:B){

        cat('########### the',itermc,'round...###########\n')

        outcomemc <- rep(NA,num)
        for(i in 1:num){
          outcomemc[i] <- rbinom(1,1,prmc[i])
        }
        
        tempmc <- c()
        for(k in 1:length(theta)){
          dataframc <- data.frame(outcomemc,covar,eex,(dum1+theta[k]*dum2))
          # glmmc.null <- glm(outcomemc~covar+dataframc[,4]+dataframc[,5],data=dataframc,family="binomial")
          glmmc.null <- glm(outcomemc~.,data=dataframc,family="binomial")
          betamc.null <- matrix(c(glmmc.null$coefficients,0))
          # ge.design <- matrix(cbind(rep(1,num),covar,eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e),num,len+4)
          ge.design <- as.matrix(cbind(rep(1,num),covar,eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e))
          
          # Score function
          scormc <- sum(outcomemc*(dum1e+theta[k]*dum2e)-exp(ge.design%*%betamc.null)*(dum1e+theta[k]*dum2e)/(1+exp(ge.design%*%betamc.null)))
        
          # Informatic metric
          fmmc <- matrix(0,nrow=len+4,ncol=len+4)
          for(fnumi in 1:(len+4)){
            for(fnumj in 1:(len+4)){
              fmmc[fnumi,fnumj] <- sum((exp(ge.design%*%betamc.null)/(1+exp(ge.design%*%betamc.null))^2)*ge.design[,fnumi]*ge.design[,fnumj])
            }
          }
          fminvmc <- solve(fmmc)
          tempmc[k] <- (scormc^2)*fminvmc[6,6]
        }
        mestmc[itermc] <- max(tempmc)
      }
      
      ## compute p-value

      temp <- c()
      for(kk in 1:length(theta)){
        datafra <- data.frame(outcome,covar,eex,(dum1+theta[kk]*dum2))
        # glm.null <- glm(outcome~covar+datafra[,4]+datafra[,5],data=datafra,family="binomial")
        glm.null <- glm(outcome~.,data=datafra,family="binomial")
        beta.null <- matrix(c(glm.null$coefficients,0))
        # ge.design <- matrix(cbind(rep(1,num),covar,eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e),num,len+4)
        ge.design <- as.matrix(cbind(rep(1,num),covar,eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e))

        # Score function
        scor <- sum(outcome*(dum1e+theta[kk]*dum2e)-exp(ge.design%*%beta.null)*(dum1e+theta[kk]*dum2e)/(1+exp(ge.design%*%beta.null)))
        
        # Informatic metric
        fm <- matrix(0,nrow=len+4,ncol=len+4)
        for(fnumi in 1:(len+4)){
          for(fnumj in 1:(len+4)){
            fm[fnumi,fnumj] <- sum((exp(ge.design%*%beta.null)/(1+exp(ge.design%*%beta.null))^2)*ge.design[,fnumi]*ge.design[,fnumj])
          }
        }
        fminv <- solve(fm)
        temp[kk] <- (scor^2)*fminv[6,6]
      }
      meststat <- max(temp)
      
      p.value <- sum(abs(mestmc)>=meststat)/B
      return(p.value)
    }
    
    ## test joint effects
    
    else{
      
      ## permutation
      
      mestmc <- c()
      for (itermc in 1:B){

        cat('########### the',itermc,'round...###########\n')

        outcomemc <- rep(NA,num)
        for(i in 1:num){
          outcomemc[i] <- rbinom(1,1,prmc[i])
        }
        
        dataxynullmc <- data.frame(outcomemc,covar)
        # glmmc.null <- glm(outcomemc~covar,data=dataxynullmc,family="binomial")
        glmmc.null <- glm(outcomemc~.,data=dataxynullmc,family="binomial")
        betamc.null <- matrix(c(glmmc.null$coefficients,0,0,0))

        tempmc <- c()
        
        for(k in 1:length(theta)){
          # ge.design <- matrix(cbind(rep(1,num),covar,eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e),num,len+4)
          ge.design <- as.matrix(cbind(rep(1,num),covar,eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e))
          
          # Score function
          # sc.design <- matrix(cbind(eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e),num,3)
          sc.design <- as.matrix(cbind(eex,dum1+theta[k]*dum2,dum1e+theta[k]*dum2e))
          scormc <- t(matrix(outcomemc))%*%sc.design-t(matrix(exp(ge.design%*%betamc.null)/(1+exp(ge.design%*%betamc.null))))%*%sc.design
        
          # Informatic metric
          fmmc <- matrix(0,nrow=len+4,ncol=len+4)
          for(fnumi in 1:(len+4)){
            for(fnumj in 1:(len+4)){
              fmmc[fnumi,fnumj] <- sum((exp(ge.design%*%betamc.null)/(1+exp(ge.design%*%betamc.null))^2)*ge.design[,fnumi]*ge.design[,fnumj])
            }
          }
          fminvmc <- solve(fmmc)
          tempmc[k] <- scormc%*%fminvmc[4:6,4:6]%*%t(scormc)
        }
        mestmc[itermc] <- max(tempmc)
      }
      
      ## compute p-value
      
      dataxynull <- data.frame(outcome,covar)
      # glm.null <- glm(outcome~covar,data=dataxynull,family="binomial")
      glm.null <- glm(outcome~.,data=dataxynull,family="binomial")
      beta.null <- matrix(c(glm.null$coefficients,0,0,0))
      
      temp <- c()
      
      for(kk in 1:length(theta)){
        # ge.design <- matrix(cbind(rep(1,num),covar,eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e),num,len+4)
        ge.design <- as.matrix(cbind(rep(1,num),covar,eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e))
        
        # Score function
        # sc.design <- matrix(cbind(eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e),num,3)
        sc.design <- as.matrix(cbind(eex,dum1+theta[kk]*dum2,dum1e+theta[kk]*dum2e))
        scor <- t(matrix(outcome))%*%sc.design-t(matrix(exp(ge.design%*%beta.null)/(1+exp(ge.design%*%beta.null))))%*%sc.design
    
        # Informatic metric
        fm <- matrix(0,nrow=len+4,ncol=len+4)
        for(fnumi in 1:(len+4)){
          for(fnumj in 1:(len+4)){
            fm[fnumi,fnumj] <- sum((exp(ge.design%*%beta.null)/(1+exp(ge.design%*%beta.null))^2)*ge.design[,fnumi]*ge.design[,fnumj])
          }
        }
        fminv <- solve(fm)
        temp[kk] <- scor%*%fminv[4:6,4:6]%*%t(scor)
      }
      meststat <- max(temp)
      
      p.value <- sum(abs(mestmc)>=meststat)/B
      return(p.value)
    }
  }
}

##-------------------------------------------------------------------

## nt-net 
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

#-----------------------------------------------------------------
## test interaction

funcmc1 <- function(x,dataymc,snp1,snp2,covar){
  num <- length(dataymc)
  len <- ncol(covar)
  
  # dummy variable
  z1 <- matrix(0,nrow=num,ncol=2)
  z1[which(snp1==0.5),2] <- 1
  z1[which(snp1==1),1] <- 1
  z2 <- matrix(0,nrow=num,ncol=2)
  z2[which(snp2==0.5),2] <- 1
  z2[which(snp2==1),1] <- 1
  
  dum1 <- z1[,1]
  dum2 <- z1[,2]
  dum3 <- z2[,1]
  dum4 <- z2[,2]
  dum13 <- z1[,1]*z2[,1]
  dum14 <- z1[,1]*z2[,2]
  dum23 <- z1[,2]*z2[,1]
  dum24 <- z1[,2]*z2[,2]
  # gg.design <- matrix(cbind(matrix(1,num,1),covar,dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24),num,len+4)
  gg.design <- as.matrix(cbind(matrix(1,num,1),covar,dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24))

  dataframc <- data.frame(dataymc,covar,(dum1+x[1]*dum2),(dum3+x[2]*dum4))
  # glmmc.null <- glm(dataymc~covar+dataframc[,4]+dataframc[,5],data=dataframc,family="binomial")
  glmmc.null <- glm(dataymc~.,data=dataframc,family="binomial")
  betamc.null <- matrix(c(glmmc.null$coefficients,0))
  
  # score function
  scormc <- sum(dataymc*(dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24)-exp(gg.design%*%betamc.null)
                *(dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24)/(1+exp(gg.design%*%betamc.null)))
  
  # Informatic metric
  fmmc <- matrix(0,nrow=len+4,ncol=len+4)
  for(fnumi in 1:(len+4)){
    for(fnumj in 1:(len+4)){
      fmmc[fnumi,fnumj] <- sum((exp(gg.design%*%betamc.null)/(1+exp(gg.design%*%betamc.null))^2)*gg.design[,fnumi]*gg.design[,fnumj])
    }
  }
  fminvmc <- solve(fmmc,tol=1e-20)
  tempmc <- (scormc^2)*fminvmc[6,6]
  
  return(tempmc)
}


func1 <- function(y, outcome, snp1, snp2, covar){
  num <- length(outcome)
  len <- ncol(covar)
  
  # dummy variable
  z1 <- matrix(0,nrow=num,ncol=2)
  z1[which(snp1==0.5),2] <- 1
  z1[which(snp1==1),1] <- 1
  z2 <- matrix(0,nrow=num,ncol=2)
  z2[which(snp2==0.5),2] <- 1
  z2[which(snp2==1),1] <- 1
  
  dum1 <- z1[,1]
  dum2 <- z1[,2]
  dum3 <- z2[,1]
  dum4 <- z2[,2]
  dum13 <- z1[,1]*z2[,1]
  dum14 <- z1[,1]*z2[,2]
  dum23 <- z1[,2]*z2[,1]
  dum24 <- z1[,2]*z2[,2]
  # gg.design <- matrix(cbind(matrix(1,num,1),covar,dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24),num,len+4)
  gg.design <- as.matrix(cbind(matrix(1,num,1),covar,dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24))
  
  datafra <- data.frame(outcome,covar,(dum1+y[1]*dum2),(dum3+y[2]*dum4))
  # glm.null <- glm(outcome~covar+datafra[,4]+datafra[,5],data=datafra,family="binomial")
  glm.null <- glm(outcome~.,data=datafra,family="binomial")
  beta.null <- matrix(c(glm.null$coefficients,0))
  
  # Score function
  scor <- sum(outcome*(dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24)-exp(gg.design%*%beta.null)
                *(dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24)/(1+exp(gg.design%*%beta.null)))
  
  # Informatic metric
  # fitemp <- matrix(NA,nrow=num,ncol=len+4)
  # for(fnum in 1:num){
  #   fitemp[fnum,] <- outcome[fnum]*gg.design[fnum,]-exp(gg.design[fnum,]%*%beta.null)/(1+exp(gg.design[fnum,]%*%beta.null))*gg.design[fnum,]
  # }
  # 
  # fm <- matrix(0,nrow=len+4,ncol=len+4)
  # for(fmnum in 1:num){
  #   fm <- fm+matrix(fitemp[fmnum,])%*%t(matrix(fitemp[fmnum,]))
  # }
  
  fm <- matrix(0,nrow=len+4,ncol=len+4)
  for(fnumi in 1:(len+4)){
    for(fnumj in 1:(len+4)){
      fm[fnumi,fnumj] <- sum((exp(gg.design%*%beta.null)/(1+exp(gg.design%*%beta.null))^2)*gg.design[,fnumi]*gg.design[,fnumj])
    }
  }
  fminv <- solve(fm,tol=1e-20)
  temp <- (scor^2)*fminv[6,6]
  
  return(temp)
}

#---------------------------------------------------------------
## test joint effect

funcmc2 <- function(x, dataymc, snp1, snp2, covar){
  num <- length(dataymc)
  len <- ncol(covar)
  
  # dummy variable
  z1 <- matrix(0,nrow=num,ncol=2)
  z1[which(snp1==0.5),2] <- 1
  z1[which(snp1==1),1] <- 1
  z2 <- matrix(0,nrow=num,ncol=2)
  z2[which(snp2==0.5),2] <- 1
  z2[which(snp2==1),1] <- 1
  
  dum1 <- z1[,1]
  dum2 <- z1[,2]
  dum3 <- z2[,1]
  dum4 <- z2[,2]
  dum13 <- z1[,1]*z2[,1]
  dum14 <- z1[,1]*z2[,2]
  dum23 <- z1[,2]*z2[,1]
  dum24 <- z1[,2]*z2[,2]
  # gg.design <- matrix(cbind(matrix(1,num,1),covar,dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24),num,len+4)
  gg.design <- as.matrix(cbind(matrix(1,num,1),covar,dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24))

  dataxynullmc <- data.frame(dataymc,covar)
  # glmmc.null <- glm(dataymc~covar,data=dataxynullmc,family="binomial")
  glmmc.null <- glm(dataymc~.,data=dataxynullmc,family="binomial")
  betamc.null <- matrix(c(glmmc.null$coefficients,0,0,0))
  
  # Score function
  # sc.design <- matrix(cbind(dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24),num,3)
  sc.design <- as.matrix(cbind(dum1+x[1]*dum2,dum3+x[2]*dum4,dum13+x[2]*dum14+x[1]*dum23+x[1]*x[2]*dum24))
  scormc <- t(matrix(dataymc))%*%sc.design-t(matrix(exp(gg.design%*%betamc.null)/(1+exp(gg.design%*%betamc.null))))%*%sc.design
  
  # Fisher
  fmmc <- matrix(0,nrow=len+4,ncol=len+4)
  for(fnumi in 1:(len+4)){
    for(fnumj in 1:(len+4)){
      fmmc[fnumi,fnumj] <- sum((exp(gg.design%*%betamc.null)/(1+exp(gg.design%*%betamc.null))^2)*gg.design[,fnumi]*gg.design[,fnumj])
    }
  }
  fminvmc <- solve(fmmc,tol=1e-20)
  tempmc <- scormc%*%fminvmc[4:6,4:6]%*%t(scormc)

  return(tempmc)

}


func2 <- function(y, outcome, snp1, snp2, covar){
  num <- length(outcome)
  len <- ncol(covar)
  
  # dummy variable
  z1 <- matrix(0,nrow=num,ncol=2)
  z1[which(snp1==0.5),2] <- 1
  z1[which(snp1==1),1] <- 1
  z2 <- matrix(0,nrow=num,ncol=2)
  z2[which(snp2==0.5),2] <- 1
  z2[which(snp2==1),1] <- 1
  
  dum1 <- z1[,1]
  dum2 <- z1[,2]
  dum3 <- z2[,1]
  dum4 <- z2[,2]
  dum13 <- z1[,1]*z2[,1]
  dum14 <- z1[,1]*z2[,2]
  dum23 <- z1[,2]*z2[,1]
  dum24 <- z1[,2]*z2[,2]
  # gg.design <- matrix(cbind(matrix(1,num,1),covar,dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24),num,len+4)
  gg.design <- as.matrix(cbind(matrix(1,num,1),covar,dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24))
  
  dataxynull <- data.frame(outcome,covar)
  # glm.null <- glm(outcome~covar,data=dataxynull,family="binomial")
  glm.null <- glm(outcome~.,data=dataxynull,family="binomial")
  beta.null <- matrix(c(glm.null$coefficients,0,0,0))

  # Score function
  # sc.design <- matrix(cbind(dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24),num,3)
  sc.design <- as.matrix(cbind(dum1+y[1]*dum2,dum3+y[2]*dum4,dum13+y[2]*dum14+y[1]*dum23+y[1]*y[2]*dum24))
  scor <- t(matrix(outcome))%*%sc.design-t(matrix(exp(gg.design%*%beta.null)/(1+exp(gg.design%*%beta.null))))%*%sc.design
  
  # Informatic metric
  fm <- matrix(0,nrow=len+4,ncol=len+4)
  for(fnumi in 1:(len+4)){
    for(fnumj in 1:(len+4)){
      fm[fnumi,fnumj] <- sum((exp(gg.design%*%beta.null)/(1+exp(gg.design%*%beta.null))^2)*gg.design[,fnumi]*gg.design[,fnumj])
    }
  }
  fminv <- solve(fm,tol=1e-20)
  temp <- scor%*%fminv[4:6,4:6]%*%t(scor)
  
  return(temp)
  
}
