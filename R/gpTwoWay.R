###############################################################################
###############################################################################
gpTwoWay<-function(formula, data, method = c("gPB","gPQ"), alpha = 0.05, na.rm = TRUE, verbose = TRUE) 
{
  data <- model.frame(formula, data)
  fml <-as.character(formula)
  ftmp <- strsplit(fml,"~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1],"[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1] #Drop spaces
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  InterFacAFacB <- paste(FacA,":",FacB, sep = "")
  method = match.arg(method)
  nM <-10000 # MonteCarlo samples
  
  if (!is.data.frame(data)) stop("Data must be in data.frame class.")
  if(length(Factors)!=2) stop("Please correct the RHS of the formula. Formula must include two factors.")
  if(!is.factor(data[,colnames(data)==FacA])) stop(paste(FacA, "must be a factor."))
  if(!is.factor(data[,colnames(data)==FacB])) stop(paste(FacB, "must be a factor."))
  if(!is.numeric(data[,colnames(data)==y])) stop(paste(y, "must be a numeric."))
  
  if (na.rm){
    completeObs <- complete.cases(data)
    data <- data[completeObs,]
  }
  
  cnames <- colnames(data)
  
  ADat <- as.vector(data[,cnames==FacA])
  BDat <- as.vector(data[,cnames==FacB])
  yDat <- as.vector(data[,cnames==y])
  
  As <- unique(ADat) #table(ADat)
  Bs <- unique(BDat) #table(BDat)
  I <- length(As)  # number of A Factor levels
  J <- length(Bs)
  # Compute Sample means and variances and call 
  XBij <- matrix(NA,I,J) # Sample means
  S2ij <- matrix(NA,I,J) # Sample variances
  Nij <- matrix(NA,I,J) # Sample sizes
  for (i in 1:I) {
    if (1==1) { # This method fails when there are cells with no data
      Yi <- yDat[ADat==As[i]]
      Bi <- BDat[ADat==As[i]]
      XBij[i,]<- tapply(Yi, Bi, mean) 
      S2ij[i,]<- tapply(Yi, Bi, var) 
      Nij[i,]<- tapply(Yi, Bi, length) 
    }
    else {
      for (j in 1:J)	{
        Yij <- yDat[ADat==As[i] & BDat==Bs[j]] # data from ijth cell
        XBij[i,j] <- mean(Yij, na.rm=T)
        S2ij[i,j] <- var(Yij, na.rm=T)
        Nij[i,j] <- length(Yij)
      }
    }
    
  }
  colnames(XBij)=colnames(S2ij)=colnames(Nij)<- Bs
  rownames(XBij)=rownames(S2ij)=rownames(Nij)<- As 

  ### ******************* Pval code below read data column wise, and so transpose; Sample size is a vector 
  
  # Compute p-values for interaction effect 
  Tobs = InterSumSqMat(Nij,XBij,S2ij)
  
  #print("---------")
  if (method=="gPB") {
    pval.Int=pvalGPB(Nij,XBij,S2ij,Tobs,nM,method="gPB") 
    
    # print("Interation effect:  Generalized p-value by PB method")		
  }else {
    pval.Int=pvalGPB(Nij, XBij,S2ij,Tobs,nM, method="gPQ") 
    #  print("Generalized p-value by GPQ method")		
    
  }
  #print(pval.Int)
  
  ########################## Compute p-values for main effect A
  TobsA = AInterSumSqMat(Nij,XBij,S2ij) 
  if (method=="gPB") {
    pval.FacA=FacApvalGPB(Nij, XBij,S2ij,TobsA, nM, method="gPB") 
   }else {
    pval.FacA=FacApvalGPB(Nij, XBij,S2ij,TobsA, nM, method="gPQ")     
  }
  #print(pval.FacA)
  
  ########################### Compute p-values for main effect B
  TobsB = BInterSumSqMat(Nij,XBij,S2ij) 
  if (method=="gPB") {
    pval.FacB=FacBpvalGPB(Nij,XBij,S2ij,TobsB, nM, method="gPB") 
    method.name<-"Generalized p-value by PB method"	
    
  }else {
    pval.FacB=FacBpvalGPB(Nij, XBij,S2ij,TobsB,nM, method="gPQ") 
    method.name<-"Generalized p-value by GPQ method"
  }
  #print(pval.FacB)
  
  store = data.frame(matrix(NA, nrow = 3, ncol = 3))
  store$X1 = c(FacA, FacB, InterFacAFacB)
  store$X2 = c(pval.FacA,pval.FacB, pval.Int)
  store$X3 = ifelse(store$X2 > alpha, "Not reject", "Reject")
  colnames(store) = c("Factor", "  P-value", "  Result")
  
  store4<-store
  
  if (verbose){

  ncharacter <- matrix(NA,dim(store4)[1],dim(store4)[2])
  for (i in 1:dim(store4)[1]){
    
    for (j in 1:dim(store4)[2]){
      
      ncharacter[i,j] <- nchar(toString((store4[i,j])))
      
    }
    
  }
  maxentry <- sum(apply(ncharacter, MARGIN = 2, function(x) max(x, na.rm=TRUE)))
  if (maxentry<25) maxentry <- 25
  
  if (method == "gPB"){
    
    line<-paste(c("  ",method.name, " (alpha = ",alpha,")"),sep = "")
    line2<-paste(c(rep("=",round((maxentry+10-32)/2,0)),rep("=",40),rep("=",round((maxentry+10-32)/2,0))),sep = "")
    
    
  }else if (method == "gPQ"){
    line<-paste(c("  ",method.name, " (alpha = ",alpha,")"),sep = "")
    line2<-paste(c(rep("=",round((maxentry+10-33)/2,0)),rep("=",40),rep("=",round((maxentry+10-33)/2,0))),sep = "")
    
  }else stop("Please correct the method argument.")
  
  cat("\n",line,sep = "")
  cat("\n",line2,sep = "","\n")
  print(store4,row.names = FALSE)
  cat(line2,sep = "","\n\n")
  }
  
  result <- list()
  
  result$output <- store4
  result$alpha <- alpha
  result$method <- method.name
  result$data <- data
  result$formula <- formula
  invisible(result)
  
}
###############################################################################
###############################################################################
# Compute Genralized PB Test for interaction
# Xobs and S2obs matrices below
# Tobs = InterSumSqMat(nn,xb,s2)
# pvalGPB(nij, Xobs, S2obs, Tobs)
pvalGPB=function(nij, Xobs, S2obs, Tobs, M=1000,method="gPB") 
{
  # All matrices are read columnn by column
  ndat = as.vector(t(nij))
  nroot = sqrt(ndat)
  ndat1 = ndat-1 # ndat is a vector
  I = nrow(Xobs)
  J = ncol(Xobs)
  IJ = I*J
  sqr = sqrt(S2obs)
  sqrvec = as.vector(t(sqr))
  Ones = matrix(1,I, J)
  pvec = ttvec =rep(0, M)
  
  # Use a sample from readily created random numbers in calling program
  
  sdd = sqrvec/nroot
  for (i in 1:M) { 	
    Xmean = mean(Xobs)
    Zvec = sdd*rnorm(IJ, 0, 1)
    Zvec = Zvec -mean(Zvec) # Mean center for better performnec
    XB = t(matrix(Zvec, J, I))
    Uvec = rchisq(IJ, ndat1)
    S2 = S2obs * t(matrix(Uvec, J, I)) / (nij-1)
    
    Tvec = rt(IJ, ndat1)
    Tmat = t(matrix(Tvec, J, I))
    if (method=="gPB") TmGPB = InterSumSqMat(nij,XB,S2)
    else { 
      wt = nij/S2obs
      Tsn = Tmat / sqrt(wt) 
      TmGPB = InterSumSqMat(nij,Tsn,S2obs) 
    }	
    pvec[i] = TmGPB > Tobs
    ttvec[i] = TmGPB
    
  }  
  
  pval=mean(pvec)
  
  # print(pval)
  return(pval)
  
}
###############################################################################
###############################################################################
################# Compute Observed Intervction Sum of Squares #############
# wtmult = 1 if weight is determined within the function; otherwise it used to multiple regular woght
# Tobs = InterSumSq(Nij,XBij,S2ij)
# test as tmp = InterSumSqMat(nn,xb, s2)
#
InterSumSqMat = function(nij,XX, S2obs) # Need both Xrand and Xobs both when passing Monte Carlo smaples
{	
  # Re test mean for testing
  I = nrow(XX)
  J = ncol(XX)
  IJ = I*J
  wt = nij/S2obs
  Ii = diag(1,I)
  Ij = diag(1,J)
  Oneij = rep(1, IJ)
  Onei = rep(1,I)
  Onej = rep(1,J)
  # Design Matrix always observed
  X = as.matrix(cbind(Oneij, kronecker(Ii, Onej), kronecker(Onei, Ij)))
  I1 = I+1
  I2 = I+2
  wtvec = as.vector(t(wt))
  wi = rowSums(wt)  
  wj = colSums(wt)
  Iw = length(wtvec)+1
  IJS = I+J+1
  l1 =l2 = rep(0, IJS)
  L = matrix(0, IJS,2)
  l1[2:I1] = wi # As ib Xu et al we can also set this to 1
  l2[I2:IJS] = wj
  L[,1]=l1
  L[,2]=l2
  LL = t(L%*%t(L))
  SIG = diag(1/wtvec) # In Xu notastion wt is reciprocal
  SigI = diag(wtvec)
  SigIhalf = diag(sqrt(wtvec))
  XbVec = as.vector(t(XX)) 
  
  leftPart <-SigIhalf %*% X  
  midPart <- t(X) %*% SigI %*% X + LL
  midPartI <- solve(midPart) # ginv(midPart) 
  midmidPart  <- SigI %*% X %*% midPartI %*% t(X)%*%SigI
  T = t(XbVec) %*% (SigI - midmidPart) %*% XbVec
  return(T)
}
###############################################################################
###############################################################################
AInterSumSqMat = function(nij,XX, S2obs) # Need both Xrand and Xobs both when passing Monte Carlo smaples
{	
  # Re test mean for testing
  I = nrow(XX)
  J = ncol(XX)
  IJ = I*J
  wt = nij/S2obs
  Ii = diag(1,I)
  Ij = diag(1,J)
  Oneij = rep(1, IJ)
  Onei = rep(1,I)
  Onej = rep(1,J)
  # Design Matrix always observed
  X = as.matrix(cbind(Oneij, kronecker(Onei, Ij)))
  I1 = I+1
  I2 = I+2
  wtvec = as.vector(t(wt))
  wj = colSums(wt)
  Iw = length(wtvec)+1
  #
  IJS = J+1
  L = rep(0, IJS)
  L[2:IJS]=wj
  LL = t(L%*%t(L))
  SIG = diag(1/wtvec) # In Xu notastion wt is reciprocal
  SigI = diag(wtvec)
  SigIhalf = diag(sqrt(wtvec))
  XbVec = as.vector(t(XX)) 
  
  leftPart <-SigIhalf %*% X  #### S^-1/2 X # Correct
  midPart <- t(X) %*% SigI %*% X + LL
  #  midpart <<- midPart
  midPartI <- solve(midPart) # ginv(midPart) # Note ^-1 is typo in Xu et al paper
  midmidPart  <- SigI %*% X %*% midPartI %*% t(X)%*%SigI
  T = t(XbVec) %*% (SigI - midmidPart) %*% XbVec
  return(T)
}
###############################################################################
###############################################################################
BInterSumSqMat = function(nij,XX, S2obs) # Need both Xrand and Xobs both when passing Monte Carlo smaples
{	
  # Re test mean for testing
  I = nrow(XX)
  J = ncol(XX)
  IJ = I*J
  wt = nij/S2obs
  Ii = diag(1,I)
  Ij = diag(1,J)
  Oneij = rep(1, IJ)
  Onei = rep(1,I)
  Onej = rep(1,J)
  # Design Matrix always observed
  X = as.matrix(cbind(Oneij, kronecker(Ii, Onej)))
  I1 = I+1
  I2 = I+2
  wtvec = as.vector(t(wt))
  wi = rowSums(wt)  
  Iw = length(wtvec)+1
  #
  IJS = I+1
  L =rep(0, IJS)
  L[2:I1] = wi 
  LL = t(L%*%t(L))
  SIG = diag(1/wtvec) # In Xu notastion wt is reciprocal
  SigI = diag(wtvec)
  SigIhalf = diag(sqrt(wtvec))
  XbVec = as.vector(t(XX)) 
  
  leftPart <-SigIhalf %*% X  #### S^-1/2 X # Correct
  midPart <- t(X) %*% SigI %*% X + LL
  #  midpart <<- midPart
  midPartI <- solve(midPart) # ginv(midPart) # Note ^-1 is typo in Xu et al paper
  midmidPart  <- SigI %*% X %*% midPartI %*% t(X)%*%SigI
  T = t(XbVec) %*% (SigI - midmidPart) %*% XbVec
  return(T)
}
###############################################################################
###############################################################################
FacApvalGPB=function(nij, Xobs, S2obs, Tobs, M=1000,method="gPB") 
{
  # All matrices are read columnn by column
  ndat = as.vector(t(nij))
  nroot = sqrt(ndat)
  ndat1 = ndat-1 # ndat is a vector
  I = nrow(Xobs)
  J = ncol(Xobs)
  IJ = I*J
  sqr = sqrt(S2obs)
  sqrvec = as.vector(t(sqr))
  Ones = matrix(1,I, J)
  pvec = ttvec =rep(0, M)
  
  # Use a sample from readily created random numbers in calling program
  
  sdd = sqrvec/nroot
  for (i in 1:M) { 	
    Xmean = mean(Xobs)
    Zvec = sdd*rnorm(IJ, 0, 1)
    Zvec = Zvec -mean(Zvec) # Mean center for better performnec
    XB = t(matrix(Zvec, J, I))
    Uvec = rchisq(IJ, ndat1)
    S2 = S2obs * t(matrix(Uvec, J, I)) / (nij-1)
    
    Tvec = rt(IJ, ndat1)
    Tmat = t(matrix(Tvec, J, I))
    if (method=="gPB") TmGPB = AInterSumSqMat(nij,XB,S2)
    else { 
      wt = nij/S2obs
      Tsn = Tmat / sqrt(wt) 
      TmGPB = AInterSumSqMat(nij,Tsn,S2obs) 
    }	
    pvec[i] = TmGPB > Tobs
    ttvec[i] = TmGPB
    
  }  
  
  pval=mean(pvec)
  
  return(pval)
}
###############################################################################
###############################################################################
FacBpvalGPB=function(nij, Xobs, S2obs, Tobs, M=1000,method="gPB") 
{
  # All matrices are read columnn by column
  ndat = as.vector(t(nij))
  nroot = sqrt(ndat)
  ndat1 = ndat-1 # ndat is a vector
  I = nrow(Xobs)
  J = ncol(Xobs)
  IJ = I*J
  sqr = sqrt(S2obs)
  sqrvec = as.vector(t(sqr))
  Ones = matrix(1,I, J)
  pvec = ttvec =rep(0, M)
  
  # Use a sample from readily created random numbers in calling program
  
  sdd = sqrvec/nroot
  for (i in 1:M) { 	
    Xmean = mean(Xobs)
    Zvec = sdd*rnorm(IJ, 0, 1)
    Zvec = Zvec -mean(Zvec) # Mean center for better performnec
    XB = t(matrix(Zvec, J, I))
    Uvec = rchisq(IJ, ndat1)
    S2 = S2obs * t(matrix(Uvec, J, I)) / (nij-1)
    
    Tvec = rt(IJ, ndat1)
    Tmat = t(matrix(Tvec, J, I))
    if (method=="gPB") TmGPB = BInterSumSqMat(nij,XB,S2)
    else { 
      wt = nij/S2obs
      Tsn = Tmat / sqrt(wt) 
      TmGPB = BInterSumSqMat(nij,Tsn,S2obs)
    }	
    pvec[i] = TmGPB > Tobs
    ttvec[i] = TmGPB
    
  }  
  
  pval=mean(pvec)
  
  #  print(pval)
  return(pval)
  
}
###############################################################################
###############################################################################



