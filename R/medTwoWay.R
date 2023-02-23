medTwoWay<-function(formula, data, alpha = 0.05, na.rm = TRUE, verbose = TRUE){
  data_ <- model.frame(formula, data)
  fml <- as.character(formula)
  ftmp <- strsplit(fml, "~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1], "[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1]
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  data_[,FacA]<-factor(data_[,FacA])
  data_[,FacB]<-factor(data_[,FacB])
  InterFacAFacB <- paste(FacA, ":", FacB, sep = "")
  nM <- 10000
  if (!is.data.frame(data_)) 
    stop("Data must be in data.frame class.")
  if (length(Factors) != 2) 
    stop("Please correct the RHS of the formula. Formula must include two factors.")
  if (!is.factor(data_[, colnames(data_) == FacA])) 
    stop(paste(FacA, "must be a factor."))
  if (!is.factor(data_[, colnames(data_) == FacB])) 
    stop(paste(FacB, "must be a factor."))
  if (!is.numeric(data_[, colnames(data_) == y])) 
    stop(paste(y, "must be a numeric."))
  if (na.rm) {
    completeObs <- complete.cases(data_)
    data_ <- data_[completeObs, ]
  }
  J<-nlevels(data_[,FacA])
  K<-nlevels(data_[,FacB])
  p<-J*K
  n<-nrow(data_)/(J*K)
  data<-groupping(data_,FacA,FacB,y)
  medians<-matrix(unlist(lapply(data,median)),J,K,byrow = T)
  h<-matrix(length(data[[1]]),J,K,byrow = T)
  sdmed<-matrix(unlist(lapply(data, function(x) medsems(x)^2)),J,K,byrow = T)
  sumfacAmedians<-apply(medians,1,sum)
  sumfacBmedians<-apply(medians,2,sum)
  facAsd<-apply(sdmed,1,function(x) (1/sum(x)))
  facBsd<-apply(sdmed,2,function(x) (1/sum(x)))
  facAinvsd<-apply(sdmed,1,function(x) sum(1/x))
  facBinvsd<-apply(sdmed,2,function(x) sum(1/x))
  facAnuhat<-0
  facBnuhat<-0
  for (j in 1:J) {facAnuhat[j]<-(sum(sdmed[j,]))^2/sum(sdmed[j, ]^2/(h[j, ] - 1))}
  for (k in 1:K) {facBnuhat[k]<-(sum(sdmed[,k]))^2/sum(sdmed[,k]^2/(h[,k] - 1))}
  SDMED<-1/sdmed
  xmed<-matrix(,J,K,byrow =T)
  xav<-matrix(,J,K,byrow =T)
  for (j in 1:J) {
    for(k in 1:K){
      xmed[j,k]<-sum(SDMED[,k]*medians[,k]/facBinvsd[k])+sum(SDMED[j,]*medians[j,]/facAinvsd[j])-sum(SDMED*medians/sum(SDMED))
      xav[j,k]<-(1-SDMED[j,k]*(1/sum(SDMED[j,])+1/sum(SDMED[,k])-1/sum(SDMED)))^2/(h[j,k]-3)
    }
  }
  sumfacAmedianshat<-sum(facAsd*sumfacAmedians)/sum(facAsd)
  sumfacBmedianshat<-sum(facBsd*sumfacBmedians)/sum(facBsd)
  BfacA<-sum((1 - facAsd/sum(facAsd))^2/facAnuhat)
  BfacB<-sum((1 - facBsd/sum(facBsd))^2/facBnuhat)
  FacA_value<-sum(facAsd * (sumfacAmedians - sumfacAmedianshat)^2)/((J - 1)*(1 + 2*(J - 2)*BfacA/(J^2 - 1)))
  FacB_value<-sum(facBsd * (sumfacBmedians - sumfacBmedianshat)^2)/((K - 1)*(1 + 2*(K - 2)*BfacB/(K^2 - 1)))
  p.val_FacA<-pf(FacA_value,J-1,Inf,lower.tail = FALSE)
  p.val_FacB<-pf(FacB_value,K-1,Inf,lower.tail = FALSE)
  InterFacAFacB_value<-sum(SDMED*(medians-xmed)^2)
  interactdof<-(J-1)*(K-1)
  p.val_InterFacAFacB<-1-pchisq(InterFacAFacB_value,interactdof)
  store<-data.frame(matrix(NA,nrow=3,ncol = 5))
  colnames(store)<-c("Factor","df","Statistic","P_value","Result")
  store$Factor<-c(FacA,FacB,InterFacAFacB)
  store$df<-c(paste0("F(",J-1,",Inf)"),paste0("F(",K-1,",Inf)"),paste0("Chisq(",(J - 1) * (K - 1),")"))
  store$Statistic<-c(FacA_value,FacB_value,InterFacAFacB_value)
  store$P_value<-c(p.val_FacA,p.val_FacB,p.val_InterFacAFacB)
  store$Result<-ifelse(store$P_value>alpha,"Not reject","Reject")
  store4<-store
  if(verbose){
    cat("\n","  Two-way ANOVA for Medians ","(alpha = ",alpha,")",sep="")
    cat("\n", "---------------------------------------------------------------------", sep = "","\n")
    print(store4, row.names = FALSE)
    cat("---------------------------------------------------------------------", sep = "", "\n")
  }
  result<-list()
  result$output<-store4
  result$alpha<-alpha
  result$method<-"Two-way ANOVA for medians"
  result$data<-data_
  result$formula<-formula
  attr(result, "class") <- "twt"
  invisible(result)
}
