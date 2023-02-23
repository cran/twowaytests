aovTwoWay<-function(formula, data, alpha = 0.05, na.rm = TRUE, verbose = TRUE){
  data <- model.frame(formula, data)
  fml <- as.character(formula)
  ftmp <- strsplit(fml, "~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1], "[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1]
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  data[,FacA]<-factor(data[,FacA])
  data[,FacB]<-factor(data[,FacB])
  InterFacAFacB <- paste(FacA, ":", FacB, sep = "")
  nM <- 10000
  if (!is.data.frame(data)) 
    stop("Data must be in data.frame class.")
  if (length(Factors) != 2) 
    stop("Please correct the RHS of the formula. Formula must include two factors.")
  if (!is.factor(data[, colnames(data) == FacA])) 
    stop(paste(FacA, "must be a factor."))
  if (!is.factor(data[, colnames(data) == FacB])) 
    stop(paste(FacB, "must be a factor."))
  if (!is.numeric(data[, colnames(data) == y])) 
    stop(paste(y, "must be a numeric."))
  if (na.rm) {
    completeObs <- complete.cases(data)
    data <- data[completeObs, ]
  }
  J<-nlevels(data[,FacA])
  K<-nlevels(data[,FacB])
  p<-J*K
  n<-nrow(data)
  data_mean<-mean(data[,y])
  x<-0
  k<-0
  for (i in levels(data[,FacA])) {
    k<-k+1
    x[k]<-mean(data[data[,FacA]==i,y])
  }
  sssample<-K*(n/p)*sum((x-data_mean)^2)
  x<-0
  k<-0
  for (i in levels(data[,FacB])) {
    k<-k+1
    x[k]<-mean(data[data[,FacB]==i,y])
  }
  ssfactor<-J*(n/p)*sum((x-data_mean)^2)
  data_x<-groupping(data,FacA,FacB,y)
  x<-0
  for (i in c(1:(J*K))) {
    x[i]<-(var(data_x[[i]]))*(length(data_x[[i]])-1)
  }
  sserror<-sum(x)
  sstoplam<-var(data[,y])*(n-1)
  ssinteract<-sstoplam-(sssample+sserror+ssfactor)
  mssample<-sssample/(J-1)
  msfactor<-ssfactor/(K-1)
  msinteract<-ssinteract/((J-1)*(K-1))
  mserror<-sserror/((n-1)-((J-1)+(K-1)+((J-1)*(K-1))))
  fsample<-mssample/mserror
  ffactor<-msfactor/mserror
  finteract<-msinteract/mserror
  store<-data.frame(matrix(NA,nrow=4,ncol =7))
  colnames(store)<-c("Factor","data","SS","MS","F","P_value","Result")
  store$Factor<-c(FacA,FacB,InterFacAFacB,"Residuals")
  store$data<-c(J-1,K-1,(J-1)*(K-1),n-p)
  store$SS<-c(sssample,ssfactor,ssinteract,sserror)
  store$MS<-c(round(mssample,digits=4),round(msfactor,digits=4),round(msinteract,digits=4),round(mserror,digits=4))
  store$F<-c(round(fsample,digits = 3),round(ffactor,digits=4),round(finteract,digits = 3),"")
  P_value<-c(round(pf(fsample,J-1,n-((J)*(K)),lower.tail = FALSE),digits=7),
             round(pf(ffactor,K-1,n-((J)*(K)),lower.tail = FALSE),digits=7),
             round(pf(finteract,((J-1)*(K-1)),n-((J)*(K)),lower.tail = FALSE),digits=7))
  Result<-c(ifelse(P_value>alpha,"Not reject","Reject"))
  store$P_value<-c(P_value,"")
  store$Result<-c(Result,"")
  if(verbose){
    cat("\n","  Two-way ANOVA ","(alpha = ",alpha,")",sep="")
    cat("\n", "----------------------------------------------------------------------------", sep = "","\n")
    print(store, row.names = FALSE)
    cat("----------------------------------------------------------------------------", sep = "", "\n")
  }
  result<-list()
  result$output<-store
  result$alpha<-alpha
  result$method<-"Two-way ANOVA"
  result$data<-data
  result$formula<-formula
  attr(result, "class") <- "twt"
  invisible(result)
} 
