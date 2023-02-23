tmeanTwoWay<-function(formula, data, tr = 0.1, alpha = 0.05, na.rm = TRUE, verbose = TRUE){
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
  n<-nrow(data)/(J*K)
  h<-matrix((n-2*floor(tr*n)),nrow = J*K,ncol = 1)
  v<-NA
  trim_means<-NA
  data_x<-groupping(data,FacA,FacB,y)
  for(t in 1:p){
    v[t]<-(n-1)*winsorvar(data_x[[t]],tr=tr)/((h[t]-1)*h[t])
    trim_means[t]<-mean(data_x[[t]],trim=tr)
  }
  v <- diag(v, p, p)
  ij <- matrix(c(rep(1, J)), 1, J)
  ik <- matrix(c(rep(1, K)), 1, K)
  cj <- diag(1, J-1, J)
  for (i in 1:(J-1)) cj[i, i + 1] <- (-1)
  ck <- diag(1, K-1, K)
  for (i in 1:K-1) ck[i, i + 1] <- (-1)
  cj_ik_mat<-kronecker(cj,ik)
  cj_ck_mat<-kronecker(cj,ck)
  ij_ck_mat<-kronecker(ij,ck)
  FacA_value<-t(trim_means)%*%t(cj_ik_mat)%*%
    solve(cj_ik_mat%*%(v)%*%t(cj_ik_mat))%*%(cj_ik_mat)%*%trim_means
  FacB_value<-t(trim_means)%*%t(ij_ck_mat)%*%
    solve(ij_ck_mat%*%(v)%*%t(ij_ck_mat))%*%(ij_ck_mat)%*%trim_means
  InterFacAFacB_value<-t(trim_means)%*%t(cj_ck_mat)%*%
    solve(cj_ck_mat%*%(v)%*%t(cj_ck_mat))%*%(cj_ck_mat)%*%trim_means
  alphas<-c(1:999)/1000
  for (i in 1:999) {
    seda <- i
    if (FacA_value > p_value_fon(cj_ik_mat,trim_means,v,h,alphas[i])) 
      break
  }
  p.val_FacA<-seda/1000
  for (i in 1:999) {
    seda <- i
    if (FacB_value > p_value_fon(ij_ck_mat,trim_means,v,h,alphas[i])) 
      break
  }
  p.val_FacB<-seda/1000
  for (i in 1:999) {
    seda <- i
    if (InterFacAFacB_value > p_value_fon(cj_ck_mat,trim_means,v,h,alphas[i])) 
      break
  }
  p.val_InterFacAFacB<-seda/1000
  store<-data.frame(matrix(NA,nrow=3,ncol = 4))
  colnames(store)<-c("Factor","Statistic","P_value","Result")
  store$Factor<-c(FacA,FacB,InterFacAFacB)
  store$Statistic<-c(FacA_value,FacB_value,InterFacAFacB_value)
  store$P_value<-c(p.val_FacA,p.val_FacB,p.val_InterFacAFacB)
  store$Result<-ifelse(store$P_value>alpha,"Not reject","Reject")
  store4<-store
  if(verbose){
    cat("\n","  Two-way ANOVA for Trimmed Means ","(alpha = ",alpha,")",sep="")
    cat("\n", "---------------------------------------------------------------------------", sep = "","\n")
    print(store4, row.names = FALSE)
    cat("---------------------------------------------------------------------------", sep = "", "\n")
  }
  result<-list()
  result$output<-store4
  result$alpha<-alpha
  result$method<-"Two-way ANOVA for Trimmed Means"
  result$data<-data
  result$formula<-formula
  attr(result, "class") <- "twt"
  invisible(result)
}
