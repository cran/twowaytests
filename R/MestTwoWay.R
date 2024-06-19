MestTwoWay <- function(formula, data, estimator=c("mom_est", "onestep_est", "median","bisquare","hampel"), nboot = 500, distance = c("mahalanobis", "projected"), seed = 123, alpha = 0.05, na.rm = TRUE, verbose = TRUE){
  set.seed(seed)
  if(exists("estimator")==FALSE) estimator<-"mom_est"
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
  
  estimator=match.arg(estimator)
  distance=match.arg(distance)
  if(estimator=="mom_est"){
    estimator_new=mom_est
    method.name<-"Two-way ANOVA for Modified One-step M-estimator"
  }else if(estimator=="onestep_est"){
    estimator_new=onestep_est
    method.name<-"Two-way ANOVA for One-step M-estimator"
  }else if(estimator=="hampel"){
    method.name<-"Two-way ANOVA for Hampel's M-estimator"
  }else if(estimator=="bisquare"){
    method.name<-"Two-way ANOVA for Tukey's Biweight (Bisquare) M-estimator"
  }else{
    estimator_new=median
    method.name<-"Two-way ANOVA for Median"
  }
  J<-nlevels(data[,FacA])
  K<-nlevels(data[,FacB])
  p<-J*K
  contmatrixs<-contrastmatrix(J,K)
  contrastA<-contmatrixs$contrastA
  contrastB<-contmatrixs$contrastB
  contrastAB<-contmatrixs$contrastAB
  data_x<-groupping(data,FacA,FacB,y)
  length_x<-list()
  estimates<-list()
  for (i in 1:length(data_x)) {
    if(estimator=="hampel"){
      model <- rlm(data_x[[i]] ~ 1, psi = psi.hampel, maxit = 200)
      estimation<-coef(model)[1]
    }else if(estimator=="bisquare"){
      model <- rlm(data_x[[i]] ~ 1, psi = psi.bisquare, maxit = 200)
      estimation<-coef(model)[1]
    }else{
      estimation<-estimator_new(data_x[[i]])
    }
    estimates<-append(estimates,estimation)
    length_x<-append(length_x,length(data_x[[i]]))
  }
  estimates<-unlist(estimates)
  length_x<-unlist(length_x)
  boot_matrix<- matrix(NA, nrow = p, ncol = nboot)
  for (j in 1:p) {
    sample_data <- matrix(sample(data_x[[j]], size =nboot*length(data_x[[j]]), replace = TRUE), nrow = nboot)
    if(estimator=="bisquare"){
      boot_matrix[j, ] <-apply(sample_data, 1, function(row){
        model_psi<-rlm(row~1, psi = psi.bisquare, maxit = 200)
        coef(model_psi)[1]
      })
    }else if(estimator=="hampel"){
      boot_matrix[j, ] <-apply(sample_data, 1, function(row){
        model_psi<-rlm(row~1, psi = psi.hampel, maxit = 200)
        coef(model_psi)[1]
      })
    }else{
      boot_matrix[j, ] <- apply(sample_data, 1, estimator_new)
    }
    
    if (sum(is.na(boot_matrix[j,])) > 0) 
      boot_matrix[j, ][which(is.na(boot_matrix[j,]))] <- mean(boot_matrix[j, ], na.rm = TRUE)
  }
  bcontA_ <- t(contrastA) %*% boot_matrix
  tvectorA_ <-(t(contrastA) %*% as.matrix(estimates))[, 1]
  tempcenA_ <- rowMeans(bcontA_)
  smatrixA_ <- var(t(bcontA_) - tempcenA_ + tvectorA_)
  bcontA_ <- rbind(t(bcontA_), rep(0, ncol(contrastA)))
  if (distance=="mahalanobis") {
    dist_matrix <- mahalanobis(bcontA_, tvectorA_, smatrixA_)
  }else{
    dist_matrix <- projected_dist(bcontA_)
  }
  p.val_FacA <- 1 - sum(dist_matrix[nboot+1] >= dist_matrix[1:nboot])/nboot
  bcontB_ <- t(contrastB) %*% boot_matrix
  tvectorB_ <- (t(contrastB) %*% as.matrix(estimates))[, 1]
  tempcenB_ <- rowMeans(bcontB_)
  smatrixB_ <- var(t(bcontB_) - tempcenB_ + tvectorB_)
  bcontB_ <- rbind(t(bcontB_), rep(0, ncol(contrastB)))
  if (distance=="mahalanobis") {
    dist_matrix <- mahalanobis(bcontB_, tvectorB_, smatrixB_)
  }else{
    dist_matrix <- projected_dist(bcontB_)
  }
  p.val_FacB <- 1 - sum(dist_matrix[nboot+1] >= dist_matrix[1:nboot])/nboot
  bcontAB_ <- t(contrastAB) %*% boot_matrix
  tvectorAB_ <- (t(contrastAB) %*% as.matrix(estimates))[, 1]
  tempcenAB_ <- rowMeans(bcontAB_)
  smatrixAB_ <- var(t(bcontAB_) - tempcenAB_ + tvectorAB_)
  bcontAB_ <- rbind(t(bcontAB_), rep(0, ncol(contrastAB)))
  if (distance=="mahalanobis") {
    dist_matrix <- mahalanobis(bcontAB_, tvectorAB_, smatrixAB_)
  }else{
    dist_matrix <- projected_dist(bcontAB_)
  }
  p.val_InterFacAFacB <- 1 - sum(dist_matrix[nboot+1] >= dist_matrix[1:nboot])/nboot
  store<-data.frame(matrix(NA,nrow=3,ncol = 3))
  colnames(store)<-c("Factor","P.value","Result")
  store$Factor<-c(FacA,FacB,InterFacAFacB)
  store$P.value<-c(p.val_FacA,p.val_FacB,p.val_InterFacAFacB)
  store$Result<-ifelse(store$P.value>alpha,"Not reject","Reject")
  store4<-store
  if(verbose){
    cat("\n","  ",method.name," ","(alpha = ",alpha,")",sep="")
    cat("\n", "----------------------------------------------------------------------------", sep = "","\n")
    print(store4, row.names = FALSE)
    cat("----------------------------------------------------------------------------", sep = "", "\n")
  }
  
  result<-list()
  result$output <- store4
  result$alpha <- alpha
  result$method <- method.name
  result$data <- data
  result$formula <- formula
  attr(result, "class") <- "twt"
  invisible(result)
}
