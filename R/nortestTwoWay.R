#nortestTwoWay <- function(formula, data, method = c("SW", "SF", "LT", "AD", "CVM", "PT"), 
#                          alpha = 0.05, plot = c("qqplot", "histogram"), 
#                          mfrow = NULL, na.rm = TRUE, verbose = TRUE){ 

nortestTwoWay <- function(formula, data, method = c("SW", "SF", "LT", "AD", "CVM", "PT"), 
                          alpha = 0.05, plot = c("qqplot", "histogram"), 
                          na.rm = TRUE, verbose = TRUE){ 


  data <- model.frame(formula, data)
  fml <-as.character(formula)
  ftmp <- strsplit(fml,"~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1],"[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1] #Drop spaces
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  InterFacAFacB <- paste(y,"~",FacB, sep = "")
  
  dname1<-y
  dname2<-FacA
  dname3<-FacB
  
  
  if (!is.data.frame(data)) stop("Data must be in data.frame class.")
  
  if(length(Factors)!=2) stop("Please correct the RHS of the formula. Formula must include two factors.")
  
  if (na.rm){
    completeObs <- complete.cases(data)
    data <- data[completeObs,]
  }
  
  
  y_vector = data[[y]]
  FacA_vector = as.factor(data[[FacA]])
  FacB_vector = as.factor(data[[FacB]])
  
  method = match.arg(method)
  
  
  FacA_levels <- levels(FacA_vector)
  FacB_levels <- levels(FacB_vector)
  
  stor1 <- stor2 <- stor3 <- stor4 <- NULL
  
  if (method == "SW") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], shapiro.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Shapiro-Wilk Normality Test"
  }

  if (method == "SF") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], sf.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Shapiro-Francia Normality Test"
  }
  
  if (method == "LT") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], lillie.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Lilliefors (Kolmogorov-Smirnov) Normality Test"
  }
  
  if (method == "AD") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], ad.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Anderson-Darling Normality Test"
  }
  
  if (method == "CVM") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], cvm.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Cramer-von Mises Normality Test"
  }
  
  if (method == "PT") {
    for (i in FacA_levels) {
      kk <- tapply(y_vector[FacA_vector==i], FacB_vector[FacA_vector==i], pearson.test)
      stor1 <- c(stor1, unlist(lapply(kk, function(x) as.numeric(x$statistic))))
      stor2 <- c(stor2, unlist(lapply(kk, function(x) as.numeric(x$p))))
      stor3 <- c(stor3, rep(i, length(FacB_levels)))
      stor4 <- c(stor4, names(kk))
    }
    method.name = "Pearson Chi-square Normality Test"
  }
  
  store = data.frame(matrix(NA, nrow = length(FacA_levels)*length(FacB_levels), ncol = 5))
  store$X1 = stor3
  store$X2 = stor4
  store$X3 = stor1
  store$X4 = stor2
  store$X5 = ifelse(store$X4 > alpha, "Not reject", "Reject")
  colnames(store) = c(dname2, dname3, "Statistic", "p.value", "  Normality")
  
  formula.text<-
  
  if (verbose) {
    cat("\n", "",method.name, paste("(alpha = ",alpha,")",sep = ""), "\n", sep = " ")
    cat("==============================================================================", 
        "\n", sep = " ")
    formula.text<-format(formula)
    cat("  formula :", formula.text, "\n\n", sep = " ")
    
    print(store)
    cat("==============================================================================", 
        "\n\n", sep = " ")
  }
  
  
  
  if (!is.null(plot)){
    plot = match.arg(plot)
    
    par(mfrow = c(length(FacA_levels),length(FacB_levels)))
    
    if (plot == "qqplot") {
      
      for (i in 1:length(FacA_levels)) {
        for (j in 1:length(FacB_levels)){
        qqnorm(y_vector[which(FacA_vector == (FacA_levels[i]) & FacB_vector == (FacB_levels[j]))], main = paste(FacA_levels[i], "&", FacB_levels[j]))
        qqline(y_vector[which(FacA_vector == (FacA_levels[i]) & FacB_vector == (FacB_levels[j]))])
        }}
      }
    
    
    if (plot == "histogram") {
      

      for (i in 1:length(FacA_levels)) {
        for (j in 1:length(FacB_levels)){
        hist(y_vector[which(FacA_vector == (FacA_levels[i]) & FacB_vector == (FacB_levels[j]))], xlab = paste(FacA_levels[i], "&", FacB_levels[j]), freq = FALSE, main = NULL)
        x<-NULL
        rm(x)
        curve(dnorm(x, mean = mean(y_vector[which(FacA_vector == (FacA_levels[i]) & FacB_vector == (FacB_levels[j]))]), sd = sd(y_vector[which(FacA_vector == (FacA_levels[i]) & FacB_vector == (FacB_levels[j]))])), col = "red", add = TRUE)
      }}
  }
  
  }
  
  invisible(store) 
}
