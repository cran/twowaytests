
homogtestTwoWay <- function(formula, data, method = c("Levene", "Bartlett", "Fligner"), alpha = 0.05, na.rm = TRUE, verbose = TRUE){


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
  
  
  y = data[[y]]
  FacA_vector = as.factor(data[[FacA]])
  FacB_vector = as.factor(data[[FacB]])
  group = interaction(FacA_vector, FacB_vector)
  
  
  method = match.arg(method)
  

if (method == "Levene"){
out=leveneTest(y, group, center="mean")
if (verbose) {
            cat("\n", " Levene's Homogeneity Test", paste("(alpha = ",alpha,")",sep = ""), "\n", 
                sep = " ")
            cat("================================================", 
                "\n", sep = " ")
            cat("  formula :", dname1,"~",dname2,"*",dname3, "\n\n", sep = " ")
            cat("  statistic  :", out$F[1], "\n", sep = " ")
            cat("  num df     :", out$Df[1], "\n", sep = " ")
		 cat("  denum df   :", out$Df[2], "\n", sep = " ")
            cat("  p.value    :", out$P[1], "\n\n", sep = " ")
            cat(if (out$P[1] > alpha) {
                "  Result     : Variances are homogeneous."
            }
            else {
                "  Result     : Variances are not homogeneous."
            }, "\n")
            cat("================================================", 
                "\n\n", sep = " ")
        }

result <- list()
result$statistic <- out$F[1]
result$parameter <- out$Df
result$p.value <- out$P[1]
result$method <- "Levene's Homogeneity Test"
result

}

if (method == "Bartlett"){
out=bartlett.test(y, group)
if (verbose) {
            cat("\n", " Bartlett's Homogeneity Test", paste("(alpha = ",alpha,")",sep = ""), "\n", 
                sep = " ")
            cat("================================================", 
                "\n", sep = " ")
            cat("  formula :", dname1,"~",dname2,"*",dname3, "\n\n", sep = " ")
            cat("  statistic  :", out$statistic, "\n", sep = " ")
            cat("  parameter  :", out$parameter, "\n", sep = " ")
		 cat("  p.value    :", out$p.value, "\n\n", sep = " ")
            cat(if (out$p.value > alpha) {
                "  Result     : Variances are homogeneous."
            }
            else {
                "  Result     : Variances are not homogeneous."
            }, "\n")
            cat("================================================", 
                "\n\n", sep = " ")
        }

result <- list()
result$statistic <- as.numeric(out$statistic)
result$parameter <- as.numeric(out$parameter)
result$p.value <- out$p.value
result$method <- "Bartlett's Homogeneity Test"

}


if (method == "Fligner"){
out=fligner.test(y, group)
if (verbose) {
            cat("\n", " Fligner-Killeen Homogeneity Test", paste("(alpha = ",alpha,")",sep = ""), "\n", 
                sep = " ")
            cat("====================================================", 
                "\n", sep = " ")
            cat("  formula :", dname1,"~",dname2,"*",dname3, "\n\n", sep = " ")
            cat("  statistic  :", out$statistic, "\n", sep = " ")
            cat("  parameter  :", out$parameter, "\n", sep = " ")
		 cat("  p.value    :", out$p.value, "\n\n", sep = " ")
            cat(if (out$p.value > alpha) {
                "  Result     : Variances are homogeneous."
            }
            else {
                "  Result     : Variances are not homogeneous."
            }, "\n")
            cat("====================================================", 
                "\n\n", sep = " ")
        }

result <- list()
result$statistic <- as.numeric(out$statistic)
result$parameter <- as.numeric(out$parameter)
result$p.value <- out$p.value
result$method <- "Fligner-Killeen Homogeneity Test"

}

invisible(result)
}