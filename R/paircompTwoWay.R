paircompTwoWay <- function(x, ...) UseMethod("paircompTwoWay")

paircompTwoWay.default <- function(x, ...)paircompTwoWay.twt(x, ...)




paircompTwoWay.twt <- function(x, adjust.method = c("bonferroni", "holm", "hochberg", "hommel", "BH", 
                                                    "BY", "fdr", "none"), verbose = TRUE, ...){
  data<-x$data
  alpha <- x$alpha
  dp<-as.character(x$formula)
  y<-data[[dp[[2L]]]]
  group1_name<-strsplit(dp[3],split=" ")[[1]][1]
  group2_name<-strsplit(dp[3],split=" ")[[1]][3]
  id_a<-levels(data[group1_name][[1]])
  comb_a<-t(combn(id_a,2))
  comb2_a<-dim(comb_a)[1]
  id_b<-levels(data[group2_name][[1]])
  comb_b<-t(combn(id_b,2))
  comb2_b<-dim(comb_b)[1]
  
  if(x$method=="Two-way ANOVA for Modified One-step M-estimators "){
    method_num<-1
    analysis<-MestTwoWay
    mest_est<-"mom_est"
  } 
  if(x$method=="Two-way ANOVA for One-step M-estimators "){
    method_num<-1
    analysis<-MestTwoWay
    mest_est<-"onestep_est"
  }
  if(x$method=="Two-way ANOVA for Medians "){
    method_num<-1
    analysis<-MestTwoWay
    mest_est<-"median"
  }
  if(x$method=="Two-way ANOVA") analysis<-aovTwoWay; method_num<-3
  if(x$method=="Two-way ANOVA for medians") analysis<-medTwoWay; method_num<-3
  if(x$method=="Two-way ANOVA for Trimmed Means") analysis<-tmeanTwoWay; method_num<-3
  if(x$method=="Generalized p-value by PB method"){
    method_num<-2
    analysis<-gpTwoWay
    gp_method<-"gPB"
  }
  if(x$method=="Generalized p-value by GPQ method"){
    method_num<-2
    analysis<-gpTwoWay
    gp_method<-"gPQ"
  }
  adjust.method <- match.arg(adjust.method)
  
  if (adjust.method == "bonferroni") method.name = paste("Bonferroni correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "holm") method.name = paste("Holm correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "hochberg") method.name = paste("Hochberg correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "hommel") method.name = paste("Hommel correction (alpha = ",alpha,")",sep = "")
  if ((adjust.method == "fdr")|(adjust.method == "BH")) method.name = paste("Benjamini-Hochberg correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "BY") method.name = paste("Benjamini-Yekutieli correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "none") method.name = paste("No correction (alpha = ",alpha,")",sep = "")
  pval_a<-NULL
  pval_b<-NULL
  pval_ab<-NULL
  if(x$output[3,"Result"]=="Reject"){
    comb_a<-t(combn(id_a,1))
    comb2_a<-dim(comb_a)[1]
    store_ab<-data.frame(matrix(NA,nrow=comb2_a*comb2_b,ncol=5))
    k<-0
    for (i in 1:comb2_a) {
      for (j in 1:comb2_b) {
        k<-k+1
        store_ab[k,1]<-comb_a[i]
        store_ab[k,2]<-comb_b[j,1]
        store_ab[k,3]<-comb_b[j,2]
        data_sub<-data[((data[group2_name]==store_ab[k,2])|(data[group2_name]==store_ab[k,3]))&(data[group1_name]==store_ab[k,1]),]
        if(x$method=="Two-way ANOVA"){
          one_sample_test<-st.test(as.formula(sprintf("%s~%s",dp[[2L]],group2_name)),data_sub,verbose = FALSE)
        }else if((x$method=="Generalized p-value by PB method")|(x$method=="Generalized p-value by GPQ method")){
          one_sample_test<-wt.test(as.formula(sprintf("%s~%s",dp[[2L]],group2_name)),data_sub,verbose = FALSE)
        }else{
          one_sample_test<-mw.test(as.formula(sprintf("%s~%s",dp[[2L]],group2_name)),data_sub,verbose = FALSE)
        }
        store_ab[k,4]<-one_sample_test$p.value[[1]]
      }
    }
    for (m in unique(store_ab$X1)) {
      store_ab[store_ab$X1==m,]$X4<-p.adjust(store_ab[store_ab$X1==m,]$X4,method=adjust.method)
    }
    store_ab$X5 = ifelse(store_ab$X4 <= alpha, "Reject", "Not reject")
    colnames(store_ab) = c(group1_name, paste(group2_name,"(a)",sep = ""),paste(group2_name,"(b)",sep = ""), "P_value", "  No difference")
    
    if (verbose==TRUE){
      cat("\n", "",method.name,"for subgroups of each",group1_name,"level", "\n", sep = " ")
      cat("--------------------------------------------------------------------------------", "\n", sep = " ")
      print(store_ab)
      cat("--------------------------------------------------------------------------------", "\n\n", sep = " ")
    }
    invisible(store_ab)
  }else{
    if(x$output[1,"Result"]=="Reject"){
      for(i in 1:comb2_a){
        data_sub<-data[(data[group1_name]==comb_a[i,1])|(data[group1_name]==comb_a[i,2]),]
        if(method_num==1){
          analysis_<-analysis(x$formula,data_sub,estimator = mest_est,verbose = FALSE)
        }else if(method_num==2){
          analysis_<-analysis(x$formula,data_sub,method = gp_method,verbose = FALSE)
        }else{
          analysis_<-analysis(x$formula,data_sub,verbose = FALSE)
        }
        pval_a<-c(pval_a,analysis_$output[1,"P_value"])
      }
      padj_a <- p.adjust(pval_a, method = adjust.method)
      store_a = data.frame(matrix(NA, nrow = comb2_a, ncol = 4))
      
      store_a$X1=comb_a[,1]
      store_a$X2=comb_a[,2]
      store_a$X3 = padj_a
      store_a$X4 = ifelse(store_a$X3 <= alpha, "Reject", "Not reject")
      colnames(store_a) = c("Level (a)", "Level (b)", "P_value", "  No difference")
      
      if (verbose==TRUE){
        cat("\n", "",method.name,"for",group1_name, "\n", sep = " ")
        cat("-----------------------------------------------------", "\n", sep = " ")
        print(store_a)
        cat("-----------------------------------------------------", "\n\n", sep = " ")
      }
      invisible(store_a)
      
    }else{cat("\n","Pairwise comparisons could not be performed since",group1_name,"is not statistically significant.","\n")}
  }
}