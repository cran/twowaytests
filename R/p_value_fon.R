p_value_fon<-function(kron_matrix,tr_means,v,h,alpha){
  R_matrix<-v%*%t(kron_matrix)%*%solve(kron_matrix%*%v%*%t(kron_matrix))%*%kron_matrix
  A_value<-sum(diag((diag(R_matrix))^2/diag(h-1)))
  crit_value<-qchisq(1-alpha,nrow(kron_matrix))
  crit_value<-crit_value+(crit_value/(2*nrow(kron_matrix)))*A_value*(1+3*crit_value/(nrow(kron_matrix)+2))
  crit_value
}
