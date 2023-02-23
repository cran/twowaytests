contrastmatrix<-function(j,k){
  contrast_A<-matrix(0,nrow=(j*k),ncol=((j^2-j)/2))
  contrast_B<-matrix(0,nrow=(j*k),ncol=((k^2-k)/2))
  contrast_AB<-matrix(0,nrow = (j*k),ncol=(((j^2-j)/2)*((k^2-k)/2)))
  t<-0
  for (i in 1:j) {
    for (ik in 1:j) {
      if(i<ik){
        t<-t+1
        empty_mat<-matrix(0,nrow=j,ncol=k)
        empty_mat[i,]<-1
        empty_mat[ik,]<-(-1)
        contrast_A[,t]<-t(empty_mat)
      }
    }
  }
  t<-0
  for (n in 1:k) {
    for (nn in 1:k) {
      if(n<nn){
        t<-t+1
        empty_mat<-matrix(0,nrow=j,ncol=k)
        empty_mat[,n]<-1
        empty_mat[,nn]<-(-1)
        contrast_B[,t]<-t(empty_mat)
      }
    }
  }
  t<-0
  for(i in 1:j){
    for(ik in 1:j){
      if(i < ik){
        for(n in 1:k){
          for(nn in 1:k){
            if(n<nn){
              t<-t+1
              empty_mat<-matrix(0,nrow=j,ncol=k)
              empty_mat[i,n]<-1
              empty_mat[i,nn]<-(-1)
              empty_mat[ik,n]<-(-1)
              empty_mat[ik,nn]<-1
            }
            contrast_AB[,t]<-t(empty_mat)
          }}}}}
  list(contrastA=contrast_A,contrastB=contrast_B,contrastAB=contrast_AB)
}
