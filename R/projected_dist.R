projected_dist<-function(data,na.rm=TRUE){
  if(na.rm==TRUE){
    data<-data[rowSums(is.na(data))!=1,]
  }
  if(is.null(ncol(data))){
    cent<-median(data)
    dist<-abs(data-cent)
    project_dist<-dist/(optimalf(dist)$uq-optimalf(dist)$lq)
  }else{
    cent<-apply(data,2,median)
    matriks<-matrix(ncol=nrow(data),nrow=nrow(data))
    for (i in 1:nrow(data)){
      b<-data[i,]-cent
      dist<-NA
      sumsquare_b<-sum(b^2)
      if(sumsquare_b!=0){
        for (j in 1:nrow(data)){
          a<-data[j,]-cent
          tempt<-sum(a*b)*b/b^2
          dist[j]<-sqrt(sum(tempt^2))
        }
        tempt<-optimalf(dist)
        matriks[,i]<-dist/(tempt$uq-tempt$lq)
        
      }}
    project_dist<-apply(matriks,1,max,na.rm=TRUE)
  }
  
  project_dist
}
