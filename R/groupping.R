groupping<-function(data_,factor1,factor2,variable){
  data_g<-list()
  k<-0
  for (i in levels(data_[,factor1])) {
    for (j in levels(data_[,factor2])) {
      k<-k+1
      data_g[[k]]<-data_[data_[,factor1]==i&data_[,factor2]==j,variable]
    }
  }
  data_g
}
