optimalf<-function(data,na.rm=FALSE){
  if(na.rm==TRUE){
    data<-data[complete.cases(data),]
  }
  i<-floor(length(data)*0.25+5/12)
  j<-(length(data)*0.25)-i+(5/12)
  lq<-(1-j)*sort(data)[i]+j*sort(data)[i+1]
  k<-length(data)-i+1
  uq<-(1-j)*sort(data)[k]+j*sort(data)[k-1]
  list(lq=lq,uq=uq)
}
