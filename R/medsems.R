medsems<- function(data){
  x<-sort(data)
  value<-round((length(x)+1)/2-qnorm(.995)*sqrt(length(x)/4))
  value<-ifelse(value==0,1,value)
  result<-sqrt(((x[length(x)-value+1]-x[value])/(2*qnorm(.995)))^2)
  result
}
