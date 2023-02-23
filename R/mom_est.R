mom_est<-function(data,bend_con=2.24){
  data <- data[!is.na(data)]
  upper<-(data>median(data)+bend_con*mad(data))
  lower<-(data<median(data)-bend_con*mad(data))
  x<-rep(T,length(data))
  x[upper]<-F
  x[lower]<-F
  mom_est<-mean(data[x])
  mom_est
}
