onestep_est<-function(data,bend_con=1.28,median=TRUE){
  data <- data[!is.na(data)]
  if(median==TRUE){
    location=median(data)
  }
  else{
    location=mom_est(data,bend_con = bend_con)
  }
  y<-(data-location)/med_abs_dev(data)
  value<-median(data)+med_abs_dev(data)*(sum(ifelse(abs(y)<=bend_con,y,bend_con*sign(y)))/length(data[abs(y) <= bend_con]))
  value
}
