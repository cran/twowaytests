med_abs_dev<-function(data,cen=median(data),cons=1.4826,lo_med=FALSE,hi_med=FALSE){
  data <- data[!is.na(data)]
  if(lo_med+hi_med==2){
    stop("Must be least option FALSE between 'lo_med' and 'hi_med'")
  }else{
    if(lo_med+hi_med+(length(data)%%2==0)==2){
      n <-length(data)%/%2 + ifelse(hi_med==TRUE,1,0)
      value<-sort(abs(data- cen), partial = n)[n]*cons
    }else{
      value<-median(abs(data-cen))*cons
    }
    value
  }
  
}
