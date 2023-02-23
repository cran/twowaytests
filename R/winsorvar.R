winsorvar <- function(data, tr = 0.2, na.rm = FALSE){
  if(na.rm){
    data<-data[complete.cases(data),]
  }
  data<-sort(data)
  xbottom<-floor(tr*length(data))+1
  xtop<-length(data)-xbottom+1
  data[data<=data[xbottom]]<-data[xbottom]
  data[data>=data[xtop]]<-data[xtop]
  var(data)
}
