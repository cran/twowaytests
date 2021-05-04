gplotTwoWay <- function (formula, data, type = c("errorbar", "boxplot"), xlab = NULL, ylab = NULL, title = NULL, legend.title = NULL, width = NULL, option = c("se", "sd"), na.rm = TRUE){
    
  data <- model.frame(formula, data)
  fml <-as.character(formula)
  ftmp <- strsplit(fml,"~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1],"[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1] #Drop spaces
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  InterFacAFacB <- paste(y,"~",FacB, sep = "")
  

  dname1<-y
  dname2<-FacA
  dname3<-FacB
  
  
  if (!is.data.frame(data)) stop("Data must be in data.frame class.")
  if(length(Factors)!=2) stop("Please correct the RHS of the formula. Formula must include two factors.")
  if(!is.factor(data[,colnames(data)==FacA])) stop(paste(FacA, "must be a factor."))
  if(!is.factor(data[,colnames(data)==FacB])) stop(paste(FacB, "must be a factor."))
  if(!is.numeric(data[,colnames(data)==y])) stop(paste(y, "must be a numeric."))
  
  if (na.rm){
    completeObs <- complete.cases(data)
    data <- data[completeObs,]
  }
  
  
  y_vector = data[[y]]
  FacA_vector = data[[FacA]]
  FacB_vector = data[[FacB]]
  data<-as.data.frame(cbind(y_vector, FacA_vector, FacB_vector))
  data$FacA_vector <- as.factor(FacA_vector)
  data$FacB_vector <- as.factor(FacB_vector)
  type = match.arg(type)
  
  if (type == "boxplot"){
    
    
    if (is.null(width)) width <- 0.75 else width <- width
    
    
    out <- ggplot(data, aes(x=FacA_vector, y=y_vector, fill=FacB_vector))
    out<-out + geom_boxplot(width = width)
    
  }
  
  if (type == "errorbar"){
  option = match.arg(option)
  Extremity<-dose<-dev<-NULL
  colnames(data)<-c(y,FacA,FacB)
  result <- descTwoWay(formula = formula, data = data)
  df2<-data.frame(matrix(0,dim(result)[1],4))
  colnames(df2)<-c("Extremity","dose","mean","dev")
  
  df2$Extremity<-result[,1]
  df2$dose<-result[,2]
  df2$mean<-result[,4]
  
  if (option=="se") df2$dev <- result[,5]/sqrt(result[,3])
  if (option=="sd") df2$dev <- result[,5]
  

  
  if (is.null(width)) width <- 0.2 else width <- width
  
  out<- ggplot(df2, aes(x=Extremity, y=mean, fill=dose)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-dev, ymax=mean+dev), width = width, size = 0.8,position=position_dodge(.9))
  }
  
  if (is.null(ylab)) y.name<-dname1 else y.name<-ylab
  if (is.null(xlab)) x.name<-dname2 else x.name<-xlab
  if (is.null(legend.title)) legend.title<-dname3 else legend.title<-legend.title
  if (is.null(title)) title.name<-"" else title.name<-title
  
  out <- out + labs(x=x.name, y=y.name, title=title.name, fill=legend.title)
  
  return(out)
}
