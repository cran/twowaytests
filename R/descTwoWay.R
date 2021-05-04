descTwoWay<-function(formula, data) 
{
  data <- model.frame(formula, data)
  fml <-as.character(formula)
  ftmp <- strsplit(fml,"~")
  y <- as.vector(ftmp[[2]])
  Factors <- strsplit(ftmp[[3]][1],"[*]")[[1]]
  FacA <- strsplit(Factors[1], " ")[[1]][1] #Drop spaces
  FacB <- strsplit(Factors[2], " ")[[1]][2]
  InterFacAFacB <- paste(y,"~",FacB, sep = "")

  if (!is.data.frame(data)) stop("Data must be in data.frame class.")
  if(length(Factors)!=2) stop("Please correct the RHS of the formula. Formula must include two factors.")
  if(!is.factor(data[,colnames(data)==FacA])) stop(paste(FacA, "must be a factor."))
  if(!is.factor(data[,colnames(data)==FacB])) stop(paste(FacB, "must be a factor."))
  if(!is.numeric(data[,colnames(data)==y])) stop(paste(y, "must be a numeric."))
  
  FacA_levels <- levels(data[,colnames(data)==FacA])
  
  out <- out2 <- NULL
  for (i in FacA_levels){
    out <- describe(as.formula(noquote(InterFacAFacB)), data = data[data[,colnames(data)==FacA]==i,])
    Factors <- cbind(rep(i,dim(out)[1]),rownames(out))
    colnames(Factors) <- c(FacA,FacB)
    out2 <- rbind(out2, cbind(Factors,out))
    out <- NULL
  }
  rownames(out2) <- NULL
  return(out2)
  }