formula.oglmx<-function(x, ...){
  # extract the formula for an oglmx object
  # for use to apply a model name in lrtest
  if (is.null(x$formula[[2]])){
    value<-x$formula[[1]]
  } else {
    # collect the names from the terms output
    # from the mean equation, include response term
    termsFM1<-terms(x$formula[[1]])
    termsFM2<-terms(x$formula[[2]])
    rhsvarsFM2<-attr(termsFM2,"term.labels")
    updateform<-as.formula(paste("~.|",paste(rhsvarsFM2,collapse="+")))
    value<-update.formula(x$formula[[1]],updateform)
  }
  return(value)
}
