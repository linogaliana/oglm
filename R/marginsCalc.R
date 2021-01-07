

getEtas.Exp<-function(thresholds,xb_matrix,sd_matrix){
  eta1numerator<-apply(xb_matrix,2,function(x){thresholds[2]-x})
  eta0numerator<-apply(xb_matrix,2,function(x){thresholds[1]-x})
  etas1<-eta1numerator/sd_matrix
  etas0<-eta0numerator/sd_matrix
  return(list(etas1,etas0))
}

continuous.margin.mean<-function(paramvec,etas,link,std.dev){
    margineffect<- scoreMean(etas[[1]],etas[[2]],std.dev,1,link)%*%paramvec
    return(margineffect)
}

D_continuous.margin.mean_mean<-function(whichMargins,whichXest,X,paramvec,etas,link,std.dev){
  firstterm<-scoreMean(etas[[1]],etas[[2]],std.dev,1,link)
  ftmatrix<-t(sapply(whichMargins,function(x){as.numeric(whichXest %in% x)}))
  d_margineffect1<- ftmatrix*mean(firstterm)
  #  return((ProbFuncDD(etas[[1]])-ProbFuncDD(etas[[2]]))/(std.dev^2))
  secondterm<-X*as.vector(hessMean_Mean(etas[[1]],etas[[2]],std.dev,1,link))
  d_margineffect2<- (matrix(paramvec,ncol=nrow(X),nrow=length(paramvec))%*%secondterm)/length(etas[[1]])
  #return(list(d_margineffect1,d_margineffect2))
  d_margineffect<- d_margineffect1+d_margineffect2
  return(d_margineffect)
}


D_continuous.margin.mean_var<-function(Z,paramvec,etas,link,std.dev,gstd.dev){
  d_margineffect<-(matrix(paramvec,ncol=nrow(Z),nrow=length(paramvec))%*%(Z*as.vector(hessMean_Var(etas[[1]],etas[[2]],std.dev,gstd.dev,1,link))))/length(etas[[1]])
  return(d_margineffect)
}

D_continuous.margin.mean_alpha<-function(estThresh,outcomematrix,paramvec,etas,link,std.dev){
  d_margineffect<- (matrix(paramvec,ncol=length(etas[[1]]),nrow=length(paramvec))%*%hessMean_Thresh(estThresh,outcomematrix,etas[[1]],etas[[2]],std.dev,1,link))/length(etas[[1]])
  return(d_margineffect)
}

continuous.margin.sd<-function(paramvec,etas,link,std.dev,gstd.dev){
  margineffect<- scoreVar(etas[[1]],etas[[2]],std.dev,gstd.dev,1,link)%*%paramvec
  return(margineffect)
}

D_continuous.margin.var_mean<-function(X,paramvec,etas,link,std.dev,gstd.dev){
  d_margineffect<-(matrix(paramvec,ncol=nrow(X),nrow=length(paramvec))%*%(X*as.vector(hessMean_Var(etas[[1]],etas[[2]],std.dev,gstd.dev,1,link))))/length(etas[[1]])
  return(d_margineffect)
}


D_continuous.margin.var_var<-function(whichMargins,whichZest,Z,paramvec,etas,link,std.dev,gstd.dev,hstd.dev){
  firstterm<-scoreVar(etas[[1]],etas[[2]],std.dev,gstd.dev,1,link)
  ftmatrix<-t(sapply(whichMargins,function(x){as.numeric(whichZest %in% x)}))
  d_margineffect1<- -ftmatrix*mean(firstterm)
  secondterm<-Z*as.vector(hessVar_Var(etas[[1]],etas[[2]],std.dev,gstd.dev,hstd.dev,1,link))
  d_margineffect2<- (matrix(paramvec,ncol=nrow(Z),nrow=length(paramvec))%*%secondterm)/length(etas[[1]])
  d_margineffect<- d_margineffect1+d_margineffect2
  return(d_margineffect)
}

D_continuous.margin.var_alpha<-function(estThresh,outcomematrix,paramvec,etas,link,std.dev,gstd.dev){
  d_margineffect<-(matrix(paramvec,ncol=length(etas[[1]]),nrow=length(paramvec))%*%hessVar_Thresh(estThresh,outcomematrix,etas[[1]],etas[[2]],std.dev,gstd.dev,1,link))/length(etas[[1]])
  return(d_margineffect)
}

discrete.margin_meanonly<-function(beta,X,whichVars,etas,link,std.dev){
  ProbFunc<-.cdf.func(link)
  XBmarginremoved<-t(apply(X[,whichVars,drop=FALSE],1,function(x){x*beta[whichVars]}))/as.vector(std.dev)
  etas1Untreated<-as.vector(etas[[1]])+XBmarginremoved
  etas0Untreated<-as.vector(etas[[2]])+XBmarginremoved
  TreatEffects<-t(beta[whichVars]%*%t(std.dev))
  etas1Treated<-etas1Untreated-TreatEffects
  etas0Treated<-etas0Untreated-TreatEffects
  margineffect<-ProbFunc(etas1Treated)-ProbFunc(etas0Treated)-(ProbFunc(etas1Untreated)-ProbFunc(etas0Untreated))
  margineffect<-apply(margineffect,2,mean)
  attr(margineffect,"inputetas")<-list(etas1Treated=etas1Treated,etas0Treated=etas0Treated,etas1Untreated=etas1Untreated,etas0Untreated=etas0Untreated)
  return(margineffect)
}

discrete.margin_varonly<-function(delta,Z,whichVars,sdmodel,etas,link,std.dev){
  ProbFunc<-.cdf.func(link)
  ZDmarginremoved<-as.vector(Z%*%delta)-t(apply(Z[,whichVars,drop=FALSE],1,function(x){x*delta[whichVars]}))
  ZDUntreated<-ZDmarginremoved
  StdUntreated<-eval({z<-ZDUntreated;sdmodel})
  TreatEffects<-matrix(delta[whichVars],ncol=length(whichVars),nrow=nrow(ZDmarginremoved),byrow = TRUE)
  ZDTreated<-ZDmarginremoved+TreatEffects
  StdTreated<-eval({z<-ZDTreated;sdmodel})
  etas1Untreated<-as.vector(etas[[1]])*as.vector(std.dev)/StdUntreated
  etas0Untreated<-as.vector(etas[[2]])*as.vector(std.dev)/StdUntreated
  etas1Treated<-as.vector(etas[[1]])*as.vector(std.dev)/StdTreated
  etas0Treated<-as.vector(etas[[2]])*as.vector(std.dev)/StdTreated
  margineffect<-ProbFunc(etas1Treated)-ProbFunc(etas0Treated)-(ProbFunc(etas1Untreated)-ProbFunc(etas0Untreated))
  margineffect<-apply(margineffect,2,mean)
  attr(margineffect,"inputetas")<-list(etas1Treated=etas1Treated,etas0Treated=etas0Treated,etas1Untreated=etas1Untreated,etas0Untreated=etas0Untreated)
  attr(margineffect,"ZDbinaries")<-list(ZDTreated=ZDTreated,ZDUntreated=ZDUntreated)
  attr(margineffect,"StdDevs")<-list(StdTreated=StdTreated,StdUntreated=StdUntreated)
  return(margineffect)
}

discrete.margin_both<-function(beta,X,delta,Z,BothEqLocs,sdmodel,etas,link,std.dev){
  ProbFunc<-.cdf.func(link)
  XBmarginremoved<-t(apply(X[,BothEqLocs$meanandvarLOC,drop=FALSE],1,function(x){x*beta[BothEqLocs$meanandvarLOC]}))
  etas1NumUntreated<-as.vector(etas[[1]])*as.vector(std.dev)+XBmarginremoved
  etas0NumUntreated<-as.vector(etas[[2]])*as.vector(std.dev)+XBmarginremoved
  TreatEffectsMean<-matrix(beta[BothEqLocs$meanandvarLOC],ncol=length(BothEqLocs$meanandvarLOC),nrow=nrow(XBmarginremoved),byrow = TRUE)
  etas1NumTreated<-etas1NumUntreated-TreatEffectsMean
  etas0NumTreated<-etas0NumUntreated-TreatEffectsMean
  ZDmarginremoved<-as.vector(Z%*%delta)-t(apply(Z[,BothEqLocs$meanandvarLOCZ,drop=FALSE],1,function(x){x*delta[BothEqLocs$meanandvarLOCZ]}))
  ZDUntreated<-ZDmarginremoved
  StdUntreated<-eval({z<-ZDmarginremoved;sdmodel})
  TreatEffectsVar<-matrix(delta[BothEqLocs$meanandvarLOCZ],ncol=length(BothEqLocs$meanandvarLOCZ),nrow=nrow(ZDmarginremoved),byrow = TRUE)
  ZDTreated<-ZDmarginremoved+TreatEffectsVar
  StdTreated<-eval({z<-ZDTreated;sdmodel})
  etas1Untreated<-etas1NumUntreated/StdUntreated
  etas0Untreated<-etas0NumUntreated/StdUntreated
  etas1Treated<-etas1NumTreated/StdTreated
  etas0Treated<-etas0NumTreated/StdTreated
  margineffect<-ProbFunc(etas1Treated)-ProbFunc(etas0Treated)-(ProbFunc(etas1Untreated)-ProbFunc(etas0Untreated))
  margineffect<-apply(margineffect,2,mean)
  attr(margineffect,"inputetas")<-list(etas1Treated=etas1Treated,etas0Treated=etas0Treated,etas1Untreated=etas1Untreated,etas0Untreated=etas0Untreated)
  attr(margineffect,"ZDbinaries")<-list(ZDTreated=ZDTreated,ZDUntreated=ZDUntreated)
  attr(margineffect,"StdDevs")<-list(StdTreated=StdTreated,StdUntreated=StdUntreated)
  return(margineffect)
}

D_discrete.margin_meanonly.mean<-function(whichVars,whichXest,X,fouretas,link,std.dev){
  # whichVars says which elements of the beta vector correspond to binary variables
  # whichXest says which elements of the beta vector are estimated
  TreatScore<-scoreMean(fouretas[[1]],fouretas[[2]],as.vector(std.dev),1,link)
  UntreatScore<-scoreMean(fouretas[[3]],fouretas[[4]],as.vector(std.dev),1,link)
  UntreatedX<-lapply(whichVars,function(x){X[,x]<-0;X})
  TreatedX<-lapply(whichVars,function(x){X[,x]<-1;X})
  d_margineffect<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(TreatedX[[x]][,whichXest,drop=FALSE]*TreatScore[,x]-UntreatedX[[x]][,whichXest,drop=FALSE]*UntreatScore[,x],2,mean)}))
  return(d_margineffect)
}

D_discrete.margin_mean.var<-function(whichZest,Z,fouretas,link,std.dev,gstd.dev){
  TreatScore<-scoreVar(fouretas[[1]],fouretas[[2]],as.vector(std.dev),as.vector(gstd.dev),1,link)
  UntreatScore<-scoreVar(fouretas[[3]],fouretas[[4]],as.vector(std.dev),as.vector(gstd.dev),1,link)
  d_margineffect<-t(sapply(c(1:ncol(TreatScore)),function(x){apply(Z[,whichZest,drop=FALSE]*(TreatScore[,x]-UntreatScore[,x]),2,mean)}))
  return(d_margineffect)
}

D_discrete.margin_mean.alpha<-function(estThresh,outcomematrix,fouretas,std.dev,link){
  TreatScore<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[1]][,x],fouretas[[2]][,x],as.vector(std.dev),1,link),2,mean)}))
  UntreatScore<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[3]][,x],fouretas[[4]][,x],as.vector(std.dev),1,link),2,mean)}))
  d_margineffect<-TreatScore-UntreatScore
}


D_discrete.margin_var.mean<-function(whichXest,X,fouretas,link,StdDevs){
  TreatScore<-scoreMean(fouretas[[1]],fouretas[[2]],StdDevs[[1]],1,link)
  UntreatScore<-scoreMean(fouretas[[3]],fouretas[[4]],StdDevs[[2]],1,link)
  d_margineffect<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(X[,whichXest,drop=FALSE]*TreatScore[,x]-X[,whichXest,drop=FALSE]*UntreatScore[,x],2,mean)}))
  return(d_margineffect)
}

D_discrete.margin_varonly.var<-function(whichVars,whichZest,Z,fouretas,ZDinputs,link,StdDevs,gsdmodel){
  GStdTreated<-eval({z<-ZDinputs[[1]];gsdmodel})
  GStdUntreated<-eval({z<-ZDinputs[[2]];gsdmodel})
  TreatScore<-scoreVar(fouretas[[1]],fouretas[[2]],StdDevs[[1]],GStdTreated,1,link)
  UntreatScore<-scoreVar(fouretas[[3]],fouretas[[4]],StdDevs[[2]],GStdUntreated,1,link)
  UntreatedZ<-lapply(whichVars,function(x){Z[,x]<-0;Z})
  TreatedZ<-lapply(whichVars,function(x){Z[,x]<-1;Z})
  d_margineffect<-t(sapply(c(1:length(whichVars)),function(x){apply(TreatedZ[[x]][,whichZest,drop=FALSE]*TreatScore[,x]-UntreatedZ[[x]][,whichZest,drop=FALSE]*UntreatScore[,x],2,mean)}))
  return(d_margineffect)
}

D_discrete.margin_var.alpha<-function(estThresh,outcomematrix,fouretas,StdDevs,link){
  if (sum(estThresh)>1){
    TreatScore<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[1]][,x],fouretas[[2]][,x],StdDevs[[1]][,x],1,link),2,mean)}))
    UntreatScore<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[3]][,x],fouretas[[4]][,x],StdDevs[[2]][,x],1,link),2,mean)}))
  } else {
    TreatScore<-sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[1]][,x],fouretas[[2]][,x],StdDevs[[1]][,x],1,link),2,mean)})
    UntreatScore<-sapply(c(1:ncol(fouretas[[1]])),function(x){apply(scoreThresh(estThresh,outcomematrix,fouretas[[3]][,x],fouretas[[4]][,x],StdDevs[[2]][,x],1,link),2,mean)})
  }
  d_margineffect<-TreatScore-UntreatScore
}

D_discrete.margin_meanvar.mean<-function(whichXest,X,BothEqLocs,fouretas,StdDevs,link){
  TreatScore<-scoreMean(fouretas[[1]],fouretas[[2]],StdDevs[[1]],1,link)
  UntreatScore<-scoreMean(fouretas[[3]],fouretas[[4]],StdDevs[[2]],1,link)
  whichVars<-BothEqLocs$meanandvarLOC
  UntreatedX<-lapply(whichVars,function(x){X[,x]<-0;X})
  TreatedX<-lapply(whichVars,function(x){X[,x]<-1;X})
  d_margineffect<-t(sapply(c(1:ncol(fouretas[[1]])),function(x){apply(TreatedX[[x]][,whichXest,drop=FALSE]*TreatScore[,x]-UntreatedX[[x]][,whichXest,drop=FALSE]*UntreatScore[,x],2,mean)}))
  return(d_margineffect)
}

D_discrete.margin_meanvar.var<-function(whichZest,Z,BothEqLocs,fouretas,ZDinputs,link,StdDevs,gsdmodel){
  GStdTreated<-eval({z<-ZDinputs[[1]];gsdmodel})
  GStdUntreated<-eval({z<-ZDinputs[[2]];gsdmodel})
  TreatScore<-scoreVar(fouretas[[1]],fouretas[[2]],StdDevs[[1]],GStdTreated,1,link)
  UntreatScore<-scoreVar(fouretas[[3]],fouretas[[4]],StdDevs[[2]],GStdUntreated,1,link)
  whichVars<-BothEqLocs$meanandvarLOCZ
  UntreatedZ<-lapply(whichVars,function(x){Z[,x]<-0;Z})
  TreatedZ<-lapply(whichVars,function(x){Z[,x]<-1;Z})
  d_margineffect<-t(sapply(c(1:length(whichVars)),function(x){apply(TreatedZ[[x]][,whichZest,drop=FALSE]*TreatScore[,x]-UntreatedZ[[x]][,whichZest,drop=FALSE]*UntreatScore[,x],2,mean)}))
  return(d_margineffect)
}

calcMEstdErrors<-function(derivME,estHess){
  StdErrors<-sapply(c(1:nrow(derivME)),function(x){(derivME[x, ,drop=FALSE]%*%estHess%*%t(derivME[x, ,drop=FALSE]))^0.5})
}
