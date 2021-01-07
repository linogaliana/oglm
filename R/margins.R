margins.oglmx<-function(object,atmeans=TRUE,AME=FALSE,location=NULL,outcomes="All",
                        ascontinuous=FALSE,Vars=NULL){
  # given object of class oglmx return marginal effects.
  # default return gives the marginal effect evaluated at the mean of the RHS variables (atmeans=TRUE)
  # AME=TRUE calculates the marginal effect at each data point in the sample and averages output, requires
  # the data frames to be returned as output in the call to oglmx
  # location is NULL or a numeric vector of length equal to the number of RHS variables, allows the user to
  # obtain marginal effects for any possible variable combination, rather than just the mean
  # Vars allows the user to specify which variables marginal effects are sought
  # outcomes indicates which outcome the marginal effect is sought


  #object<-results.oprobhet1
  #atmeans=TRUE
  #AME=FALSE
  #location=NULL
  #outcomes="All"
  #ascontinuous=FALSE
  #dummyzero=FALSE
  #Vars=NULL

  if (outcomes=="All"){
    outcomeMat<-diag(object$NOutcomes)
    outcomeNames<-object$Outcomes
  } else {
    outcomeMat<-diag(object$NOutcomes)[match(outcomes,object$outcomes,0L),]
    if (nrow(outcomeMat)==0){
      stop("Invalid Outcome Requested")
    }
    outcomeNames<-outcomes
  }

  beta<-object$allparams$beta
  delta<-object$allparams$delta

  paramtypes<-attr(object$coefficients,"coefftypes")

  whichXest<-object$Est.Parameters[[1]]
  whichZest<-object$Est.Parameters[[2]]
  whichAlphaest<-object$Est.Parameters[[3]]

  if (is.null(Vars)){
    betalocs<-c(1:length(beta))
    if (object$Hetero){
      deltalocs<-c(1:length(delta))
    }
  } else {
    betalocs<-match(Vars,names(object$varMeans[[1]]))
    if (object$Hetero){
      deltalocs<-match(Vars,names(object$varMeans[[2]]))
    }
  }
  betalocs<-betalocs[names(object$varMeans[[1]])[betalocs]!="(Intercept)"]
  if (object$Hetero){deltalocs<-deltalocs[names(object$varMeans[[2]])[deltalocs]!="(Intercept)"]}

  thresholds<-getThresholds(outcomeMat,object$allparams$threshparam)

  if (object$Hetero){
    gsdmodel<-D(object$sdmodel,"z")
    hsdmodel<-D(gsdmodel,"z")
  }

  if (atmeans){
    X<-matrix(object$varMeans[[1]],nrow=1)
    Z<-matrix(object$varMeans[[2]],nrow=1)
    MeanValues<-X%*%beta
    Std.Values<-eval({z<-Z%*%delta;object$sdmodel})
    if (object$Hetero){
      GStd.Values<-eval({z<-Z%*%delta;gsdmodel})
      HStd.Values<-eval({z<-Z%*%delta;hsdmodel})
    }
  } else if (!is.null(location)){
    if (!object$Hetero){
      X<-matrix(location,nrow = 1)
      MeanValues<-X%*%beta
      Std.Values<-eval({z<-1%*%delta;object$sdmodel})
    } else {
      if (length(location)!=2){
        stop("Provide location as a list with two elements: a vector of variable values for the mean equation \n
             and a vector of values for the variance equation.")
      }
      X<-matrix(location[[1]],nrow=1)
      Z<-matrix(location[[2]],nrow=1)
      MeanValues<-X%*%object$allparams$beta
      Std.Values<-eval({z<-Z%*%delta;object$sdmodel})
      GStd.Values<-eval({z<-Z%*%delta;gsdmodel})
      HStd.Values<-eval({z<-Z%*%delta;hsdmodel})
    }
  }

  if (AME){
    X<-object$modelframes$X
    Z<-object$modelframes$Z
    MeanValues<-X%*%beta
    Std.Values<-eval({z<-Z%*%delta;object$sdmodel})
    if (object$Hetero){
      GStd.Values<-eval({z<-Z%*%delta;gsdmodel})
      HStd.Values<-eval({z<-Z%*%delta;hsdmodel})
    }
  }

  if (sum(paramtypes[[3]])>0){
    outcomeMat<-lapply(c(1:nrow(thresholds)),function(x){matrix(outcomeMat[x,],nrow=length(MeanValues),ncol=ncol(outcomeMat),byrow = TRUE)})
  }

  etas<-lapply(c(1:nrow(thresholds)),function(x){getEtas.Exp(thresholds[x,],MeanValues,Std.Values)})

  # need to find the locations of the beta coefficients that will be used in the continuous calculations
  if (!ascontinuous){
    betalocscon<-betalocs[!(betalocs %in% c(1:length(beta))[object$varBinary[[1]]])]
    betalocsbin<-betalocs[(betalocs %in% c(1:length(beta))[object$varBinary[[1]]])]
    # if (length(object$factorvars)>0){
    #   betafactorlocs<-lapply(c(1:length(object$factorvars)),function(x){match(paste(object$factorvars[x],attr(object$factorvars,"levels")[[x]],sep=""),names(object$varMeans[[1]]),0)})
    # } else {
    #   betafactorlocs<-NULL
    # }

    if (object$Hetero){
      deltalocscon<-deltalocs[!(deltalocs %in% c(1:length(delta))[object$varBinary[[2]]])]
      deltalocsbin<-deltalocs[(deltalocs %in% c(1:length(delta))[object$varBinary[[2]]])]
      # if (length(object$factorvars)>0){
      #   deltafactorlocs<-lapply(c(1:length(object$factorvars)),function(x){match(paste(object$factorvars[x],attr(object$factorvars,"levels")[[x]],sep=""),names(object$varMeans[[2]]),0)})
      # } else {
      #   deltafactorlocs<-NULL
      # }
    }
  } else {
    betalocscon<-betalocs
    betalocsbin<-numeric()
    if (object$Hetero){
      deltalocscon<-deltalocs
      deltalocsbin<-numeric()
    }
  }

  betacont<-beta[betalocscon]
  betabin<-beta[betalocsbin]
  whichXest<-c(1:length(whichXest))[whichXest]
  if (object$Hetero){
    deltacont<-delta[deltalocscon]
    deltabin<-delta[deltalocsbin]
    whichZest<-c(1:length(whichZest))[whichZest]
  }


  if (length(betacont)>0){
    meanmargin<-lapply(c(1:length(etas)),function(x){continuous.margin.mean(betacont,etas[[x]],object$link,Std.Values)})
    marginsmean_cont<-lapply(meanmargin,function(x){apply(x,2,mean)})

    if (sum(paramtypes[[1]])>0){
      d_meanmargin_db<-lapply(etas,function(x){D_continuous.margin.mean_mean(betalocscon,whichXest,X[,whichXest,drop=FALSE],betacont,x,object$link,Std.Values)})
    } else {d_meanmargin_db<-NULL}
    if (sum(paramtypes[[2]])>0){
      d_meanmargin_dd<-lapply(etas,function(x){D_continuous.margin.mean_var(Z[,whichZest,drop=FALSE],betacont,x,object$link,Std.Values,GStd.Values)})
    } else {d_meanmargin_dd<-NULL}
    if (sum(paramtypes[[3]])>0){
      d_meanmargin_da<-lapply(c(1:length(etas)),function(x){D_continuous.margin.mean_alpha(whichAlphaest,outcomeMat[[x]],betacont,etas[[x]],object$link,Std.Values)})
    } else {d_meanmargin_da<-NULL}
    d_meanmargin_cont<-list(d_meanmargin_db,d_meanmargin_dd,d_meanmargin_da)
    d_meanmargin_cont<-d_meanmargin_cont[!sapply(d_meanmargin_cont,is.null)]
    d_meanmargin_cont<-lapply(c(1:length(etas)),function(x){do.call(cbind,lapply(d_meanmargin_cont,function(y){y[[x]]}))})

  } else {
    marginsmean_cont<-NULL
    d_meanmargin_cont<-NULL
  }

  if (object$Hetero){
    if (length(deltacont)>0){
      sdmargin<-lapply(c(1:length(etas)),function(x){continuous.margin.sd(deltacont,etas[[x]],object$link,Std.Values,GStd.Values)})
      marginsvar_cont<-lapply(sdmargin,function(x){apply(x,2,mean)})

      if (sum(paramtypes[[1]])>0){
        d_varmargin_db<-lapply(etas,function(x){D_continuous.margin.var_mean(X[,whichXest,drop=FALSE],deltacont,x,object$link,Std.Values,GStd.Values)})
      } else {d_varmargin_db<-NULL}
      if (sum(paramtypes[[2]])>0){
        d_varmargin_dd<-lapply(etas,function(x){D_continuous.margin.var_var(deltalocscon,whichZest,Z[,whichZest,drop=FALSE],deltacont,x,object$link,Std.Values,GStd.Values,HStd.Values)})
      } else {d_varmargin_dd<-NULL}
      if (sum(paramtypes[[3]])>0){
        d_varmargin_da<-lapply(c(1:length(etas)),function(x){D_continuous.margin.var_alpha(whichAlphaest,outcomeMat[[x]],deltacont,etas[[x]],object$link,Std.Values,GStd.Values)})
      } else {d_varmargin_da<-NULL}
      d_varmargin_cont<-list(d_varmargin_db,d_varmargin_dd,d_varmargin_da)
      d_varmargin_cont<-d_varmargin_cont[!sapply(d_varmargin_cont,is.null)]
      d_varmargin_cont<-lapply(c(1:length(etas)),function(x){do.call(cbind,lapply(d_varmargin_cont,function(y){y[[x]]}))})

    } else {
      marginsvar_cont<-NULL
      d_varmargin_cont<-NULL
    }

    # need to collect together margins and terms for the calculation of standard errors that occur in both the mean
    # and variance equations

    # check that there are variables to be collected that are in both the mean and variance equation
    BothEq<-object$BothEq[(object$BothEq$meanandvarLOC %in% betalocscon),,drop=FALSE]
    if (nrow(BothEq)>0){
      Xrows<-match(BothEq$meanandvarLOC,betalocscon)
      Zrows<-match(BothEq$meanandvarLOCZ,deltalocscon)
      marginsmeanvar_cont<-lapply(c(1:length(etas)),function(x){marginsmean_cont[[x]][Xrows]+marginsvar_cont[[x]][Zrows]})
      marginsmean_cont<-lapply(c(1:length(etas)),function(x){marginsmean_cont[[x]][-Xrows]})
      marginsvar_cont<-lapply(c(1:length(etas)),function(x){marginsvar_cont[[x]][-Zrows]})
      d_meanvarmargin_cont<-lapply(c(1:length(etas)),function(x){d_meanmargin_cont[[x]][Xrows,,drop=FALSE]+d_varmargin_cont[[x]][Zrows,,drop=FALSE]})
      d_meanmargin_cont<-lapply(c(1:length(etas)),function(x){d_meanmargin_cont[[x]][-Xrows,,drop=FALSE]})
      d_varmargin_cont<-lapply(c(1:length(etas)),function(x){d_varmargin_cont[[x]][-Zrows,,drop=FALSE]})

    } else {
      d_meanvarmargin_cont<-NULL
      marginsmeanvar_cont<-NULL
    }
  } else {
    marginsvar_cont<-NULL
    d_varmargin_cont<-NULL
    d_meanvarmargin_cont<-NULL
    marginsmeanvar_cont<-NULL
  }

  # now for binary variables, need to separate variables into three categories:
  #     1. Only mean
  #     2. Only variance
  #     3. Both mean variance
  if (object$Hetero){
    BothEq<-object$BothEq[(object$BothEq$meanandvarLOC %in% betalocsbin),,drop=FALSE]
    if (nrow(BothEq)>0){
      # take the locations that are in BothEq out of the vectors to be collected
      betalocsbin<-betalocsbin[!(betalocsbin %in% BothEq$meanandvarLOC)]
      deltalocsbin<-deltalocsbin[!(deltalocsbin %in% BothEq$meanandvarLOCZ)]
    }
  }

  if (length(betalocsbin)>0){
    marginsmean_Disc<-lapply(etas,function(x){discrete.margin_meanonly(beta,X,betalocsbin,x,object$link,Std.Values)})
    if (sum(paramtypes[[1]])>0){
      d_meanmargin_bin_db<-lapply(marginsmean_Disc,function(x){D_discrete.margin_meanonly.mean(betalocsbin,whichXest,X,attr(x,"inputetas"),object$link,Std.Values)})
    } else {d_meanmargin_bin_db<-NULL}
    if (sum(paramtypes[[2]])>0){
      d_meanmargin_bin_dd<-lapply(marginsmean_Disc,function(x){D_discrete.margin_mean.var(whichZest,Z,attr(x,"inputetas"),object$link,Std.Values,GStd.Values)})
    } else {d_meanmargin_bin_dd<-NULL}
    if (sum(paramtypes[[3]])>0){
      d_meanmargin_bin_da<-lapply(c(1:length(etas)),function(x){D_discrete.margin_mean.alpha(whichAlphaest,outcomeMat[[x]],attr(marginsmean_Disc[[x]],"inputetas"),Std.Values,object$link)})
    } else {d_meanmargin_bin_da<-NULL}
    d_meanmargin_bin<-list(d_meanmargin_bin_db,d_meanmargin_bin_dd,d_meanmargin_bin_da)
    d_meanmargin_bin<-d_meanmargin_bin[!sapply(d_meanmargin_bin,is.null)]
    d_meanmargin_bin<-lapply(c(1:length(etas)),function(x){do.call(cbind,lapply(d_meanmargin_bin,function(y){y[[x]]}))})
  } else {
    marginsmean_Disc<-NULL
    d_meanmargin_bin<-NULL
  }

  if (object$Hetero){
    if (length(deltalocsbin)>0){
      marginsvar_Disc<-lapply(etas,function(x){discrete.margin_varonly(delta,Z,deltalocsbin,object$sdmodel,x,object$link,Std.Values)})
      if (sum(paramtypes[[1]])>0){
        d_varmargin_bin_db<-lapply(marginsvar_Disc,function(x){D_discrete.margin_var.mean(whichXest,X,attr(x,"inputetas"),object$link,attr(x,"StdDevs"))})
      } else {d_varmargin_bin_db<-NULL}
      if (sum(paramtypes[[2]])>0){
        d_varmargin_bin_dd<-lapply(marginsvar_Disc,function(x){D_discrete.margin_varonly.var(deltalocsbin,whichZest,Z,attr(x,"inputetas"),attr(x,"ZDbinaries"),object$link,attr(x,"StdDevs"),gsdmodel)})
      } else {d_varmargin_bin_dd<-NULL}
      if (sum(paramtypes[[3]])>0){
        d_varmargin_bin_da<-lapply(c(1:length(etas)),function(x){D_discrete.margin_var.alpha(whichAlphaest,outcomeMat[[x]],attr(marginsvar_Disc[[x]],"inputetas"),attr(marginsvar_Disc[[x]],"StdDevs"),object$link)})
      } else {d_varmargin_bin_da<-NULL}
      d_varmargin_bin<-list(d_varmargin_bin_db,d_varmargin_bin_dd,d_varmargin_bin_da)
      d_varmargin_bin<-d_varmargin_bin[!sapply(d_varmargin_bin,is.null)]
      d_varmargin_bin<-lapply(c(1:length(etas)),function(x){do.call(cbind,lapply(d_varmargin_bin,function(y){y[[x]]}))})
    } else {
      marginsvar_Disc<-NULL
      d_varmargin_bin<-NULL
    }

    if (nrow(BothEq)>0){
      marginsmeanvar_Disc<-lapply(etas,function(x){discrete.margin_both(beta,X,delta,Z,BothEq,object$sdmodel,x,object$link,Std.Values)})
      if (sum(paramtypes[[1]])>0){
        d_meanvarmargin_bin_db<-lapply(marginsmeanvar_Disc,function(x){D_discrete.margin_meanvar.mean(whichXest,X,BothEq,attr(x,"inputetas"),attr(x,"StdDevs"),object$link)})
      } else {d_meanvarmargin_bin_db<-NULL}
      if (sum(paramtypes[[2]])>0){
        d_meanvarmargin_bin_dd<-lapply(marginsmeanvar_Disc,function(x){D_discrete.margin_meanvar.var(whichZest,Z,BothEq,attr(x,"inputetas"),attr(x,"ZDbinaries"),object$link,attr(x,"StdDevs"),gsdmodel)})
      } else {d_meanvarmargin_bin_dd<-NULL}
      if (sum(paramtypes[[3]])>0){
        d_meanvarmargin_bin_da<-lapply(c(1:length(etas)),function(x){D_discrete.margin_var.alpha(whichAlphaest,outcomeMat[[x]],attr(marginsmeanvar_Disc[[x]],"inputetas"),attr(marginsmeanvar_Disc[[x]],"StdDevs"),object$link)})
      } else {d_meanvarmargin_bin_da<-NULL}
      d_meanvarmargin_bin<-list(d_meanvarmargin_bin_db,d_meanvarmargin_bin_dd,d_meanvarmargin_bin_da)
      d_meanvarmargin_bin<-d_meanvarmargin_bin[!sapply(d_meanvarmargin_bin,is.null)]
      d_meanvarmargin_bin<-lapply(c(1:length(etas)),function(x){do.call(cbind,lapply(d_meanvarmargin_bin,function(y){y[[x]]}))})
    } else {
      marginsmeanvar_Disc<-NULL
      d_meanvarmargin_bin<-NULL
    }

  } else {
    marginsvar_Disc<-NULL
    d_varmargin_bin<-NULL
    marginsmeanvar_Disc<-NULL
    d_meanvarmargin_bin<-NULL
  }
  allmargins<-list(marginsmeanvar_Disc,marginsmeanvar_cont,marginsmean_Disc,marginsmean_cont,marginsvar_Disc,marginsvar_cont)
  keeplist<-sapply(allmargins,is.null)
  allmargins<-lapply(c(1:length(outcomeNames)),function(x){do.call(c,lapply(allmargins[!keeplist],function(y){y[[x]]}))})
  alld_margins<-list(d_meanvarmargin_bin,d_meanvarmargin_cont,d_meanmargin_bin,d_meanmargin_cont,d_varmargin_bin,d_varmargin_cont)[!keeplist]
  alld_margins<-lapply(c(1:length(outcomeNames)),function(x){do.call(rbind,lapply(alld_margins,function(y){y[[x]]}))})

  vcovmat<-vcov.oglmx(object)
  StdErrors<-lapply(alld_margins,function(x){calcMEstdErrors(x,vcovmat)})
  TValues<-lapply(c(1:length(StdErrors)),function(x){allmargins[[x]]/StdErrors[[x]]})
  PValues<-lapply(TValues,function(x){2*pnorm(-abs(x))})

  Xnames<-names(object$varMeans[[1]])
  marginnames<-character()
  if (object$Hetero){
    Znames<-names(object$varMeans[[2]])
    marginnames<-c(marginnames,BothEq$meanandvarNAME)
    marginnames<-c(marginnames,Xnames[betalocscon[betalocscon %in% object$BothEq$meanandvarLOC]])
    marginnames<-c(marginnames,Xnames[betalocsbin])
    marginnames<-c(marginnames,Xnames[betalocscon[!(betalocscon %in% object$BothEq$meanandvarLOC)]])
    marginnames<-c(marginnames,Znames[deltalocsbin])
    marginnames<-c(marginnames,Znames[deltalocscon[!(deltalocscon %in% object$BothEq$meanandvarLOCZ)]])
  } else {
    marginnames<-c(marginnames,Xnames[betalocsbin])
    marginnames<-c(marginnames,Xnames[betalocscon])
  }

  output<-lapply(c(1:length(StdErrors)),function(x){cbind("Marg. Eff"=allmargins[[x]],"Std. error"=StdErrors[[x]],
                                                          "t value"=TValues[[x]],"Pr(>|t|)"=PValues[[x]])})
  output<-lapply(output,function(x){rownames(x)<-marginnames; x})
  names(output)<-outcomeNames
  class(output)<-c("margins.oglmx")
  output
}
