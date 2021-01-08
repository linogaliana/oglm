(params, X, Z, outcomeMat, w,
  link,
  whichparametersmean,
  whichparametersscale,
  whichparametersthresh,
  sdmodel, gsdmodel, hsdmodel,
  analhessian = TRUE)


oglm:::eval_llk(
  params = as.numeric(mod1$coefficients),
  X = mod1$modelframes$X,
  Z = mod1$modelframes$Z,
  outcomeMat = mod1$modelframes$outcomeMatrix,
  whichparametersmean = attr(mod1$coefficients, "coefftypes")[[1]],
  whichparametersscale = attr(mod1$coefficients, "coefftypes")[[2]],
  whichparametersthresh = attr(mod1$coefficients, "coefftypes")[[3]],
  sdmodel = mod1$sdmodel,
  analhessian = TRUE
)
