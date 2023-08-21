#nflods-cross validation of lambda1 and lambda2
Sigmainvband.cv_CP <- function(y, nlam=10,flmin=0.01, lam1list=NULL, lam2list=NULL,folds=NULL,nfolds=5){
  n=nrow(y)
  S=cov(y)*(n-1)/n
  if(is.null(folds)){
    folds=makefolds(n, nfolds)
    nfolds=length(folds)
  }
  if (is.null(lam1list) &&  (is.null(lam2list)) ){
    lambda1max=2*q*sqrt(log(m*q)/n)
    lambda2max=2*sqrt(log(q)/n)
    lam1list=sort(lambda1max*exp(seq(0, log(flmin), length= nlam)))
    lam2list=sort(lambda2max*exp(seq(0, log(flmin), length= nlam)))
  }
  Clam = expand.grid(lam1list, lam2list)
  Lam = apply(Clam, MARGIN = 1, list)
  errs=matrix(0, nlam*nlam, nfolds)
  for (i in seq(nfolds)) {
    y.tr = y[-folds[[1]], ]
    S.tr = crossprod(y.tr)/(dim(y.tr)[1])
    clusterExport(cl,  c("Lam", "S.tr", "m", "q"), envir = environment())
    P  = parLapply(cl, 1:length(Lam), function(i)
     {
      lambda1 = as.double(Lam[[i]][[1]][1])
      lambda2 = as.double(Lam[[i]][[1]][2])
      Sigmainvband_CC(S.tr, m, q, lambda1, lambda2)
     })
    y.te=y[folds[[i]], ]
    S.te=crossprod(y.te)/(dim(y.te)[1])
    for(j in seq(nlam*nlam)){
      errs[j,i]=nlikelihood(S.te, P[[j]])
    }
  }
  me=rowMeans(errs)
  ibest=which.min(me)
  for(i in seq(nlam)){
    for(j in seq(nlam)){
      if(ibest == (j-1)*nlam+i){
        ibest1=i
        ibest2=j
      }
    }
  }
  lam1best = lam1list[ibest1]
  lam2best = lam2list[ibest2]
  S.cv = Sigmainvband_CC(S, m, q, lam1best, lam2best)
  list(errs=errs, me=me, lam1best = lam1list[ibest1], lam2best = lam2list[ibest2],
       ibest = c(ibest1,ibest2),S.cv=S.cv)
}


