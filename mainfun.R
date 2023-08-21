elliproj=function(y,tau){
  # This function performs the ellipsoid projection
  R=y
  for(l in 1:(m-1)){
    R_l=R[((m-l)*q+1):(m*q),1:(l*q)]
    e_l=lower.tri(matrix(0,l,l),diag =TRUE)
    ind_l=kronecker(e_l, matrix(rep(TRUE,q^2),q,q), FUN = "&")
    tmpnorm=sqrt(crossprod(R_l[ind_l]))
    R_l[ind_l]=max(1-tau/tmpnorm,0)*R_l[ind_l]
    R[((m-l)*q+1):(m*q),1:(l*q)]=R_l
  }
  return(R)
}
bT_update<- function(S,hatOmega,rho,u,bgamma){
  res=diag(m*q)
  #give Phi_j first then bT
  for(j in 2:m){
    A=2*hatOmega[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]
    B=S[1:((j-1)*q),1:((j-1)*q)]
    C=2*hatOmega[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]%*%S[((j-1)*q+1):(j*q),1:((j-1)*q)] +
      rho*bgamma[((j-1)*q+1):(j*q),1:((j-1)*q)]-u[((j-1)*q+1):(j*q),1:((j-1)*q)]
    res[((j-1)*q+1):(j*q),1:((j-1)*q)]=BS(A,B,rho,C)
  }
  return(res)
}

# S,hatOmega,init_bT,lambda1
bTadmm=function(S,hatOmega,init_bT,lambda1,tol=1e-4,itermax=1e5){
  #This function solve the following estimation problem
  #\argmin_T \{tr\{T^T\hat \Omega T S\}  + \lambda_1 \sum_{\ell=1}^{m-1}\|T_{g_{\ell}}\|_F \}
  #using an ADMM algorithm with changing rho
  r=(m*q)^2
  #Default parameter in ADMM 
  tolabs=tol
  tolrel=tol
  #changing rho
  rho=2
  mu=10
  inc=2
  dec=2
  #Initialize the result
  bT=init_bT
  bgamma=init_bT
  #dual variable
  u=matrix(0,m*q,m*q)
  for(i in seq(itermax)){
    #Primal&Dual Updates
    bT_new=bT_update(S=S,hatOmega=hatOmega,rho=rho,u=u,bgamma=bgamma)
    bgamma_new=elliproj(y=bT_new+u/rho,tau = lambda1/rho)
    u=u+rho*(bT_new-bgamma_new)
    #check convergence 
    #primal residual
    pres=sqrt(crossprod(as.vector(bT_new-bgamma_new))) 
    # dual residual 
    dres=rho*sqrt(crossprod(as.vector(bgamma_new-bgamma)))
    #primal tolerance
    peps=tolabs*sqrt(r)+tolrel*max(sqrt(crossprod(as.vector(bT_new))), 
                                   sqrt(crossprod(as.vector(bgamma_new))))
    #dual tolerance
    deps=tolabs*sqrt(r) + tolrel*sqrt(crossprod(as.vector(u))) 
    
    if( dres <= deps & pres <=peps) 
    { 
      break 
    } 
    else{ 
      # if not, update estimates and rho 
      bT=bT_new 
      bgamma=bgamma_new
      #update rho if needed
      if(pres>mu*dres){
        rho=rho*inc
      }
      else if(dres>mu*pres){
        rho=rho/dec 
      }
    }
  }
  if(i==itermax){
    cat("ADMM fails to converge",fill=TRUE)
  }
  return (bgamma_new) 
}
#glasso for Omega
hatomega_update<-function(S,bT,lambda2){
  bT=diag(m*q)-bT
  res=matrix(0,m*q,m*q)
  res[1:q,1:q]=glasso(s=S[1:q,1:q],rho=lambda2)$wi
  for(j in 2:m){
    S0=S[((j-1)*q+1):(j*q),1:((j-1)*q)]%*%t(bT[((j-1)*q+1):(j*q),1:((j-1)*q)])
    S1=S[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]-S0-t(S0)+bT[((j-1)*q+1):(j*q),1:((j-1)*q)]%*%
       S[1:((j-1)*q),1:((j-1)*q)]%*%t(bT[((j-1)*q+1):(j*q),1:((j-1)*q)])
    res[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]=glasso(s=S1,rho=lambda2)$wi
  }
  return(res)
}
Sigmainvband<-function(S,lambda1,lambda2,eps=1e-4,itermax=1000){
  T0=diag(m*q)
  #find the optimal Omega with T fixed at T0
  Omega0=hatomega_update(S,T0,lambda2)
  #find a candidate T with Omega fixed at Omega0
  T1=2*diag(m*q)-bTadmm(S,Omega0,T0,lambda1)
  #find the optimal Omega with T fixed at T1
  Omega1=hatomega_update(S,T1,lambda2)
  #alternative search for optimal T and Omega
  itercount=0
  while(sum(abs(T1-T0))/sum(abs(T0))>eps || sum(abs(Omega1-Omega0))/sum(abs(Omega0))>eps){
    T0=T1
    Omega0=Omega1
    Omega1=hatomega_update(S,T0,lambda2)
    T1=2*diag(m*q)-bTadmm(S,Omega1,T0,lambda1)
    itercount=itercount+1
    if(itercount>itermax)
    { 
      cat("Something is wrong ")
      break
    }
    #else{
    #cat("We have finished",itercount,"iterations",fill = TRUE)
    #}
  }
  Sigma.inv=t(T1)%*%Omega1%*%T1
  return(Sigma.inv)
}

#nflods-cross validation of lambda1 and lambda2
Sigmainvband.cv<-function(y, nlam=10,flmin=0.01, lam1list=NULL, lam2list=NULL,folds=NULL,nfolds=5){
  n=nrow(y)
  S=cov(y)*(n-1)/n
  if(is.null(folds)){
    folds=makefolds(n, nfolds)
    nfolds=length(folds)
  }
  if (is.null(lam1list) &&  (is.null(lam2list)) ){
    lambda1max=2.5*sqrt(q^2*log(m*q)/n)
    lambda2max=2*sqrt(log(q)/n)
    lam1list=sort(lambda1max*exp(seq(0, log(flmin), length= nlam)))
    lam2list=sort(lambda2max*exp(seq(0, log(flmin), length= nlam)))
  }
  P=array(0,c(m*q,m*q,nlam*nlam))
  errs=matrix(0, nlam*nlam, nfolds)
  for(k in seq(nfolds)){
    y.tr=y[-folds[[k]], ]
    S.tr=crossprod(y.tr)/(dim(y.tr)[1])
    for(i in seq(nlam)){
      for(j in seq(nlam)){
        lambda1 = lam1list[i]
        lambda2 = lam2list[j]
        P[ , ,(i-1)*nlam+j]=Sigmainvband(S, lambda1, lambda2)
      }
    }
    y.te=y[folds[[k]], ]
    S.te=crossprod(y.te)/(dim(y.te)[1])
    for(j in seq(nlam*nlam)){
      errs[j,k]=BIC(S.te, P[, , j])
    }
  }
  me=rowMeans(errs)
  ibest=which.min(me)
  for(i in seq(nlam)){
    for(j in seq(nlam)){
      if(ibest == (i-1)*nlam+j){
        ibest1=i
        ibest2=j
      }
    }
  }
  S.cv=Sigmainvband(S,lambda1 = lam1list[ibest1],lambda2 = lam2list[ibest2])
  list(errs=errs, me=me, lam1.best = lam1list[ibest1], lam2.best = lam2list[ibest2],
       ibest = c(ibest1,ibest2),S.cv=S.cv)
}




