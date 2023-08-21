library (MASS)
library(glasso)
library(plotrix)
library(graphics)
library(Matrix)
library(varband)
library(QZ)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
source("otherfun.R")
sourceCpp("mainfun.cpp") 
#Give bT,bOmega,bSigma and Y=(y1,y2,...,yn)

q=10
m=50
p=m*q
n=50
k=1
Phi0=matrix(0,m,m)
for(i in 2:m){
  for(j in 1:m){
    if(j<i){ 
      Phi0[i,j]=0.5*ifelse((i-j)<=k,1,0)
      
    }
  }
}
Phi1=ar1(0.5,q)
bT0=diag(m*q)-kronecker(Phi0,Phi1)
Omega = SparseOmega(q)
bOmega=kronecker(diag(m),Omega)
bSigma.inv=t(bT0)%*%bOmega%*%bT0
bSigma=ginv(bSigma.inv)
Y=mvrnorm(n,mu=rep(0,m*q),bSigma)
S=cov(Y)*(n-1)/n


lambda1 <- 8*q*sqrt(log(m*q)/n)

lambda2 <- 0.75*sqrt(log(q)/n)

system.time(A <- Sigmainvband_CC(S, m, q, lambda1, lambda2))

S.tpn(A, bSigma.inv)

opnloss(bSigma.inv, A)

write.table(A,file ="A", sep ="\t", quote =F, eol = "\n",col.names = F) 

write.table(bSigma.inv,file ="bS", sep ="\t", quote =F, eol = "\n",col.names = F) 
