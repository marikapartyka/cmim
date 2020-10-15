#fast version

# Loading required packages

if(!require(mlbench)){
  install.packages("mlbench")
  library(mlbench)
}

if(!require(mlbench)){
  install.packages("infotheo")
  library(infotheo)
}

# Implementation

CMIMselection<-function(X,y, kmax){
  
  stopifnot(all(!is.na(X)) & nrow(X)==length(y) & ncol(X)>1 & kmax<=ncol(X) & kmax==floor(kmax))
  N<-ncol(X)
  S<-rep(0,N)
  score<-rep(0,N)
  nu<-rep(0,N)
  ps<-rep(0,N)
  
  for(i in (1:N))
  {
    score[i]<-mutinformation(y, X[,i])
  }
  
  for(k in (1:kmax))
  {
    ps[k]<-0
    for(j in (1:N))
    {
      while((score[j]>ps[k]) & (S[j]<k-1) )
      {
        S[j]<-S[j]+1
        score[j]<-min(score[j], condinformation(y, X[,j], X[,nu[S[j]]]))
      }
      if(score[j]>ps[k])
      {
        ps[k]<-score[j]
        nu[k]<-j
      }
      
    }
    
  }
  return(list(S = nu[1:kmax], score = ps[1:kmax]))

} 


data("BreastCancer")
rows_without_NA <- complete.cases(BreastCancer)
X <- BreastCancer[obserwacje_bez_NA , -c(1, 11)]
y <- BreastCancer[obserwacje_bez_NA , 11]

CMIM_sel <- CMIMselection(X=X, y=y, kmax=9)
CMIM_sel$S
CMIM_sel$score

