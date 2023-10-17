# libraries

# install.packages("mvtnorm")
library(mvtnorm)

#' create time row
#' 
#' Generates p time points equally spaced within [0,1].
#'
#' @param p, integer, Number of time points.
#'
#' @return timeRow, numeric vector of length p, p time points.
createTimeRow<-function(p){
  timeRow<-seq(from=0,to=1,length.out=p+1)[-1]
  return(timeRow)
}

#' create e matrix
#'
#' Generates n samples over p time points equally spaced within [0,1],
#' for a zero-mean Gaussian process,
#' with covariance function covariance:[0,1]×[0,1]↦ℝ,covariance(t1,t2)=e^(-|t1-t2|).
#' 
#' @references
#' ntguardian. (2018, January 26). R function for simulating Gaussian processes: R-bloggers. R. https://www.r-bloggers.com/2018/01/r-function-for-simulating-gaussian-processes/ 
#'
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return eMatrix, n×p matrix, n samples over p time points.
createeMatrix<-function(n,p){
  timeRow<-createTimeRow(p=p)
  covariance<-function(t1,t2){exp(x=-abs(x=t1-t2))}
  sigma<-sapply(X=timeRow,FUN=function(t1){sapply(X=timeRow,FUN=function(t2){covariance(t1,t2)})})
  eMatrix<-rmvnorm(n=n,sigma=sigma)
  return(eMatrix)
}

#' create c column
#'
#' Generates n samples,
#' from a Bernoulli distribution with probability parameter rho=1/10.
#'
#' @param n, integer, Number of samples.
#'
#' @return cColumn, numeric vector of length n, n samples.
createcColumn<-function(n){
  rho<-1/10
  cColumn<-rbinom(n=n,size=1,prob=rho)
  return(cColumn)
}

#' create sigma column
#'
#' Generates n samples,
#' from {-1,1}, where both -1 and 1 have probability 1/2.
#'
#' @param n, integer, Number of samples.
#'
#' @return sigmaColumn, numeric vector of length n, n samples.
createsigmaColumn<-function(n){
  sigmaColumn<-sample(x=c(1,-1),size=n,replace=TRUE)
  return(sigmaColumn)
}

#' create eTilde matrix
#'
#' Generates n samples over p time points equally spaced within [0,1],
#' for a zero-mean Gaussian process,
#' with covariance function covariance:[0,1]×[0,1]↦ℝ,covariance(t1,t2)=6e^(-|t1-t2|^0.1).
#' 
#' @references
#' ntguardian. (2018, January 26). R function for simulating Gaussian processes: R-bloggers. R. https://www.r-bloggers.com/2018/01/r-function-for-simulating-gaussian-processes/ 
#'
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return eTildeMatrix, n×p matrix, n samples over p time points.
createeTildeMatrix<-function(n,p){
  timeRow<-createTimeRow(p=p)
  covariance<-function(t1,t2){6*exp(x=-abs(x=t1-t2)**0.1)}
  sigma<-sapply(X=timeRow,FUN=function(t1){sapply(X=timeRow,FUN=function(t2){covariance(t1,t2)})})
  eTildeMatrix<-rmvnorm(n=n,sigma=sigma)
  return(eTildeMatrix)
}

#' create Model 1 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 1, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel1Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  samples<-t(x=4*timeRow+t(x=eMatrix))
  outlierIndices<-NULL
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 2 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 2, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel2Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  cColumn<-createcColumn(n=n)
  sigmaColumn<-createsigmaColumn(n=n)
  samples<-t(x=4*timeRow+t(x=eMatrix))+6*cColumn*sigmaColumn
  outlierIndices<-which(x=cColumn==1)
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 3 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 3, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel3Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  TColumnNumeric<-runif(n=n)
  TMatrixBoolean<-sapply(X=timeRow,FUN=function(time){TColumnNumeric<=time})
  cColumn<-createcColumn(n=n)
  sigmaColumn<-createsigmaColumn(n=n)
  samples<-t(x=4*timeRow+t(x=eMatrix))+TMatrixBoolean*6*cColumn*sigmaColumn
  outlierIndices<-which(x=cColumn==1)
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 4 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 4, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel4Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  l=0.08
  TColumnNumeric<-runif(n=n)
  TMatrixBoolean<-sapply(X=timeRow,FUN=function(time){(TColumnNumeric<=time)&(time<=(TColumnNumeric+l))})
  cColumn<-createcColumn(n=n)
  sigmaColumn<-createsigmaColumn(n=n)
  samples<-t(x=4*timeRow+t(x=eMatrix))+TMatrixBoolean*6*cColumn*sigmaColumn
  outlierIndices<-which(x=(0<apply(X=TMatrixBoolean,MARGIN=1,FUN=sum)&cColumn==1))
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 5 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 5, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel5Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  cColumn<-createcColumn(n=n)
  eMatrix<-createeMatrix(n=n,p=p)
  eTildeMatrix<-createeTildeMatrix(n=n,p=p)
  samples<-t(x=4*timeRow+t(x=(1-cColumn)*eMatrix+cColumn*eTildeMatrix))
  outlierIndices<-which(x=cColumn==1)
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 6 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 6, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel6Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  cColumn<-createcColumn(n=n)
  samples<-t(x=4*timeRow+t(x=eMatrix)+t(x=sapply(X=timeRow,FUN=function(time){cColumn*(0.5*sin(x=40*pi*time))})))
  outlierIndices<-which(x=cColumn==1)
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create Model 7 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for Model 7, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
#' @return outlierIndices, numeric vector, Outlier indices.
createModel7Samples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  cColumn<-createcColumn(n=n)
  eMatrix<-createeMatrix(n=n,p=p)
  samples<-sapply(X=timeRow,FUN=function(time){2*sin(x=15*pi*time+2*cColumn)})+eMatrix
  outlierIndices<-which(x=cColumn==1)
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}

#' create base model of Model 7 samples
#' 
#' Generates n samples over p time points equally spaced within [0,1],
#' for the base model of Model 7, and
#' returns outlier indices.
#' 
#' @param n, integer, Number of samples.
#' @param p, integer, Number of time points.
#'
#' @return samples, n×p matrix, n samples over p time points.
createModel7BaseSamples<-function(n,p){
  timeRow<-createTimeRow(p=p)
  eMatrix<-createeMatrix(n=n,p=p)
  samples<-t(x=2*sin(x=15*pi*timeRow)+t(x=eMatrix))
  outlierIndices<-NULL
  return(list(samples=samples,
              outlierIndices=outlierIndices))
}