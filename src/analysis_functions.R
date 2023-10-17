# libraries

# install.packages("devtools")
# library(devtools)
# install_github("hhuang90/TVD")
library(TVD)

#' create pointwise median matrix
#' 
#' Determines the pointwise median pointwiseMedian(t_i),
#' of n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#'
#' @return pointwiseMedianMatrix, k×p matrix, pointwiseMedian(t_i)=pointwiseMedianMatrix[l,i] for all l∈{1,...,k}.
createPointwiseMedianMatrix<-function(dataMatrix,k=1){
  p<-ncol(x=dataMatrix)
  pointwiseMedianMatrix<-matrix(data=apply(X=dataMatrix,MARGIN=2,FUN=median),nrow=k,ncol=p,byrow=TRUE)
  return(pointwiseMedianMatrix)
}

#' create fTilde(s,Δ) matrix for s=t-Δ
#' 
#' Determines fTilde_l(s_i,Δ), for s_i=t_i-Δ,
#' for k curves f_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i], and
#' for n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#'
#' @return fTildetMinusDeltaMatrix, k×p matrix, fTilde_l(s_i,Δ)=fTildetMinusDeltaMatrix[l,i-1] for s_i=t_i-Δ.
createfTildetMinusDeltaMatrix<-function(dataMatrix,fMatrix=dataMatrix){
  p<-ncol(x=dataMatrix)
  k<-nrow(x=fMatrix)
  pointwiseMediansPlusDeltaRow<-createPointwiseMedianMatrix(dataMatrix=dataMatrix)[,-1]
  fsMatrix<-fMatrix[,-p]
  fsPlusDeltaMatrix<-fMatrix[,-1]
  fTildetMinusDeltaMatrix<-cbind(t(x=t(x=fsMatrix-fsPlusDeltaMatrix)+pointwiseMediansPlusDeltaRow),NA)
  return(fTildetMinusDeltaMatrix)
}

#' create R_f(t) matrix
#' 
#' Estimates R_f(t) by R_f(t_i),
#' for curve f,
#' over p time points t_i,
#' with f(t_i)=fRow[i], 
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fRow, numeric vector of length p, Curve over p time points.
#'
#' @return RftMatrix, n×p matrix, R_f(t_i)=RftMatrix[j,i] for X_j.
createRftMatrix<-function(dataMatrix,fRow){
  RftMatrix<-t(x=t(x=dataMatrix)<=fRow)
  return(RftMatrix)
}

#' create E(R_f(t)) matrix
#' 
#' Estimates E(R_f(t)) by E(R_{f_l}(t_i)),
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @details If f_l is the pointwise median of X(t) for all l∈{1,...,k}, then returns 1/2.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param fIsPointwiseMedian, boolean, Set TRUE if f_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#'
#' @return ERftMatrix, numeric, 1/2, if fIsPointwiseMedian==TRUE.
#' @return ERftMatrix, k×p matrix, E(R_{f_l}(t_i))=ERftMatrix[l,i], if fIsPointwiseMedian==FALSE.
createERftMatrix<-function(dataMatrix,fMatrix=dataMatrix,fIsPointwiseMedian=FALSE){
  if(fIsPointwiseMedian){
    ERftMatrix<-1/2
  }else{
    ERftMatrix<-t(x=apply(X=fMatrix,MARGIN=1,FUN=function(fRow){apply(X=createRftMatrix(dataMatrix=dataMatrix,fRow=fRow),MARGIN=2,FUN=mean)}))
  }
  return(ERftMatrix)
}

#' create Var(R_f(t)) matrix
#' 
#' Estimates Var(R_f(t)) by Var(R_{f_l}(t_i)),
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param fIsPointwiseMedian, boolean, Set TRUE if f_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#'
#' @return VarRftMatrix, k×p matrix, Var(R_{f_l}(t_i))=VarRftMatrix[l,i].
createVarRftMatrix<-function(dataMatrix,fMatrix=dataMatrix,fIsPointwiseMedian=FALSE){
  ERftMatrix<-createERftMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix,fIsPointwiseMedian=fIsPointwiseMedian)
  VarRftMatrix<-ERftMatrix*(1-ERftMatrix)
  return(VarRftMatrix)
}

#' create E(E(R_f(t)|R_g(t-Δ))^2) matrix
#' 
#' Estimates E(E(R_f(t)|R_g(t-Δ))^2) by E(E(R_{f_l}(t_i)|R_{g_l}(t_{i-1}))^2),
#' for k curves f_l and k curves g_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i] and g_l(t_i)=gMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param gMatrix, k×p matrix, k curves over p time points.
#' @param gIsPointwiseMedian, boolean, Set TRUE if g_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#' @param isTilde, boolean, Set TRUE if fMatrix==createPointwiseMedianMatrix(dataMatrix=dataMatrix,k=k) and gMatrix==createfTildetMinusDeltaMatrix(dataMatrix=dataMatrix), FALSE otherwise.
#'
#' @return EERftGivenRgtMinusDeltaSquaredMatrix, k×p matrix, E(E(R_{f_l}(t_i)|R_{g_l}(t_{i-1}))^2)=EERftGivenRgtMinusDeltaSquaredMatrix[l,i].
createEERftGivenRgtMinusDeltaSquaredMatrix<-function(dataMatrix,fMatrix=dataMatrix,gMatrix=fMatrix,gIsPointwiseMedian=FALSE,isTilde=FALSE){
  p<-ncol(dataMatrix)
  k<-nrow(fMatrix)
  
  # E(E(R_f(t)|R_g(t-Δ))^2)
  # =E(R_f(t)|R_g(t-Δ)=1)^2*P(R_g(t-Δ)=1)+E(R_f(t)|R_g(t-Δ)=0)^2*P(R_g(t-Δ)=0)
  # =(1P(R_f(t)=1|R_g(t-Δ)=1)+0P(R_f(t)=0|R_g(t-Δ)=1))^2*P(R_g(t-Δ)=1)+(1P(R_f(t)=1|R_g(t-Δ)=0)+0P(R_f(t)=0|R_g(t-Δ)=0))^2*P(R_g(t-Δ)=0)
  # =P(R_f(t)=1|R_g(t-Δ)=1)^2*P(R_g(t-Δ)=1)+P(R_f(t)=1|R_g(t-Δ)=0)^2*P(R_g(t-Δ)=0)
  # =P(R_f(t)=1 and R_g(t-Δ)=1)^2/P(R_g(t-Δ)=1)+P(R_f(t)=1 and R_g(t-Δ)=0)^2/P(R_g(t-Δ)=0)
  EERftGivenRgtMinusDeltaSquaredMatrix<-matrix(data=NA,nrow=k,ncol=p)
  for(l in 1:k){
    flRow<-fMatrix[l,]
    glRow<-gMatrix[l,]
    dataMatrixfl<-dataMatrix
    dataMatrixgl<-dataMatrix
    if(isTilde){
      dataMatrixfl[l,]<-flRow
      dataMatrixgl[l,]<-glRow
    }
    PRflt1AndRgltMinusDelta1Row<-apply(X=createRftMatrix(dataMatrix=dataMatrixfl,fRow=flRow)&cbind(NA,createRftMatrix(dataMatrix=dataMatrixgl,fRow=glRow)[,-p]),MARGIN=2,FUN=mean)
    PRflt1AndRgltMinusDelta0Row<-apply(X=createRftMatrix(dataMatrix=dataMatrixfl,fRow=flRow)&cbind(NA,!createRftMatrix(dataMatrix=dataMatrixgl,fRow=glRow)[,-p]),MARGIN=2,FUN=mean)
    PRgltMinusDelta1Row<-cbind(NA,data=createERftMatrix(dataMatrix=dataMatrixgl,fMatrix=t(x=as.matrix(x=glRow)),fIsPointwiseMedian=gIsPointwiseMedian)[,-p,drop=FALSE])
    PRgltMinusDelta0Row<-1-PRgltMinusDelta1Row
    for(i in 2:p){
      PRflti1AndRgltiMinusDelta1<-PRflt1AndRgltMinusDelta1Row[i]
      PRflti1AndRgltiMinusDelta0<-PRflt1AndRgltMinusDelta0Row[i]
      PRgltiMinusDelta1<-PRgltMinusDelta1Row[i]
      PRgltiMinusDelta0<-PRgltMinusDelta0Row[i]
      if(PRgltiMinusDelta1==0){
        EERftGivenRgtMinusDeltaSquaredMatrix[l,i]<-PRflti1AndRgltiMinusDelta0**2/PRgltiMinusDelta0
      }else if(PRgltiMinusDelta0==0){
        EERftGivenRgtMinusDeltaSquaredMatrix[l,i]<-PRflti1AndRgltiMinusDelta1**2/PRgltiMinusDelta1
      }else{
        EERftGivenRgtMinusDeltaSquaredMatrix[l,i]<-PRflti1AndRgltiMinusDelta1**2/PRgltiMinusDelta1+PRflti1AndRgltiMinusDelta0**2/PRgltiMinusDelta0
      }
    }
  }
  
  return(EERftGivenRgtMinusDeltaSquaredMatrix)
}

#' create Var(E(R_f(t)|R_g(t-Δ))) matrix
#' 
#' Estimates Var(E(R_f(t)|R_g(t-Δ))) by Var(E(R_{f_l}(t_i)|R_{g_l}(t_{i-1}))),
#' for k curves f_l and k curves g_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i] and g_l(t_i)=gMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param gMatrix, k×p matrix, k curves over p time points.
#' @param fIsPointwiseMedian, boolean, Set TRUE if f_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#' @param gIsPointwiseMedian, boolean, Set TRUE if g_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#' @param isTilde, boolean, Set TRUE if fMatrix==createPointwiseMedianMatrix(dataMatrix=dataMatrix,k=k) and gMatrix==createfTildetMinusDeltaMatrix(dataMatrix=dataMatrix), FALSE otherwise.
#'
#' @return VarERftGivenRgtMinusDeltaMatrix, k×p matrix, Var(E(R_{f_l}(t_i)|R_{g_l}(t_{i-1})))=VarERftGivenRgtMinusDeltaMatrix[l,i].
createVarERftGivenRgtMinusDeltaMatrix<-function(dataMatrix,fMatrix=dataMatrix,gMatrix=fMatrix,fIsPointwiseMedian=FALSE,gIsPointwiseMedian=FALSE,isTilde=FALSE){
  EERftGivenRgtMinusDeltaSquaredMatrix<-createEERftGivenRgtMinusDeltaSquaredMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix,gMatrix=gMatrix,gIsPointwiseMedian=gIsPointwiseMedian,isTilde=isTilde)
  
  ERftMatrix<-createERftMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix,fIsPointwiseMedian=fIsPointwiseMedian)
  
  # Var(E(R_f(t)|R_g(t)))
  # =E(E(R_f(t)|R_g(t))^2)-E(E(R_f(t)|R_g(t)))^2
  # =E(E(R_f(t)|R_g(t))^2)-E(R_f(t))^2 (by law of total expectation)
  VarERftGivenRgtMinusDeltaMatrix<-EERftGivenRgtMinusDeltaSquaredMatrix-ERftMatrix**2
  
  return(VarERftGivenRgtMinusDeltaMatrix)
}

#' create E(Var(R_f(t)|R_g(t))) matrix
#' 
#' Estimates E(Var(R_f(t)|R_g(t))) by E(Var(R_{f_l}(t_i)|R_{g_l}(t_{i-1}))),
#' for k curves f_l and k curves g_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i] and g_l(t_i)=gMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param gMatrix, k×p matrix, k curves over p time points.
#' @param fIsPointwiseMedian, boolean, Set TRUE if f_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#' @param gIsPointwiseMedian, boolean, Set TRUE if g_l is the pointwise median of X(t) for all l∈{1,...,k}, FALSE otherwise.
#' @param isTilde, boolean, Set TRUE if fMatrix==createPointwiseMedianMatrix(dataMatrix=dataMatrix,k=k) and gMatrix==createfTildetMinusDeltaMatrix(dataMatrix=dataMatrix), FALSE otherwise.
#'
#' @return EVarRftGivenRgtMinusDeltaMatrix, k×p matrix, E(Var(R_{f_l}(t_i)|R_{g_l}(t_{i-1})))=EVarRftGivenRgtMinusDeltaMatrix[l,i].
createEVarRftGivenRftMinusDeltaMatrix<-function(dataMatrix,fMatrix=dataMatrix,fIsPointwiseMedian=FALSE,gIsPointwiseMedian=FALSE,isTilde=FALSE){
  
  # E(R_f(t)**2)
  # =1**2*P(R_f(t)=1)+0**2*P(R_f(t)=0)
  # =P(R_f(t)=1)
  # =1P(R_f(t)=1)+0P(R_f(t)=0)
  # =E(R_f(t))
  ERftSquaredMatrix<-createERftMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix,fIsPointwiseMedian=fIsPointwiseMedian)
  
  EERftGivenRftMinusDeltaSquaredMatrix<-createEERftGivenRgtMinusDeltaSquaredMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix,gMatrix=fMatrix,gIsPointwiseMedian=gIsPointwiseMedian,isTilde=FALSE)
  
  # E(Var(R_f(t)|R_f(t-Δ)))
  # =E(E(R_f(t)^2|R_f(t-Δ))-E(R_f(t)|R_f(t-Δ))^2)
  # =E(E(R_f(t)^2|R_f(t-Δ)))-E(E(R_f(t)|R_f(t-Δ))^2) (by linearity of expectation)
  # =E(R_f(t)^2)-E(E(R_f(t)|R_f(t-Δ))^2) (by law of total expectation)
  EVarRftGivenRftMinusDeltaMatrix<-ERftSquaredMatrix-EERftGivenRftMinusDeltaSquaredMatrix
  
  return(EVarRftGivenRftMinusDeltaMatrix)
}

#' create v(t,Δ) matrix
#' 
#' Determines v(t_i,Δ),
#' for k curves f_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, l∈{1,...,k}.
#' 
#' @details I couldn't find the value of v(t,Δ) for constant functions, and I set it to be 0.
#' 
#' @param fMatrix, k×p matrix, k curves over p time points.
#'
#' @return vtMatrix, k×p matrix, v(t,Δ)=vtMatrix[l,i] for f=fMatrix[l,].
createvtMatrix<-function(fMatrix){
  p<-ncol(fMatrix)
  vtMatrix<-t(x=apply(X=fMatrix,MARGIN=1,FUN=function(fRow){if(sum(abs(x=diff(x=fRow)))!=0){c(NA,abs(x=diff(x=fRow))/sum(abs(x=diff(x=fRow))))}else{c(NA,rep(x=0,times=p-1))}}))
  return(vtMatrix)
}

#' create S_fTilde(t,Δ) matrix
#' 
#' Estimates S_fTilde(t,Δ) by S_{fTilde_l}(t_i,Δ),
#' for k curves f_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#'
#' @return SfTildetMatrix, k×p matrix, S_{fTilde_l}(t_i,Δ)=SfTildetColumn[l,i].
createSfTildetMatrix<-function(dataMatrix,fMatrix=dataMatrix){
  k<-nrow(x=fMatrix)
  fTildeMatrix<-createPointwiseMedianMatrix(dataMatrix=dataMatrix,k=k)
  fTildetMinusDeltaMatrix<-createfTildetMinusDeltaMatrix(dataMatrix=dataMatrix)
  VarERfTildetGivenRfTildetMinusDeltaMatrix<-createVarERftGivenRgtMinusDeltaMatrix(dataMatrix=dataMatrix,fMatrix=fTildeMatrix,gMatrix=fTildetMinusDeltaMatrix,fIsPointwiseMedian=TRUE,isTilde=TRUE)
  VarRfTildetMatrix<-createVarRftMatrix(dataMatrix=dataMatrix,fMatrix=fTildeMatrix,fIsPointwiseMedian=TRUE)
  SfTildetMatrix<-VarERfTildetGivenRfTildetMinusDeltaMatrix/VarRfTildetMatrix
  return(SfTildetMatrix)
}

#' create MSS(f,Δ) column
#' 
#' Estimates MSS(f,Δ) by MSS(f_l,Δ),
#' for k curves f_l,
#' over p equally spaced time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' where t_i-Δ=t_{i-1},
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#'
#' @details
#' createMSSfColumn(dataMatrix=dataMatrix) is the same as: 
#' TVDMSS(data=t(x=dataMatrix),nCurve=nrow(x=dataMatrix),nPoint=ncol(x=DataMatrix))$MSS,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, p curves over p time points.
#'
#' @return MSSfColumn, numeric vector of length k, MSS(f_l,Δ)=MSSfColumn[l].
createMSSfColumn<-function(dataMatrix,fMatrix=dataMatrix){
  MSSfColumn<-apply(X=createvtMatrix(fMatrix=fMatrix)[,-1]*createSfTildetMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix)[,-1],MARGIN=1,FUN=sum)
  return(MSSfColumn)
}

#' create w(t)
#' 
#' Determines w(t),
#' for p time points t_i,
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#'
#' @return wt, numeric, w(t_i)=wt.
createwt<-function(dataMatrix){
  p<-ncol(x=dataMatrix)
  wt<-1/p
  return(wt)
}

#' create TVD(f) column
#' 
#' Estimates TVD(f) by TVD(f_l),
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @details
#' createTVDfColumn(dataMatrix=dataMatrix) is the same as: 
#' TVDMSS(data=t(x=dataMatrix),nCurve=nrow(x=dataMatrix),nPoint=ncol(x=dataMatrix))$TVD,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#'
#' @return TVDfColumn, numeric vector of length k, TVD(f_l)=TVDfColumn[l].
createTVDfColumn<-function(dataMatrix,fMatrix=dataMatrix){
  TVDfColumn<-apply(X=createwt(dataMatrix=dataMatrix)*createVarRftMatrix(dataMatrix=dataMatrix,fMatrix=fMatrix),MARGIN=1,FUN=sum)
  return(TVDfColumn)
}

#' create MBD(f) column
#' 
#' Estimates MBD(f) by MBD(f_l),
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @references
#' López-Pintado, S., & Romo, J. (2009). On the concept of depth for functional data. Journal of the American Statistical Association, 104(486), 718–734. https://doi.org/10.1198/jasa.2009.0108 
#' Sun, Y., & Genton, M. G. (2011). Functional boxplots. Journal of Computational and Graphical Statistics, 20(2), 316–334. https://doi.org/10.1198/jcgs.2011.09224 
#' 
#' @details
#' createMBdfColumn(dataMatrix=dataMatrix) is the same as: 
#' fbplot(fit=dataMatrix)$depth
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#'
#' @return MBDfColumn, numberic vector of length k, MBD(f_l)=MBDfColumn[l].
createMBDfColumn<-function(dataMatrix,fMatrix=dataMatrix){
  n<-nrow(x=dataMatrix)
  p<-ncol(x=dataMatrix)
  k<-nrow(x=fMatrix)
  MBDfColumn<-rep(x=0,times=k)
  for (l1 in 1:(n-1)){
    for (l2 in (l1+1):n){
      upperBorderOfBandRow<-apply(X=dataMatrix[c(l1,l2),],MARGIN=2,FUN=max)
      lowerBorderOfBandRow<-apply(X=dataMatrix[c(l1,l2),],MARGIN=2,FUN=min)
      lambdaAfl1fl2<-apply(X=t(x=lowerBorderOfBandRow<=t(x=fMatrix))&t(x=t(x=fMatrix)<=upperBorderOfBandRow),MARGIN=1,FUN=sum)
      lambdarAfl1fl2<-lambdaAfl1fl2/p
      MBDfColumn<-MBDfColumn+lambdarAfl1fl2
    }
  }
  MBDfColumn<-1/choose(n=n,k=2)*MBDfColumn
  return(MBDfColumn)
}

#' create ED(f) column
#' 
#' Estimates ED(f) by ED(f_l),
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=dataMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @references
#' Narisetty, N. N., & Nair, V. N. (2016). Extremal depth for functional data and applications. Journal of the American Statistical Association, 111(516), 1705–1714. https://doi.org/10.1080/01621459.2015.1110033 
#' 
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#'
#' @return EDfColumn, numeric vector of length k, ED(f_l)=EDfColumn[l].
createEDfColumn<-function(dataMatrix,fMatrix=dataMatrix){
  n<-nrow(x=dataMatrix)
  k<-nrow(x=fMatrix)
  DgtMatrix<-round(x=t(x=apply(X=fMatrix,MARGIN=1,FUN=function(g){1-abs(x=apply(X=t(x=t(x=dataMatrix)<g)-t(x=g<t(x=dataMatrix)),MARGIN=2,FUN=sum))/n})),digits=3)
  rangeOfDgt<-round(x=(1:n)/n,digits=3)
  PhigrMatrix<-sapply(X=rangeOfDgt,FUN=function(r){apply(X=t(x=t(x=DgtMatrix)<=r),MARGIN=1,FUN=sum)})
  fiLeqgMatrix<-matrix(data=NA,nrow=n,ncol=k)
  for(fiIndex in 1:n){
    for(gIndex in 1:k){
      fiLeqg<-NULL
      depthIndex<-1
      while(is.null(x=fiLeqg)&(depthIndex<=n)){
        if(PhigrMatrix[gIndex,depthIndex]<PhigrMatrix[fiIndex,depthIndex]){
          fiLeqg<-TRUE
        }else if(PhigrMatrix[fiIndex,depthIndex]<PhigrMatrix[gIndex,depthIndex]){
          fiLeqg<-FALSE
        }
        depthIndex<-depthIndex+1
      }
      if(is.null(x=fiLeqg)){
        fiLeqg<-TRUE
      }
      fiLeqgMatrix[fiIndex,gIndex]<-fiLeqg
    }
  }
  EDfColumn<-apply(X=fiLeqgMatrix,MARGIN=2,FUN=sum)/n
  return(EDfColumn)
}

#' detect outlier indices using boxplot
#' 
#' Detects outlier indices using a boxplot with depth and factor for determining whiskers equal to factor,
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @details
#' detectOutlierIndicesBoxplot(x=dataMatrix=dataMatrix,depth=depth,factor=factor) is the same as: 
#' .boxDetect(depth=match.fun(FUN=paste("create",depth,"fColumn",sep=""))(dataMatrix=dataMatrix,fMatrix=fMatrix),boxFactor=factor) in increasing order by depth value,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, n×p matrix, n curves over p time points.
#' @param depth, character, "MSS", "TVD", "MBD", or "ED".
#' @param factor, integer, Factor for determining whiskers.
#'
#' @return outlierIndicesInIncreasingOrderByDepth, integer vector, Outlier indices in increasing order by depth.
detectOutlierIndicesBoxplot<-function(dataMatrix,fMatrix=dataMatrix,depth,factor){
  depthValuesColumn<-match.fun(FUN=paste("create",depth,"fColumn",sep=""))(dataMatrix=dataMatrix,fMatrix=fMatrix)
  fiveNumberSummary<-fivenum(x=depthValuesColumn)
  lowerHingeOfDepthValues<-fiveNumberSummary[2]
  IQROfDepthValues<-fiveNumberSummary[4]-fiveNumberSummary[2]
  lowerBoundOfNonOutliers<-lowerHingeOfDepthValues-factor*IQROfDepthValues
  outlierIndices<-which(depthValuesColumn<lowerBoundOfNonOutliers)
  numberOfOutliers<-length(x=outlierIndices)
  if(0<numberOfOutliers){
    indicesInIncreasingOrderByDepth<-order(depthValuesColumn)
    outlierIndicesInIncreasingOrderByDepth<-indicesInIncreasingOrderByDepth[sort(x=match(x=outlierIndices,table=indicesInIncreasingOrderByDepth))]
  }else{
    outlierIndicesInIncreasingOrderByDepth<-NULL
  }
  return(outlierIndicesInIncreasingOrderByDepth)
}