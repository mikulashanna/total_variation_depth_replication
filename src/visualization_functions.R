#' create ylim
#' 
#' Creates ylim for a plot.
#' 
#' @param dataMatrix, n×p matrix, n curves p time points to plot.
#' 
#' @return ylim, numeric vector of length 2, ylim for a plot.
createylim<-function(dataMatrix){
  yMinimum<-min(dataMatrix)
  yMaximum<-max(dataMatrix)
  yLowerLimit<-yMinimum-(yMaximum-yMinimum)/2
  yUpperLimit<-yMaximum+(yMaximum-yMinimum)/2
  ylim<-c(yLowerLimit,yUpperLimit)
}

#' create Figure 1
#' 
#' Creates Figure 1 with the given outlier.
#' 
#' @param nonOutlierMatrix, (n-1)×p matrix, n-1 non-outliers over p time points.
#' @param outlierMatrix, 1×p matrix, Outlier over p time points.
#' @param color, character, Color of outlier.
#' @param ylim, numeric vector of length 2, ylim for the plot.
#' @param title, character, Title of plot.
createFigure1<-function(nonOutlierMatrix,outlierMatrix,color,ylim,title){
  p<-ncol(x=nonOutlierMatrix)
  timeRow<-createTimeRow(p=p)
  dataMatrix<-rbind(nonOutlierMatrix,outlierMatrix)
  matplot(x=timeRow,y=t(x=nonOutlierMatrix),type="b",lty=1,pch=1,col="lightgrey",xlab="t",ylab="X(t)",ylim=ylim,main=title)
  points(x=timeRow,y=outlierMatrix,type="b",pch=1,col=color)
}

#' pointwise variation depth of outlier in Figure 1
#' 
#' Generates 
#' Var(R_pointwiseMedian(t)) for the pointwise median, and
#' Var(R_f(t)),
#' Var(E(R_f(t)|R_f(t-Δ))),
#' E(Var(R_f(t)|R_f(t-Δ))), and
#' values for outlier f in Figure 1,
#' where Δ is the difference between two consecutive time points.
#' 
#' @param nonOutlierMatrix, (n-1)×p matrix, n-1 non-outliers over p time points.
#' @param outlieMatrix, 1×p matrix, Outlier over p time points.
#' @param legend, character, Legend of table.
#' 
#' @return legend, character, Legend of table.
#' @return results, Table with pointwise variation depth values.
pointwiseVariationDepthOfOutlierFigure1<-function(nonOutlierMatrix,outlierMatrix,legend){
  p<-ncol(x=nonOutlierMatrix)
  timeRow<-createTimeRow(p=p)
  dataMatrix<-rbind(nonOutlierMatrix,outlierMatrix)
  pointwiseMedianMatrix<-createPointwiseMedianMatrix(dataMatrix=dataMatrix)
  VarRpointwiseMediant<-createVarRftMatrix(dataMatrix=dataMatrix,fMatrix=pointwiseMedianMatrix)
  VarRft<-createVarRftMatrix(dataMatrix=dataMatrix,fMatrix=outlierMatrix)
  VarERftRftMinusDelta<-createVarERftGivenRgtMinusDeltaMatrix(dataMatrix=dataMatrix,fMatrix=outlierMatrix)
  EVarRftRftMinusDelta<-createEVarRftGivenRftMinusDeltaMatrix(dataMatrix=dataMatrix,fMatrix=outlierMatrix)
  results<-round(x=rbind(t=timeRow,"Var(R_pointwiseMedian(t))"=c(VarRpointwiseMediant),"Var(R_f(t))"=c(VarRft),"Var(E(R_f(t)|R_f(t-Δ)))"=c(VarERftRftMinusDelta),"E(Var(R_f(t)|R_f(t-Δ)))"=c(EVarRftRftMinusDelta)),digits=3)
  return(list(legend=legend,results=results))
}

#' create functional boxplot
#' 
#' Detects outlier indices,
#' pointwise minimum of central region,
#' pointwise maximum of central region, and
#' functional median index,
#' using a functional boxplot with depth and factor for determining whiskers equal to factor,
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @references
#' Sun, Y., & Genton, M. G. (2011). Functional boxplots. Journal of Computational and Graphical Statistics, 20(2), 316–334. https://doi.org/10.1198/jcgs.2011.09224 
#' 
#' @details
#' Plots functional boxplot if plot==TRUE.
#' 
#' createFunctionalBoxplot(dataMatrix=dataMatrix,depth=depth,factor=1.5,indicesToRemove=indicesToRemove)$outlierIndicesInIncreasingOrderByDepth is the same as: 
#' .fbplotDetect(data=t(x=dataMatrix),depth=match.fun(FUN=paste("create",depth,"fColumn",sep=""))(dataMatrix=dataMatrix),shapeOutlier=indicesToRemove) in increasing order by depth value,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#' 
#' createFunctionalBoxplot(dataMatrix=dataMatrix,depth=depth,factor=factor,plot=TRUE) plots the same as
#' fbplot(fit=t(x=dataMatrix),depth=match.fun(FUN=paste("create",depth,"fColumn",sep=""))(dataMatrix=dataMatrix),plot=TRUE,factor=factor).
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param depth, character, "MSS", "TVD", "MBD", or "ED".
#' @param factor, integer, Factor for determining whiskers.
#' @param indicesToRemove, integer vector, Indices of curves to exclude from outlier detection.
#' @param plot, boolean, Set TRUE if plot, FALSE otherwise.
#' @param plotTitle, character, Title of plot.
#' 
#' @return pointwiseMinimumOfCentralRegionRow, numeric vector of length p, Pointwise minimum of central region.
#' @return pointwiseMaximumOfCentralRegionRow, numeric vector of length p, Pointwise maximum of central region.
#' @return outlierIndicesInIncreasingOrderByDepth, integer vector, Outlier indices in increasing order by depth.
#' @return functionalMedianIndex, integer, Functional median index.
createFunctionalBoxplot<-function(dataMatrix,fMatrix=dataMatrix,depth,factor,indicesToRemove=NULL,plot=FALSE,plotTitle=NULL){
  k<-nrow(x=fMatrix)
  originalIndices<-1:k
  originalfMatrix<-fMatrix
  depthValuesColumn<-match.fun(FUN=paste("create",depth,"fColumn",sep=""))(dataMatrix=dataMatrix,fMatrix=fMatrix)
  if(!is.null(x=indicesToRemove)){
    originalIndices<-originalIndices[-indicesToRemove]
    fMatrix<-fMatrix[-indicesToRemove,]
    depthValuesColumn<-depthValuesColumn[-indicesToRemove]
  }
  sizeOfCentralRegion<-ceiling(x=k/2)
  indicesInDecreasingOrderByDepth<-order(depthValuesColumn,decreasing=TRUE)
  centralRegionMatrix<-fMatrix[indicesInDecreasingOrderByDepth[1:sizeOfCentralRegion],]
  pointwiseMinimumOfCentralRegionRow<-apply(X=centralRegionMatrix,MARGIN=2,FUN=min)
  pointwiseMaximumOfCentralRegionRow<-apply(X=centralRegionMatrix,MARGIN=2,FUN=max)
  IQRRow<-pointwiseMaximumOfCentralRegionRow-pointwiseMinimumOfCentralRegionRow
  lowerFenceOfNonOutliersRow<-pointwiseMinimumOfCentralRegionRow-factor*IQRRow
  upperFenceOfNonOutliersRow<-pointwiseMaximumOfCentralRegionRow+factor*IQRRow
  outlierIndices<-which(x=0<apply(X=t(x=t(x=fMatrix)<lowerFenceOfNonOutliersRow)|t(x=upperFenceOfNonOutliersRow<t(x=fMatrix)),MARGIN=1,FUN=sum))
  numberOfOutliers<-length(x=outlierIndices)
  if(0<numberOfOutliers){
    indicesInIncreasingOrderByDepth<-order(depthValuesColumn)
    originalOutlierIndicesInIncreasingOrderByDepth<-originalIndices[indicesInIncreasingOrderByDepth[sort(x=match(x=outlierIndices,table=indicesInIncreasingOrderByDepth))]]
  }else{
    originalOutlierIndicesInIncreasingOrderByDepth<-NULL
  }
  originalFunctionalMedianIndex<-originalIndices[indicesInDecreasingOrderByDepth[1]]
  if(plot){
    p<-ncol(x=dataMatrix)
    timeRow<-createTimeRow(p=p)
    
    # plot pointwise minimum and pointwise maximum of non-outliers
    if(0<numberOfOutliers){
      nonOutlierMatrix<-fMatrix[-outlierIndices,]
    }else{
      nonOutlierMatrix<-fMatrix
    }
    pointwiseMinimumOfNonOutliersRow<-apply(X=nonOutlierMatrix,MARGIN=2,FUN=min)
    pointwiseMaximumOfNonOutliersRow<-apply(X=nonOutlierMatrix,MARGIN=2,FUN=max)
    ylim<-createylim(dataMatrix=fMatrix)
    matplot(x=timeRow,y=cbind(pointwiseMinimumOfNonOutliersRow,pointwiseMaximumOfNonOutliersRow),type="l",lty=1,col="black",xlab="t",ylab="X(t)",ylim=ylim,main=plotTitle)
    
    # color central region
    x=c(timeRow,rev(x=timeRow))
    y=c(pointwiseMinimumOfCentralRegionRow,rev(x=pointwiseMaximumOfCentralRegionRow))
    polygon(x=x,y=y,col="magenta",border="black")
    
    # plot outliers
    outlierMatrix<-fMatrix[outlierIndices,,drop=FALSE]
    matlines(x=timeRow,y=t(x=outlierMatrix),lty=2,col="red")
    
    # plot functional median
    functionalMedianRow<-originalfMatrix[originalFunctionalMedianIndex,]
    lines(x=timeRow,y=functionalMedianRow)
  }
  return(list(pointwiseMinimumOfCentralRegionRow=pointwiseMinimumOfCentralRegionRow,
              pointwiseMaximumOfCentralRegionRow=pointwiseMaximumOfCentralRegionRow,
              outlierIndicesInIncreasingOrderByDepth=originalOutlierIndicesInIncreasingOrderByDepth,
              functionalMedianIndex=originalFunctionalMedianIndex))
}

#' create functional boxplot using MSS+TVD
#' 
#' Detects shape and magnitude outlier indices,
#' pointwise minimum of central region,
#' pointwise maximum of central region, and
#' functional median index,
#' using a functional boxplot using MSS+TVD,
#' for k curves f_l,
#' over p time points t_i,
#' with f_l(t_i)=fMatrix[l,i],
#' based on n observations at p time points X_j(t_i)=dataMatrix[j,i],
#' for i∈{1,...,p}, j∈{1,...,n}, l∈{1,...,k}.
#' 
#' @details
#' Plots functional boxplot if plot==TRUE.
#' 
#' createFunctionalBoxplotMSSTVD(dataMatrix=dataMatrix)$shapeOutlierIndicesInIncreasingOrderByDepth is the same as:
#' detectOutlier(data=t(x=dataMatrix),nCurve=nrow(x=dataMatrix),nPoint=ncol(x=dataMatrix),empFactor=3)$sout,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#' 
#' createFunctionalBoxplotMSSTVD(dataMatrix=dataMatrix)$magnitudeOutlierIndicesInIncreasingOrderByDepth is the same as:
#' detectOutlier(data=t(x=dataMatrix),nCurve=nrow(x=dataMatrix),nPoint=ncol(x=dataMatrix),empFactor=3)$mout,
#' at https://github.com/hhuang90/TVD/blob/master/R/TVD.R.
#'
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points.
#' @param plot, boolean, Set TRUE if plot, FALSE otherwise.
#' @param plotTitle, character, Title of plot.
#'
#' @return pointwiseMinimumOfCentralRegionRow, numeric vector of length p, Pointwise minimum of central region.
#' @return pointwiseMaximumOfCentralRegionRow, numeric vector of length p, Pointwise maximum of central region.
#' @return shapeOutlierIndicesInIncreasingOrderByDepth, integer vector, Shape outlier indices in increasing order by depth.
#' @return magnitudeOutlierIndicesInIncreasingOrderByDepth, integer vector, Magnitude outlier indices in increasing order by depth.
#' @return functionalMedianIndex, integer, Functional median index.
createFunctionalBoxplotMSSTVD<-function(dataMatrix,fMatrix=dataMatrix,plot=FALSE,plotTitle=NULL){
  MSSFactor<-3
  TVDFactor<-1.5
  shapeOutlierIndicesInIncreasingOrderByMSS<-detectOutlierIndicesBoxplot(dataMatrix=dataMatrix,fMatrix=fMatrix,depth="MSS",factor=MSSFactor)
  functionalBoxplot<-createFunctionalBoxplot(dataMatrix=dataMatrix,fMatrix=fMatrix,depth="TVD",factor=TVDFactor,indicesToRemove=shapeOutlierIndicesInIncreasingOrderByMSS,plot=plot,plotTitle=plotTitle)
  pointwiseMinimumOfCentralRegionRow<-functionalBoxplot$pointwiseMinimumOfCentralRegionRow
  pointwiseMaximumOfCentralRegionRow<-functionalBoxplot$pointwiseMaximumOfCentralRegionRow
  magnitudeOutlierIndicesInIncreasingOrderByTVD<-functionalBoxplot$outlierIndicesInIncreasingOrderByDepth
  functionalMedianIndex<-functionalBoxplot$functionalMedianIndex
  return(list(pointwiseMinimumOfCentralRegionRow=pointwiseMinimumOfCentralRegionRow,
              pointwiseMaximumOfCentralRegionRow=pointwiseMaximumOfCentralRegionRow,
              shapeOutlierIndicesInIncreasingOrderByMSS=shapeOutlierIndicesInIncreasingOrderByMSS,
              magnitudeOutlierIndicesInIncreasingOrderByTVD=magnitudeOutlierIndicesInIncreasingOrderByTVD,
              functionalMedianIndex=functionalMedianIndex))
}

#' create Figure 4
#' 
#' Creates Figure 4.
#' 
#' @param dataMatrix, n×p matrix, n observations at p time points.
#' @param fMatrix, k×p matrix, k curves over p time points. 
#' @param depth, character, "MBD", "ED", or "MSSTVD".
createFigure4<-function(dataMatrix,fMatrix=dataMatrix,depth){
  p<-ncol(x=dataMatrix)
  timeRow<-createTimeRow(p=p)
  
  # plot samples
  ylim=createylim(dataMatrix=dataMatrix)
  matplot(x=timeRow,y=t(x=dataMatrix),type="l",lty=1,pch=1,col="lightgrey",xlab="",ylab="",xaxt="n",yaxt="n",ylim=ylim)
  title(xlab="t",line=0.5,cex.lab=0.5)
  title(ylab="X(t)",line=1,cex.lab=0.5)
  axis(side=1,line=0,cex.axis=0.5,padj=-3,tck=-0.05)
  axis(side=2,line=0,cex.axis=0.5,padj=2,tck=-0.05)
  
  # color central region
  if(depth=="MSSTVD"){
    functionalBoxplot<-createFunctionalBoxplotMSSTVD(dataMatrix=dataMatrix)
  }else{
    functionalBoxplot<-createFunctionalBoxplot(dataMatrix=dataMatrix,depth=depth,factor=1.5)
  }
  pointwiseMinimumOfCentralRegionRow<-functionalBoxplot$pointwiseMinimumOfCentralRegionRow
  pointwiseMaximumOfCentralRegionRow<-functionalBoxplot$pointwiseMaximumOfCentralRegionRow
  x=c(timeRow,rev(x=timeRow))
  y=c(pointwiseMinimumOfCentralRegionRow,rev(x=pointwiseMaximumOfCentralRegionRow))
  polygon(x=x,y=y,col="magenta",border=FALSE)
  
  # plot functional median
  functionalMedianIndex<-functionalBoxplot$functionalMedianIndex
  functionalMedianRow<-dataMatrix[functionalMedianIndex,]
  lines(x=timeRow,y=functionalMedianRow)
}