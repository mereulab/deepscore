#' This function plots the correlation matrix between two datasets
#'
#' This function generates a correlation plot
#'
#' @param x.data The Seurat Object that contains the first data
#' @param x.assay The assay in the Seurat Object that contains the first data
#' @param y.data The Seurat Object that contains the second data
#' @param y.assay The assay in the Seurat Object that contains the second data
#' @param downsample The number of cells to downsample in each class
#' @param features The features to take into account in the correlation matrix
#'
#' @import corrplot
#' @import viridisLite
#' @import Seurat

ds_correlation_plot<-function(x.data,x.assay,y.data,y.assay,features,palette="turbo",method="circle"){

  DefaultAssay(x.data)<-x.assay
  DefaultAssay(y.data)<-y.assay

  x.small <- subset(x.data,cells=colnames(x.data))
  y.small <- subset(y.data,cells=colnames(x.small))

  DefaultAssay(x.small)<-x.assay
  ref_x<-x.small
  DefaultAssay(y.small)<-y.assay
  ref_y<-y.small

  if(missing(features)){
    common_features<-intersect(rownames(ref_x),rownames(ref_y))
    ref_x<-subset(ref_x,features = common_features)
    ref_y<-subset(ref_y,features = common_features)
  }else{
    ref_x<-subset(ref_x,features = features)
    ref_y<-subset(ref_y,features = features)
    common_features<-intersect(rownames(ref_x),rownames(ref_y))
    ref_x<-subset(ref_x,features = common_features)
    ref_y<-subset(ref_y,features = common_features)
  }

  ref_x_data<-AverageExpression(
    ref_x,
    assays = x.assay)
  ref_x_data<-scale(ref_x_data[[1]])

  ref_y_data<-AverageExpression(
    ref_y,
    assays = y.assay)
  ref_y_data<-scale(ref_y_data[[1]])

  cormat<-cor(ref_x_data,ref_y_data)

  return(corrplot(cormat, type="upper",method = method,tl.col="black",col=get(palette)(100)))
}










