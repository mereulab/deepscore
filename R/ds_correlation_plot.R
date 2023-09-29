#' This function plots the correlation matrix between two datasets
#'
#' This function generates a correlation plot
#'
#' @param x.data The Seurat Object that contains the first data
#' @param x.assay The assay in the Seurat Object that contains the first data
#' @param y.data The Seurat Object that contains the second data
#' @param y.assay The assay in the Seurat Object that contains the second data
#' @param downsample If a downsample is applied or not to the data. If `TRUE`,`n_downsample' is needed.(Default=TRUE)
#' @param n_downsample The number of cells to downsample in each class
#' @param markers The markers to take into account in the correlation matrix
#'
#' @import corrplot
#' @import viridisLite
#' @import Seurat

ds_correlation_plot<-function(x.data,x.assay,y.data,y.assay,downsample,n_downsample,features){
  if(downsample){
    x.small <- subset(x.data,downsample=n_downsample)
    y.small <- subset(y.data,downsample=n_downsample)
  }
  DefaultAssay(x.small)<-x.assay
  ref_x<-x.small
  DefaultAssay(y.small)<-y.assay
  ref_y<-y.small

  ref_x<-subset(ref_x,features = features)
  ref_y<-subset(ref_y,features=features)
  common_features<-intersect(rownames(ref_x),rownames(ref_y))
  ref_x<-subset(ref_x,features = common_features)
  ref_y<-subset(ref_y,features = common_features)

  ref_x_data<-AverageExpression(
    ref_x,
    assays = x.assay)
  ref_x_data<-scale(ref_x_data[[1]])

  ref_y_data<-AverageExpression(
    ref_y,
    assays = y.assay)
  ref_y_data<-scale(ref_y_data[[1]])

  plot<-cormat<-cor(ref_y_data,ref_x_data)

  return(corrplot(plot, type="upper",tl.col="black",col=mako(100)))
}










