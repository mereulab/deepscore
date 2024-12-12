#' Plot performance as radarchart/spider plot
#'
#' This function plots the performance calculated by `ds_performance()`
#'
#' @param data Output of `ds_performance()`
#' @param parameters The parameters that will be plotted
#' @param legend Plot with legend (Default=TRUE)
#' @param colors Vectors of colors for the plot
#' @param plwd Line width of the plot
#' @param xlim Customize x limits
#' @param ylim Customize y limits
#' @param title Title of the plot
#'
#' @import fmsb
#' @return A radarchart plot
#' @export
#' @examples
#' #not run
#' ds_plot_performance(data = performance$Performance,
#'   parameters = c("Sensitivity","Specificity","Accuracy","Classification Rate"),
#'   colors = colors_border1,legend = FALSE)

ds_plot_performance<-function(data,parameters,legend=TRUE,colors,plwd=1,xlim=c(-1,3.5),ylim=NULL,title=NULL){
  n<-length(parameters)
  colors<-colors[rownames(data)]
  radarchart(rbind(rep(1,n) , rep(0,n) , data[,parameters]) , axistype=1 ,
              #custom polygon
              pcol=colors  , plwd=3 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="white", cglwd=0.8,
              #custom labels
              vlcex=1,
              xlim = xlim,ylim=ylim,title = title)
  if(legend){
    legend(x=2, y=1, legend = rownames(data),
           bty = "n", pch=20 , col=colors , text.col = "grey1",
           cex=1.2, pt.cex=3)
  }
}
