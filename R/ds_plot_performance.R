#' Plot performance as radarchart/spider plot
#'
#' This function plots the performance calculated by `ds_performance()`
#'
#' @param data Output of `ds_performance()`
#' @param parameters The parameters that will be plotted
#' @param legend Plot with legend (Default=TRUE)
#' @param colors Vectors of colors for the plot
#' @param plwd Line width of the plot
#'
#' @return A radarchart plot
#' @export
#' @examples
#' #not run
#' ds_plot_performance(data = performance$Performance,
#'   parameters = c("Sensitivity","Specificity","Accuracy","Classification Rate"),
#'   colors = colors_border1,legend = FALSE)

ds_plot_performance<-function(data,parameters,legend=TRUE,colors,plwd=1){
  n<-length(parameters)
  radarchart( rbind(rep(1,n) , rep(0,n) , data[,parameters]) , axistype=1 ,
              #custom polygon
              pcol=colors  , plwd=3 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="white", cglwd=0.8,
              #custom labels
              vlcex=1)
  if(legend){
    legend(x=2, y=1, legend = rownames(package_ct[[1]]),
           bty = "n", pch=20 , col=colors_border1 , text.col = "grey1",
           cex=1.2, pt.cex=3)
  }
}
