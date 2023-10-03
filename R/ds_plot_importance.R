#' This function plots the importance of the features
#'
#' This function plots the n top importance features decreasing or increasing.
#'
#' @param data The output of `ds_feature_importance`.
#' @param n The number of features to be plotted (Default=30)
#' @param Decrease If the importance are plotted decreasing or increasing.
#' @param general_alteration If plot the absolut scores, if `TRUE`, the plot would not take into account if the
#' alteration means an increase or decrease. (Default=TRUE)
#' @param main The title of the plot
#'
#' @import ggplot2
#'
#' @return A feature importance plot
#'
#' @export
#'
#' @examples
#' ds_plot_importance(data,n=50,Decrease=TRUE,general_alteration=FALSE,main="Feature Importance Plot")
#'

ds_plot_importance<-function(data,n=30,Decrease=TRUE,general_alteration=FALSE,main=""){
  rownames(data)<-NULL
  data<-data[order(data$Score, decreasing=Decrease),]
  if(general_alteration){
    data$Score<-abs(data$Score)
  }
  data<-data[1:n,]
  Selected_features<-data$Feature
  plot<-ggplot(data, aes(x=reorder(Feature,Score), y=Score, fill=Score))+
    geom_bar(stat = "identity",position="dodge")+coord_flip()+
    ylab("Variable Importance Score") + xlab("") +
    ggtitle(main)+
    scale_fill_gradientn(colors = c("yellow","orange","purple"))+theme_minimal()
  output_list<-list()
  output_list[[1]]<-Selected_features
  output_list[[2]]<-plot
  names(output_list)<-c("Selected","Feature_plot")
  return(output_list)
}
