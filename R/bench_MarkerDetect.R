#' Benchmarking of the marker detection ability
#'
#' This function calculates the marker detection rate per case
#'
#' @param input_list List combining several outputs of ds_feature_importance
#' @param ref_list List of reference markers
#' @param xlab The text for the x axis name
#' @param ylab The text for the y axis name
#' @param title The text for the title
#' @param legend_title The text fot the legend title
#' @param round_by Integer indicating the number of decimal places
#'
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#'
#' @return A heatmap with the marker detection rate per case
#'
#' @export
#'
#'

bench_MarkerDetect<-function(threshold,input_list,ref_list,xlab=NULL,ylab=NULL,title=NULL,legend_title=NULL,round_by=3){
  for (i in 1:length(input_list)) {
    Filter<-c()
    for (s in 1:nrow(input_list[[i]])) {
      if(input_list[[i]][s,]$Score >= threshold){
        Filter<-c(Filter,"Y")
      }
      else{
        Filter<-c(Filter,"N")
      }
    }
    input_list[[i]]$Filter<-Filter
  }

  input_list_Y<-list()
  for (i in 1:length(input_list)) {
    input_list_Y[[i]]<-input_list[[i]][which(input_list[[i]]$Filter=="Y"),]
  }

  names(input_list_Y)<-names(input_list)

  DM<-c()
  for (i in 1:length(input_list_Y)) {
    DM_vec<-c()
    for (y in 1:length(ref_list)) {
      DM_vec<-c(DM_vec,length(ref_list[[y]][which(ref_list[[y]] %in% input_list_Y[[i]]$Feature)])/length(ref_list[[y]]))
    }
    DM<-rbind(DM,DM_vec)
  }
  DM<-as.data.frame(DM)
  colnames(DM)<-names(ref_list)
  rownames(DM)<-names(input_list_Y)
  DM$Case<-rownames(DM)

  print("Number of detected markers per case:")
  for (i in 1:length(ref_list)) {
    print(paste0(names(ref_list)[i]," - ",length(ref_list[[i]])))
  }

  data_long<-melt(DM)
  data_long$variable<-factor(data_long$variable,levels = names(ref_list))
  data_long$Case<-factor(data_long$Case,levels = names(input_list_Y))

  plot<-ggplot(data_long, aes(x = variable, y = Case, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, round_by)), color = "black", size = 4) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(x = xlab, y = ylab, fill = legend_title,title = title)

  return(list("Marker_Detection"=DM,"plot"=plot))
}
