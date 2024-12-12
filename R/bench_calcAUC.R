#' Benchmarking: Calculate of the AUC of different cases
#'
#' This function calculates the AUC of different cases
#'
#' @param out_list The list of prediction probability matrix of each case
#' @param cluster The ground-truth
#' @param cell_types Character vector of cell types of the data
#'
#' @import pROC
#' @import dplyr
#'
#' @return A dataset of AUC per cell type of each case
#' @export
#'

bench_calcAUC<-function(out_list,cluster,cell_types){

  AUC_list<-list()

  cell_types<-as.character(cell_types)

  for (i in 1:length(out_list)) {
    AUC_i<-c()
    for (y in 1:length(cell_types)) {
      class_to_predict <- cell_types[y]
      print(class_to_predict)
      selected_class <- ifelse(cluster == class_to_predict, 1, 0)
      print(length(selected_class))
      predicted_probabilities<-out_list[[i]][,class_to_predict]
      print(length(predicted_probabilities))
      AUC_y<-auc(selected_class~predicted_probabilities)
      AUC_i<-c(AUC_i,AUC_y)
    }
    AUC_list[[i]]<-AUC_i
  }
  names(AUC_list)<-names(out_list)
  x<-as.data.frame(do.call(rbind, AUC_list))
  colnames(x)<-cell_types
  rownames(x)<-names(out_list)
  x$Method<-rownames(x)
  return(x)
}
