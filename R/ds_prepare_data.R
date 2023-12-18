#' This function prepares the needed inputs to train the dnn model
#'
#' This function generates a list with the needed inputs to generate the dnn model.
#' @param ref.data The Seurat Object that contains the reference data
#' @param ref.assay The assay in the Seurat Object that contains the reference data
#' @param query.data The Seurat Object that contains the query data
#' @param query.assay The assay in the Seurat Object that contains the query data
#' @param markers Pre-defined marker genes. If missing, `markers_data` is required, otherwise, the next parameters must be empty.
#' @param markers_data The output of the R function `FindAllMarkers` provided by `"Seurat"`. If not provided, the function will automatically calculate it.
#'
#' @import matchSCore2
#' @import Seurat
#' @import Signac
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import tensorflow
#' @return The list containing all input data needed to train the dnn model.
#' @export
#' @examples
#' #not run
#' no_markers_no_top_markers<-ds_prepare_data(ref.data,ref.assay,query.data,query.assay,markers_data,top_markers=F)
#' no_markers_top_markers<-ds_prepare_data(ref.data,ref.assay,query.data,query.assay,markers_data,top_markers=T,ntop=100)
#' markers<-ds_prepare_data(ref.data,ref.assay,query.data,query.assay,markers)
#'

ds_prepare_data<-function(ref.data,ref.assay,query.data,query.assay,markers,markers_data){
  DefaultAssay(ref.data)<-ref.assay
  DefaultAssay(query.data)<-query.assay

  ref.cluster<-ref.data@active.ident
  query.cluster<-query.data@active.ident

  if(missing(markers)){

    if(missing(markers_data)){
      markers_data<-FindAllMarkers(ref.data,assay = ref.assay,only.pos = T)
      markers<-markers_data$gene
      markers<-unique(markers)
    }

    else{
      markers<-markers_data$gene
      markers<-unique(markers)
    }
  }

  gg <- rownames(ref.data)

  sel_gg <- intersect(markers,gg)

  query_gg<-rownames(query.data)
  sel_gg<-intersect(query_gg,sel_gg)

  ref.data <- ref.data[sel_gg,]
  ref.data <- GetAssayData(ref.data, slot = "data", assay = ref.assay)
  ref.data<-as.matrix(ref.data)

  query.data<-query.data[sel_gg,]
  query.data <- GetAssayData(query.data, slot = "data", assay = query.assay)
  query.data<-as.matrix(query.data)

  out<-list(ref.data=ref.data,query.data=query.data,ref.cluster=ref.cluster,query.cluster=query.cluster,markers=sel_gg)

  return(out)
}
