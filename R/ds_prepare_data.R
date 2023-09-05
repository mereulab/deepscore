#' This function prepares the needed inputs to train the dnn model
#'
#' This function generates a list with the needed inputs to generate the dnn model.
#' @param ref.data The Seurat Object that contains the reference data
#' @param ref.assay The assay in the Seurat Object that contains the reference data
#' @param query.data The Seurat Object that contains the query data
#' @param query.assay The assay in the Seurat Object that contains the query data
#' @param markers Pre-defined marker genes. If missing, `markers_data` is required, otherwise, the next parameters must be empty.
#' @param markers_data The output of the R function `FindAllMarkers` provided by `"Seurat"`. If not provided, the function will automatically calculate it.
#' @param top_markers A parameter controlling the use of all features in `markers_data` or only the n top markers. If set to `TRUE`, an integer should be provided in the `ntop` parameter (default=FALSE).
#' @param ntop A parameter determining the number of top markers to be selected (default=100).
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

ds_prepare_data<-function(ref.data,ref.assay,query.data,query.assay,markers,markers_data,top_markers=F,ntop=100){
  DefaultAssay(ref.data)<-ref.assay

  ref.cluster<-ref.data@active.ident

  if(missing(markers)){

    if(missing(markers_data)){
      markers_data<-FindAllMarkers(ref.data,assay = ref.assay,only.pos = T)
    }

    if(top_markers){
      markers <- ms_top_markers(levels(markers_data$cluster),markers = markers_data,ntop = 100)
      markers<- unlist(markers)
      markers<-unique(markers)
    }
    else{
      markers<-markers_data$gene
      markers<-unique(markers)
    }
  }

  else{
    markers<-markers
  }

  ref.data <- ref.data@assays[[ref.assay]]@data
  genes <- intersect(rownames(ref.data),markers)

  gg <- rownames(ref.data)

  sel_gg <- intersect(markers,gg)

  ref.data<-as.matrix(ref.data)

  if(missing(query.data)){
    query.data<-ref.data@assays[[query.assay]]@data
    query.cluster<-ref.cluster
  }
  else{
    query.cluster<-Idents(query.data)
    query.data<-query.data@assays[[query.assay]]@data
  }

  query.data<-as.matrix(query.data)

  query.data<-query.data[sel_gg,]

  out<-list(ref.data=ref.data,query.data=query.data,ref.cluster=ref.cluster,query.cluster=query.cluster,markers=sel_gg)

  return(out)
}
