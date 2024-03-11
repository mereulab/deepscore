#' Heatmap comparing two cell type annotations
#'
#' This function generates a heatmap which compares the predicted cell types with another cell type annotation
#'
#' @param out The predicted cell types
#' @param cluster A cell type annotation
#' @param size Heatmap axis size
#'
#' @import ggplot2
#' @import hrbrthemes
#'
#' @return A heatmap comparing out and cluster
#' @export
#' @examples
#' #not run
#' heatmap<-ds_identity_heatmap(out,cluster)

ds_identity_heatmap <- function(out,cluster,size=14) {
  confusion<-table(cluster,out)
  confusion<-as.matrix(confusion)

  real<-as.data.frame(table(cluster))

  t_c<-as.data.frame(table(out))
  ord_c<-t_c[order(t_c$Freq,decreasing=TRUE),]
  ord_c<-ord_c$out
  t_r<-as.data.frame(table(cluster))
  ord_r<-t_r[order(t_r$Freq,decreasing=TRUE),]
  ord_r<-ord_r$cluster

  confusion<-confusion[,ord_c]
  confusion<-confusion[ord_r,]

  matrix.sort <- function(matrix) {
    if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
    row.max <- apply(matrix,1,which.max)
    if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
    matrix[names(sort(row.max)),]
  }

  m<-t(confusion)
  m<-matrix.sort(m)
  confusion<-t(m)

  ord_r<-rev(rownames(confusion))
  confusion<-confusion[ord_r,]
  if("unclassified" %in% unique(out)){
    ord_c<-colnames(confusion[,which(colnames(confusion) != "unclassified")])
    ord_c<-c(ord_c,"unclassified")
    confusion<-confusion[,ord_c]
  }

  confusion<-as.data.frame(confusion)


  colnames(confusion)<-c("Cluster","Predicted","Shared")
  order<-confusion$Cluster[1:nrow(real)]
  real<-real[match(order, real$cluster),]

  confusion$Proportion<-confusion$Shared/real$Freq
  confusion <- confusion[order(confusion$Shared,decreasing=TRUE),]

  if(nrow(confusion[is.na(confusion$Proportion),])!=0){
    confusion[is.na(confusion$Proportion),]$Proportion<-0
  }

  gg<- confusion %>%
    ggplot(aes(x=Predicted,y=Cluster,fill=Proportion)) +
    geom_tile()+
    scale_fill_distiller(palette = "RdPu") + theme_ipsum(axis_title_size = size,base_family = "Arial Narrow")  +
    theme(axis.text.x = element_text(angle = 90,size = size,hjust=0,vjust=0.2),
          axis.text.y = element_text(size=14)) +
    geom_text(aes(label=round(Shared,2)),size=3)+
    labs(title = paste0(""))+scale_x_discrete(position = "top")
  return(gg)
}

