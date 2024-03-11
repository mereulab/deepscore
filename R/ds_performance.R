#' Performance of the model
#'
#' This function calculates the performance of the model by generating a confusion table between the predicted and real data
#'
#' @param out The predicted cell types
#' @param cluster A cell type annotation
#'
#' @return A list containing a performance dataframe and the confusion table
#' @import pracma
#' @export
#' @examples
#' #not run
#' ds_performance(out,cluster)

ds_performance<-function(out,cluster){
  cluster<-as.factor(cluster)
  out<-as.factor(out)
  L<-append(levels(cluster),levels(out))
  L<-unique(L)
  out<-factor(out,levels = L)
  confusion<-table(cluster,out)
  ct_dt<-data.frame(row.names = levels(cluster))
  count_cluster<-as.data.frame(table(cluster))
  Un<-levels(out)[which(levels(out) == "unclassified")]

  Sensitivity<-c()
  Specificity<-c()
  Accuracy<-c()
  Classification_rate<-c()
  Error_rate<-c()
  for (i in 1:nrow(confusion)) {
    for (j in 1:ncol(confusion)) {
      if(colnames(confusion)[j] == rownames(confusion)[i]){
        se<-confusion[i,j]/sum(confusion[i,])
        Sensitivity<-append(Sensitivity,se)
        sp<-sum(count_cluster$Freq[-i])/(sum(count_cluster$Freq[-i])+sum(confusion[-i,j]))
        Specificity<-append(Specificity,sp)
        Accuracy<-append(Accuracy,(se+sp)/2)
        if(isempty(Un)==FALSE){
          Error_rate<-append(Error_rate,sum(confusion[,which(colnames(confusion)!=Un)][i,-j])/count_cluster$Freq[i])
        }
        else{
          Error_rate<-append(Error_rate,sum(confusion[i,-j])/count_cluster$Freq[i])
        }
      }
    }
    if(nrow(confusion)==ncol(confusion) || isempty(Un)){
      Classification_rate<-rep(1,nrow(confusion))
    }
    else{
      Classification_rate<-append(Classification_rate,1-(confusion[i,Un]/count_cluster$Freq[i]))
    }
  }
  ct_dt$Sensitivity<-Sensitivity
  ct_dt$Specificity<-Specificity
  ct_dt$Accuracy<-Accuracy
  ct_dt$`Classification Rate`<-Classification_rate
  ct_dt$`Error Rate`<-Error_rate
  ct_dt<-round(ct_dt,3)
  out<-list("Confusion_table"=confusion, "Performance"=ct_dt)
  return(out)
}
