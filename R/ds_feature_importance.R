#' This function calculates the feature importance of the deep learning model
#'
#' This function generates a dataframe with the feature importance score of the input features calculated
#' based on the "Permutation Feature Importance" approach.
#'
#' @param X The data with the feature information
#' @param Y The data with the output information
#' @param model The deep learning model
#' @param metric The evaluation metric used to calculate the feature importance, "loss" or "accuracy". (Default="accuracy")
#' @param features The manually selected features whose importance will be calculated. If not provided,
#' the function will use all the features. You can control the number of used features by using the parameter
#' `max_n`.
#' @param  max_n Number of features selected from X if `features` not provided.
#'
#' @import keras
#' @import tensorflow
#'
#' @return A dataframe with the importance score of each feature
#'
#' @export
#' @examples
#' ds_feature_importance(X=x_test,Y=y_test,metric="accuracy",max_n=500)
#'

ds_feature_importance<-function(X,Y,model,metric="accuracy",features,max_n=1000){
  test_metrics<- model %>% evaluate(x=X,y=Y)
  importance<-c()
  if(missing(features)){
    if(missing(max_n)){
      features<-colnames(X)
    }
    else{
    features<-colnames(X)[1:max_n]
    }
  }
  if(metric=="loss"){
    baseline<-test_metrics[1]
    num_permutations<-length(features)
    for (i in 1:num_permutations) {
      xi = rep(0,nrow(X))
      Xi = cbind(xi, X[,-which(colnames(X) == features[i])])
      permuted_metrics <- model %>% evaluate(x = Xi, y = Y,verbose = 0)
      importance <- c(importance,(((permuted_metrics[1]-baseline)/baseline))*100)
      print(paste0(i,' --- ',features[i]))
    }
  }
  if(metric=="accuracy"){
    baseline<-test_metrics[2]
    num_permutations<-length(features)
    for (i in 1:num_permutations) {
      xi = rep(0,nrow(X))
      Xi = cbind(xi, X[,-which(colnames(X) == features[i])])
      permuted_metrics <- model %>% evaluate(x = Xi, y = Y,verbose = 0)
      importance <- c(importance,(((permuted_metrics[2]-baseline)/baseline))*100)
      print(paste0(i,' --- ',features[i]))
    }
  }
  PFI_final = data.frame("Feature"=features,"Score" = importance)
  PFI_final$Score<-abs(PFI_final$Score)
  PFI_final <- PFI_final[order(PFI_final$Score, decreasing = TRUE),]
  return(PFI_final)
}
