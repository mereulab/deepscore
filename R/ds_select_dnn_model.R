#' This function trains, fits and evaluates n deep neural network (DNN) models, and select the best one according to the selected metric
#'
#' This function creates n DNN model from the reference dataset (scRNA-seq) by using the keras package,
#' and selects the best model based on the user selected metric.
#' @param n_mod Number of model iterations
#' @param metric Which metric use for the model selection, `loss` or `accuracy`(Default="accuracy")
#' @param out The output of the ds_split_data_encoder function. The input data used to train the DNN model.
#' @param hnodes The number of the nodes of the intermediate layers of the encoder. The first and last layers are not included.
#' @param epochs Keras parameter for the layer_dense function (Default= 30).
#' @param batch_size Keras parameter for the layer_dense function (Default= 32).
#' @param activation Keras parameter for the layer_dense function (Default= "relu").
#' @param add_dropout Keras parameter for the layer_dense function (Default= TRUE).
#' @param pct_dropout Keras parameter for the layer_dense function (Default= 0.2).
#' @param name_mod name of the model. This is mandatory if you want to include the model in a stacked model.
#' @param lr Keras parameter that acts as a tuning parameter in the model that determines the step size at each iteration while
#' moving toward a minimum of a loss function (Default=0.001)
#' @param weight_reg Keras parameter for the layer weight regularizers (Default=TRUE)
#' @param l1 Keras parameter for the layer weight regularizers (Default=0)
#' @param l2 Keras parameter for the layer weight regularizers (Default=0)
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param earlystopping Keras parameter to stop training when a monitored metric has stopped improving (Default=TRUE)
#' @param lr_scheduler Keras parameter to modulate how the learning rate of your optimizer changes over time (Default=FALSE)
#' @param schedule Keras parameter to modulate how the learning rate of your optimizer changes over time (Default=NULL)
#' @param patience Keras parameter for the number of epochs with no improvement after which training will be stopped (Default=2)
#' @param ... Other keras parameters.
#' @import keras
#' @import glue
#' @import dplyr
#' @return The DNN model.
#' @export
#' @examples
#' #not run
#' out <- ds_split_data_encoder(features = features,clus = cluster,prop = 0.8,verbose = T)
#' model <- ds_select_dnn_model(n_mod=6,metric="accuracy",out = out,hnodes = c(1000),verbose = T,epochs = 10,batch_size = 32)

ds_select_dnn_model<-function(n_mod,metric="accuracy",out,hnodes, epochs=10,
                          batch_size=32, activation="relu", add_dropout=TRUE,
                          pct_dropout=0.2,name_mod="mod", lr=0.001,
                          weight_reg=TRUE, l1=0, l2=0,verbose = TRUE, earlystopping=TRUE, lr_scheduler=FALSE,
                          schedule=NULL, patience=2, ...){
  n_models<-list()
  loss_v<-c()
  accuracy_v<-c()
  for (n in 1:n_mod) {
    n_models[[n]]<-ds_dnn_model_1(out = out1,hnodes = hnodes,epochs = epochs,
                   batch_size = batch_size,activation = activation,add_dropout = add_dropout,
                   pct_dropout = pct_dropout,name_mod = paste0(name_mod,"_",n),lr=lr,
                   weight_reg = weight_reg,l1=l1,l2=l2,verbose=verbose,earlystopping = earlystopping,
                   lr_scheduler = lr_scheduler,schedule = schedule,patience = 2)
    quality_metrics<- n_models[[n]] %>% evaluate(x=out1$test_x,y=out1$test_y)
    quality_metrics<-unname(quality_metrics)
    loss_v<-c(loss_v,quality_metrics[1])
    accuracy_v<-c(accuracy_v,quality_metrics[2])
  }

  if(metric=="accuracy"){
    sel_idx<-which.max(accuracy_v)
    print(paste0("loss: ",loss_v[sel_idx]))
    print(paste0("accuracy: ",accuracy_v[sel_idx]))
    return(n_models[[sel_idx]])
  }

  if(metric=="loss"){
    sel_idx<-which.min(loss_v)
    print(paste0("loss: ",loss_v[sel_idx]))
    print(paste0("accuracy: ",accuracy_v[sel_idx]))
    return(n_models[[sel_idx]])
  }
}

