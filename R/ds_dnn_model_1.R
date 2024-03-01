#' This function trains, fits and evaluates a deep neural network (DNN) model.
#'
#' This function creates a DNN model from the reference dataset (scRNA-seq) by using the keras package.
#' @param out The output of the ds_split_data_encoder function. The input data used to train the DNN model.
#' @param hnodes The number of the nodes of the intermediate layers of the encoder. The first and last layers are not included.
#' @param epochs Keras parameter for the layer_dense function (Default= 30).
#' @param set.seed Set a seed (Default=`TRUE`)
#' @param seed Set a seed and a generator
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
#' @import tensorflow
#' @return The DNN model.
#' @export
#' @examples
#' #not run
#' enc <- ds_encoder(data = scaled,genes = gg,dims = 2000,verbose = T,hnodes = c(5000))
#' features <- ds_get_features(enc = enc,data = scaled,genes=gg)
#' out <- ds_split_data_encoder(features = features,clus = cluster,prop = 0.8,verbose = T)
#' model <- ds_dnn_model(out = out,hnodes = c(1000),verbose = T,epochs = 10,batch_size = 32)

ds_dnn_model_1 <- function(out,hnodes, epochs=10,set.seed=TRUE,seed,
                         batch_size=32, activation="relu", add_dropout=TRUE,
                         pct_dropout=0.2,name_mod="mod", lr=0.001,
                         weight_reg=TRUE, l1=0, l2=0,verbose = TRUE, earlystopping=TRUE, lr_scheduler=FALSE,
                         schedule=NULL, patience=2, ...){

  library(keras)
  tensorflow::set_random_seed(seed)
  train_x <- out$train_x
  train_y <- out$train_y

  test_x <- out$test_x
  test_y <- out$test_y
  Nclass<-ncol(out$train_y)

  common_features<-colnames(train_x)

  ds <- list(
    model = NULL,
    common_features = common_features,
    epochs = epochs,
    batch_size = batch_size,
    activation = activation,
    add_dropout = add_dropout,
    pct_dropout = pct_dropout,
    lr = lr,
    w_norm = NULL
  )

  n_features <- length(common_features)

  if (weight_reg) {
    ds$w_norm <- constraint_maxnorm(max_value = 2)
  }

  model <- keras_model_sequential() %>%
    layer_dense(units = hnodes[1], activation = activation, input_shape = ncol(train_x),name = paste(name_mod,"dense_1",sep="_"))

  if(add_dropout==TRUE){
    model <- model %>%
      layer_dropout(rate = pct_dropout,name = paste(name_mod,"dropout_1",sep="_"))
  }

  N<-length(hnodes)
  for (i in 1:N) {
    model <- model %>%
      layer_dense(units=hnodes[i], activation=activation, input_shape = c(n_features),
                  kernel_constraint=ds$w_norm, name = paste(name_mod,"dense",i+1,sep="_"),
                  kernel_regularizer=regularizer_l1_l2(l1=l1, l2=l2))
    if (add_dropout) {
      model <- layer_dropout(model, rate=pct_dropout,name = paste(name_mod,"dropout",i+1,sep="_"))
    }
    # N_features as input shape is only for the first layer
    n_features <- NULL
  }

  model <- model %>%
    layer_dense(Nclass, activation="softmax", name=paste(name_mod,"out",sep="_"))

  model %>% compile(
    optimizer = optimizer_adam(learning_rate=lr),
    loss = 'categorical_crossentropy',
    metrics = c('accuracy')
  )

  callbacks <- c(callback_tensorboard(log_dir = "./deepscore_logs"),
                 callback_progbar_logger())

  if (earlystopping) {
    callbacks <- c(callbacks, callback_early_stopping(
      monitor = "val_loss",
      patience = patience,
      restore_best_weights = TRUE
    ))
  }

  if (lr_scheduler) {
    callbacks <- c(callbacks, callback_learning_rate_scheduler(schedule))
  }

  history <- model %>% fit(
    x = train_x,
    y = train_y,
    epochs = ds$epochs,
    batch_size = ds$batch_size,
    validation_split = .2,
    callbacks=callbacks,
    verbose = 2
  )

  if(verbose==TRUE){
    test_performance <- model %>% evaluate(test_x,
                                           test_y)
    print("Model evaluation:")
    print(test_performance)
  }



  return(model)

}

