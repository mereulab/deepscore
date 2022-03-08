

DeepScore <- function(hidden_nodes, common_features, n_labels, epochs=10, 
                 batch_size=32, activation="relu", dropout=TRUE, 
                 dropout_rate=0.2, batchnorm=TRUE, lr=0.001,
                 weight_reg=TRUE, l1=0, l2=0) {
        
    ds <- list(
        model = NULL,
        common_features = common_features,
        epochs = epochs,
        batch_size = batch_size,
        activation = activation,
        dropout = dropout,
        dropout_rate = dropout_rate,
        lr = lr,
        w_norm = NULL
        )
    
    n_features <- length(common_features)

    if (weight_reg) {
        ds$w_norm <- constraint_maxnorm(max_value = 2)
    }

    model <- keras_model_sequential()

    if (batchnorm) {
        model <- layer_batch_normalization(model)
    }
    
    for (i in seq_along(hidden_nodes)) {
        model <- model %>% 
                    layer_dense(hidden_nodes[i], activation=activation, input_shape = c(n_features),
                                kernel_constraint=ds$w_norm, name=glue("dense{hidden_nodes[i]}"),
                                kernel_regularizer=regularizer_l1_l2(l1=l1, l2=l2))
        if (dropout) {
            model <- layer_dropout(model, rate=dropout_rate)
        }
        if (batchnorm) {
            model <- layer_batch_normalization(model)
        }
        
        # N_features as input shape is only for the first layer
        n_features <- NULL
    }
    
    model <- model %>% 
                layer_dense(n_labels, activation="softmax", name="output")
    
    model %>% compile(
        optimizer = optimizer_adam(learning_rate=lr), 
        loss = 'categorical_crossentropy',
        metrics = c('accuracy')
    )
    
    #print(summary(model))
    
    ds$model <- model

    class(ds) <- append(class(ds), "DeepScore")
    return(ds)
}




set_reference <- function(ds, reference, labels, assay="RNA",
                          test_prop = 0.2) {

    reference <- reference[ds$common_features,]
    reference <- GetAssayData(reference, slot = "data", assay = assay)
    
    labels <- as.factor(unname(unlist(labels)))
    classes <- levels(labels)
    levels(labels) <- seq(1:length(levels(labels)))

    sizes <- unlist(lapply(levels(labels), function(x) length(labels[which(labels %in% x)])))
    names(sizes) <- levels(labels)

    coln <- colnames(reference)
    rown <- rownames(reference)
                           
    # Random sampling for training data but making sure we take from all clusters
    train.idx <- unlist(sapply(levels(labels), function(x) 
        sample(x = which(labels %in% x), 
        size = sizes[names(sizes) == x] * (1 - test_prop)))
        )
    train.data <- reference[, train.idx] %>% Matrix::t()
    train_y <- labels[train.idx]

    test.data <- reference[, -train.idx] %>% Matrix::t()
    test_y <- labels[-train.idx]
                       
    # Shuffling the training data
    model.train <- data.frame(train_y, train.data, check.names = F)
    model.train <- model.train[sample(seq(1:nrow(model.train))),]
    
    ds$split_data <- list(
        train_x = as.matrix(model.train[,-1]),
        train_y = to_categorical(model.train[,1]),
        test_x = test.data,
        test_y = to_categorical(test_y),
        classes = classes
        )
                
    return(ds)
}
                        
                        

train <- function(ds, plot = TRUE, earlystopping=TRUE, lr_scheduler=FALSE, 
              schedule=NULL, patience=2) {
    
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
    
    history <- ds$model %>% fit(
        x = ds$split_data$train_x,
        y = ds$split_data$train_y,
        epochs = ds$epochs,
        batch_size = ds$batch_size,
        validation_split = .2,
        callbacks=callbacks,
        verbose = 2
    )
    
    test_performance <- ds$model %>% evaluate(ds$split_data$test_x,
                                               ds$split_data$test_y)
    print("Model evaluation on unseen data:")
    print(test_performance)
    
    if (plot) {
        print(plot(history))
    }
    
    return(ds)
}
                        
                        
                        
annotate <- function(ds, query, assay="RNA") {
    
    query.data <- GetAssayData(query, slot = "data", assay = assay)[ds$common_features,]
    query.data <- Matrix::t(query.data)

    probs <- ds$model %>% predict(query.data)
    probs <- data.frame(probs, check.names = F)
    names(probs) <- c("unclassified", ds$classes)
    predicted <- apply(probs, 1, function(x) {
        ifelse(max(x) > 0.5, 
               names(probs)[which(x==max(x))],
               "unclassified"
              )
    })

    query$ds_annotation <- factor(predicted)

    return(query)
}
