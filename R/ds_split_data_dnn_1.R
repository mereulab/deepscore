#' This function split the reference data in training and test set to use it as input for the ds_dnn_model function.
#'
#' @param scale.data scale.data A scaled/normalized matrix of gene expressions like in the `scale.data`
#' of the Seurat object. Rows are genes and columns are cells from the reference
#' dataset.
#' @param clus A named factor with the cell type annotations from the reference data.
#' @param genes The set of genes you want to use.
#' @param prop The proportion of cells for the training dataset (default=0.5).
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param scale_the_data Set if data has to be normalized. (default=TRUE)
#' @param train.samples Set train samples
#' @return A list with the following data: train_x= training_data,train_y= categorical cell type variables,
#' test_x =test_data,test_y=categorical celltype from the test data, classes=original cell types, train_samples=Set of train samples).
#' @export
#' @examples
#' #not run
#  out <- ds_split_data(scale.data = scaled,clus = cluster,genes = sel_gg,prop = 0.6)
#  model <- ds_dnn_model(out = out,hnodes = c(100),verbose = T,epochs = 30,batch_size = 32)

ds_split_data_dnn_1 <- function (scale.data, clus, genes, prop = NULL,
                                 verbose = TRUE, scale_the_data=FALSE,train.samples=NULL)
{
  classes <- levels(clus)
  levels(clus) <- seq(1:length(levels(clus)))
  features <- unique(genes)
  if (verbose)
    message("Splitting the reference into train and test datasets...")
  total <- 6
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  Sys.sleep(0.1)
  progress <- 1
  setTxtProgressBar(pb, progress)
  start.time <- Sys.time()
  sizes <- unlist(lapply(levels(clus), function(x) length(clus[which(clus %in%
                                                                       x)])))
  names(sizes) <- levels(clus)
  progress <- 2
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  if (is.null(prop)) {
    prop <- 0.5
  }
  if(is.null(train.samples)){
    train.sample <- unlist(sapply(levels(clus), function(x) sample(x = which(clus %in%
                                                                                x), size = sizes[names(sizes) == x] * prop)))
    s <- colnames(scale.data)[train.sample]
  }
  else{
    s <- train.samples
  }
  train.sample <- which(colnames(scale.data) %in% s)
  test.sample <- which(!colnames(scale.data) %in% s)
  train.data <- scale.data[rownames(scale.data) %in% features, train.sample]
  progress <- 3
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  test.data <- scale.data[rownames(scale.data) %in% features, test.sample]
  progress <- 4
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  out.train <- clus[which(names(clus) %in% s)]
  out.test <- clus[which(!names(clus) %in% s)]

  if(scale_the_data){
    var.train <- apply(train.data, 1, function(x) (x - min(x))/(max(x) - min(x)))
    var.test <- apply(test.data, 1, function(x) (x - min(x))/(max(x) - min(x)))
  }
  else{
    var.train <- t(scale.data[, train.sample])
    var.test <- t(scale.data[, test.sample])
  }
  progress <- 5
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)

  model.train <- data.frame(out.train, var.train,check.names = F)
  n <- nrow(model.train)
  shuffle_idx<-sample(x=seq(1:n),size = n)
  progress <- 6
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)

  # model.test <- data.frame(out.test, var.test,check.names = F)
  # n <- nrow(model.test)
  # model.test <- model.test[sample(x=seq(1:n),size = n),]
  library(keras)
  train_x <- as.matrix(model.train[,-1])[shuffle_idx,]
  train_y <- to_categorical(model.train[,1])[shuffle_idx,]
  test_x <- var.test
  test_y <- to_categorical(out.test)

  return(list(train_x=train_x,train_y=train_y,test_x =test_x,test_y=test_y,classes=classes,train_samples=s))
}
