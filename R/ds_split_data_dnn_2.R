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
#' @return A list with the following data: train_x= training_data,train_y= categorical cell type variables,
#' test_x =test_data,test_y=categorical celltype from the test data, classes=original cell types).
#' @export
#' @examples
#' #not run
#  out <- ds_split_data(scale.data = scaled,clus = cluster,genes = sel_gg,prop = 0.6)
#  model <- ds_dnn_model(out = out,hnodes = c(100),verbose = T,epochs = 30,batch_size = 32)


ds_split_data_dnn_2 <- function (scale.data, clus, genes, prop = 0.5,
                        verbose = TRUE)
{
  reference <- scale.data[genes,]

  clus <- as.factor(unname(unlist(clus)))
  classes <- levels(clus)
  levels(clus) <- seq(1:length(levels(clus)))

  sizes <- unlist(lapply(levels(clus), function(x) length(clus[which(clus %in% x)])))
  names(sizes) <- levels(clus)

  coln <- colnames(reference)
  rown <- rownames(reference)

  # Random sampling for training data but making sure we take from all clusters
  train.idx <- unlist(sapply(levels(clus), function(x)
    sample(x = which(clus %in% x),
           size = sizes[names(sizes) == x] * (1 - prop)))
  )
  train.data <- reference[, train.idx] %>% Matrix::t()
  train_y <- clus[train.idx]

  test.data <- reference[, -train.idx] %>% Matrix::t()
  test_y <- clus[-train.idx]

  # Shuffling the training data
  model.train <- data.frame(train_y, train.data, check.names = F)
  model.train <- model.train[sample(seq(1:nrow(model.train))),]

  split.data<- list(
    train_x = as.matrix(model.train[,-1]),
    train_y = to_categorical(model.train[,1]),
    test_x = test.data,
    test_y = to_categorical(test_y),
    classes = classes
  )

  return(split.data)
}
