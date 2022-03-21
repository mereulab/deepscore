#!/usr/bin/env Rscript

library(Seurat)
library(glue)
library(reticulate)
library(sceasy)

file_path <- toString(commandArgs(TRUE)[1])
ass <- toString(commandArgs(TRUE)[2])
conda_env <- toString(commandArgs(TRUE)[3])
drop_predictions <- toString(commandArgs(TRUE)[4])

use_condaenv(conda_env)
sc <- import("scanpy", convert = FALSE)

if (drop_predictions == "TRUE") drop_predictions <- TRUE
if (drop_predictions == "FALSE") drop_predictions <- FALSE
if (drop_predictions == "NA") drop_predictions <- TRUE
print(glue("Drop predictions is: {drop_predictions}"))

if (ass == "NA") ass <- "RNA"
print(glue("Selected assay is: {ass}"))

if (file.exists(file_path)) {
  
  if (grepl(".rds", file_path, fixed = TRUE)) {
    print("Converting from rds to h5ad")
    out_path <- gsub(".rds", glue("_{ass}.h5ad"), file_path)
    seu <- readRDS(file_path)
    
  } else if (grepl(".RData", file_path, fixed = TRUE)) {
    print("Converting from RData to h5ad")
    out_path <- gsub(".RData", glue("_{ass}.h5ad"), file_path)
    assign('seu', get(load(file_path)))
    
  } else {
    print("Please format your file as .rds or .RData")
  }
  
  # Drop metadata values starting with prediction
  if (drop_predictions) {
    for (metacol in colnames(seu@meta.data)) {
      if (substr(metacol, 1, 10) == "prediction") {
        seu[[metacol]] <- NULL
      }
    }
  }
  
  # Delete repeated metadata variables for the features in case of ATAC
  for (assay in Assays(seu)) {
    seu@assays[[assay]]@meta.features <- seu@assays[[assay]]@meta.features[, 
                          !duplicated(colnames(seu@assays[[assay]]@meta.features))]
  }

  #seu <- DietSeurat(seu, assays = ass, dimreducs = names(seu@reductions), 
  #                  graphs = names(seu@graphs))
  convertFormat(seu, from="seurat", to="anndata", main_layer="counts", assay=ass,
                drop_single_values=FALSE, outFile=out_path)
} else {
  print('File not found')
}
