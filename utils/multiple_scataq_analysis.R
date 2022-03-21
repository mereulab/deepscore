### A collection of function to perform downstream analysis on several
### scATACseq samples at the same time.

get_peaks <- function(folder) {
  
  peak_path <- glue("{folder}/peaks.bed")
  peaks <- read.table(peak_path, col.names = c("chr", "start", "end"))
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}



get_peaks_from_object <- function(folder_file, assay="peaks") {
    
    folder = folder_file[1]
    object_file = folder_file[2]
  
    assign('seu', get(load(glue("{folder}/{object_file}"))))
    DefaultAssay(seu) <- assay
    return(granges(seu))
}



build_combined_peak_seurat <- function(folder_file, combined.peaks, 
                                       assay="peaks", drop_predictions = F) {
  
    folder <- folder_file[1]
    object_file <- folder_file[2]
    file_path <- glue("{folder}/{object_file}")
    print(file_path)
    
    if (file.exists(file_path)) {
  
        if (grepl(".rds", file_path, fixed = TRUE)) {
            
            extension <- ".rds"
            seu <- readRDS(file_path)

        } else if (grepl(".RData", file_path, fixed = TRUE)) {
            
            extension <- ".RData"
            assign('seu', get(load(file_path)))

        } else {
            print("Please format your file as .rds or .RData")
        }
    }

    DefaultAssay(seu) <- assay

    frag_file <- list.files(path=toString(folder), pattern = "fragments.tsv.gz$")
    frags <- CreateFragmentObject(path = glue("{folder}/{frag_file}"), 
                                  cells = rownames(seu@meta.data))

    counts <- FeatureMatrix(fragments = frags, features = combined.peaks,
                          cells = rownames(seu@meta.data))
    print("Feature Matrix has been created")
    seu[[assay]] <- CreateChromatinAssay(counts, fragments = frags)
        
    # Drop metadata values starting with prediction
    if (drop_predictions) {
      for (metacol in colnames(seu@meta.data)) {
        if (substr(metacol, 1, 10) == "prediction") {
          seu[[metacol]] <- NULL
        }
      }
    }    

    out_file <- gsub(extension, glue("_combpeaks{extension}"), object_file)
        
    if (extension == ".rds") {
        saveRDS(seu, file = glue("{folder}/{out_file}"))
    } else if (extension == ".RData") {
        save(seu, file = glue("{folder}/{out_file}"))
    }
    
    return(seu)
}