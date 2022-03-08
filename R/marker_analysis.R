

top_markers <- function(markers, ntop=10) {

    c_names <- levels(markers$cluster)

    top <- lapply(c_names, function(x) markers$gene[markers$cluster == x][1:ntop])
    top <- lapply(top, function(x) x[!is.na(x)])
    names(top) <- c_names
    return(top)
}
                  

find_common_variable_genes <- function(ref, sample, ref.assay="RNA",
                                       sample.assay="RNA", target_n_genes=4000) {

    DefaultAssay(ref) <- ref.assay
    DefaultAssay(sample) <- sample.assay

    n <- 10000
    common <- c()
    while (length(common) < target_n_genes) {
        
        print(glue('Looking for {n} HVG'))
        ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = n)
        sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = n)
        
        ref_genes <- VariableFeatures(ref)
        sample_genes <- VariableFeatures(sample)
        
        common <- intersect(ref_genes, sample_genes)
        print(glue('Found {length(common)} genes in common'))
        n <- n + 1000
    }

    return(common)
}
                  

find_refmarkers_in_variable_genes <- function(ref.markers, sample, sample.assay="RNA", target_n_genes=2000) {

    DefaultAssay(sample) <- sample.assay
    n_cat <- length(levels(ref.markers$cluster))
    n <- as.integer(round(target_n_genes*2/n_cat))

    common <- c()
    
    prev.n.variable <- 0
    max.reached <- FALSE
    
    while ((length(common) < target_n_genes) & (!max.reached)) {
        
        print(glue('Looking for {n} marker genes from each label'))
        ref_genes <- top_markers(ref.markers, ntop = n)
        ref_genes <- unname(unlist(ref_genes))

        sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = length(ref_genes))

        common <- intersect(ref_genes, VariableFeatures(sample))
        print(glue('Found {length(common)} genes in common'))
        n <- n + 100

        if (length(VariableFeatures(sample)) == prev.n.variable) {
            print('No more DEG found, returning last result')
            max.reached <- TRUE
        }
        
        prev.n.variable <- length(VariableFeatures(sample))
    }

    return(common)
}
                  
                  
find_refmarkers_in_genes <- function(ref.markers, sample, sample.assay="RNA", target_n_genes=2000) {

    DefaultAssay(sample) <- sample.assay
    n_cat <- length(levels(ref.markers$cluster))
    n <- as.integer(round(target_n_genes/n_cat))

    common <- c()
    max.reached <- FALSE
    
    while ((length(common) < target_n_genes) & (!max.reached)) {
        
        print(glue('Looking for {n} marker genes from each label'))
        ref_genes <- top_markers(ref.markers, ntop = n)
        ref_genes <- unname(unlist(ref_genes))

        common <- intersect(ref_genes, rownames(endocrine))
        print(glue('Found {length(common)} genes in common'))
        n <- n + 100
        
        if (length(ref.markers$gene) == length(ref_genes)) {
            print('No more markers to go through, returning last result')
            max.reached <- TRUE
        }
    }

    return(common)
}
                  
                  
summary_ggplot <- function(data, ylab, xlab) {
    
    df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
    df <- melt(cbind(pred_id = rownames(df), df), id = "pred_id")
    df$pred_id <- factor(df$pred_id, levels = levels(df$variable))

    gg <- ggplot(df, aes_string(x = "pred_id", y = "variable", fill = "value")) + 
        labs(x = xlab, y = ylab) + geom_tile(aes_string(fill = "value")) +
        geom_text(aes(label = round(value, 2))) + scale_fill_viridis_c() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16),
            axis.title = element_text(size = 16)
        )
    
    return(gg)
}

                  
                  
matchSCore2 <- function(ref.markers, query.markers,
                        xlab="Predicted labels", ylab="Reference labels") {
    
    score <- 0
    lab.ref <- seq(1:length(ref.markers))
    lab.obs <- seq(1:length(query.markers))

    JI_matrix <- vector()
    
    for (ref.cluster in ref.markers) {
        
        JI <- vector()
        
        for (q.cluster in query.markers) {
        
            I <- length(intersect(ref.cluster, q.cluster))
            J <- I / (length(ref.cluster) + length(q.cluster) - I)
            JI <- append(JI, J, length(JI))
        }
        
        score <- sum(score, max(JI))
        JI_matrix <- rbind(JI_matrix, JI)
    }
    
    score <- score / length(ref.markers)
    colnames(JI_matrix) <- names(query.markers)
    rownames(JI_matrix) <- names(ref.markers)

    gg <- summary_ggplot(data = JI_matrix, ylab, xlab)

    return(
        list(
            matchScore = score,
            labels = apply(JI_matrix, 2, which.max),
            max_JI = apply(JI_matrix, 2, max),
            JI.mat = JI_matrix,
            ggplot = gg
        )
    )
}
                  
                  
                  
