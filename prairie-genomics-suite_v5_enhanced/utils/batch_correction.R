# Batch Correction Utilities for Prairie Genomics Suite v5
# Based on Emory methodology using ComBat-seq and limma methods

#' Extract batch information from sample names
#' @param sample_names Character vector of sample names
#' @param batch_info List containing detected batch patterns
#' @return Named vector of batch assignments
extract_batch_vector <- function(sample_names, batch_info) {
  
  if (is.null(batch_info) || length(batch_info) == 0) {
    return(NULL)
  }
  
  # Try each detected batch pattern
  for (pattern_name in names(batch_info)) {
    pattern_info <- batch_info[[pattern_name]]
    
    # Create batch vector for all samples
    batch_vector <- rep(NA, length(sample_names))
    names(batch_vector) <- sample_names
    
    # Extract batch information using regex
    batch_regex <- paste0(pattern_name, "([0-9]+)")
    
    for (i in seq_along(sample_names)) {
      sample <- sample_names[i]
      
      # Try to extract batch number
      if (grepl(batch_regex, sample, ignore.case = TRUE)) {
        batch_id <- gsub(paste0(".*", batch_regex, ".*"), "\\1", sample, ignore.case = TRUE)
        batch_vector[i] <- paste0(pattern_name, batch_id)
      } else {
        # If no batch pattern found, assign to default batch
        batch_vector[i] <- paste0(pattern_name, "1")
      }
    }
    
    # Check if we have meaningful batch separation
    unique_batches <- unique(batch_vector)
    if (length(unique_batches) > 1 && length(unique_batches) <= length(sample_names) / 2) {
      return(factor(batch_vector))
    }
  }
  
  return(NULL)
}

#' Apply ComBat-seq batch correction
#' @param count_data Matrix of count data (genes x samples)
#' @param batch_vector Factor vector indicating batch membership
#' @param condition_vector Factor vector indicating biological conditions
#' @param shrink Logical, whether to apply shrinkage
#' @return Batch-corrected count matrix
apply_combat_seq <- function(count_data, batch_vector, condition_vector, shrink = TRUE) {
  
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("sva package is required for ComBat-seq batch correction")
  }
  
  tryCatch({
    # Ensure we have valid batch and condition information
    if (length(batch_vector) != ncol(count_data)) {
      stop("Batch vector length must match number of samples")
    }
    
    if (length(condition_vector) != ncol(count_data)) {
      stop("Condition vector length must match number of samples")
    }
    
    # Apply ComBat-seq
    corrected_counts <- sva::ComBat_seq(
      counts = count_data,
      batch = batch_vector,
      group = condition_vector,
      shrink = shrink,
      shrink.disp = TRUE,
      gene.subset.n = min(10000, nrow(count_data))
    )
    
    # Ensure counts are non-negative integers
    corrected_counts <- round(corrected_counts)
    corrected_counts[corrected_counts < 0] <- 0
    
    return(corrected_counts)
    
  }, error = function(e) {
    stop(paste("ComBat-seq correction failed:", e$message))
  })
}

#' Apply limma removeBatchEffect
#' @param count_data Matrix of count data (genes x samples)
#' @param batch_vector Factor vector indicating batch membership  
#' @param condition_vector Factor vector indicating biological conditions
#' @return Batch-corrected count matrix (pseudo-counts)
apply_limma_batch_correction <- function(count_data, batch_vector, condition_vector) {
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package is required for batch correction")
  }
  
  tryCatch({
    # Convert to log2 scale for limma
    log_counts <- log2(count_data + 1)
    
    # Create design matrix for biological conditions
    design <- model.matrix(~ condition_vector)
    
    # Apply batch correction
    corrected_log <- limma::removeBatchEffect(
      log_counts,
      batch = batch_vector,
      design = design
    )
    
    # Convert back to count scale (approximate)
    corrected_counts <- 2^corrected_log - 1
    corrected_counts <- round(corrected_counts)
    corrected_counts[corrected_counts < 0] <- 0
    
    return(corrected_counts)
    
  }, error = function(e) {
    stop(paste("Limma batch correction failed:", e$message))
  })
}

#' Detect potential batch effects in data
#' @param count_data Matrix of count data (genes x samples)
#' @param sample_names Character vector of sample names
#' @param batch_vector Factor vector indicating batch membership (optional)
#' @return List with batch effect diagnostics
detect_batch_effects <- function(count_data, sample_names, batch_vector = NULL) {
  
  diagnostics <- list(
    batch_detected = FALSE,
    methods_used = c(),
    results = list()
  )
  
  # Method 1: Look for batch patterns in sample names
  batch_patterns <- c(
    "Batch", "batch", "Run", "run", "Seq", "seq",
    "Lane", "lane", "Chip", "chip", "Plate", "plate",
    "Rep", "rep", "Replicate", "replicate"
  )
  
  name_based_batches <- list()
  for (pattern in batch_patterns) {
    matches <- grep(pattern, sample_names, value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) {
      # Extract batch identifiers
      batch_ids <- gsub(paste0(".*", pattern, "([0-9]+).*"), "\\1", matches, ignore.case = TRUE)
      unique_ids <- unique(batch_ids)
      
      if (length(unique_ids) > 1 && length(unique_ids) <= length(sample_names) / 2) {
        name_based_batches[[pattern]] <- list(
          samples = matches,
          batch_ids = batch_ids,
          n_batches = length(unique_ids)
        )
      }
    }
  }
  
  if (length(name_based_batches) > 0) {
    diagnostics$batch_detected <- TRUE
    diagnostics$methods_used <- c(diagnostics$methods_used, "name_pattern")
    diagnostics$results$name_patterns <- name_based_batches
  }
  
  # Method 2: PCA-based batch detection (if we have batch information)
  if (!is.null(batch_vector) && length(unique(batch_vector)) > 1) {
    
    tryCatch({
      # Normalize data for PCA
      log_counts <- log2(count_data + 1)
      
      # Select most variable genes for PCA
      gene_vars <- apply(log_counts, 1, var)
      top_var_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, nrow(log_counts))]
      
      # Perform PCA
      pca_result <- prcomp(t(log_counts[top_var_genes, ]), center = TRUE, scale. = TRUE)
      
      # Test for association between PCs and batch
      pc_batch_associations <- c()
      for (i in 1:min(5, ncol(pca_result$x))) {
        pc_values <- pca_result$x[, i]
        
        # ANOVA test for association
        aov_result <- aov(pc_values ~ batch_vector)
        p_value <- summary(aov_result)[[1]][1, "Pr(>F)"]
        
        pc_batch_associations <- c(pc_batch_associations, p_value)
      }
      
      # If any of the first 5 PCs are significantly associated with batch
      if (any(pc_batch_associations < 0.05, na.rm = TRUE)) {
        diagnostics$batch_detected <- TRUE
        diagnostics$methods_used <- c(diagnostics$methods_used, "pca_association")
        diagnostics$results$pca_batch_association <- list(
          pc_p_values = pc_batch_associations,
          significant_pcs = which(pc_batch_associations < 0.05),
          variance_explained = summary(pca_result)$importance[2, 1:length(pc_batch_associations)]
        )
      }
      
    }, error = function(e) {
      # PCA analysis failed, continue without it
      diagnostics$results$pca_error <- e$message
    })
  }
  
  # Method 3: Hierarchical clustering batch detection
  if (!is.null(batch_vector) && length(unique(batch_vector)) > 1) {
    
    tryCatch({
      # Calculate sample distances
      log_counts <- log2(count_data + 1)
      sample_dists <- dist(t(log_counts))
      
      # Hierarchical clustering
      hclust_result <- hclust(sample_dists, method = "ward.D2")
      
      # Cut tree to get clusters
      n_clusters <- length(unique(batch_vector))
      cluster_assignments <- cutree(hclust_result, k = n_clusters)
      
      # Calculate how well clusters match batches
      # Use adjusted rand index or similar metric
      if (requireNamespace("cluster", quietly = TRUE)) {
        batch_cluster_agreement <- cluster::adjustedRandIndex(batch_vector, cluster_assignments)
        
        if (batch_cluster_agreement > 0.3) {  # Threshold for considering batch effect
          diagnostics$batch_detected <- TRUE
          diagnostics$methods_used <- c(diagnostics$methods_used, "hierarchical_clustering")
          diagnostics$results$clustering_batch_agreement <- batch_cluster_agreement
        }
      }
      
    }, error = function(e) {
      # Clustering analysis failed, continue without it
      diagnostics$results$clustering_error <- e$message
    })
  }
  
  return(diagnostics)
}

#' Visualize batch effects before and after correction
#' @param original_data Original count matrix
#' @param corrected_data Batch-corrected count matrix
#' @param batch_vector Batch assignments
#' @param condition_vector Biological conditions
#' @return List of plots for batch effect visualization
visualize_batch_effects <- function(original_data, corrected_data = NULL, 
                                   batch_vector, condition_vector) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for batch effect visualization")
  }
  
  plots <- list()
  
  # PCA plot of original data
  plots$pca_original <- create_batch_pca_plot(
    count_data = original_data,
    batch_vector = batch_vector,
    condition_vector = condition_vector,
    title = "PCA - Original Data"
  )
  
  # PCA plot of corrected data (if provided)
  if (!is.null(corrected_data)) {
    plots$pca_corrected <- create_batch_pca_plot(
      count_data = corrected_data,
      batch_vector = batch_vector,
      condition_vector = condition_vector,
      title = "PCA - Batch Corrected Data"
    )
  }
  
  return(plots)
}

#' Create PCA plot colored by batch and shaped by condition
#' @param count_data Count matrix
#' @param batch_vector Batch assignments
#' @param condition_vector Biological conditions
#' @param title Plot title
#' @return ggplot object
create_batch_pca_plot <- function(count_data, batch_vector, condition_vector, title) {
  
  # Log transform for PCA  
  log_counts <- log2(count_data + 1)
  
  # Select most variable genes
  gene_vars <- apply(log_counts, 1, var)
  top_var_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, nrow(log_counts))]
  
  # Perform PCA
  pca_result <- prcomp(t(log_counts[top_var_genes, ]), center = TRUE, scale. = TRUE)
  
  # Create data frame for plotting
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Batch = factor(batch_vector),
    Condition = factor(condition_vector),
    Sample = colnames(count_data)
  )
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Create plot
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(color = Batch, shape = Condition), size = 3, alpha = 0.8) +
    ggplot2::scale_color_brewer(type = "qual", palette = "Set1") +
    ggplot2::scale_shape_manual(values = c(16, 17, 18, 15, 3, 4, 8, 10)[1:length(unique(condition_vector))]) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  return(p)
}