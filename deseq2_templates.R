# Prairie Genomics Suite - DESeq2 R Analysis Templates
# Modularized R scripts based on provided examples

# ================================
# PACKAGE LOADING AND SETUP
# ================================

load_required_packages <- function() {
  required_packages <- c(
    "DESeq2", "dplyr", "pheatmap", "ggplot2", "ComplexHeatmap", 
    "sva", "EnhancedVolcano", "ggVennDiagram", "plotly", "car", 
    "rgl", "MASS", "RColorBrewer", "ggrepel", "tidyverse"
  )
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg %in% c("DESeq2", "sva", "EnhancedVolcano", "ComplexHeatmap")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, ask = FALSE)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  return(TRUE)
}

# ================================
# DATA PREPARATION FUNCTIONS
# ================================

prepare_deseq_data <- function(count_matrix, metadata, design_formula = "~ Condition") {
  # Ensure count data is integer matrix
  count_matrix <- round(as.matrix(count_matrix))
  
  # Ensure sample names match between count matrix and metadata
  common_samples <- intersect(colnames(count_matrix), rownames(metadata))
  count_matrix <- count_matrix[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = as.formula(design_formula)
  )
  
  # Filter low-expression genes (following provided examples)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Additional filtering: genes expressed > 10 in at least n samples
  min_samples <- max(7, ncol(dds) * 0.1)  # At least 10% of samples or 7, whichever is higher
  keep <- rowSums(counts(dds) > 10) >= min_samples
  dds <- dds[keep,]
  
  return(dds)
}

apply_batch_correction <- function(count_matrix, metadata, batch_column = "Batch", condition_column = "Condition") {
  if (!batch_column %in% colnames(metadata)) {
    return(count_matrix)  # No batch correction if batch column not present
  }
  
  batch <- metadata[[batch_column]]
  condition <- metadata[[condition_column]]
  
  # Apply ComBat-seq batch correction
  adjusted_counts <- ComBat_seq(
    as.matrix(count_matrix), 
    batch = batch, 
    group = condition,
    shrink = TRUE, 
    shrink.disp = TRUE, 
    gene.subset.n = 10000
  )
  
  return(adjusted_counts)
}

# ================================
# CORE DESEQ2 ANALYSIS
# ================================

run_deseq2_analysis <- function(dds, reference_condition = NULL) {
  # Set reference level if specified
  if (!is.null(reference_condition)) {
    dds$Condition <- relevel(dds$Condition, ref = reference_condition)
  }
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  return(dds)
}

extract_results <- function(dds, contrast = NULL, alpha = 0.05) {
  if (is.null(contrast)) {
    res <- results(dds, alpha = alpha)
  } else {
    res <- results(dds, contrast = contrast, alpha = alpha)
  }
  
  # Order by adjusted p-value
  res_ordered <- res[order(res$padj),]
  
  # Merge with normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  res_data <- merge(
    as.data.frame(res_ordered), 
    as.data.frame(normalized_counts), 
    by = 'row.names', 
    sort = FALSE
  )
  names(res_data)[1] <- 'gene'
  
  return(res_data)
}

get_all_pairwise_contrasts <- function(dds) {
  conditions <- unique(dds$Condition)
  contrasts <- list()
  
  for (i in 1:(length(conditions)-1)) {
    for (j in (i+1):length(conditions)) {
      contrast_name <- paste0(conditions[j], "_vs_", conditions[i])
      contrast_vector <- c("Condition", as.character(conditions[j]), as.character(conditions[i]))
      contrasts[[contrast_name]] <- contrast_vector
    }
  }
  
  return(contrasts)
}

# ================================
# PCA ANALYSIS
# ================================

perform_pca_analysis <- function(dds, intgroup = "Condition", ntop = 500) {
  # Variance stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  
  # Standard 2D PCA
  pca_plot <- plotPCA(vsd, intgroup = intgroup, returnData = FALSE)
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  
  return(list(
    vsd = vsd,
    pca_plot = pca_plot,
    pca_data = pca_data
  ))
}

create_3d_pca_plot <- function(vsd, intgroup = "Condition", ntop = 500, confidence = 0.95) {
  # Calculate PCA manually for 3D
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  group <- colData(vsd)[[intgroup]]
  d <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    group = group,
    name = colnames(vsd)
  )
  
  # Create base 3D plot
  p <- plot_ly(
    data = d,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~group,
    mode = "markers",
    type = "scatter3d",
    marker = list(size = 8)
  ) %>%
    layout(
      title = "3D PCA Plot",
      scene = list(
        xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100, 1), "% variance")),
        yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100, 1), "% variance")),
        zaxis = list(title = paste0("PC3: ", round(percentVar[3] * 100, 1), "% variance"))
      )
    )
  
  return(p)
}

# ================================
# VISUALIZATION FUNCTIONS
# ================================

create_enhanced_volcano_plot <- function(results_data, title = "Volcano Plot", 
                                       fc_cutoff = 1.0, p_cutoff = 0.05,
                                       label_genes = NULL) {
  df <- results_data
  df$diffexpressed <- 'NO'
  df$diffexpressed[df$log2FoldChange > fc_cutoff & df$padj < p_cutoff] <- 'UP'
  df$diffexpressed[df$log2FoldChange < -fc_cutoff & df$padj < p_cutoff] <- 'DOWN'
  
  # Add gene labels if specified
  if (!is.null(label_genes)) {
    df$delabel <- ifelse(df$gene %in% label_genes, as.character(df$gene), NA)
  } else {
    df$delabel <- NA
  }
  
  volcano_plot <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), 
                                       col = diffexpressed, label = delabel)) + 
    geom_point(size = 2.5) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = 'gray', linetype = 'dashed') +
    geom_hline(yintercept = -log10(p_cutoff), col = 'gray', linetype = 'dashed') +
    scale_color_manual(values = c("blue", "gray", "red")) +
    labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
    ggtitle(title) +
    theme_bw() +
    guides(col = guide_legend(override.aes = aes(label = '')))
  
  if (!is.null(label_genes)) {
    volcano_plot <- volcano_plot + 
      geom_text_repel(max.overlaps = Inf, color = "black", box.padding = 0.5, 
                     size = 4, segment.curvature = 0.1, segment.color = "grey")
  }
  
  return(volcano_plot)
}

create_heatmap <- function(dds, results_data, top_n = 50, sample_subset = NULL) {
  # Get top differentially expressed genes
  top_genes <- head(results_data$gene[order(results_data$padj)], top_n)
  
  # Get variance stabilized data
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top_genes, ]
  
  # Subset samples if specified
  if (!is.null(sample_subset)) {
    mat <- mat[, sample_subset]
  }
  
  # Create annotation data frame
  df <- as.data.frame(colData(dds)[colnames(mat), "Condition", drop = FALSE])
  colnames(df) <- "Condition"
  
  # Create heatmap
  pheatmap(
    mat,
    annotation_col = df,
    scale = "row",
    fontsize_row = 7,
    fontsize = 8,
    cellheight = 8,
    cellwidth = 8,
    treeheight_col = 20,
    treeheight_row = 20,
    main = paste("Top", top_n, "Differentially Expressed Genes")
  )
}

# ================================
# VENN DIAGRAM ANALYSIS
# ================================

create_venn_diagram <- function(contrast_results, fc_cutoff = 0.5, p_cutoff = 0.05, direction = "both") {
  # Extract DEGs for each contrast
  deg_lists <- list()
  
  for (contrast_name in names(contrast_results)) {
    df <- contrast_results[[contrast_name]]
    
    if (direction == "up") {
      degs <- df$gene[df$log2FoldChange > fc_cutoff & df$padj < p_cutoff]
    } else if (direction == "down") {
      degs <- df$gene[df$log2FoldChange < -fc_cutoff & df$padj < p_cutoff]
    } else {  # both
      degs <- df$gene[abs(df$log2FoldChange) > fc_cutoff & df$padj < p_cutoff]
    }
    
    deg_lists[[contrast_name]] <- degs[!is.na(degs)]
  }
  
  # Create Venn diagram
  venn_plot <- ggVennDiagram(deg_lists, label_alpha = 0, color = "black") +
    scale_fill_gradient(low = "pink", high = "coral") +
    theme(legend.position = "none")
  
  return(list(plot = venn_plot, gene_lists = deg_lists))
}

# ================================
# EXPORT FUNCTIONS
# ================================

export_results <- function(results_list, output_dir = ".") {
  for (contrast_name in names(results_list)) {
    filename <- file.path(output_dir, paste0(contrast_name, "_DESeq2_results.csv"))
    write.csv(results_list[[contrast_name]], file = filename, row.names = FALSE)
  }
  return(TRUE)
}

# ================================
# ADVANCED VISUALIZATION FUNCTIONS
# ================================

load_visualization_packages <- function() {
  viz_packages <- c(
    "ComplexHeatmap", "circlize", "RColorBrewer", "viridis",
    "survival", "survminer", "ggpubr", "gridExtra", "cowplot",
    "igraph", "visNetwork", "networkD3", "GGally"
  )
  
  for (pkg in viz_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg %in% c("ComplexHeatmap")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, ask = FALSE)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  return(TRUE)
}

create_publication_heatmap <- function(expression_data, metadata, 
                                     top_genes = 50, 
                                     condition_column = "Condition",
                                     cluster_rows = TRUE, 
                                     cluster_cols = TRUE,
                                     scale_data = TRUE,
                                     annotation_columns = NULL,
                                     color_scheme = "RdYlBu") {
  
  load_visualization_packages()
  
  # Prepare data
  if (is.character(top_genes)) {
    # If top_genes is a vector of gene names
    selected_genes <- intersect(top_genes, rownames(expression_data))
  } else {
    # If top_genes is a number, select most variable genes
    gene_vars <- apply(expression_data, 1, var, na.rm = TRUE)
    selected_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_genes, nrow(expression_data))]
  }
  
  # Subset expression data
  mat <- as.matrix(expression_data[selected_genes, ])
  
  # Scale data if requested
  if (scale_data) {
    mat <- t(scale(t(mat)))
  }
  
  # Prepare annotations
  col_annotations <- data.frame(row.names = colnames(mat))
  col_annotations[[condition_column]] <- metadata[[condition_column]][match(colnames(mat), rownames(metadata))]
  
  if (!is.null(annotation_columns)) {
    for (col in annotation_columns) {
      if (col %in% colnames(metadata)) {
        col_annotations[[col]] <- metadata[[col]][match(colnames(mat), rownames(metadata))]
      }
    }
  }
  
  # Create color schemes
  condition_colors <- RColorBrewer::brewer.pal(
    min(length(unique(col_annotations[[condition_column]])), 11), 
    "Set3"
  )
  names(condition_colors) <- unique(col_annotations[[condition_column]])
  
  anno_colors <- list()
  anno_colors[[condition_column]] <- condition_colors
  
  # Create ComplexHeatmap
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = ifelse(scale_data, "Z-score", "Expression"),
    col = circlize::colorRamp2(
      c(-2, 0, 2), 
      RColorBrewer::brewer.pal(3, color_scheme)
    ),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      df = col_annotations,
      col = anno_colors,
      annotation_height = unit(0.5, "cm")
    ),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = nrow(mat) <= 100,
    show_column_names = ncol(mat) <= 50,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )
  
  return(ht)
}

create_survival_analysis <- function(survival_data, expression_data = NULL, 
                                    gene_of_interest = NULL,
                                    time_column = "time",
                                    event_column = "event",
                                    group_column = "group",
                                    split_method = "median") {
  
  load_visualization_packages()
  
  # Prepare survival data
  surv_df <- survival_data
  
  # If gene expression provided, create gene-based groups
  if (!is.null(expression_data) && !is.null(gene_of_interest)) {
    if (gene_of_interest %in% rownames(expression_data)) {
      gene_expr <- expression_data[gene_of_interest, ]
      
      if (split_method == "median") {
        threshold <- median(gene_expr, na.rm = TRUE)
        surv_df$gene_group <- ifelse(gene_expr > threshold, "High", "Low")
      } else if (split_method == "tertiles") {
        tertiles <- quantile(gene_expr, c(0.33, 0.67), na.rm = TRUE)
        surv_df$gene_group <- cut(gene_expr, 
                                 breaks = c(-Inf, tertiles[1], tertiles[2], Inf),
                                 labels = c("Low", "Medium", "High"))
      }
      group_column <- "gene_group"
    }
  }
  
  # Create survival object
  surv_formula <- as.formula(paste("Surv(", time_column, ",", event_column, ") ~", group_column))
  fit <- survival::survfit(surv_formula, data = surv_df)
  
  # Create survival plot
  p <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = RColorBrewer::brewer.pal(length(unique(surv_df[[group_column]])), "Set1"),
    title = paste("Survival Analysis:", 
                  ifelse(!is.null(gene_of_interest), gene_of_interest, "Groups")),
    xlab = "Time",
    ylab = "Survival Probability",
    legend.title = group_column,
    legend.labs = levels(as.factor(surv_df[[group_column]])),
    font.main = c(14, "bold"),
    font.x = c(12, "plain"),
    font.y = c(12, "plain"),
    font.tickslab = c(10, "plain")
  )
  
  # Cox proportional hazards model
  cox_formula <- as.formula(paste("Surv(", time_column, ",", event_column, ") ~", group_column))
  cox_fit <- survival::coxph(cox_formula, data = surv_df)
  
  return(list(
    plot = p,
    cox_model = cox_fit,
    survival_fit = fit,
    summary = summary(cox_fit)
  ))
}

create_pathway_network <- function(pathway_results, 
                                 gene_expression = NULL,
                                 top_pathways = 20,
                                 layout = "fr",
                                 node_size_by = "gene_count",
                                 color_by = "pvalue") {
  
  load_visualization_packages()
  
  # Prepare pathway data
  if (is.data.frame(pathway_results)) {
    pathways <- head(pathway_results[order(pathway_results[[color_by]]), ], top_pathways)
  } else {
    pathways <- pathway_results
  }
  
  # Create network nodes
  nodes <- data.frame(
    id = pathways$ID,
    label = pathways$Description,
    size = pathways[[node_size_by]],
    color = -log10(pathways[[color_by]]),
    stringsAsFactors = FALSE
  )
  
  # Create edges based on gene overlap (simplified)
  edges <- data.frame(from = character(), to = character(), weight = numeric())
  
  if ("geneID" %in% colnames(pathways)) {
    for (i in 1:(nrow(pathways)-1)) {
      for (j in (i+1):nrow(pathways)) {
        genes_i <- strsplit(pathways$geneID[i], "/")[[1]]
        genes_j <- strsplit(pathways$geneID[j], "/")[[1]]
        overlap <- length(intersect(genes_i, genes_j))
        if (overlap > 0) {
          edges <- rbind(edges, data.frame(
            from = pathways$ID[i],
            to = pathways$ID[j],
            weight = overlap
          ))
        }
      }
    }
  }
  
  # Create igraph object
  if (nrow(edges) > 0) {
    g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  } else {
    g <- igraph::graph_from_data_frame(data.frame(from = character(), to = character()), 
                                       vertices = nodes, directed = FALSE)
  }
  
  # Create visualization
  if (nrow(edges) > 0) {
    # Network plot with igraph
    coords <- igraph::layout_with_fr(g)
    
    p <- ggplot() +
      geom_segment(
        data = data.frame(
          x = coords[edges$from, 1],
          y = coords[edges$from, 2],
          xend = coords[edges$to, 1],
          yend = coords[edges$to, 2],
          weight = edges$weight
        ),
        aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
        color = "gray50"
      ) +
      geom_point(
        data = data.frame(
          x = coords[, 1],
          y = coords[, 2],
          nodes
        ),
        aes(x = x, y = y, size = size, color = color),
        alpha = 0.8
      ) +
      scale_color_viridis_c(name = paste("-log10(", color_by, ")")) +
      scale_size_continuous(name = node_size_by, range = c(3, 12)) +
      theme_void() +
      labs(title = "Pathway Network Analysis")
  } else {
    # Simple scatter plot if no edges
    p <- ggplot(nodes, aes(x = 1:nrow(nodes), y = color, size = size, color = color)) +
      geom_point(alpha = 0.8) +
      scale_color_viridis_c(name = paste("-log10(", color_by, ")")) +
      scale_size_continuous(name = node_size_by, range = c(3, 12)) +
      theme_minimal() +
      labs(title = "Pathway Significance", x = "Pathway Rank", y = paste("-log10(", color_by, ")"))
  }
  
  return(list(
    plot = p,
    network = g,
    nodes = nodes,
    edges = edges
  ))
}

# ================================
# MAIN WORKFLOW FUNCTION
# ================================

run_complete_deseq2_workflow <- function(count_matrix, metadata, 
                                       reference_condition = NULL,
                                       apply_batch_correction = FALSE,
                                       batch_column = "Batch",
                                       output_dir = ".",
                                       create_plots = TRUE) {
  
  # Load required packages
  load_required_packages()
  
  # Apply batch correction if requested
  if (apply_batch_correction && batch_column %in% colnames(metadata)) {
    count_matrix <- apply_batch_correction(count_matrix, metadata, batch_column)
  }
  
  # Prepare DESeq2 data
  dds <- prepare_deseq_data(count_matrix, metadata)
  
  # Run DESeq2 analysis
  dds <- run_deseq2_analysis(dds, reference_condition)
  
  # Get all pairwise contrasts
  contrasts <- get_all_pairwise_contrasts(dds)
  
  # Extract results for all contrasts
  results_list <- list()
  for (contrast_name in names(contrasts)) {
    results_list[[contrast_name]] <- extract_results(dds, contrasts[[contrast_name]])
  }
  
  # PCA analysis
  pca_results <- perform_pca_analysis(dds)
  
  # Export results
  export_results(results_list, output_dir)
  
  return(list(
    dds = dds,
    results = results_list,
    pca_results = pca_results,
    contrasts = contrasts
  ))
}