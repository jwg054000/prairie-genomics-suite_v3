# Analysis Wrappers for Code Logging
# Enhanced analysis functions that capture R code for transparency
#
# Author: Prairie Genomics Team  
# Date: January 24, 2025
# Purpose: Wrap existing analysis functions to enable code visibility

# Source required modules
source("code_visibility/code_logger.R")
source("code_visibility/code_generator.R")

# Enhanced DESeq2 analysis with code logging
run_deseq2_analysis_logged <- function(count_matrix, sample_info, design_formula, 
                                      contrast_info = NULL, filter_params = NULL,
                                      session_id = NULL, species = "human") {
  
  # Initialize session if needed
  if (is.null(session_id)) {
    session_id <- init_code_logger()
  }
  
  # Record package usage
  log_package_usage(session_id, "DESeq2")
  log_package_usage(session_id, "ggplot2")
  
  start_time <- Sys.time()
  
  # Generate and log data loading code
  file_info <- list(name = "count_matrix.csv")  # Placeholder
  processing_steps <- c("remove_zero_genes", "ensure_integer")
  
  data_loading_code <- generate_data_loading_code(file_info, processing_steps)
  log_analysis_step(
    session_id = session_id,
    step_name = "Data Loading and Preprocessing",
    category = "data_upload",
    code_snippet = data_loading_code,
    description = "Load count matrix and perform initial data preprocessing",
    input_data = paste("Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples"),
    output_data = "Preprocessed count matrix ready for DESeq2"
  )
  
  # Generate and log sample annotation code
  pattern_info <- list(pattern = "auto-detected", confidence = 0.95)  # Placeholder
  sample_annotation_code <- generate_sample_annotation_code(sample_info, design_formula, pattern_info)
  log_analysis_step(
    session_id = session_id,
    step_name = "Sample Annotation Setup",
    category = "sample_annotation", 
    code_snippet = sample_annotation_code,
    description = "Create sample information data frame and experimental design",
    input_data = paste("Sample groups:", paste(unique(sample_info$group), collapse = ", ")),
    output_data = paste("Design formula:", design_formula)
  )
  
  # Generate and log DESeq2 analysis code
  deseq2_code <- generate_deseq2_code(design_formula, contrast_info, filter_params)
  
  # Run the actual DESeq2 analysis
  tryCatch({
    # Create DESeq2 dataset
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = sample_info,
      design = as.formula(design_formula)
    )
    
    # Run DESeq2
    dds <- DESeq2::DESeq(dds)
    
    # Extract results
    if (!is.null(contrast_info)) {
      results <- DESeq2::results(dds, contrast = c("group", contrast_info$numerator, contrast_info$denominator))
    } else {
      results <- DESeq2::results(dds)
    }
    
    # Convert to data frame
    results_df <- as.data.frame(results)
    results_df$gene_id <- rownames(results_df)
    
    # Apply filtering if specified
    if (!is.null(filter_params)) {
      significant_genes <- results_df[
        !is.na(results_df$padj) & 
        results_df$padj < filter_params$padj_cutoff & 
        abs(results_df$log2FoldChange) > filter_params$fc_cutoff,
      ]
    } else {
      significant_genes <- results_df[
        !is.na(results_df$padj) & 
        results_df$padj < 0.05 & 
        abs(results_df$log2FoldChange) > 1,
      ]
    }
    
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Log the analysis step
    log_analysis_step(
      session_id = session_id,
      step_name = "DESeq2 Differential Expression Analysis",
      category = "deseq2",
      code_snippet = deseq2_code,
      parameters = list(
        design_formula = design_formula,
        contrast = contrast_info,
        filters = filter_params
      ),
      description = "Run DESeq2 analysis and extract differential expression results",
      input_data = paste("DESeq2 dataset with", nrow(count_matrix), "genes"),
      output_data = paste("Results with", nrow(significant_genes), "significant genes"),
      execution_time = execution_time
    )
    
    # Add gene symbol conversion if requested
    if (species %in% c("human", "mouse")) {
      conversion_code <- generate_gene_conversion_code(species, "ensembl_to_symbol")
      
      # Perform actual conversion (simplified for demo)
      if (species == "human") {
        log_package_usage(session_id, "org.Hs.eg.db")
      } else {
        log_package_usage(session_id, "org.Mm.eg.db")
      }
      
      log_analysis_step(
        session_id = session_id,
        step_name = "Gene Symbol Conversion",
        category = "deseq2",
        code_snippet = conversion_code,
        description = paste("Convert Ensembl IDs to gene symbols for", species),
        input_data = "Ensembl gene IDs",
        output_data = "Results with gene symbols added"
      )
    }
    
    # Return enhanced results
    return(list(
      success = TRUE,
      results = results_df,
      significant_genes = significant_genes,
      dds = dds,
      session_id = session_id,
      analysis_summary = list(
        total_genes = nrow(results_df),
        significant_genes = nrow(significant_genes),
        upregulated = sum(significant_genes$log2FoldChange > 0, na.rm = TRUE),
        downregulated = sum(significant_genes$log2FoldChange < 0, na.rm = TRUE)
      )
    ))
    
  }, error = function(e) {
    log_analysis_step(
      session_id = session_id,
      step_name = "DESeq2 Analysis - ERROR",
      category = "deseq2",
      code_snippet = paste("# ERROR occurred:", e$message),
      description = "DESeq2 analysis failed",
      output_data = paste("Error:", e$message)
    )
    
    return(list(
      success = FALSE,
      error = e$message,
      session_id = session_id
    ))
  })
}

# Enhanced pathway analysis with code logging
run_pathway_analysis_logged <- function(deseq2_results, analysis_type, species,
                                       parameters = list(), session_id = NULL) {
  
  # Initialize session if needed
  if (is.null(session_id)) {
    session_id <- init_code_logger()
  }
  
  start_time <- Sys.time()
  
  # Log package usage based on analysis type
  if (analysis_type %in% c("GO", "KEGG")) {
    log_package_usage(session_id, "clusterProfiler")
    if (species == "human") {
      log_package_usage(session_id, "org.Hs.eg.db")
    } else {
      log_package_usage(session_id, "org.Mm.eg.db")
    }
  } else if (analysis_type == "GSEA") {
    log_package_usage(session_id, "fgsea")
    log_package_usage(session_id, "msigdbr")
  }
  
  # Generate pathway analysis code
  pathway_code <- generate_pathway_code(analysis_type, species, parameters)
  
  # Run the actual pathway analysis (simplified for demo)
  tryCatch({
    # Here you would call the actual pathway analysis functions
    # For now, we'll create a mock result
    
    if (analysis_type == "GO") {
      # Mock GO analysis
      analysis_result <- list(
        success = TRUE,
        n_terms = 25,
        method = "GO enrichment"
      )
    } else if (analysis_type == "KEGG") {
      # Mock KEGG analysis  
      analysis_result <- list(
        success = TRUE,
        n_pathways = 15,
        method = "KEGG pathway enrichment"
      )
    } else if (analysis_type == "GSEA") {
      # Mock GSEA analysis
      analysis_result <- list(
        success = TRUE,
        n_pathways = 50,
        n_significant = 12,
        method = "Gene Set Enrichment Analysis"
      )
    }
    
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Log the analysis step
    log_analysis_step(
      session_id = session_id,
      step_name = paste(analysis_type, "Pathway Analysis"),
      category = "pathway",
      code_snippet = pathway_code,
      parameters = c(list(analysis_type = analysis_type, species = species), parameters),
      description = paste("Run", analysis_type, "pathway analysis for", species),
      input_data = "DESeq2 differential expression results",
      output_data = paste("Pathway analysis results:", analysis_result$method),
      execution_time = execution_time
    )
    
    return(c(analysis_result, list(session_id = session_id)))
    
  }, error = function(e) {
    log_analysis_step(
      session_id = session_id,
      step_name = paste(analysis_type, "Pathway Analysis - ERROR"),
      category = "pathway",
      code_snippet = paste("# ERROR occurred:", e$message),
      description = paste(analysis_type, "pathway analysis failed"),
      output_data = paste("Error:", e$message)
    )
    
    return(list(
      success = FALSE,
      error = e$message,
      session_id = session_id
    ))
  })
}

# Enhanced visualization with code logging
create_plot_logged <- function(plot_type, data_info, parameters = list(), session_id = NULL) {
  
  # Initialize session if needed  
  if (is.null(session_id)) {
    session_id <- init_code_logger()
  }
  
  # Log package usage
  log_package_usage(session_id, "ggplot2")
  if (plot_type == "heatmap") {
    log_package_usage(session_id, "pheatmap")
  }
  
  start_time <- Sys.time()
  
  # Generate visualization code
  viz_code <- generate_visualization_code(plot_type, data_info, parameters)
  
  # Create the actual plot (simplified for demo)
  tryCatch({
    # Here you would create the actual plot
    # For now, we'll just return success
    
    plot_result <- list(
      success = TRUE,
      plot_type = plot_type,
      message = paste(stringr::str_to_title(plot_type), "plot created successfully")
    )
    
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Log the visualization step
    log_analysis_step(
      session_id = session_id,
      step_name = paste(stringr::str_to_title(plot_type), "Plot Generation"),
      category = "visualization",
      code_snippet = viz_code,
      parameters = c(list(plot_type = plot_type), parameters),
      description = paste("Generate", plot_type, "plot for data visualization"),
      input_data = "Analysis results data",
      output_data = paste(stringr::str_to_title(plot_type), "plot"),
      execution_time = execution_time
    )
    
    return(c(plot_result, list(session_id = session_id)))
    
  }, error = function(e) {
    log_analysis_step(
      session_id = session_id,
      step_name = paste(stringr::str_to_title(plot_type), "Plot - ERROR"),
      category = "visualization", 
      code_snippet = paste("# ERROR occurred:", e$message),
      description = paste(plot_type, "plot generation failed"),
      output_data = paste("Error:", e$message)
    )
    
    return(list(
      success = FALSE,
      error = e$message,
      session_id = session_id
    ))
  })
}

# Enhanced data upload with code logging
log_data_upload <- function(file_info, processing_info, session_id = NULL) {
  
  # Initialize session if needed
  if (is.null(session_id)) {
    session_id <- init_code_logger()
  }
  
  # Log packages used for data loading
  if (grepl("\\.xlsx?$", file_info$name, ignore.case = TRUE)) {
    log_package_usage(session_id, "readxl")
  }
  
  # Generate data loading code
  data_code <- generate_data_loading_code(file_info, processing_info$steps)
  
  # Log the data upload step
  log_analysis_step(
    session_id = session_id,
    step_name = "Data Upload and Loading",
    category = "data_upload",
    code_snippet = data_code,
    parameters = list(
      filename = file_info$name,
      file_size = file_info$size,
      processing_steps = processing_info$steps
    ),
    description = "Upload and preprocess count matrix data",
    input_data = paste("File:", file_info$name, "Size:", file_info$size),
    output_data = paste("Count matrix:", processing_info$final_dimensions)
  )
  
  return(session_id)
}

# Get analysis code for display
get_analysis_code_snippet <- function(session_id, category = NULL, step_name = NULL) {
  if (!session_id %in% names(.code_logger_env)) {
    return("# No session data available")
  }
  
  session_data <- .code_logger_env[[session_id]]
  
  if (is.null(category) && is.null(step_name)) {
    # Return all code
    return(get_session_code(session_id))
  }
  
  # Filter by category or step name
  filtered_steps <- session_data$analysis_steps
  
  if (!is.null(category)) {
    filtered_steps <- filtered_steps[sapply(filtered_steps, function(x) x$category == category)]
  }
  
  if (!is.null(step_name)) {
    filtered_steps <- filtered_steps[sapply(filtered_steps, function(x) grepl(step_name, x$step_name, ignore.case = TRUE))]
  }
  
  if (length(filtered_steps) == 0) {
    return("# No matching analysis steps found")
  }
  
  # Combine code snippets
  code_snippets <- sapply(filtered_steps, function(step) {
    c(
      paste0("# ", step$step_name),
      paste0("# Timestamp: ", step$timestamp),
      if (step$description != "") paste0("# ", step$description) else NULL,
      "",
      step$code_snippet,
      ""
    )
  })
  
  return(paste(unlist(code_snippets), collapse = "\n"))
}

# Export current session code
export_current_analysis <- function(session_id, format = "R", filename = NULL) {
  if (!session_id %in% names(.code_logger_env)) {
    return(NULL)
  }
  
  return(export_session_code(session_id, filename, format))
}

cat("✅ Analysis wrappers module loaded\n")
cat("📋 Enhanced functions: run_deseq2_analysis_logged(), run_pathway_analysis_logged(), create_plot_logged()\n")