# Code Logger Module for Prairie Genomics Suite
# Captures R code for each analysis step to enhance scientific rigor and reproducibility
#
# Author: Prairie Genomics Team
# Date: January 24, 2025
# Purpose: Enable users to view, understand, and reproduce analysis steps

# Global code logging environment
.code_logger_env <- new.env()

# Initialize code logging system
init_code_logger <- function(session_id = NULL, user_name = "user") {
  if (is.null(session_id)) {
    session_id <- paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Create session data structure
  session_data <- list(
    session_id = session_id,
    user_name = user_name,
    start_time = Sys.time(),
    analysis_steps = list(),
    packages_used = character(),
    R_version = R.version.string,
    package_versions = list(),
    current_step = 0
  )
  
  # Store in global environment
  .code_logger_env[[session_id]] <- session_data
  
  cat("🔧 Code logger initialized for session:", session_id, "\n")
  return(session_id)
}

# Log an analysis step with code and metadata
log_analysis_step <- function(session_id, step_name, category, 
                             code_snippet, parameters = list(), 
                             description = "", input_data = "", 
                             output_data = "", execution_time = NULL) {
  
  if (!session_id %in% names(.code_logger_env)) {
    warning("Session not found. Initializing new session.")
    session_id <- init_code_logger(session_id)
  }
  
  session_data <- .code_logger_env[[session_id]]
  session_data$current_step <- session_data$current_step + 1
  
  # Create step record
  step_record <- list(
    step_number = session_data$current_step,
    step_name = step_name,
    category = category,  # "data_upload", "sample_annotation", "deseq2", "pathway", "visualization"
    timestamp = Sys.time(),
    code_snippet = code_snippet,
    parameters = parameters,
    description = description,
    input_data = input_data,
    output_data = output_data,
    execution_time = execution_time
  )
  
  # Add to session steps
  session_data$analysis_steps[[length(session_data$analysis_steps) + 1]] <- step_record
  
  # Update session data
  .code_logger_env[[session_id]] <- session_data
  
  cat("📝 Logged step", session_data$current_step, ":", step_name, "\n")
  return(step_record)
}

# Record package usage
log_package_usage <- function(session_id, package_name, version = NULL) {
  if (!session_id %in% names(.code_logger_env)) {
    return(NULL)
  }
  
  session_data <- .code_logger_env[[session_id]]
  
  # Add package to list if not already present
  if (!package_name %in% session_data$packages_used) {
    session_data$packages_used <- c(session_data$packages_used, package_name)
  }
  
  # Record version if provided or can be determined
  if (is.null(version)) {
    version <- tryCatch({
      as.character(packageVersion(package_name))
    }, error = function(e) "unknown")
  }
  
  session_data$package_versions[[package_name]] <- version
  .code_logger_env[[session_id]] <- session_data
}

# Generate data loading code snippet
generate_data_loading_code <- function(file_info, processing_steps = list()) {
  code_lines <- c(
    "# Data Loading and Preprocessing",
    "# ==============================",
    ""
  )
  
  # File loading based on type
  if (grepl("\\.csv$", file_info$name, ignore.case = TRUE)) {
    code_lines <- c(code_lines, 
      paste0("# Load count matrix from CSV file"),
      paste0("count_matrix <- read.csv('", file_info$name, "', row.names = 1, check.names = FALSE)"),
      ""
    )
  } else if (grepl("\\.xlsx?$", file_info$name, ignore.case = TRUE)) {
    code_lines <- c(code_lines,
      "library(readxl)",
      paste0("# Load count matrix from Excel file"),
      paste0("count_matrix <- as.data.frame(read_excel('", file_info$name, "', sheet = 1))"),
      "rownames(count_matrix) <- count_matrix[,1]",
      "count_matrix <- count_matrix[,-1]",
      ""
    )
  } else {
    code_lines <- c(code_lines,
      paste0("# Load count matrix from ", file_info$name),
      paste0("count_matrix <- read.delim('", file_info$name, "', row.names = 1, check.names = FALSE)"),
      ""
    )
  }
  
  # Add preprocessing steps
  if (length(processing_steps) > 0) {
    code_lines <- c(code_lines, "# Data preprocessing:")
    
    if ("remove_zero_genes" %in% processing_steps) {
      code_lines <- c(code_lines,
        "# Remove genes with zero counts across all samples",
        "count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]"
      )
    }
    
    if ("ensure_integer" %in% processing_steps) {
      code_lines <- c(code_lines,
        "# Ensure count data is integer (required for DESeq2)",
        "count_matrix <- round(count_matrix)"
      )
    }
    
    code_lines <- c(code_lines, "")
  }
  
  # Add data summary
  code_lines <- c(code_lines,
    "# Data summary",
    "cat('Loaded count matrix with', nrow(count_matrix), 'genes and', ncol(count_matrix), 'samples\\n')",
    "head(count_matrix[,1:min(5, ncol(count_matrix))])",
    ""
  )
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate sample annotation code
generate_sample_annotation_code <- function(sample_info, design_formula, pattern_info = NULL) {
  code_lines <- c(
    "# Sample Annotation and Design Setup",
    "# ==================================",
    ""
  )
  
  # If pattern was detected automatically
  if (!is.null(pattern_info)) {
    code_lines <- c(code_lines,
      "# Extract sample groups from sample names using pattern matching",
      paste0("# Pattern detected: ", pattern_info$pattern),
      paste0("# Confidence: ", round(pattern_info$confidence * 100, 1), "%"),
      ""
    )
  }
  
  # Create sample info data frame
  sample_groups <- unique(sample_info$group)
  code_lines <- c(code_lines,
    "# Create sample information data frame",
    "sample_info <- data.frame(",
    "  sample_id = c(",
    paste0("    ", paste(paste0("'", sample_info$sample_id, "'"), collapse = ",\n    ")),
    "  ),",
    "  group = c(",
    paste0("    ", paste(paste0("'", sample_info$group, "'"), collapse = ",\n    ")),
    "  ),",
    "  stringsAsFactors = FALSE",
    ")",
    "rownames(sample_info) <- sample_info$sample_id",
    ""
  )
  
  # Design formula
  code_lines <- c(code_lines,
    "# Set up experimental design for DESeq2",
    paste0("design_formula <- ~ ", design_formula),
    "",
    "# Ensure sample order matches between count matrix and sample info",
    "sample_info <- sample_info[colnames(count_matrix), ]",
    ""
  )
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate DESeq2 analysis code
generate_deseq2_code <- function(design_formula, contrast_info, filter_params) {
  code_lines <- c(
    "# DESeq2 Differential Expression Analysis",
    "# ======================================",
    "",
    "library(DESeq2)",
    ""
  )
  
  # Create DESeq2 dataset
  code_lines <- c(code_lines,
    "# Create DESeq2 dataset object",
    paste0("dds <- DESeqDataSetFromMatrix("),
    "  countData = count_matrix,",
    "  colData = sample_info,",
    paste0("  design = ", design_formula),
    ")",
    ""
  )
  
  # Run DESeq2 analysis
  code_lines <- c(code_lines,
    "# Run DESeq2 analysis with default parameters (recommended)",
    "dds <- DESeq(dds)",
    ""
  )
  
  # Extract results
  if (!is.null(contrast_info)) {
    code_lines <- c(code_lines,
      "# Extract results for specified contrast",
      paste0("results <- results(dds, contrast = c('group', '", contrast_info$numerator, "', '", contrast_info$denominator, "'))"),
      ""
    )
  } else {
    code_lines <- c(code_lines,
      "# Extract results (using last coefficient)",
      "results <- results(dds)",
      ""
    )
  }
  
  # Convert to data frame and add gene symbols
  code_lines <- c(code_lines,
    "# Convert to data frame for easier manipulation",
    "results_df <- as.data.frame(results)",
    "results_df$gene_id <- rownames(results_df)",
    ""
  )
  
  # Add filtering
  if (!is.null(filter_params)) {
    code_lines <- c(code_lines,
      "# Filter for significant genes",
      paste0("significant_genes <- results_df["),
      paste0("  !is.na(results_df$padj) & "),
      paste0("  results_df$padj < ", filter_params$padj_cutoff, " & "),
      paste0("  abs(results_df$log2FoldChange) > ", filter_params$fc_cutoff),
      paste0(", ]"),
      "",
      "cat('Found', nrow(significant_genes), 'significant genes\\n')",
      ""
    )
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Generate pathway analysis code
generate_pathway_code <- function(analysis_type, species, parameters) {
  code_lines <- c(
    paste0("# ", analysis_type, " Pathway Analysis"),
    paste0("# ", paste(rep("=", nchar(analysis_type) + 17), collapse = "")),
    ""
  )
  
  # Load required packages
  if (analysis_type %in% c("GO", "KEGG")) {
    code_lines <- c(code_lines,
      "library(clusterProfiler)",
      "library(org.Hs.eg.db)  # or org.Mm.eg.db for mouse",
      ""
    )
  } else if (analysis_type == "GSEA") {
    code_lines <- c(code_lines,
      "library(fgsea)",
      "library(msigdbr)",
      ""
    )
  }
  
  # Gene preparation
  if (analysis_type %in% c("GO", "KEGG")) {
    code_lines <- c(code_lines,
      "# Prepare gene list for over-representation analysis",
      "significant_genes <- results_df[",
      paste0("  !is.na(results_df$padj) & results_df$padj < ", parameters$padj_cutoff, " &"),
      paste0("  abs(results_df$log2FoldChange) > ", parameters$fc_cutoff),
      ", ]",
      "",
      "# Convert Ensembl IDs to Entrez IDs",
      "gene_list <- mapIds(org.Hs.eg.db, keys = significant_genes$gene_id,",
      "                    column = 'ENTREZID', keytype = 'ENSEMBL')",
      "gene_list <- gene_list[!is.na(gene_list)]",
      ""
    )
  } else if (analysis_type == "GSEA") {
    code_lines <- c(code_lines,
      "# Prepare ranked gene list for GSEA",
      "# Remove genes with NA values",
      "clean_results <- results_df[",
      "  !is.na(results_df$log2FoldChange) & !is.na(results_df$pvalue), ]",
      "",
      "# Create ranking statistic (signed p-value method)",
      "gene_ranks <- sign(clean_results$log2FoldChange) * (-log10(clean_results$pvalue))",
      "names(gene_ranks) <- clean_results$gene_id",
      "",
      "# Sort in decreasing order",
      "gene_ranks <- sort(gene_ranks, decreasing = TRUE)",
      ""
    )
  }
  
  # Analysis-specific code
  if (analysis_type == "GO") {
    code_lines <- c(code_lines,
      "# Run GO enrichment analysis",
      "go_results <- enrichGO(",
      "  gene = gene_list,",
      "  OrgDb = org.Hs.eg.db,",
      paste0("  ont = '", parameters$ontology, "',"),
      "  pAdjustMethod = 'BH',",
      "  pvalueCutoff = 0.05,",
      "  qvalueCutoff = 0.2,",
      "  readable = TRUE",
      ")",
      ""
    )
  } else if (analysis_type == "KEGG") {
    code_lines <- c(code_lines,
      "# Run KEGG pathway analysis",
      "kegg_results <- enrichKEGG(",
      "  gene = gene_list,",
      paste0("  organism = '", if(species == "human") "hsa" else "mmu", "',"),
      "  pvalueCutoff = 0.05,",
      "  pAdjustMethod = 'BH'",
      ")",
      ""
    )
  } else if (analysis_type == "GSEA") {
    code_lines <- c(code_lines,
      "# Get gene sets from MSigDB",
      paste0("gene_sets <- msigdbr(species = '", if(species == "human") "Homo sapiens" else "Mus musculus", "',"),
      paste0("                     collection = '", parameters$collection, "')"),
      "",
      "# Convert to list format for fgsea",
      "pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)",
      "",
      "# Run GSEA analysis",
      "gsea_results <- fgsea(",
      "  pathways = pathways,",
      "  stats = gene_ranks,",
      "  minSize = 15,",
      "  maxSize = 500",
      ")",
      ""
    )
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Get complete analysis script for a session
get_session_code <- function(session_id, include_comments = TRUE, include_sessioninfo = TRUE) {
  if (!session_id %in% names(.code_logger_env)) {
    return("# No session data found")
  }
  
  session_data <- .code_logger_env[[session_id]]
  
  # Header
  code_lines <- c(
    "# Prairie Genomics Suite - Complete Analysis Script",
    "# =================================================",
    paste0("# Generated: ", Sys.time()),
    paste0("# Session ID: ", session_id),
    paste0("# User: ", session_data$user_name),
    paste0("# R Version: ", session_data$R_version),
    "",
    "# This script reproduces the complete analysis performed in the",
    "# Prairie Genomics Suite web application.",
    ""
  )
  
  # Package loading
  if (length(session_data$packages_used) > 0) {
    code_lines <- c(code_lines,
      "# Load required packages",
      "# =====================",
      ""
    )
    
    for (pkg in session_data$packages_used) {
      version <- session_data$package_versions[[pkg]]
      version_comment <- if (!is.null(version) && version != "unknown") {
        paste0("  # Version: ", version)
      } else {
        ""
      }
      
      code_lines <- c(code_lines,
        paste0("library(", pkg, ")", version_comment)
      )
    }
    code_lines <- c(code_lines, "")
  }
  
  # Add each analysis step
  for (i in seq_along(session_data$analysis_steps)) {
    step <- session_data$analysis_steps[[i]]
    
    if (include_comments) {
      code_lines <- c(code_lines,
        paste0("# Step ", step$step_number, ": ", step$step_name),
        paste0("# Category: ", step$category),
        paste0("# Timestamp: ", step$timestamp),
        if (step$description != "") paste0("# Description: ", step$description) else NULL,
        ""
      )
    }
    
    # Add the actual code
    code_lines <- c(code_lines, step$code_snippet, "")
  }
  
  # Session info
  if (include_sessioninfo) {
    code_lines <- c(code_lines,
      "# Session Information",
      "# ==================",
      "sessionInfo()",
      ""
    )
  }
  
  return(paste(code_lines, collapse = "\n"))
}

# Export session code to file
export_session_code <- function(session_id, filename = NULL, format = "R") {
  if (is.null(filename)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("prairie_genomics_analysis_", timestamp, ".", tolower(format))
  }
  
  code_content <- get_session_code(session_id)
  
  if (format == "R") {
    writeLines(code_content, filename)
  } else if (format == "Rmd") {
    # Convert to R Markdown format
    rmd_content <- convert_to_rmarkdown(code_content, session_id)
    writeLines(rmd_content, filename)
  }
  
  cat("📄 Analysis code exported to:", filename, "\n")
  return(filename)
}

# Clear session data
clear_session <- function(session_id) {
  if (session_id %in% names(.code_logger_env)) {
    rm(list = session_id, envir = .code_logger_env)
    cat("🗑️ Session", session_id, "cleared\n")
  }
}

# Get all active sessions
get_active_sessions <- function() {
  return(names(.code_logger_env))
}

# Get session summary
get_session_summary <- function(session_id) {
  if (!session_id %in% names(.code_logger_env)) {
    return(NULL)
  }
  
  session_data <- .code_logger_env[[session_id]]
  
  summary <- list(
    session_id = session_data$session_id,
    user_name = session_data$user_name,
    start_time = session_data$start_time,
    num_steps = length(session_data$analysis_steps),
    packages_used = session_data$packages_used,
    categories = unique(sapply(session_data$analysis_steps, function(x) x$category)),
    last_activity = if (length(session_data$analysis_steps) > 0) {
      max(sapply(session_data$analysis_steps, function(x) x$timestamp))
    } else {
      session_data$start_time
    }
  )
  
  return(summary)
}

cat("✅ Code logger module loaded\n")
cat("📋 Functions available: init_code_logger(), log_analysis_step(), get_session_code(), export_session_code()\n")