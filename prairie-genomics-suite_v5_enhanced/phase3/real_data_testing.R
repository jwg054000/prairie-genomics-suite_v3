# Phase 3 - Real Data Testing & Validation System
# Upload and validate real RNA-seq datasets for ground truth testing
# 
# This module provides:
# - Secure data upload interface for expression matrices
# - Data format validation and standardization
# - Ground truth collection and comparison
# - System validation against expert assessments
# - Template creation from proven workflows

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(readr)
library(readxl)

# ===========================================
# DATA UPLOAD CONFIGURATION
# ===========================================

# Define supported file formats and size limits
UPLOAD_CONFIG <- list(
  max_file_size_mb = 100,
  supported_formats = list(
    expression_matrix = c(".csv", ".tsv", ".txt", ".xlsx", ".xls", ".RData", ".rds"),
    sample_metadata = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
    additional_files = c(".txt", ".md", ".pdf")
  ),
  required_columns = list(
    expression_matrix = NULL,  # Will be validated dynamically
    sample_metadata = c("Sample_ID")  # Minimum required column
  )
)

# ===========================================
# DATA VALIDATION RULES
# ===========================================

# Define comprehensive validation rules for real data
DATA_VALIDATION_RULES <- list(
  
  expression_matrix = list(
    
    format_detection = list(
      name = "Matrix Format Detection",
      description = "Automatically detect if genes are rows or columns",
      validation_function = "detect_matrix_orientation",
      importance = "critical",
      auto_fix = TRUE,
      explanation = "RNA-seq data can be formatted as genes×samples or samples×genes. We'll automatically detect and standardize the format."
    ),
    
    gene_id_validation = list(
      name = "Gene ID Format Validation", 
      description = "Validate and standardize gene identifiers",
      validation_function = "validate_gene_ids",
      importance = "high",
      auto_fix = TRUE,
      explanation = "We support Ensembl IDs, Gene Symbols, and RefSeq IDs. Mixed formats will be flagged for review."
    ),
    
    numeric_validation = list(
      name = "Numeric Data Validation",
      description = "Ensure all expression values are numeric",
      validation_function = "validate_numeric_data",
      importance = "critical",
      auto_fix = FALSE,
      explanation = "Expression data must be numeric. Non-numeric values suggest formatting issues."
    ),
    
    sample_alignment = list(
      name = "Sample-Metadata Alignment",
      description = "Verify sample names match between expression and metadata",
      validation_function = "validate_sample_alignment",
      importance = "critical", 
      auto_fix = FALSE,
      explanation = "Sample names in the expression matrix must exactly match those in the metadata file."
    ),
    
    data_distribution = list(
      name = "Expression Distribution Analysis",
      description = "Analyze expression data distribution patterns",
      validation_function = "analyze_expression_distribution",
      importance = "medium",
      auto_fix = FALSE,
      explanation = "RNA-seq count data should follow expected statistical distributions. Unusual patterns may indicate processing issues."
    )
  ),
  
  sample_metadata = list(
    
    required_columns = list(
      name = "Required Column Validation",
      description = "Check for essential metadata columns",
      validation_function = "validate_required_columns",
      importance = "critical",
      auto_fix = FALSE,
      explanation = "Sample_ID is required. Condition column is highly recommended for differential analysis."
    ),
    
    data_types = list(
      name = "Metadata Data Types",
      description = "Validate appropriate data types for each column",
      validation_function = "validate_metadata_types",
      importance = "medium",
      auto_fix = TRUE,
      explanation = "Categorical variables should be factors, numeric variables should be numbers."
    ),
    
    experimental_design = list(
      name = "Experimental Design Assessment",
      description = "Analyze experimental design for common issues",
      validation_function = "assess_experimental_design",
      importance = "high",
      auto_fix = FALSE,
      explanation = "We'll check for balanced designs, adequate replication, and potential confounding."
    )
  )
)

# ===========================================
# GROUND TRUTH COLLECTION FRAMEWORK
# ===========================================

# Define structure for collecting expert knowledge
GROUND_TRUTH_SCHEMA <- list(
  
  experiment_assessment = list(
    experiment_type = list(
      question = "What type of experiment is this?",
      type = "select",
      options = c("cell_line_treatment", "tissue_comparison", "time_course", "patient_samples", "dose_response", "other"),
      required = TRUE,
      explanation = "This helps us validate our auto-detection algorithm against your expert assessment."
    ),
    
    confidence_level = list(
      question = "How confident are you in this experiment classification?",
      type = "slider",
      range = c(1, 10),
      default = 8,
      required = TRUE,
      explanation = "This helps us understand when our auto-detection should defer to expert judgment."
    ),
    
    design_challenges = list(
      question = "What challenges did you encounter with this experimental design?",
      type = "textarea",
      required = FALSE,
      explanation = "This helps us improve our error prevention system by learning from real issues."
    )
  ),
  
  analysis_parameters = list(
    
    significance_threshold = list(
      question = "What adjusted p-value threshold did you use?",
      type = "numeric",
      default = 0.05,
      min = 0.001,
      max = 0.1,
      required = TRUE,
      explanation = "We'll compare this to our intelligent parameter recommendations."
    ),
    
    fold_change_threshold = list(
      question = "What fold change threshold did you use?",
      type = "numeric", 
      default = 2.0,
      min = 1.1,
      max = 10.0,
      required = TRUE,
      explanation = "This helps us validate our effect size recommendations."
    ),
    
    parameter_rationale = list(
      question = "Why did you choose these specific parameters?",
      type = "textarea",
      required = FALSE,
      explanation = "Your reasoning helps us improve our educational explanations."
    )
  ),
  
  biological_context = list(
    
    positive_controls = list(
      question = "List genes that should respond to your treatment (positive controls)",
      type = "textarea",
      placeholder = "e.g., TP53, MYC, CDKN1A (one per line or comma-separated)",
      required = FALSE,
      explanation = "These help validate our positive control detection system."
    ),
    
    negative_controls = list(
      question = "List genes that should NOT change (negative controls)",
      type = "textarea", 
      placeholder = "e.g., ACTB, GAPDH, RPL13A (one per line or comma-separated)",
      required = FALSE,
      explanation = "These help validate our housekeeping gene stability checks."
    ),
    
    expected_pathways = list(
      question = "What biological pathways do you expect to be affected?",
      type = "textarea",
      placeholder = "e.g., Cell cycle, Apoptosis, DNA repair",
      required = FALSE,
      explanation = "This helps us test our pathway coherence analysis."
    ),
    
    known_issues = list(
      question = "What issues did you encounter and how did you solve them?",
      type = "textarea",
      required = FALSE,
      explanation = "This helps us improve our troubleshooting guidance and error prevention."
    )
  ),
  
  results_validation = list(
    
    final_gene_count = list(
      question = "How many genes were significantly differentially expressed in your final analysis?",
      type = "numeric",
      min = 0,
      required = FALSE,
      explanation = "We'll compare this to our analysis results to validate accuracy."
    ),
    
    key_findings = list(
      question = "What were your key biological findings from this experiment?",
      type = "textarea",
      required = FALSE,
      explanation = "This helps us validate that our system reaches similar conclusions."
    ),
    
    publication_status = list(
      question = "What is the publication status of this data?",
      type = "select",
      options = c("published", "submitted", "in_preparation", "internal_use", "preliminary"),
      required = FALSE,
      explanation = "This helps us understand the quality and maturity of the analysis."
    )
  )
)

# ===========================================
# DATA UPLOAD INTERFACE
# ===========================================

# Create real data upload UI
create_real_data_upload_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Upload header
    div(
      class = "upload-header",
      style = "background: linear-gradient(135deg, #3b82f6 0%, #1d4ed8 100%); 
               color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
      
      div(
        h2("📊 Real Data Integration & Testing", style = "margin: 0;"),
        p("Upload your RNA-seq data to validate our scientific guardrails", style = "margin: 5px 0 0 0; opacity: 0.9;")
      )
    ),
    
    # Upload workflow
    tabsetPanel(
      id = ns("upload_tabs"),
      type = "pills",
      
      # Data Upload Tab
      tabPanel(
        "📁 Data Upload",
        value = "upload",
        
        div(style = "padding: 20px 0;",
          
          # Instructions
          div(
            class = "alert alert-info",
            style = "margin-bottom: 20px;",
            icon("info-circle"),
            strong("Upload Instructions: "),
            "Start with your most confident/successful experiment. We'll validate our system against your expert knowledge."
          ),
          
          fluidRow(
            # Expression Matrix Upload
            column(6,
              div(
                class = "upload-card",
                style = "background: white; padding: 20px; border-radius: 12px; 
                         border: 2px dashed #d1d5db; margin-bottom: 20px;",
                
                h4("📊 Expression Matrix", style = "margin-top: 0;"),
                p("Upload your gene expression count matrix", style = "color: #64748b;"),
                
                fileInput(
                  ns("expression_file"),
                  "Choose Expression File",
                  accept = paste(UPLOAD_CONFIG$supported_formats$expression_matrix, collapse = ","),
                  width = "100%"
                ),
                
                div(
                  style = "font-size: 12px; color: #64748b; margin-top: 10px;",
                  "Supported formats: CSV, TSV, Excel, RData. Max size: 100MB"
                )
              )
            ),
            
            # Sample Metadata Upload
            column(6,
              div(
                class = "upload-card",
                style = "background: white; padding: 20px; border-radius: 12px; 
                         border: 2px dashed #d1d5db; margin-bottom: 20px;",
                
                h4("📋 Sample Metadata", style = "margin-top: 0;"),
                p("Upload sample information and experimental design", style = "color: #64748b;"),
                
                fileInput(
                  ns("metadata_file"),
                  "Choose Metadata File",
                  accept = paste(UPLOAD_CONFIG$supported_formats$sample_metadata, collapse = ","),
                  width = "100%"
                ),
                
                div(
                  style = "font-size: 12px; color: #64748b; margin-top: 10px;",
                  "Required: Sample_ID column. Recommended: Condition, Batch columns"
                )
              )
            )
          ),
          
          # Upload Status
          div(
            id = ns("upload_status"),
            style = "margin-top: 20px;",
            uiOutput(ns("upload_status_display"))
          )
        )
      ),
      
      # Data Preview Tab
      tabPanel(
        "👁️ Data Preview",
        value = "preview",
        
        div(style = "padding: 20px 0;",
          uiOutput(ns("data_preview_display"))
        )
      ),
      
      # Ground Truth Tab
      tabPanel(
        "🎯 Ground Truth Collection",
        value = "ground_truth",
        
        div(style = "padding: 20px 0;",
          h4("Share Your Expert Knowledge", style = "margin-bottom: 20px;"),
          p("Help us validate our system by sharing your analysis decisions and biological insights.", 
            style = "color: #64748b; margin-bottom: 30px;"),
          
          uiOutput(ns("ground_truth_form"))
        )
      ),
      
      # Validation Results Tab
      tabPanel(
        "✅ System Validation",
        value = "validation",
        
        div(style = "padding: 20px 0;",
          uiOutput(ns("validation_results_display"))
        )
      )
    )
  )
}

# ===========================================
# DATA VALIDATION ENGINE
# ===========================================

# Main data validation function
validate_uploaded_data <- function(expression_data, sample_metadata) {
  
  validation_results <- list()
  
  # Validate expression matrix
  if (!is.null(expression_data)) {
    validation_results$expression_matrix <- validate_expression_matrix(expression_data)
  }
  
  # Validate sample metadata
  if (!is.null(sample_metadata)) {
    validation_results$sample_metadata <- validate_sample_metadata(sample_metadata)
  }
  
  # Cross-validate alignment between files
  if (!is.null(expression_data) && !is.null(sample_metadata)) {
    validation_results$data_alignment <- validate_data_alignment(expression_data, sample_metadata)
  }
  
  # Overall validation summary
  validation_results$overall_status <- calculate_validation_status(validation_results)
  
  return(validation_results)
}

# Validate expression matrix
validate_expression_matrix <- function(expression_data) {
  
  results <- list()
  
  # Detect matrix orientation
  orientation_result <- detect_matrix_orientation(expression_data)
  results$orientation <- list(
    rule = DATA_VALIDATION_RULES$expression_matrix$format_detection,
    result = orientation_result,
    status = if (orientation_result$confidence > 0.8) "pass" else "warning",
    timestamp = Sys.time()
  )
  
  # Validate gene IDs
  gene_id_result <- validate_gene_ids(rownames(expression_data))
  results$gene_ids <- list(
    rule = DATA_VALIDATION_RULES$expression_matrix$gene_id_validation,
    result = gene_id_result,
    status = if (gene_id_result$format_confidence > 0.7) "pass" else "warning",
    timestamp = Sys.time()
  )
  
  # Validate numeric data
  numeric_result <- validate_numeric_data(expression_data)
  results$numeric_data <- list(
    rule = DATA_VALIDATION_RULES$expression_matrix$numeric_validation,
    result = numeric_result,
    status = if (numeric_result$all_numeric) "pass" else "fail",
    timestamp = Sys.time()
  )
  
  return(results)
}

# Detect matrix orientation (genes as rows vs columns)
detect_matrix_orientation <- function(data_matrix) {
  
  n_rows <- nrow(data_matrix)
  n_cols <- ncol(data_matrix)
  
  # Heuristics for detecting orientation
  scores <- list()
  
  # 1. Typical gene count vs sample count ratios
  if (n_rows > n_cols * 5) {
    scores$ratio_evidence <- "genes_as_rows"
  } else if (n_cols > n_rows * 5) {
    scores$ratio_evidence <- "genes_as_columns"
  } else {
    scores$ratio_evidence <- "ambiguous"
  }
  
  # 2. Gene ID pattern detection in row names
  row_gene_score <- detect_gene_id_patterns(rownames(data_matrix))
  col_gene_score <- detect_gene_id_patterns(colnames(data_matrix))
  
  if (row_gene_score > col_gene_score * 2) {
    scores$id_evidence <- "genes_as_rows"
  } else if (col_gene_score > row_gene_score * 2) {
    scores$id_evidence <- "genes_as_columns"
  } else {
    scores$id_evidence <- "ambiguous"
  }
  
  # 3. Sample name pattern detection
  row_sample_score <- detect_sample_name_patterns(rownames(data_matrix))
  col_sample_score <- detect_sample_name_patterns(colnames(data_matrix))
  
  if (col_sample_score > row_sample_score) {
    scores$sample_evidence <- "genes_as_rows"
  } else if (row_sample_score > col_sample_score) {
    scores$sample_evidence <- "genes_as_columns"
  } else {
    scores$sample_evidence <- "ambiguous"
  }
  
  # Calculate confidence
  evidence_for_rows <- sum(c(
    scores$ratio_evidence == "genes_as_rows",
    scores$id_evidence == "genes_as_rows", 
    scores$sample_evidence == "genes_as_rows"
  ))
  
  evidence_for_cols <- sum(c(
    scores$ratio_evidence == "genes_as_columns",
    scores$id_evidence == "genes_as_columns",
    scores$sample_evidence == "genes_as_columns"
  ))
  
  if (evidence_for_rows > evidence_for_cols) {
    predicted_orientation <- "genes_as_rows"
    confidence <- evidence_for_rows / 3
  } else if (evidence_for_cols > evidence_for_rows) {
    predicted_orientation <- "genes_as_columns"
    confidence <- evidence_for_cols / 3
  } else {
    predicted_orientation <- "genes_as_rows"  # Default assumption
    confidence <- 0.5
  }
  
  return(list(
    predicted_orientation = predicted_orientation,
    confidence = confidence,
    evidence = scores,
    dimensions = list(rows = n_rows, cols = n_cols),
    recommendation = if (confidence > 0.8) {
      paste("High confidence:", predicted_orientation)
    } else {
      "Please verify orientation - evidence is ambiguous"
    }
  ))
}

# Detect gene ID patterns
detect_gene_id_patterns <- function(identifiers) {
  
  if (is.null(identifiers) || length(identifiers) == 0) return(0)
  
  # Remove NA and empty strings
  clean_ids <- identifiers[!is.na(identifiers) & identifiers != ""]
  if (length(clean_ids) == 0) return(0)
  
  # Count different ID pattern matches
  patterns <- list(
    ensembl = "^ENS[A-Z]*[0-9]{11}",        # ENSEMBL IDs
    entrez = "^[0-9]+$",                     # Entrez IDs (numeric)
    symbol = "^[A-Z][A-Z0-9-]*[0-9]*$",     # Gene symbols
    refseq = "^[NX][MR]_[0-9]+",            # RefSeq IDs
    uniprot = "^[A-Z][0-9][A-Z0-9]{3}[0-9]" # UniProt IDs
  )
  
  pattern_scores <- sapply(patterns, function(pattern) {
    matches <- sum(grepl(pattern, clean_ids, ignore.case = FALSE))
    return(matches / length(clean_ids))
  })
  
  # Return the maximum score (best pattern match rate)
  return(max(pattern_scores))
}

# Detect sample name patterns  
detect_sample_name_patterns <- function(identifiers) {
  
  if (is.null(identifiers) || length(identifiers) == 0) return(0)
  
  clean_ids <- identifiers[!is.na(identifiers) & identifiers != ""]
  if (length(clean_ids) == 0) return(0)
  
  # Common sample naming patterns
  sample_patterns <- list(
    sample_prefix = "^[Ss]ample[_-]?[0-9]",
    replicate_pattern = "[Rr]ep[_-]?[0-9]",
    condition_pattern = "(control|ctrl|treat|treatment|case|patient)",
    numeric_suffix = "[_-][0-9]+$",
    batch_pattern = "[Bb]atch[_-]?[0-9]"
  )
  
  pattern_scores <- sapply(sample_patterns, function(pattern) {
    matches <- sum(grepl(pattern, clean_ids, ignore.case = TRUE))
    return(matches / length(clean_ids))
  })
  
  # Return the maximum score
  return(max(pattern_scores))
}

# Validate gene IDs
validate_gene_ids <- function(gene_ids) {
  
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    return(list(
      format_detected = "none",
      format_confidence = 0,
      issues = "No gene IDs provided"
    ))
  }
  
  # Detect the primary gene ID format
  pattern_scores <- c(
    ensembl = sum(grepl("^ENS[A-Z]*[0-9]{11}", gene_ids)) / length(gene_ids),
    symbol = sum(grepl("^[A-Z][A-Z0-9-]*[0-9]*$", gene_ids)) / length(gene_ids),
    entrez = sum(grepl("^[0-9]+$", gene_ids)) / length(gene_ids),
    refseq = sum(grepl("^[NX][MR]_[0-9]+", gene_ids)) / length(gene_ids)
  )
  
  primary_format <- names(which.max(pattern_scores))
  format_confidence <- max(pattern_scores)
  
  # Check for mixed formats
  mixed_format <- sum(pattern_scores > 0.1) > 1
  
  # Identify potential issues
  issues <- c()
  if (format_confidence < 0.5) {
    issues <- c(issues, "Low confidence in gene ID format detection")
  }
  if (mixed_format) {
    issues <- c(issues, "Mixed gene ID formats detected")
  }
  if (any(duplicated(gene_ids))) {
    issues <- c(issues, "Duplicate gene IDs found")
  }
  
  return(list(
    format_detected = primary_format,
    format_confidence = format_confidence,
    mixed_format = mixed_format,
    issues = if (length(issues) > 0) issues else "No issues detected",
    pattern_scores = pattern_scores
  ))
}

# Validate numeric data
validate_numeric_data <- function(data_matrix) {
  
  # Check if all data is numeric
  numeric_columns <- sapply(data_matrix, is.numeric)
  all_numeric <- all(numeric_columns)
  
  # Count non-numeric values
  non_numeric_count <- 0
  if (!all_numeric) {
    non_numeric_cols <- names(data_matrix)[!numeric_columns]
    non_numeric_count <- sum(sapply(data_matrix[non_numeric_cols], function(x) sum(!is.na(x) & !is.numeric(x))))
  }
  
  # Check for missing values
  missing_values <- sum(is.na(data_matrix))
  missing_percentage <- missing_values / (nrow(data_matrix) * ncol(data_matrix)) * 100
  
  # Check for negative values (unusual in count data)
  if (all_numeric) {
    negative_values <- sum(data_matrix < 0, na.rm = TRUE)
    negative_percentage <- negative_values / sum(!is.na(data_matrix)) * 100
  } else {
    negative_values <- 0
    negative_percentage <- 0
  }
  
  return(list(
    all_numeric = all_numeric,
    non_numeric_count = non_numeric_count,
    missing_values = missing_values,
    missing_percentage = missing_percentage,
    negative_values = negative_values,
    negative_percentage = negative_percentage,
    issues = generate_numeric_issues(all_numeric, missing_percentage, negative_percentage)
  ))
}

# Generate issues list for numeric validation
generate_numeric_issues <- function(all_numeric, missing_pct, negative_pct) {
  issues <- c()
  
  if (!all_numeric) {
    issues <- c(issues, "Non-numeric values detected in expression data")
  }
  if (missing_pct > 10) {
    issues <- c(issues, paste("High percentage of missing values:", round(missing_pct, 1), "%"))
  }
  if (negative_pct > 1) {
    issues <- c(issues, paste("Negative values detected:", round(negative_pct, 1), "% of data"))
  }
  
  if (length(issues) == 0) {
    return("No issues detected")
  } else {
    return(issues)
  }
}

# Validate sample metadata
validate_sample_metadata <- function(sample_metadata) {
  
  results <- list()
  
  # Check required columns
  required_cols <- c("Sample_ID", "Condition")
  missing_cols <- required_cols[!required_cols %in% colnames(sample_metadata)]
  
  results$required_columns <- list(
    rule = DATA_VALIDATION_RULES$sample_metadata$column_validation,
    result = list(
      missing_columns = missing_cols,
      has_required = length(missing_cols) == 0,
      issues = if (length(missing_cols) > 0) paste("Missing required columns:", paste(missing_cols, collapse = ", ")) else "No issues detected"
    ),
    status = if (length(missing_cols) == 0) "pass" else "fail",
    timestamp = Sys.time()
  )
  
  # Check sample ID uniqueness
  duplicate_samples <- sample_metadata$Sample_ID[duplicated(sample_metadata$Sample_ID)]
  
  results$sample_uniqueness <- list(
    rule = DATA_VALIDATION_RULES$sample_metadata$uniqueness_validation,
    result = list(
      duplicate_samples = duplicate_samples,
      has_duplicates = length(duplicate_samples) > 0,
      issues = if (length(duplicate_samples) > 0) paste("Duplicate sample IDs:", paste(duplicate_samples, collapse = ", ")) else "No issues detected"
    ),
    status = if (length(duplicate_samples) == 0) "pass" else "fail",
    timestamp = Sys.time()
  )
  
  # Check condition distribution
  condition_counts <- table(sample_metadata$Condition)
  min_replicates <- min(condition_counts)
  
  results$condition_balance <- list(
    rule = DATA_VALIDATION_RULES$sample_metadata$balance_validation,
    result = list(
      condition_counts = as.list(condition_counts),
      min_replicates = min_replicates,
      is_balanced = min_replicates >= 2,
      issues = if (min_replicates < 2) "Some conditions have fewer than 2 replicates" else "No issues detected"
    ),
    status = if (min_replicates >= 2) "pass" else "warning",
    timestamp = Sys.time()
  )
  
  return(results)
}

# Validate data alignment between expression and metadata
validate_data_alignment <- function(expression_data, sample_metadata) {
  
  results <- list()
  
  # Check sample ID overlap
  expr_samples <- colnames(expression_data)
  meta_samples <- sample_metadata$Sample_ID
  
  common_samples <- intersect(expr_samples, meta_samples)
  expr_only <- setdiff(expr_samples, meta_samples)
  meta_only <- setdiff(meta_samples, expr_samples)
  
  overlap_percentage <- length(common_samples) / max(length(expr_samples), length(meta_samples)) * 100
  
  results$sample_overlap <- list(
    rule = DATA_VALIDATION_RULES$data_alignment$sample_matching,
    result = list(
      common_samples = common_samples,
      expr_only_samples = expr_only,
      meta_only_samples = meta_only,
      overlap_percentage = overlap_percentage,
      total_common = length(common_samples),
      issues = if (overlap_percentage < 50) {
        "Low overlap between expression and metadata sample IDs"
      } else if (length(expr_only) > 0 || length(meta_only) > 0) {
        "Some samples present in only one file"
      } else {
        "No issues detected"
      }
    ),
    status = if (overlap_percentage >= 80) "pass" else if (overlap_percentage >= 50) "warning" else "fail",
    timestamp = Sys.time()
  )
  
  return(results)
}

# Calculate overall validation status
calculate_validation_status <- function(validation_results) {
  
  # Extract all validation statuses
  all_statuses <- c()
  all_issues <- c()
  
  # Check expression matrix validation
  if (!is.null(validation_results$expression_matrix)) {
    for (check_name in names(validation_results$expression_matrix)) {
      check <- validation_results$expression_matrix[[check_name]]
      if (!is.null(check$status)) {
        all_statuses <- c(all_statuses, check$status)
        
        # Extract issues if status is warning or fail
        if (check$status %in% c("warning", "fail") && !is.null(check$result$issues)) {
          all_issues <- c(all_issues, check$result$issues)
        }
      }
    }
  }
  
  # Check sample metadata validation
  if (!is.null(validation_results$sample_metadata)) {
    for (check_name in names(validation_results$sample_metadata)) {
      check <- validation_results$sample_metadata[[check_name]]
      if (!is.null(check$status)) {
        all_statuses <- c(all_statuses, check$status)
        
        if (check$status %in% c("warning", "fail") && !is.null(check$result$issues)) {
          all_issues <- c(all_issues, check$result$issues)
        }
      }
    }
  }
  
  # Check data alignment validation
  if (!is.null(validation_results$data_alignment)) {
    for (check_name in names(validation_results$data_alignment)) {
      check <- validation_results$data_alignment[[check_name]]
      if (!is.null(check$status)) {
        all_statuses <- c(all_statuses, check$status)
        
        if (check$status %in% c("warning", "fail") && !is.null(check$result$issues)) {
          all_issues <- c(all_issues, check$result$issues)
        }
      }
    }
  }
  
  # Determine overall status
  if (any(all_statuses == "fail")) {
    overall_status <- "fail"
  } else if (any(all_statuses == "warning")) {
    overall_status <- "warning"
  } else {
    overall_status <- "pass"
  }
  
  return(list(
    status = overall_status,
    issues = all_issues,
    total_checks = length(all_statuses),
    passed_checks = sum(all_statuses == "pass"),
    warning_checks = sum(all_statuses == "warning"),
    failed_checks = sum(all_statuses == "fail")
  ))
}

# Export functions
list(
  create_real_data_upload_ui = create_real_data_upload_ui,
  validate_uploaded_data = validate_uploaded_data,
  calculate_validation_status = calculate_validation_status,
  detect_matrix_orientation = detect_matrix_orientation,
  validate_gene_ids = validate_gene_ids,
  UPLOAD_CONFIG = UPLOAD_CONFIG,
  DATA_VALIDATION_RULES = DATA_VALIDATION_RULES,
  GROUND_TRUTH_SCHEMA = GROUND_TRUTH_SCHEMA
)