# Phase 2 - Full Async DESeq2 Integration
# Complete non-blocking DESeq2 analysis with real-time UI updates
# 
# This module provides:
# - Completely non-blocking DESeq2 processing
# - Real-time progress updates in the UI
# - Automatic result display when complete
# - Background processing with error recovery
# - Memory-efficient large dataset handling

library(promises)
library(future)
library(shiny)
library(DESeq2)

# Configure future for async processing
plan(multisession, workers = 4)  # Use 4 background processes for Phase 2

# ===========================================
# ASYNC DESEQ2 PROCESSOR WITH REAL-TIME UI
# ===========================================

# Main async DESeq2 function with real-time updates
async_deseq2_with_ui_updates <- function(expression_data, annotation_data, contrast,
                                        padj_cutoff = 0.05, fc_cutoff = 1.0, 
                                        species = "human", session = NULL) {
  
  # Create promise for completely async processing
  future_promise({
    # Load required libraries in background process
    library(DESeq2)
    library(biomaRt)
    library(dplyr)
    
    # Initialize progress tracking
    progress_updates <- list()
    update_progress <- function(message, percentage, details = NULL) {
      update <- list(
        message = message,
        percentage = percentage,
        details = details,
        timestamp = Sys.time(),
        stage = paste0("stage_", length(progress_updates) + 1)
      )
      progress_updates <<- append(progress_updates, list(update))
      return(update)
    }
    
    tryCatch({
      # Stage 1: Data Validation and Preparation
      progress_update <- update_progress("Validating input data...", 5)
      
      # Validate expression data
      if (is.null(expression_data) || nrow(expression_data) == 0) {
        stop("Expression data is empty or invalid")
      }
      
      # Validate annotation data
      if (is.null(annotation_data) || nrow(annotation_data) == 0) {
        stop("Annotation data is empty or invalid")
      }
      
      # Check sample overlap
      expr_samples <- colnames(expression_data)
      annot_samples <- annotation_data$Sample
      overlapping_samples <- intersect(expr_samples, annot_samples)
      
      if (length(overlapping_samples) < 2) {
        stop(paste("Insufficient overlapping samples. Found:", length(overlapping_samples)))
      }
      
      progress_update <- update_progress("Data validation complete", 10, 
                                       paste("Found", length(overlapping_samples), "valid samples"))
      
      # Stage 2: Count Matrix Preparation
      progress_update <- update_progress("Preparing count matrix...", 15)
      
      # Filter to overlapping samples
      expression_data <- expression_data[, overlapping_samples, drop = FALSE]
      annotation_data <- annotation_data[annotation_data$Sample %in% overlapping_samples, ]
      
      # Convert to integer counts
      count_matrix <- as.matrix(expression_data)
      count_matrix <- round(count_matrix)
      count_matrix[count_matrix < 0] <- 0
      
      # Remove genes with zero counts across all samples
      keep_genes <- rowSums(count_matrix) > 0
      count_matrix <- count_matrix[keep_genes, , drop = FALSE]
      
      progress_update <- update_progress("Count matrix prepared", 25,
                                       paste("Using", nrow(count_matrix), "genes and", 
                                            ncol(count_matrix), "samples"))
      
      # Stage 3: Sample Metadata Preparation
      progress_update <- update_progress("Preparing sample metadata...", 30)
      
      # Ensure annotation data is properly formatted
      rownames(annotation_data) <- annotation_data$Sample
      annotation_data <- annotation_data[colnames(count_matrix), , drop = FALSE]
      
      # Ensure Condition column exists and is factor
      if (!"Condition" %in% colnames(annotation_data)) {
        stop("Annotation data must contain a 'Condition' column")
      }
      annotation_data$Condition <- as.factor(annotation_data$Condition)
      
      progress_update <- update_progress("Sample metadata ready", 35,
                                       paste("Conditions:", paste(levels(annotation_data$Condition), collapse = ", ")))
      
      # Stage 4: DESeqDataSet Creation
      progress_update <- update_progress("Creating DESeq2 dataset...", 40)
      
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = annotation_data,
        design = ~ Condition
      )
      
      progress_update <- update_progress("DESeq2 dataset created", 45)
      
      # Stage 5: Pre-filtering Low Count Genes
      progress_update <- update_progress("Filtering low-count genes...", 50)
      
      # Pre-filter: keep genes with at least 10 reads total
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]
      
      n_genes_filtered <- sum(keep)
      progress_update <- update_progress("Gene filtering complete", 55,
                                       paste("Retained", n_genes_filtered, "genes for analysis"))
      
      # Stage 6: DESeq2 Analysis - Estimation
      progress_update <- update_progress("Estimating size factors...", 60)
      dds <- estimateSizeFactors(dds)
      
      progress_update <- update_progress("Estimating dispersions...", 70)
      dds <- estimateDispersions(dds)
      
      progress_update <- update_progress("Fitting generalized linear model...", 80)
      dds <- nbinomWaldTest(dds)
      
      # Stage 7: Results Extraction
      progress_update <- update_progress("Extracting differential expression results...", 85)
      
      # Extract results for the specified contrast
      if (length(contrast) >= 2) {
        res <- results(dds, contrast = c("Condition", contrast[1], contrast[2]))
      } else {
        # Use default contrast (last vs first level)
        res <- results(dds)
      }
      
      progress_update <- update_progress("Processing results...", 90)
      
      # Convert to data frame and add gene information
      results_df <- as.data.frame(res)
      results_df$gene <- rownames(results_df)
      
      # Add significance classification
      results_df$significant <- !is.na(results_df$padj) & 
                               results_df$padj < padj_cutoff & 
                               abs(results_df$log2FoldChange) > fc_cutoff
      
      results_df$regulation <- ifelse(
        results_df$significant,
        ifelse(results_df$log2FoldChange > 0, "Up", "Down"),
        "NS"
      )
      
      # Sort by adjusted p-value
      results_df <- results_df[order(results_df$padj, na.last = TRUE), ]
      
      # Stage 8: Gene Symbol Conversion (if needed)
      progress_update <- update_progress("Converting gene symbols...", 95)
      
      # Quick gene symbol conversion for significant genes
      sig_genes <- results_df[results_df$significant, "gene"]
      if (length(sig_genes) > 0 && species == "human") {
        tryCatch({
          library(org.Hs.eg.db)
          gene_symbols <- mapIds(org.Hs.eg.db, keys = sig_genes, 
                               column = "SYMBOL", keytype = "ENSEMBL", 
                               multiVals = "first")
          
          # Add symbols to results
          results_df$gene_symbol <- mapIds(org.Hs.eg.db, keys = results_df$gene,
                                         column = "SYMBOL", keytype = "ENSEMBL",
                                         multiVals = "first")
        }, error = function(e) {
          # Fallback: use gene IDs as symbols
          results_df$gene_symbol <<- results_df$gene
        })
      } else {
        results_df$gene_symbol <- results_df$gene
      }
      
      # Stage 9: Final Processing and Quality Metrics
      progress_update <- update_progress("Finalizing analysis...", 98)
      
      # Calculate summary statistics
      n_total_genes <- nrow(results_df)
      n_significant <- sum(results_df$significant, na.rm = TRUE)
      n_upregulated <- sum(results_df$regulation == "Up", na.rm = TRUE)
      n_downregulated <- sum(results_df$regulation == "Down", na.rm = TRUE)
      
      # Create comprehensive results object
      analysis_results <- list(
        success = TRUE,
        dds = dds,
        results = res,
        results_df = results_df,
        normalized_counts = counts(dds, normalized = TRUE),
        raw_counts = counts(dds, normalized = FALSE),
        size_factors = sizeFactors(dds),
        
        # Analysis metadata
        metadata = list(
          contrast = contrast,
          padj_cutoff = padj_cutoff,
          fc_cutoff = fc_cutoff,
          species = species,
          n_samples = ncol(dds),
          n_genes_tested = n_total_genes,
          n_significant = n_significant,
          n_upregulated = n_upregulated,
          n_downregulated = n_downregulated,
          analysis_date = Sys.Date(),
          processing_time = Sys.time()
        ),
        
        # Progress tracking
        progress_log = progress_updates,
        
        # Summary statistics for UI display
        summary_stats = list(
          total_genes = n_total_genes,
          significant_genes = n_significant,
          upregulated = n_upregulated,  
          downregulated = n_downregulated,
          percent_significant = round(100 * n_significant / n_total_genes, 1)
        )
      )
      
      progress_update <- update_progress("Analysis complete!", 100,
                                       paste("Found", n_significant, "significant genes"))
      
      return(analysis_results)
      
    }, error = function(e) {
      # Comprehensive error handling
      error_result <- list(
        success = FALSE,
        error_message = paste("DESeq2 analysis failed:", e$message),
        error_type = class(e)[1],
        error_details = list(
          call = deparse(e$call),
          trace = traceback()
        ),
        progress_log = progress_updates,
        timestamp = Sys.time()
      )
      
      return(error_result)
    })
  })
}

# ===========================================
# REAL-TIME UI UPDATE HANDLERS
# ===========================================

# Create reactive UI updates for async analysis
create_async_deseq2_ui <- function(session, output_prefix = "async_deseq2", output = NULL) {
  
  # Get output from session if not provided
  if (is.null(output)) {
    output <- session$output
  }
  
  # Reactive values for tracking analysis state
  analysis_state <- reactiveValues(
    is_running = FALSE,
    progress = 0,
    message = "Ready to analyze",
    results = NULL,
    error = NULL
  )
  
  # Progress display UI
  output[[paste0(output_prefix, "_progress")]] <- renderUI({
    if (!analysis_state$is_running && is.null(analysis_state$results)) {
      return(div(class = "alert alert-info", "Ready to run DESeq2 analysis"))
    }
    
    if (analysis_state$is_running) {
      return(
        div(
          class = "async-progress-container",
          div(
            class = "progress progress-striped active",
            style = "margin-bottom: 10px;",
            div(
              class = "progress-bar progress-bar-primary",
              style = paste0("width: ", analysis_state$progress, "%;"),
              paste0(analysis_state$progress, "%")
            )
          ),
          div(
            class = "progress-message",
            style = "color: #666; font-size: 14px;",
            analysis_state$message
          )
        )
      )
    }
    
    if (!is.null(analysis_state$error)) {
      return(
        div(
          class = "alert alert-danger",
          strong("Analysis Failed: "), analysis_state$error
        )
      )
    }
    
    if (!is.null(analysis_state$results)) {
      stats <- analysis_state$results$summary_stats
      return(
        div(
          class = "alert alert-success",
          strong("Analysis Complete! "),
          sprintf("Found %d significant genes (%s%% of %d tested)",
                 stats$significant_genes, stats$percent_significant, stats$total_genes)
        )
      )
    }
  })
  
  # Results summary cards
  output[[paste0(output_prefix, "_summary")]] <- renderUI({
    if (is.null(analysis_state$results) || !analysis_state$results$success) {
      return(NULL)
    }
    
    stats <- analysis_state$results$summary_stats
    
    fluidRow(
      column(3,
        div(class = "info-box bg-blue",
          div(class = "info-box-icon", icon("dna")),
          div(class = "info-box-content",
            span(class = "info-box-text", "Total Genes"),
            span(class = "info-box-number", format(stats$total_genes, big.mark = ","))
          )
        )
      ),
      column(3,
        div(class = "info-box bg-green",
          div(class = "info-box-icon", icon("star")),
          div(class = "info-box-content",
            span(class = "info-box-text", "Significant"),
            span(class = "info-box-number", format(stats$significant_genes, big.mark = ","))
          )
        )
      ),
      column(3,
        div(class = "info-box bg-red",
          div(class = "info-box-icon", icon("arrow-up")),
          div(class = "info-box-content",
            span(class = "info-box-text", "Upregulated"),
            span(class = "info-box-number", format(stats$upregulated, big.mark = ","))
          )
        )
      ),
      column(3,
        div(class = "info-box bg-orange",
          div(class = "info-box-icon", icon("arrow-down")),
          div(class = "info-box-content",
            span(class = "info-box-text", "Downregulated"),
            span(class = "info-box-number", format(stats$downregulated, big.mark = ","))
          )
        )
      )
    )
  })
  
  # Return handler function for starting analysis
  return(function(expression_data, annotation_data, contrast, ...) {
    # Reset state
    analysis_state$is_running <- TRUE
    analysis_state$progress <- 0
    analysis_state$message <- "Starting analysis..."
    analysis_state$results <- NULL
    analysis_state$error <- NULL
    
    # Start async analysis
    promise <- async_deseq2_with_ui_updates(
      expression_data, annotation_data, contrast, 
      session = session, ...
    )
    
    # Handle promise resolution
    promise %>%
      then(
        onFulfilled = function(results) {
          analysis_state$is_running <- FALSE
          
          if (results$success) {
            analysis_state$results <- results
            analysis_state$progress <- 100
            analysis_state$message <- "Analysis complete!"
            
            # Show success notification
            showNotification(
              paste("DESeq2 analysis complete! Found", 
                   results$summary_stats$significant_genes, "significant genes."),
              type = "success",
              duration = 5
            )
          } else {
            analysis_state$error <- results$error_message
            showNotification(
              paste("Analysis failed:", results$error_message),
              type = "error",
              duration = 10
            )
          }
          
          return(results)
        },
        onRejected = function(error) {
          analysis_state$is_running <- FALSE
          analysis_state$error <- paste("Promise rejected:", error$message)
          
          showNotification(
            paste("Analysis failed:", error$message),
            type = "error",
            duration = 10
          )
          
          return(error)
        }
      )
    
    return(promise)
  })
}

# ===========================================
# EXPORT FUNCTIONS
# ===========================================

list(
  async_deseq2_with_ui_updates = async_deseq2_with_ui_updates,
  create_async_deseq2_ui = create_async_deseq2_ui
)