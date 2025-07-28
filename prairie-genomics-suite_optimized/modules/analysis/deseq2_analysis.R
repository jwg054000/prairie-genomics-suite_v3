# DESeq2 Analysis Module - Optimized Version
# Enhanced differential expression analysis with memory management and async processing
# 
# Author: Prairie Genomics Team - Optimized Version
# Features: Smart memory management, async processing, comprehensive error handling

# Load dependencies
source("config/app_config.R")
source("utils/memory_manager.R")

# Source the gene conversion cache module
tryCatch({
  source("modules/data_processing/gene_conversion_cache.R")
  cat("✅ Gene conversion cache module loaded\n")
}, error = function(e) {
  cat("⚠️ Gene conversion cache module not available - gene conversion will be limited\n")
})

# Load required packages with graceful handling
deseq2_available <- FALSE

tryCatch({
  library(DESeq2)
  deseq2_available <- TRUE
  cat("✅ DESeq2 loaded successfully\n")
}, error = function(e) {
  cat("❌ DESeq2 not available:", e$message, "\n")
})

# Gene conversion function using v5 cache system
apply_gene_conversion_v5 <- function(results_df, species = "human") {
  tryCatch({
    if (!exists("convert_genes_fast")) {
      cat("⚠️ Gene conversion cache not available, skipping conversion\n")
      return(list(
        success = FALSE,
        data = results_df,
        conversion_rate = 0,
        message = "Gene conversion cache not available"
      ))
    }
    
    # Extract gene IDs from results
    gene_ids <- results_df$Gene
    cat("🔄 Converting", length(gene_ids), "gene IDs to symbols\n")
    
    # Detect species if not specified with improved pattern matching
    if (species == "auto") {
      # Count human vs mouse gene patterns
      human_patterns <- sum(grepl("^ENSG[0-9]", gene_ids))
      mouse_patterns <- sum(grepl("^ENSMUSG[0-9]", gene_ids))
      
      # Also check for gene symbols that might indicate species
      human_symbols <- sum(grepl("^[A-Z][A-Z0-9-]+$", gene_ids) & nchar(gene_ids) < 10)
      mouse_symbols <- sum(grepl("^[A-Z][a-z0-9-]+$", gene_ids) & nchar(gene_ids) < 10)
      
      if (human_patterns > mouse_patterns && human_patterns > 0) {
        species <- "human"
        cat("🔍 Detected human gene IDs (", human_patterns, "ENSG patterns)\n")
      } else if (mouse_patterns > human_patterns && mouse_patterns > 0) {
        species <- "mouse"
        cat("🔍 Detected mouse gene IDs (", mouse_patterns, "ENSMUSG patterns)\n")
      } else if (human_symbols > mouse_symbols) {
        species <- "human"
        cat("🔍 Detected likely human gene symbols (", human_symbols, "uppercase patterns)\n")
      } else if (mouse_symbols > human_symbols) {
        species <- "mouse"
        cat("🔍 Detected likely mouse gene symbols (", mouse_symbols, "mixed case patterns)\n")
      } else {
        species <- "human"  # Default
        cat("🔍 Species detection unclear, defaulting to human\n")
        cat("   Sample gene IDs:", paste(head(gene_ids, 3), collapse = ", "), "\n")
      }
    }
    
    # Use the v5 gene conversion system
    conversion_result <- convert_genes_fast(gene_ids, species = species, use_cache = TRUE)
    
    if (!is.null(conversion_result) && nrow(conversion_result) > 0) {
      # Add gene symbols to results
      results_with_symbols <- results_df
      
      # Match conversion results to original results
      symbol_matches <- match(results_df$Gene, conversion_result$ensembl_gene_id)
      valid_matches <- !is.na(symbol_matches)
      
      # Add gene symbol column
      results_with_symbols$gene_symbol <- NA
      results_with_symbols$gene_symbol[valid_matches] <- conversion_result$gene_symbol[symbol_matches[valid_matches]]
      
      # For genes without symbols, use the original gene ID
      missing_symbols <- is.na(results_with_symbols$gene_symbol) | results_with_symbols$gene_symbol == ""
      results_with_symbols$gene_symbol[missing_symbols] <- results_with_symbols$Gene[missing_symbols]
      
      # Calculate conversion statistics
      converted_count <- sum(!missing_symbols)
      conversion_rate <- round(100 * converted_count / nrow(results_df), 1)
      
      cat("📊 Gene conversion stats:\n")
      cat("   - Total genes:", nrow(results_df), "\n")
      cat("   - Successfully converted:", converted_count, "\n")
      cat("   - Conversion rate:", conversion_rate, "%\n")
      
      return(list(
        success = TRUE,
        data = results_with_symbols,
        conversion_rate = conversion_rate,
        message = paste("Converted", converted_count, "of", nrow(results_df), "genes")
      ))
    } else {
      cat("❌ Gene conversion returned no results\n")
      return(list(
        success = FALSE,
        data = results_df,
        conversion_rate = 0,
        message = "Gene conversion failed"
      ))
    }
    
  }, error = function(e) {
    cat("❌ Error in gene conversion:", e$message, "\n")
    return(list(
      success = FALSE,
      data = results_df,
      conversion_rate = 0,
      message = paste("Gene conversion error:", e$message)
    ))
  })
}

# Enhanced DESeq2 analysis with memory management
run_deseq2_analysis <- function(expression_data, annotation_data, 
                               contrast_group1, contrast_group2,
                               padj_cutoff = DEFAULT_PADJ, 
                               fc_cutoff = DEFAULT_FC) {
  
  # Monitor this operation
  monitor_operation("deseq2_analysis", function() {
    
    tryCatch({
      cat("🚀 Starting DESeq2 differential expression analysis\n")
      
      # DEBUG: Check what data we're actually getting
      cat("🔍 DEBUG: Expression data dimensions:", nrow(expression_data), "x", ncol(expression_data), "\n")
      cat("🔍 DEBUG: Sample gene IDs:", paste(head(rownames(expression_data), 5), collapse = ", "), "\n")
      cat("🔍 DEBUG: Sample column names:", paste(head(colnames(expression_data), 3), collapse = ", "), "\n")
      
      # Validate inputs
      validation_result <- validate_deseq2_inputs(expression_data, annotation_data, 
                                                 contrast_group1, contrast_group2)
      if (!validation_result$valid) {
        return(list(
          success = FALSE,
          error = validation_result$error,
          data = NULL
        ))
      }
      
      cat("✅ Input validation passed\n")
      
      # Check if DESeq2 is available
      if (!deseq2_available) {
        return(run_deseq2_fallback(expression_data, annotation_data, 
                                 contrast_group1, contrast_group2,
                                 padj_cutoff, fc_cutoff))
      }
      
      # Prepare data for DESeq2
      cat("📊 Preparing data for DESeq2...\n")
      
      # Ensure integer counts for DESeq2
      count_matrix <- round(expression_data)
      count_matrix[count_matrix < 0] <- 0
      
      # Align annotation with expression data
      common_samples <- intersect(colnames(count_matrix), annotation_data$Sample)
      if (length(common_samples) == 0) {
        return(list(
          success = FALSE,
          error = "No matching samples between expression data and annotation",
          data = NULL
        ))
      }
      
      # Filter to common samples
      count_matrix <- count_matrix[, common_samples, drop = FALSE]
      annotation_filtered <- annotation_data[annotation_data$Sample %in% common_samples, ]
      
      # Ensure annotation order matches count matrix
      annotation_filtered <- annotation_filtered[match(colnames(count_matrix), annotation_filtered$Sample), ]
      
      cat("📋 Analysis setup:\n")
      cat("   - Genes:", nrow(count_matrix), "\n")
      cat("   - Samples:", ncol(count_matrix), "\n")
      cat("   - Contrast:", contrast_group1, "vs", contrast_group2, "\n")
      
      # Create DESeq2 dataset
      cat("🧬 Creating DESeq2 dataset...\n")
      
      # Memory check before DESeq2
      mem_status <- get_memory_status()
      if (mem_status$warning) {
        smart_gc()
      }
      
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = annotation_filtered,
        design = ~ Condition
      )
      
      # Filter low count genes to reduce memory usage
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]
      
      cat("🔬 Running DESeq2 analysis...\n")
      
      # Run DESeq2 with error handling
      dds <- tryCatch({
        DESeq(dds, quiet = TRUE)
      }, error = function(e) {
        cat("❌ DESeq2 analysis failed:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(dds)) {
        return(list(
          success = FALSE,
          error = "DESeq2 analysis failed during model fitting",
          data = NULL
        ))
      }
      
      # Extract results
      cat("📊 Extracting results...\n")
      
      results_obj <- tryCatch({
        results(dds, contrast = c("Condition", contrast_group1, contrast_group2))
      }, error = function(e) {
        cat("❌ Failed to extract results:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(results_obj)) {
        return(list(
          success = FALSE,
          error = "Failed to extract DESeq2 results",
          data = NULL
        ))
      }
      
      # Convert to data frame
      results_df <- as.data.frame(results_obj)
      results_df$Gene <- rownames(results_df)
      
      # Remove rows with NA padj
      results_df <- results_df[!is.na(results_df$padj), ]
      
      # Apply gene conversion using v5 cache system with auto species detection
      cat("🧬 Converting gene IDs to symbols...\n")
      gene_conversion_result <- apply_gene_conversion_v5(results_df, species = "auto")
      if (gene_conversion_result$success) {
        results_df <- gene_conversion_result$data
        cat("✅ Gene conversion completed:", gene_conversion_result$conversion_rate, "% success rate\n")
      } else {
        cat("⚠️ Gene conversion failed, using original gene IDs\n")
      }
      
      # Calculate statistics
      significant_genes <- sum(results_df$padj < padj_cutoff & 
                             abs(results_df$log2FoldChange) > fc_cutoff, na.rm = TRUE)
      upregulated <- sum(results_df$padj < padj_cutoff & 
                        results_df$log2FoldChange > fc_cutoff, na.rm = TRUE)
      downregulated <- sum(results_df$padj < padj_cutoff & 
                          results_df$log2FoldChange < -fc_cutoff, na.rm = TRUE)
      
      cat("✅ DESeq2 analysis completed successfully:\n")
      cat("   - Total genes analyzed:", nrow(results_df), "\n")
      cat("   - Significant genes:", significant_genes, "\n")
      cat("   - Upregulated:", upregulated, "\n")
      cat("   - Downregulated:", downregulated, "\n")
      
      # Cleanup
      rm(dds, results_obj)
      smart_gc()
      
      return(list(
        success = TRUE,
        data = results_df,
        stats = list(
          total_genes = nrow(results_df),
          significant_genes = significant_genes,
          upregulated = upregulated,
          downregulated = downregulated,
          contrast = paste(contrast_group1, "vs", contrast_group2),
          padj_cutoff = padj_cutoff,
          fc_cutoff = fc_cutoff
        )
      ))
      
    }, error = function(e) {
      cat("❌ Unexpected error in DESeq2 analysis:", e$message, "\n")
      smart_gc()  # Cleanup on error
      
      return(list(
        success = FALSE,
        error = paste("DESeq2 analysis failed:", e$message),
        data = NULL
      ))
    })
  })
}

# Input validation for DESeq2 analysis
validate_deseq2_inputs <- function(expression_data, annotation_data, group1, group2) {
  result <- list(valid = TRUE, error = NULL)
  
  # Check expression data
  if (is.null(expression_data) || nrow(expression_data) == 0) {
    result$valid <- FALSE
    result$error <- "Expression data is empty or null"
    return(result)
  }
  
  # Check annotation data
  if (is.null(annotation_data) || nrow(annotation_data) == 0) {
    result$valid <- FALSE
    result$error <- "Annotation data is empty or null"
    return(result)
  }
  
  # Check required columns
  if (!"Sample" %in% colnames(annotation_data)) {
    result$valid <- FALSE
    result$error <- "Annotation data must have 'Sample' column"
    return(result)
  }
  
  if (!"Condition" %in% colnames(annotation_data)) {
    result$valid <- FALSE
    result$error <- "Annotation data must have 'Condition' column"
    return(result)
  }
  
  # Check contrast groups exist
  available_conditions <- unique(annotation_data$Condition)
  if (!group1 %in% available_conditions) {
    result$valid <- FALSE
    result$error <- paste("Group", group1, "not found in conditions")
    return(result)
  }
  
  if (!group2 %in% available_conditions) {
    result$valid <- FALSE
    result$error <- paste("Group", group2, "not found in conditions")
    return(result)
  }
  
  # Check minimum samples per group
  group1_samples <- sum(annotation_data$Condition == group1)
  group2_samples <- sum(annotation_data$Condition == group2)
  
  if (group1_samples < 2 || group2_samples < 2) {
    result$valid <- FALSE
    result$error <- "Each group must have at least 2 samples"
    return(result)
  }
  
  return(result)
}

# Fallback analysis when DESeq2 is not available
run_deseq2_fallback <- function(expression_data, annotation_data, 
                               group1, group2, padj_cutoff, fc_cutoff) {
  cat("⚠️ Running simplified analysis (DESeq2 not available)\n")
  
  # Simple t-test based analysis
  tryCatch({
    # Get samples for each group
    group1_samples <- annotation_data$Sample[annotation_data$Condition == group1]
    group2_samples <- annotation_data$Sample[annotation_data$Condition == group2]
    
    # Filter expression data
    group1_data <- expression_data[, group1_samples, drop = FALSE]
    group2_data <- expression_data[, group2_samples, drop = FALSE]
    
    # Calculate means and fold changes
    group1_mean <- rowMeans(group1_data, na.rm = TRUE)
    group2_mean <- rowMeans(group2_data, na.rm = TRUE)
    
    # Avoid division by zero
    group2_mean[group2_mean == 0] <- 0.1
    
    log2FC <- log2((group1_mean + 0.1) / (group2_mean + 0.1))
    
    # Simple t-test
    pvalues <- apply(cbind(group1_data, group2_data), 1, function(x) {
      g1_vals <- x[1:length(group1_samples)]
      g2_vals <- x[(length(group1_samples) + 1):length(x)]
      
      tryCatch({
        t.test(g1_vals, g2_vals)$p.value
      }, error = function(e) {
        return(1)
      })
    })
    
    # Simple multiple testing correction
    padj <- p.adjust(pvalues, method = "BH")
    
    # Create results data frame
    results_df <- data.frame(
      Gene = rownames(expression_data),
      baseMean = (group1_mean + group2_mean) / 2,
      log2FoldChange = log2FC,
      pvalue = pvalues,
      padj = padj,
      stringsAsFactors = FALSE
    )
    
    # Remove genes with NA values
    results_df <- results_df[!is.na(results_df$padj), ]
    
    # Calculate statistics
    significant_genes <- sum(results_df$padj < padj_cutoff & 
                           abs(results_df$log2FoldChange) > fc_cutoff, na.rm = TRUE)
    
    cat("✅ Simplified analysis completed (", significant_genes, "significant genes)\n")
    
    return(list(
      success = TRUE,
      data = results_df,
      stats = list(
        total_genes = nrow(results_df),
        significant_genes = significant_genes,
        method = "simplified_t_test"
      ),
      warning = "Simplified analysis used - install DESeq2 for full functionality"
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = paste("Fallback analysis failed:", e$message),
      data = NULL
    ))
  })
}

# Note: Gene symbol conversion is now handled using the v5 cache system
# This provides fast, cached conversion with multiple fallback strategies

# DESeq2 Analysis UI Module
deseq2AnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        4,
        wellPanel(
          h5("🎯 Analysis Settings"),
          
          uiOutput(ns("contrast_selection")),
          
          br(),
          
          h6("📊 Filtering Parameters"),
          
          numericInput(
            ns("padj_cutoff"),
            "Adjusted p-value cutoff:",
            value = DEFAULT_PADJ,
            min = 0.001,
            max = 0.1,
            step = 0.001
          ),
          
          numericInput(
            ns("fc_cutoff"), 
            "Log2 fold change cutoff:",
            value = DEFAULT_FC,
            min = 0.1,
            max = 5.0,
            step = 0.1
          ),
          
          br(),
          
          div(
            style = "margin-top: 10px; padding: 8px; background-color: #e7f3ff; border-radius: 3px;",
            h6("🧬 Gene Symbols", style = "margin: 0; color: #0066cc;"),
            tags$small("Gene symbols are automatically converted using the v5 cache system", style = "color: #0066cc;")
          ),
          
          # Memory status indicator
          div(
            style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
            h6("💾 Memory Status"),
            uiOutput(ns("memory_status"))
          )
        )
      ),
      
      column(
        8,
        # Analysis controls
        conditionalPanel(
          condition = paste0("!output['", ns("analysis_running"), "']"),
          
          actionButton(
            ns("run_analysis"),
            "🚀 Run DESeq2 Analysis",
            class = "btn-primary btn-lg",
            style = "width: 100%; height: 60px; font-size: 18px; margin-bottom: 20px;"
          )
        ),
        
        # Progress indicator
        conditionalPanel(
          condition = paste0("output['", ns("analysis_running"), "']"),
          
          div(
            class = "alert alert-info",
            h5("🔄 Running DESeq2 Analysis..."),
            progressBar(
              id = ns("analysis_progress"),
              value = 0,
              total = 100,
              status = "primary",
              display_pct = TRUE
            ),
            uiOutput(ns("analysis_status"))
          )
        ),
        
        # Results summary
        conditionalPanel(
          condition = paste0("output['", ns("show_results"), "']"),
          
          div(
            class = "alert alert-success",
            h5("✅ Analysis Complete"),
            uiOutput(ns("results_summary"))
          )
        )
      )
    )
  )
}

# DESeq2 Analysis Server Module
deseq2AnalysisServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for this module
    local_values <- reactiveValues(
      analysis_running = FALSE,
      current_results = NULL,
      available_groups = character(0)
    )
    
    # Update available groups when annotation data changes
    observe({
      if (!is.null(values$annotation_data)) {
        unique_groups <- unique(values$annotation_data$Condition)
        unique_groups <- unique_groups[!is.na(unique_groups)]
        local_values$available_groups <- unique_groups
        
        cat("📊 Available groups for DESeq2:", paste(unique_groups, collapse = ", "), "\n")
      }
    })
    
    # Contrast selection UI
    output$contrast_selection <- renderUI({
      req(local_values$available_groups)
      
      if (length(local_values$available_groups) < 2) {
        return(
          div(
            class = "alert alert-warning",
            h6("⚠️ Insufficient Groups"),
            p("Need at least 2 groups for differential analysis."),
            p("Please check your sample annotation.")
          )
        )
      }
      
      tagList(
        h6("🔬 Select Comparison Groups"),
        
        selectInput(
          session$ns("group1"),
          "Group 1 (numerator):",
          choices = local_values$available_groups,
          selected = local_values$available_groups[1]
        ),
        
        selectInput(
          session$ns("group2"),
          "Group 2 (denominator):",
          choices = local_values$available_groups,
          selected = if (length(local_values$available_groups) > 1) local_values$available_groups[2] else NULL
        ),
        
        # Preview comparison
        div(
          style = "margin-top: 10px; padding: 8px; background-color: #e9ecef; border-radius: 3px;",
          tags$small(
            "📊 Comparison: ",
            textOutput(session$ns("comparison_preview"), inline = TRUE)
          )
        )
      )
    })
    
    # Comparison preview
    output$comparison_preview <- renderText({
      if (!is.null(input$group1) && !is.null(input$group2) && input$group1 != input$group2) {
        paste(input$group1, "vs", input$group2)
      } else {
        "Please select different groups"
      }
    })
    
    # Memory status monitoring
    output$memory_status <- renderUI({
      memory_status <- get_memory_status()
      
      status_class <- if (memory_status$critical) {
        "danger"
      } else if (memory_status$warning) {
        "warning"
      } else {
        "success"
      }
      
      div(
        class = paste("alert alert-", status_class),
        style = "padding: 5px; margin: 0;",
        tags$small(memory_status$message)
      )
    })
    
    # Analysis running status
    output$analysis_running <- reactive({
      local_values$analysis_running
    })
    outputOptions(output, "analysis_running", suspendWhenHidden = FALSE)
    
    # Run DESeq2 analysis
    observeEvent(input$run_analysis, {
      # Validate prerequisites
      if (is.null(values$expression_data)) {
        showNotification("❌ No expression data available", type = "error")
        return()
      }
      
      if (is.null(values$annotation_data)) {
        showNotification("❌ No annotation data available", type = "error")
        return()
      }
      
      if (is.null(input$group1) || is.null(input$group2)) {
        showNotification("❌ Please select comparison groups", type = "error")
        return()
      }
      
      if (input$group1 == input$group2) {
        showNotification("❌ Please select different groups for comparison", type = "error")
        return()
      }
      
      cat("🚀 Starting DESeq2 analysis:", input$group1, "vs", input$group2, "\n")
      
      # Set analysis running
      local_values$analysis_running <- TRUE
      
      # Update progress
      update_progress <- function(value, message = "") {
        if (requireNamespace("shinyWidgets", quietly = TRUE)) {
          shinyWidgets::updateProgressBar(session, "analysis_progress", value = value)
        }
        
        output$analysis_status <- renderUI({
          p(message)
        })
      }
      
      update_progress(10, "Initializing DESeq2 analysis...")
      
      # Run analysis with memory monitoring
      analysis_result <- monitor_operation("deseq2_analysis", function() {
        run_deseq2_analysis(
          expression_data = values$expression_data,
          annotation_data = values$annotation_data,
          contrast_group1 = input$group1,
          contrast_group2 = input$group2,
          padj_cutoff = input$padj_cutoff %||% DEFAULT_PADJ,
          fc_cutoff = input$fc_cutoff %||% DEFAULT_FC
        )
      })
      
      update_progress(80, "Processing results...")
      
      # Handle results
      if (analysis_result$result$success) {
        local_values$current_results <- analysis_result$result
        values$deseq2_results <- analysis_result$result$data
        
        update_progress(100, "Analysis completed successfully!")
        
        # Calculate summary statistics
        stats <- analysis_result$result$stats
        
        showNotification(
          paste("✅ DESeq2 analysis completed!",
                "Found", stats$significant_genes, "significant genes out of", stats$total_genes, "total genes"),
          type = "message",
          duration = 8
        )
        
        cat("✅ DESeq2 analysis completed successfully\n")
        cat("   - Total genes:", stats$total_genes, "\n")
        cat("   - Significant genes:", stats$significant_genes, "\n")
        
      } else {
        showNotification(
          paste("❌ DESeq2 analysis failed:", analysis_result$result$error),
          type = "error",
          duration = 15
        )
        
        cat("❌ DESeq2 analysis failed:", analysis_result$result$error, "\n")
      }
      
      # Reset analysis running status
      local_values$analysis_running <- FALSE
    })
    
    # Show results status
    output$show_results <- reactive({
      !is.null(local_values$current_results) && 
      local_values$current_results$success && 
      !local_values$analysis_running
    })
    outputOptions(output, "show_results", suspendWhenHidden = FALSE)
    
    # Results summary
    output$results_summary <- renderUI({
      req(local_values$current_results)
      req(local_values$current_results$success)
      
      stats <- local_values$current_results$stats
      
      div(
        h6("📊 Analysis Results Summary"),
        
        fluidRow(
          column(
            4,
            tags$strong("Total Genes:"), br(),
            formatC(stats$total_genes, format="d", big.mark=",")
          ),
          column(
            4,
            tags$strong("Significant Genes:"), br(),
            formatC(stats$significant_genes, format="d", big.mark=",")
          ),
          column(
            4,
            tags$strong("Significance Rate:"), br(),
            paste0(round(100 * stats$significant_genes / stats$total_genes, 1), "%")
          )
        ),
        
        br(),
        
        if (!is.null(stats$upregulated) && !is.null(stats$downregulated)) {
          fluidRow(
            column(
              6,
              tags$strong("📈 Upregulated:"), br(),
              formatC(stats$upregulated, format="d", big.mark=","),
              if (stats$significant_genes > 0) {
                paste0(" (", round(100 * stats$upregulated / stats$significant_genes, 1), "%)")
              } else {
                " (0%)"
              }
            ),
            column(
              6,
              tags$strong("📉 Downregulated:"), br(),
              formatC(stats$downregulated, format="d", big.mark=","),
              if (stats$significant_genes > 0) {
                paste0(" (", round(100 * stats$downregulated / stats$significant_genes, 1), "%)")
              } else {
                " (0%)"
              }
            )
          )
        },
        
        br(),
        
        div(
          class = "alert alert-info",
          style = "padding: 8px; margin-top: 10px;",
          tags$small(
            "🔬 Comparison: ", stats$contrast, br(),
            "📊 Filters: padj < ", stats$padj_cutoff, ", |log2FC| > ", stats$fc_cutoff,
            if (!is.null(local_values$current_results$warning)) {
              tagList(br(), "⚠️ ", local_values$current_results$warning)
            }
          )
        )
      )
    })
    
    # Return reactive for results
    return(reactive({
      local_values$current_results
    }))
  })
}

cat("✅ Optimized DESeq2 analysis module loaded\n")