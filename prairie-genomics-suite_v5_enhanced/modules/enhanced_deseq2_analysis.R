# Enhanced DESeq2 Analysis Module for Prairie Genomics Suite v5
# Based on Emory methodology with multi-group support and batch correction
# Implements advanced statistical analysis with multiple contrasts

# Load required libraries
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  DESeq2, dplyr, sva, MASS, ggplot2, 
  pheatmap, RColorBrewer, plotly, DT
)

# UI Function
enhancedDESeq2AnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        4,
        h4("🚀 Enhanced DESeq2 Analysis"),
        
        wellPanel(
          h5("📊 Analysis Configuration"),
          
          # Multi-group analysis options
          uiOutput(ns("analysis_type_selection")),
          
          conditionalPanel(
            condition = paste0("input['", ns("analysis_type"), "'] == 'multiple_pairwise'"),
            
            h6("🔄 Pairwise Comparisons"),
            uiOutput(ns("pairwise_selection")),
            
            checkboxInput(
              ns("run_all_pairwise"),
              "Run all possible pairwise comparisons",
              value = TRUE
            )
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("analysis_type"), "'] == 'custom_contrast'"),
            
            h6("🎯 Custom Contrast"),
            uiOutput(ns("custom_contrast_ui"))
          ),
          
          br(),
          
          h5("⚙️ Statistical Parameters"),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("padj_cutoff"),
                "Adj. p-value cutoff:",
                value = 0.05,
                min = 0.001,
                max = 0.1,
                step = 0.001
              )
            ),
            column(
              6,
              numericInput(
                ns("fc_cutoff"),
                "Log₂FC cutoff:",
                value = 1.0,
                min = 0.1,
                max = 5.0,
                step = 0.1
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("min_count_threshold"),
                "Min count threshold:",
                value = 10,
                min = 1,
                max = 100
              )
            ),
            column(
              6,
              numericInput(
                ns("min_samples_expressed"),
                "Min samples expressed:",
                value = 3,
                min = 1,
                max = 20
              )
            )
          )
        ),
        
        wellPanel(
          h5("🧪 Advanced Options"),
          
          selectInput(
            ns("fit_type"),
            "Fit type:",
            choices = list(
              "Parametric" = "parametric",
              "Local" = "local",
              "Mean" = "mean"
            ),
            selected = "parametric"
          ),
          
          selectInput(
            ns("test_type"),
            "Statistical test:",
            choices = list(
              "Wald" = "Wald",
              "Likelihood Ratio Test" = "LRT"
            ),
            selected = "Wald"
          ),
          
          checkboxInput(
            ns("independent_filtering"),
            "Independent filtering",
            value = TRUE
          ),
          
          checkboxInput(
            ns("lfc_shrinkage"),
            "Apply LFC shrinkage (recommended)",
            value = TRUE
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("lfc_shrinkage"), "']"),
            selectInput(
              ns("shrinkage_type"),
              "Shrinkage type:",
              choices = list(
                "apeglm" = "apeglm",
                "ashr" = "ashr", 
                "normal" = "normal"
              ),
              selected = "normal"
            )
          )
        ),
        
        # Batch correction options
        wellPanel(
          h5("🔬 Batch Effect Correction"),
          
          conditionalPanel(
            condition = paste0("output['", ns("batch_detected"), "']"),
            
            div(
              class = "alert alert-warning",
              h6("⚠️ Batch Effects Detected"),
              p("Consider applying batch correction methods.")
            ),
            
            checkboxInput(
              ns("apply_batch_correction"),
              "Apply batch correction",
              value = FALSE
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("apply_batch_correction"), "']"),
              
              selectInput(
                ns("batch_method"),
                "Batch correction method:",
                choices = list(
                  "ComBat-seq (recommended)" = "combat_seq",
                  "Limma removeBatchEffect" = "limma_batch",
                  "Include batch in design" = "design_batch"
                ),
                selected = "combat_seq"
              ),
              
              uiOutput(ns("batch_variable_selection"))
            )
          ),
          
          conditionalPanel(
            condition = paste0("!output['", ns("batch_detected"), "']"),
            p("No batch effects detected", style = "color: green; font-style: italic;")
          )
        )
      ),
      
      column(
        8,
        h4("📈 Analysis Progress & Results"),
        
        # Analysis control buttons
        fluidRow(
          column(
            6,
            conditionalPanel(
              condition = paste0("!output['", ns("analysis_running"), "']"),
              actionButton(
                ns("run_analysis"),
                "🚀 Run Enhanced DESeq2 Analysis",
                class = "btn-primary btn-lg",
                style = "width: 100%; height: 80px; font-size: 18px;"
              )
            ),
            
            conditionalPanel(
              condition = paste0("output['", ns("analysis_running"), "']"),
              div(
                style = "text-align: center; padding: 20px;",
                h5("🔄 Running Analysis..."),
                div(
                  style = "background-color: #3c8dbc; color: white; padding: 10px; border-radius: 5px; margin: 10px 0;",
                  id = ns("analysis_progress_text"),
                  "Initializing analysis pipeline..."
                )
              )
            )
          ),
          
          column(
            6,
            conditionalPanel(
              condition = paste0("output['", ns("has_results"), "']"),
              actionButton(
                ns("export_results"),
                "📥 Export All Results",
                class = "btn-success btn-lg",
                style = "width: 100%; height: 80px; font-size: 16px;"
              )
            )
          )
        ),
        
        br(),
        
        # Results summary
        conditionalPanel(
          condition = paste0("output['", ns("show_results"), "']"),
          
          wellPanel(
            h4("📊 Analysis Results Summary"),
            
            # Multi-comparison results tabs
            conditionalPanel(
              condition = paste0("input['", ns("analysis_type"), "'] == 'multiple_pairwise'"),
              
              tabsetPanel(
                id = ns("results_tabs"),
                
                # Summary tab
                tabPanel(
                  "📈 Overview",
                  br(),
                  fluidRow(
                    column(6, uiOutput(ns("total_comparisons_box"))),
                    column(6, uiOutput(ns("significant_genes_summary_box")))
                  ),
                  br(),
                  h5("Comparison Results:"),
                  DT::dataTableOutput(ns("comparison_summary_table"))
                ),
                
                # Individual comparison results
                tabPanel(
                  "🔍 Detailed Results",
                  br(),
                  fluidRow(
                    column(
                      4,
                      selectInput(
                        ns("selected_comparison"),
                        "Select comparison:",
                        choices = NULL
                      )
                    ),
                    column(
                      4,
                      downloadButton(
                        ns("download_selected_comparison"),
                        "📥 Download Selected",
                        class = "btn-info"
                      )
                    ),
                    column(
                      4,
                      actionButton(
                        ns("create_volcano_plot"),
                        "🌋 Create Volcano Plot",
                        class = "btn-warning"
                      )
                    )
                  ),
                  br(),
                  DT::dataTableOutput(ns("detailed_results_table"))
                )
              )
            ),
            
            # Single comparison results
            conditionalPanel(
              condition = paste0("input['", ns("analysis_type"), "'] != 'multiple_pairwise'"),
              
              fluidRow(
                column(3, valueBoxOutput(ns("total_genes_box"), width = 12)),
                column(3, valueBoxOutput(ns("significant_genes_box"), width = 12)),
                column(3, valueBoxOutput(ns("upregulated_box"), width = 12)),
                column(3, valueBoxOutput(ns("downregulated_box"), width = 12))
              ),
              
              br(),
              
              h5("Top Significant Genes:"),
              DT::dataTableOutput(ns("top_genes_table"))
            )
          )
        )
      )
    ),
    
    # Status messages
    div(id = ns("status_messages"))
  )
}

# Server Function
enhancedDESeq2Analysis <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    dds_object = NULL,
    batch_corrected_counts = NULL,
    analysis_results = list(),
    comparison_results = list(),
    available_comparisons = NULL,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    status_messages = list()
  )
  
  # Update available groups and comparisons when annotation changes
  observe({
    req(values$annotation_data)
    
    conditions <- unique(values$annotation_data$Condition)
    
    if (length(conditions) >= 2) {
      # Generate all possible pairwise comparisons
      comparisons <- combn(conditions, 2, simplify = FALSE)
      comparison_names <- sapply(comparisons, function(x) paste(x[1], "vs", x[2]))
      names(comparisons) <- comparison_names
      
      local_values$available_comparisons <- comparisons
      
      # Update comparison selection UI
      updateSelectInput(session, "selected_comparison",
                       choices = comparison_names,
                       selected = comparison_names[1])
    }
  })
  
  # Main analysis function
  observeEvent(input$run_analysis, {
    req(values$expression_data, values$annotation_data)
    
    local_values$analysis_running <- TRUE
    update_progress("Starting enhanced DESeq2 analysis...")
    
    tryCatch({
      # Step 1: Data preparation and validation
      update_progress("Step 1/6: Preparing and validating data...")
      analysis_data <- prepare_analysis_data(
        expression_data = values$expression_data,
        annotation_data = values$annotation_data,
        min_count = input$min_count_threshold,
        min_samples = input$min_samples_expressed
      )
      
      if (is.null(analysis_data)) {
        stop("Data preparation failed")
      }
      
      # Step 2: Batch effect correction (if enabled)
      if (input$apply_batch_correction && !is.null(values$batch_data)) {
        update_progress("Step 2/6: Applying batch effect correction...")
        analysis_data <- apply_batch_correction(
          count_matrix = analysis_data$counts,
          col_data = analysis_data$col_data,
          batch_info = values$batch_data,
          method = input$batch_method
        )
      } else {
        update_progress("Step 2/6: Skipping batch correction...")
      }
      
      # Step 3: Create DESeq2 object
      update_progress("Step 3/6: Creating DESeq2 dataset...")
      dds <- create_deseq_object(
        count_data = analysis_data$counts,
        col_data = analysis_data$col_data,
        design = ~ Condition
      )
      
      local_values$dds_object <- dds
      
      # Step 4: Run DESeq2 analysis
      update_progress("Step 4/6: Running DESeq2 statistical analysis...")
      dds <- run_deseq_analysis(
        dds = dds,
        fit_type = input$fit_type,
        test = input$test_type
      )
      
      local_values$dds_object <- dds
      
      # Step 5: Extract results based on analysis type
      update_progress("Step 5/6: Extracting and processing results...")
      
      if (input$analysis_type == "multiple_pairwise" && input$run_all_pairwise) {
        # Run all pairwise comparisons
        results <- run_multiple_comparisons(
          dds = dds,
          comparisons = local_values$available_comparisons,
          padj_cutoff = input$padj_cutoff,
          fc_cutoff = input$fc_cutoff,
          lfc_shrinkage = input$lfc_shrinkage,
          shrinkage_type = input$shrinkage_type,
          independent_filtering = input$independent_filtering
        )
        
        local_values$comparison_results <- results
        
      } else {
        # Single comparison or custom contrast
        contrast_info <- get_contrast_info(input$analysis_type, dds, local_values$available_comparisons)
        
        result <- extract_single_result(
          dds = dds,
          contrast = contrast_info$contrast,
          comparison_name = contrast_info$name,
          padj_cutoff = input$padj_cutoff,
          fc_cutoff = input$fc_cutoff,
          lfc_shrinkage = input$lfc_shrinkage,
          shrinkage_type = input$shrinkage_type,
          independent_filtering = input$independent_filtering
        )
        
        local_values$analysis_results <- result
        values$deseq2_results <- result$results_df
      }
      
      # Step 6: Finalize and save results
      update_progress("Step 6/6: Finalizing results...")
      
      # Save filtered expression data
      values$filtered_data <- counts(dds, normalized = TRUE)
      
      local_values$analysis_complete <- TRUE
      add_status_message("✅ Enhanced DESeq2 analysis completed successfully!", "success")
      
    }, error = function(e) {
      add_status_message(paste("❌ Analysis failed:", e$message), "danger")
      logger::log_error("DESeq2 analysis error: {e$message}")
    }, finally = {
      local_values$analysis_running <- FALSE
    })
  })
  
  # Data preparation function
  prepare_analysis_data <- function(expression_data, annotation_data, min_count = 10, min_samples = 3) {
    tryCatch({
      # Prepare count matrix
      count_matrix <- as.matrix(expression_data)
      count_matrix <- round(count_matrix)
      count_matrix[count_matrix < 0] <- 0
      
      # Prepare sample data
      sample_data <- annotation_data
      rownames(sample_data) <- sample_data$Sample
      sample_data <- sample_data[colnames(count_matrix), , drop = FALSE]
      
      # Ensure factors are properly set
      sample_data$Condition <- factor(sample_data$Condition)
      
      # Filter low-count genes (Emory methodology)
      # Remove genes with total counts < min_count
      keep1 <- rowSums(count_matrix) >= min_count
      
      # Remove genes not expressed in at least min_samples samples
      keep2 <- rowSums(count_matrix > 0) >= min_samples
      
      keep <- keep1 & keep2
      
      filtered_matrix <- count_matrix[keep, ]
      
      add_status_message(
        paste("Filtered", sum(!keep), "low-count genes.", nrow(filtered_matrix), "genes retained for analysis."),
        "info"
      )
      
      return(list(
        counts = filtered_matrix,
        col_data = sample_data,
        original_counts = count_matrix
      ))
      
    }, error = function(e) {
      add_status_message(paste("Data preparation error:", e$message), "danger")
      return(NULL)
    })
  }
  
  # Batch correction function
  apply_batch_correction <- function(count_matrix, col_data, batch_info, method = "combat_seq") {
    
    if (!requireNamespace("sva", quietly = TRUE)) {
      add_status_message("sva package not available for batch correction", "warning")
      return(list(counts = count_matrix, col_data = col_data))
    }
    
    tryCatch({
      # Extract batch information
      batch_vector <- extract_batch_vector(colnames(count_matrix), batch_info)
      
      if (is.null(batch_vector)) {
        add_status_message("Could not extract batch information", "warning")
        return(list(counts = count_matrix, col_data = col_data))
      }
      
      # Apply batch correction based on method
      if (method == "combat_seq") {
        # ComBat-seq method (recommended for RNA-seq)
        corrected_counts <- sva::ComBat_seq(
          counts = count_matrix,
          batch = batch_vector,
          group = col_data$Condition
        )
        
        add_status_message("Applied ComBat-seq batch correction", "success")
        
      } else if (method == "limma_batch") {
        # Limma removeBatchEffect (on log-transformed data)
        if (!requireNamespace("limma", quietly = TRUE)) {
          stop("limma package not available")
        }
        
        log_counts <- log2(count_matrix + 1)
        design <- model.matrix(~ col_data$Condition)
        
        corrected_log <- limma::removeBatchEffect(
          log_counts,
          batch = batch_vector,
          design = design
        )
        
        # Convert back to counts (approximate)
        corrected_counts <- round(2^corrected_log - 1)
        corrected_counts[corrected_counts < 0] <- 0
        
        add_status_message("Applied limma batch correction", "success")
        
      } else if (method == "design_batch") {
        # Include batch in design matrix (handled later in DESeq2)
        col_data$Batch <- batch_vector
        corrected_counts <- count_matrix
        
        add_status_message("Batch will be included in statistical model", "info")
      }
      
      return(list(
        counts = corrected_counts,
        col_data = col_data
      ))
      
    }, error = function(e) {
      add_status_message(paste("Batch correction failed:", e$message), "warning")
      return(list(counts = count_matrix, col_data = col_data))
    })
  }
  
  # Create DESeq2 object
  create_deseq_object <- function(count_data, col_data, design) {
    tryCatch({
      # Include batch in design if present
      if ("Batch" %in% colnames(col_data)) {
        design <- ~ Batch + Condition
        col_data$Batch <- factor(col_data$Batch)
      }
      
      dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        design = design
      )
      
      return(dds)
      
    }, error = function(e) {
      stop(paste("Failed to create DESeq2 object:", e$message))
    })
  }
  
  # Run DESeq2 analysis
  run_deseq_analysis <- function(dds, fit_type = "parametric", test = "Wald") {
    tryCatch({
      dds <- DESeq(
        dds,
        fitType = fit_type,
        test = test,
        quiet = FALSE
      )
      
      return(dds)
      
    }, error = function(e) {
      stop(paste("DESeq2 analysis failed:", e$message))
    })
  }
  
  # Run multiple pairwise comparisons (Emory style)
  run_multiple_comparisons <- function(dds, comparisons, padj_cutoff, fc_cutoff, 
                                     lfc_shrinkage, shrinkage_type, independent_filtering) {
    
    results_list <- list()
    
    for (i in seq_along(comparisons)) {
      comparison_name <- names(comparisons)[i]
      contrast <- c("Condition", comparisons[[i]][1], comparisons[[i]][2])
      
      update_progress(paste("Processing comparison", i, "of", length(comparisons), ":", comparison_name))
      
      # Extract results
      res <- results(
        dds,
        contrast = contrast,
        independentFiltering = independent_filtering,
        alpha = padj_cutoff
      )
      
      # Apply LFC shrinkage if requested
      if (lfc_shrinkage) {
        tryCatch({
          res_lfc <- lfcShrink(
            dds,
            contrast = contrast,
            res = res,
            type = shrinkage_type
          )
          res <- res_lfc
        }, error = function(e) {
          add_status_message(paste("LFC shrinkage failed for", comparison_name, ":", e$message), "warning")
        })
      }
      
      # Convert to data frame and add metadata
      results_df <- process_results_dataframe(res, padj_cutoff, fc_cutoff)
      
      # Add normalized counts
      normalized_counts <- counts(dds, normalized = TRUE)
      results_with_counts <- merge(
        as.data.frame(results_df),
        as.data.frame(normalized_counts),
        by = 'row.names',
        sort = FALSE
      )
      names(results_with_counts)[1] <- 'gene'
      
      results_list[[comparison_name]] <- list(
        results_object = res,
        results_df = results_df,
        results_with_counts = results_with_counts,
        comparison = comparisons[[i]],
        summary = list(
          total_genes = nrow(results_df),
          significant_genes = sum(results_df$significant, na.rm = TRUE),
          upregulated = sum(results_df$regulation == "Up", na.rm = TRUE),
          downregulated = sum(results_df$regulation == "Down", na.rm = TRUE)
        )
      )
    }
    
    return(results_list)
  }
  
  # Extract single comparison result
  extract_single_result <- function(dds, contrast, comparison_name, padj_cutoff, fc_cutoff,
                                   lfc_shrinkage, shrinkage_type, independent_filtering) {
    
    # Extract results
    res <- results(
      dds,
      contrast = contrast,
      independentFiltering = independent_filtering,
      alpha = padj_cutoff
    )
    
    # Apply LFC shrinkage if requested
    if (lfc_shrinkage) {
      tryCatch({
        res_lfc <- lfcShrink(
          dds,
          contrast = contrast,
          res = res,
          type = shrinkage_type
        )
        res <- res_lfc
      }, error = function(e) {
        add_status_message(paste("LFC shrinkage failed:", e$message), "warning")
      })
    }
    
    # Process results
    results_df <- process_results_dataframe(res, padj_cutoff, fc_cutoff)
    
    return(list(
      results_object = res,
      results_df = results_df,
      comparison_name = comparison_name
    ))
  }
  
  # Process results dataframe (common formatting)
  process_results_dataframe <- function(res, padj_cutoff, fc_cutoff) {
    
    results_df <- as.data.frame(res)
    results_df$gene <- rownames(results_df)
    
    # Add significance flags
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
    
    return(results_df)
  }
  
  # Helper functions
  extract_batch_vector <- function(sample_names, batch_info) {
    # Implementation would extract batch information from sample names
    # Based on detected patterns in batch_info
    return(NULL) # Placeholder
  }
  
  get_contrast_info <- function(analysis_type, dds, available_comparisons) {
    # Implementation would return appropriate contrast based on analysis type
    # For now, return first available comparison
    if (length(available_comparisons) > 0) {
      return(list(
        contrast = c("Condition", available_comparisons[[1]][1], available_comparisons[[1]][2]),
        name = names(available_comparisons)[1]
      ))
    }
    return(NULL)
  }
  
  update_progress <- function(message) {
    session$sendCustomMessage("updateProgress", list(
      id = ns("analysis_progress_text"),
      message = message
    ))
  }
  
  # Helper function to add status messages
  add_status_message <- function(message, type = "info") {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    
    alert_class <- switch(
      type,
      "success" = "alert-success",
      "info" = "alert-info",
      "warning" = "alert-warning",
      "danger" = "alert-danger",
      "alert-info"
    )
    
    new_message <- list(
      message = message,
      type = alert_class,
      timestamp = timestamp
    )
    
    local_values$status_messages <- append(local_values$status_messages, list(new_message))
    
    if (length(local_values$status_messages) > 5) {
      local_values$status_messages <- tail(local_values$status_messages, 5)
    }
  }
  
  # UI Outputs
  
  # Analysis type selection
  output$analysis_type_selection <- renderUI({
    req(values$annotation_data)
    
    n_conditions <- length(unique(values$annotation_data$Condition))
    
    choices <- list(
      "🔄 Multiple pairwise comparisons" = "multiple_pairwise"
    )
    
    if (n_conditions == 2) {
      choices[["🎯 Single comparison"]] <- "single_comparison"
    }
    
    choices[["⚙️ Custom contrast"]] <- "custom_contrast"
    
    selectInput(
      ns("analysis_type"),
      "Analysis type:",
      choices = choices,
      selected = "multiple_pairwise"
    )
  })
  
  # Pairwise selection UI
  output$pairwise_selection <- renderUI({
    req(local_values$available_comparisons)
    
    comparison_names <- names(local_values$available_comparisons)
    
    tagList(
      p(paste("Available comparisons:", length(comparison_names))),
      tags$ul(
        lapply(comparison_names, function(name) {
          tags$li(name)
        })
      )
    )
  })
  
  # Batch variable selection
  output$batch_variable_selection <- renderUI({
    req(values$batch_data)
    
    if (length(values$batch_data) > 0) {
      selectInput(
        ns("batch_variable"),
        "Batch variable:",
        choices = names(values$batch_data),
        selected = names(values$batch_data)[1]
      )
    }
  })
  
  # Comparison summary table
  output$comparison_summary_table <- DT::renderDataTable({
    req(local_values$comparison_results)
    
    summary_data <- data.frame(
      Comparison = names(local_values$comparison_results),
      Total_Genes = sapply(local_values$comparison_results, function(x) x$summary$total_genes),
      Significant = sapply(local_values$comparison_results, function(x) x$summary$significant_genes),
      Upregulated = sapply(local_values$comparison_results, function(x) x$summary$upregulated),
      Downregulated = sapply(local_values$comparison_results, function(x) x$summary$downregulated),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      summary_data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'tip'
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleInterval(c(100, 500), c("white", "#fff2cc", "#d4edda"))
      )
  })
  
  # Detailed results table
  output$detailed_results_table <- DT::renderDataTable({
    req(input$selected_comparison, local_values$comparison_results)
    
    selected_result <- local_values$comparison_results[[input$selected_comparison]]
    req(selected_result)
    
    # Show top 50 most significant genes
    display_data <- head(selected_result$results_df, 50)
    
    # Format for display
    display_data <- display_data %>%
      dplyr::select(gene, log2FoldChange, pvalue, padj, regulation) %>%
      dplyr::mutate(
        log2FoldChange = round(log2FoldChange, 3),
        pvalue = formatC(pvalue, format = "e", digits = 2),
        padj = formatC(padj, format = "e", digits = 2)
      )
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'tip'
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "regulation",
        backgroundColor = DT::styleEqual(
          c("Up", "Down", "NS"),
          c("#ffcccc", "#ccccff", "#f0f0f0")
        )
      )
  })
  
  # Reactive outputs for UI control
  output$analysis_running <- reactive({
    local_values$analysis_running
  })
  outputOptions(output, "analysis_running", suspendWhenHidden = FALSE)
  
  output$has_results <- reactive({
    local_values$analysis_complete && 
    (length(local_values$comparison_results) > 0 || !is.null(local_values$analysis_results))
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  output$show_results <- reactive({
    local_values$analysis_complete
  })
  outputOptions(output, "show_results", suspendWhenHidden = FALSE)
  
  output$batch_detected <- reactive({
    !is.null(values$batch_data) && length(values$batch_data) > 0
  })
  outputOptions(output, "batch_detected", suspendWhenHidden = FALSE)
  
  # Value boxes for summary statistics
  output$total_comparisons_box <- renderValueBox({
    n_comparisons <- length(local_values$comparison_results)
    
    valueBox(
      value = n_comparisons,
      subtitle = "Comparisons Analyzed",
      icon = icon("chart-line"),
      color = "blue",
      width = 12
    )
  })
  
  output$significant_genes_summary_box <- renderValueBox({
    total_sig <- if (length(local_values$comparison_results) > 0) {
      sum(sapply(local_values$comparison_results, function(x) x$summary$significant_genes))
    } else 0
    
    valueBox(
      value = total_sig,
      subtitle = "Total Significant Genes",
      icon = icon("star"),
      color = "yellow",
      width = 12
    )
  })
  
  # Return analysis results
  return(reactive({ 
    if (length(local_values$comparison_results) > 0) {
      local_values$comparison_results
    } else {
      local_values$analysis_results
    }
  }))
}