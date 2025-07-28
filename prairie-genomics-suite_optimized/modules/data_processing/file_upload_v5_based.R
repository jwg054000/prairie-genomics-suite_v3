# Data Upload Module - Based on Working v5 Implementation
# Simplified, focused approach without built-in gene conversion
# 
# Author: Prairie Genomics Team - v5 Working Approach
# Features: Proper duplicate handling, clean data loading, no gene conversion

# Load dependencies
source("config/app_config.R")
source("utils/memory_manager.R")

# Load required packages
required_packages <- c("dplyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg)
      library(pkg, character.only = TRUE)
      cat("✅", pkg, "installed and loaded\n")
    }, error = function(e) {
      cat("⚠️", pkg, "not available - some features may be limited\n")
    })
  } else {
    library(pkg, character.only = TRUE)
    cat("✅", pkg, "loaded\n")
  }
}

# Memory monitoring function
monitor_memory <- function() {
  tryCatch({
    gc_info <- gc()
    memory_mb <- round(sum(gc_info[, 2]), 1)
    
    list(
      used_mb = memory_mb,
      warning = memory_mb > MEMORY_WARNING_MB,
      critical = memory_mb > MEMORY_CRITICAL_MB,
      message = if (memory_mb > MEMORY_CRITICAL_MB) {
        "Critical memory usage! Consider restarting the session."
      } else if (memory_mb > MEMORY_WARNING_MB) {
        "High memory usage detected. Performance may be affected."
      } else {
        paste("Memory usage:", memory_mb, "MB")
      }
    )
  }, error = function(e) {
    list(used_mb = 0, warning = FALSE, critical = FALSE, message = "Memory monitoring unavailable")
  })
}

# Smart file validation
validate_expression_data <- function(data) {
  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    info = list()
  )
  
  tryCatch({
    # Basic structure checks
    if (!is.data.frame(data) && !is.matrix(data)) {
      validation_result$valid <- FALSE
      validation_result$errors <- c(validation_result$errors, "Data must be a data frame or matrix")
      return(validation_result)
    }
    
    # Dimensions check
    n_rows <- nrow(data)
    n_cols <- ncol(data)
    
    if (n_rows < 100) {
      validation_result$warnings <- c(validation_result$warnings, 
                                    paste("Only", n_rows, "genes detected. Consider if this is expected."))
    }
    
    if (n_cols < get_config("analysis", "deseq2")$min_samples) {
      validation_result$valid <- FALSE
      validation_result$errors <- c(validation_result$errors, 
                                  paste("Need at least", get_config("analysis", "deseq2")$min_samples, "samples"))
      return(validation_result)
    }
    
    # Check for gene names
    if (is.character(data[[1]]) || is.factor(data[[1]])) {
      validation_result$info$has_gene_names <- TRUE
      validation_result$info$gene_column <- 1
      numeric_cols <- 2:n_cols
    } else {
      validation_result$info$has_gene_names <- FALSE
      validation_result$info$gene_column <- NULL
      numeric_cols <- 1:n_cols
    }
    
    # Check if data is numeric
    numeric_data <- data[, numeric_cols, drop = FALSE]
    non_numeric_cols <- !sapply(numeric_data, function(x) is.numeric(x) || all(is.na(x)))
    
    if (any(non_numeric_cols)) {
      problematic_cols <- names(numeric_data)[non_numeric_cols]
      validation_result$valid <- FALSE
      validation_result$errors <- c(validation_result$errors, 
                                  paste("Non-numeric data found in columns:", paste(problematic_cols, collapse = ", ")))
      return(validation_result)
    }
    
    # Check for negative values
    if (any(numeric_data < 0, na.rm = TRUE)) {
      validation_result$warnings <- c(validation_result$warnings, 
                                    "Negative values detected. RNA-seq data should be non-negative.")
    }
    
    # Check for excessive zeros
    zero_rate <- sum(numeric_data == 0, na.rm = TRUE) / (nrow(numeric_data) * ncol(numeric_data))
    if (zero_rate > 0.8) {
      validation_result$warnings <- c(validation_result$warnings, 
                                    paste("High zero rate detected:", round(zero_rate * 100, 1), "%"))
    }
    
    # Memory estimation
    estimated_memory_mb <- object.size(data) / (1024^2)
    validation_result$info$estimated_memory_mb <- round(estimated_memory_mb, 1)
    
    if (estimated_memory_mb > get_config("memory", "memory_warning_mb")) {
      validation_result$warnings <- c(validation_result$warnings, 
                                    "Large dataset detected. Processing may be slow.")
    }
    
    # Set processing info
    validation_result$info$dimensions <- c(n_rows, n_cols)
    validation_result$info$processing_method <- if (n_rows > CHUNK_SIZE) "chunked" else "standard"
    
  }, error = function(e) {
    validation_result$valid <- FALSE
    validation_result$errors <- c(validation_result$errors, paste("Validation error:", e$message))
  })
  
  return(validation_result)
}

# Function to handle duplicate gene IDs properly (from v5 working implementation)
handle_duplicate_genes_v5 <- function(raw_data, gene_names_column = 1) {
  cat("🔍 Checking for duplicate gene IDs (v5 method)...\n")
  
  # Extract gene names
  if (is.character(gene_names_column)) {
    gene_names <- raw_data[[gene_names_column]]
  } else {
    gene_names <- raw_data[[1]]  # First column
  }
  
  # Remove gene names column from data
  if (is.character(gene_names_column)) {
    data_only <- raw_data[, !names(raw_data) %in% gene_names_column, drop = FALSE]
  } else {
    data_only <- raw_data[, -1, drop = FALSE]
  }
  
  original_count <- length(gene_names)
  cat("📊 Original genes:", original_count, "\n")
  
  # Check for duplicates
  duplicated_genes <- duplicated(gene_names)
  n_duplicates <- sum(duplicated_genes)
  
  if (n_duplicates > 0) {
    cat("⚠️ Found", n_duplicates, "duplicate gene IDs\n")
    
    # Strategy: Aggregate duplicates by summing counts (v5 method)
    cat("🔧 Aggregating duplicate genes by summing counts...\n")
    
    # Combine gene names with data
    combined_data <- cbind(gene_id = gene_names, data_only)
    
    # Convert expression columns to numeric
    expr_cols <- names(data_only)
    for (col in expr_cols) {
      combined_data[[col]] <- as.numeric(as.character(combined_data[[col]]))
    }
    
    # Aggregate by gene_id (sum duplicates)
    aggregated <- combined_data %>%
      group_by(gene_id) %>%
      summarise(across(all_of(expr_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    
    # Extract gene names and data
    final_gene_names <- aggregated$gene_id
    final_data <- aggregated[, expr_cols, drop = FALSE]
    
    cat("✅ After aggregation:", length(final_gene_names), "unique genes\n")
    cat("📉 Removed", original_count - length(final_gene_names), "duplicate entries\n")
    
  } else {
    cat("✅ No duplicate gene IDs found\n")
    final_gene_names <- gene_names
    final_data <- data_only
  }
  
  # Set rownames
  final_data <- as.data.frame(final_data)
  rownames(final_data) <- final_gene_names
  
  # Additional cleaning - remove genes with all zeros
  gene_sums <- rowSums(final_data, na.rm = TRUE)
  non_zero_genes <- gene_sums > 0
  
  if (sum(!non_zero_genes) > 0) {
    cat("🧹 Removing", sum(!non_zero_genes), "genes with zero counts\n")
    final_data <- final_data[non_zero_genes, , drop = FALSE]
  }
  
  cat("🎯 Final dataset:", nrow(final_data), "genes x", ncol(final_data), "samples\n")
  
  return(list(
    data = final_data,
    original_genes = original_count,
    final_genes = nrow(final_data),
    duplicates_removed = n_duplicates,
    zero_genes_removed = sum(!non_zero_genes)
  ))
}

# Optimized data processing with v5 duplicate handling
process_expression_data_v5 <- function(raw_data, chunk_size = CHUNK_SIZE) {
  tryCatch({
    cat("🔄 Starting v5-based data processing...\n")
    
    # Validate first
    validation <- validate_expression_data(raw_data)
    if (!validation$valid) {
      return(list(
        success = FALSE,
        error = paste("Validation failed:", paste(validation$errors, collapse = "; ")),
        data = NULL
      ))
    }
    
    # Monitor memory before processing
    mem_before <- monitor_memory()
    cat("📊 Memory before processing:", mem_before$message, "\n")
    
    # Extract dimensions
    n_rows <- nrow(raw_data)
    n_cols <- ncol(raw_data)
    
    # Handle gene names and duplicates using v5 method
    if (validation$info$has_gene_names) {
      cat("📋 Processing with gene names using v5 duplicate handling...\n")
      duplicate_result <- handle_duplicate_genes_v5(raw_data, 1)
      result_matrix <- duplicate_result$data
      duplicate_stats <- duplicate_result
    } else {
      cat("📋 Processing without gene names...\n")
      gene_names <- paste0("Gene_", seq_len(n_rows))
      expression_data <- raw_data
      result_matrix <- as.matrix(expression_data)
      rownames(result_matrix) <- gene_names
      
      duplicate_stats <- list(
        original_genes = n_rows,
        final_genes = n_rows,
        duplicates_removed = 0,
        zero_genes_removed = 0
      )
    }
    
    # Ensure numeric data
    if (!is.numeric(result_matrix)) {
      tryCatch({
        result_matrix <- apply(result_matrix, c(1,2), as.numeric)
      }, error = function(e) {
        cat("⚠️ Warning: Some data may not be numeric\n")
      })
    }
    
    # Final validation
    if (is.null(result_matrix) || nrow(result_matrix) == 0) {
      return(list(
        success = FALSE,
        error = "No valid genes remaining after processing",
        data = NULL
      ))
    }
    
    # Memory cleanup
    rm(raw_data)
    gc()
    
    # Monitor memory after processing
    mem_after <- monitor_memory()
    cat("📊 Memory after processing:", mem_after$message, "\n")
    
    # Calculate processing stats
    genes_retained <- nrow(result_matrix)
    retention_rate <- round((genes_retained / n_rows) * 100, 1)
    
    cat("✅ v5-based processing completed successfully:\n")
    cat("   - Original genes:", formatC(n_rows, format="d", big.mark=","), "\n")
    if (duplicate_stats$duplicates_removed > 0) {
      cat("   - Duplicates aggregated:", duplicate_stats$duplicates_removed, "\n")
    }
    if (duplicate_stats$zero_genes_removed > 0) {
      cat("   - Zero-count genes removed:", duplicate_stats$zero_genes_removed, "\n")
    }
    cat("   - Final genes retained:", formatC(genes_retained, format="d", big.mark=","), 
        paste0("(", retention_rate, "%)"), "\n")
    cat("   - Samples:", ncol(result_matrix), "\n")
    
    return(list(
      success = TRUE,
      data = result_matrix,
      stats = list(
        genes_input = n_rows,
        genes_retained = genes_retained,
        retention_rate = retention_rate,
        samples = ncol(result_matrix),
        processing_method = validation$info$processing_method
      ),
      duplicate_stats = duplicate_stats,
      warnings = validation$warnings
    ))
    
  }, error = function(e) {
    cat("❌ Error in process_expression_data_v5:", e$message, "\n")
    return(list(
      success = FALSE,
      error = paste("Data processing failed:", e$message),
      data = NULL
    ))
  })
}

# Data Upload UI Module (simplified, no gene conversion options)
dataUploadUI_v5 <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        fileInput(
          ns("file"),
          "Choose Expression Data File",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          width = "100%"
        ),
        
        conditionalPanel(
          condition = paste0("output['", ns("file_uploaded"), "']"),
          
          radioButtons(
            ns("sep"),
            "Separator:",
            choices = list("Comma (,)" = ",", "Tab" = "\t", "Semicolon (;)" = ";"),
            selected = ",",
            inline = TRUE
          ),
          
          radioButtons(
            ns("header"),
            "Header:",
            choices = list("Yes" = TRUE, "No" = FALSE),
            selected = TRUE,
            inline = TRUE
          ),
          
          br(),
          
          div(
            class = "alert alert-info",
            style = "padding: 10px; margin-bottom: 15px;",
            h6("ℹ️ Note: Gene ID Conversion", style = "margin-top: 0;"),
            p("Gene conversion (Ensembl → Symbols) is handled separately in the analysis modules.", style = "margin: 0; font-size: 12px;"),
            p("This module focuses on data loading and duplicate gene handling.", style = "margin: 0; font-size: 12px;")
          ),
          
          actionButton(
            ns("process_file"),
            "🚀 Process Data",
            class = "btn-primary btn-lg",
            style = "width: 100%;"
          )
        )
      ),
      
      column(
        6,
        h5("📋 File Requirements"),
        tags$ul(
          tags$li("CSV, TSV, or Excel format"),
          tags$li("Genes as rows, samples as columns"),
          tags$li("First column: Gene IDs (optional)"),
          tags$li("Numeric expression values"),
          tags$li(paste("Maximum size:", get_config("memory", "max_file_size_mb"), "MB"))
        ),
        
        # Memory status
        div(
          id = ns("memory_status"),
          style = "margin-top: 15px; padding: 10px; border-radius: 5px; background-color: #f8f9fa;",
          h6("💾 Memory Status"),
          uiOutput(ns("memory_info"))
        )
      )
    ),
    
    # File preview
    conditionalPanel(
      condition = paste0("output['", ns("show_preview"), "']"),
      
      br(),
      
      fluidRow(
        box(
          title = "📄 File Preview", 
          status = "info", 
          solidHeader = TRUE, 
          width = 12,
          collapsible = TRUE,
          
          DT::dataTableOutput(ns("file_preview"))
        )
      )
    ),
    
    # Processing status
    conditionalPanel(
      condition = paste0("output['", ns("show_processing"), "']"),
      
      br(),
      
      div(
        class = "alert alert-info",
        h5("🔄 Processing Data..."),
        progressBar(
          id = ns("processing_progress"),
          value = 0,
          total = 100,
          status = "primary",
          display_pct = TRUE
        ),
        uiOutput(ns("processing_status"))
      )
    ),
    
    # Results summary
    conditionalPanel(
      condition = paste0("output['", ns("show_results"), "']"),
      
      br(),
      
      div(
        class = "alert alert-success",
        h5("✅ Data Processing Complete"),
        uiOutput(ns("processing_summary"))
      ),
      
      # Processed data preview
      conditionalPanel(
        condition = paste0("output['", ns("show_processed_preview"), "']"),
        
        br(),
        
        fluidRow(
          box(
            title = "📊 Processed Data Preview", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            collapsible = TRUE,
            
            p("This shows your processed data after duplicate handling and cleaning. Gene conversion (if needed) will be handled in the analysis modules."),
            
            DT::dataTableOutput(ns("processed_preview"))
          )
        )
      )
    )
  )
}

# Data Upload Server Module (simplified, no gene conversion)
dataUploadServer_v5 <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for this module
    local_values <- reactiveValues(
      raw_data = NULL,
      processed_data = NULL,
      processing = FALSE,
      validation_result = NULL,
      duplicate_stats = NULL
    )
    
    # Monitor file upload
    output$file_uploaded <- reactive({
      !is.null(input$file)
    })
    outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)
    
    # Memory monitoring
    output$memory_info <- renderUI({
      memory_status <- monitor_memory()
      
      status_class <- if (memory_status$critical) {
        "danger"
      } else if (memory_status$warning) {
        "warning"
      } else {
        "success"
      }
      
      div(
        class = paste("alert alert-", status_class),
        p(memory_status$message, style = "margin: 0; font-size: 12px;")
      )
    })
    
    # Read uploaded file
    observe({
      req(input$file)
      
      tryCatch({
        cat("📁 Reading uploaded file:", input$file$name, "\n")
        
        # Check file size
        file_size_mb <- file.size(input$file$datapath) / (1024^2)
        if (file_size_mb > get_config("memory", "max_file_size_mb")) {
          showNotification(
            paste("File too large:", round(file_size_mb, 1), "MB. Maximum allowed:", 
                  get_config("memory", "max_file_size_mb"), "MB"),
            type = "error"
          )
          return()
        }
        
        # Determine file type and separator
        file_ext <- tolower(tools::file_ext(input$file$datapath))
        separator <- input$sep %||% ","
        has_header <- as.logical(input$header %||% TRUE)
        
        # Read file with robust error handling
        if (file_ext %in% c("csv", "txt", "tsv", "")) {
          local_values$raw_data <- tryCatch({
            if (requireNamespace("data.table", quietly = TRUE)) {
              dt_result <- data.table::fread(
                input$file$datapath,
                sep = separator,
                header = has_header,
                stringsAsFactors = FALSE,
                data.table = FALSE,
                showProgress = FALSE
              )
              dt_result
            } else {
              read.csv(
                input$file$datapath,
                sep = separator,
                header = has_header,
                stringsAsFactors = FALSE,
                check.names = FALSE,
                row.names = NULL
              )
            }
          }, error = function(e) {
            read.table(
              input$file$datapath,
              sep = separator,
              header = has_header,
              stringsAsFactors = FALSE,
              check.names = FALSE,
              fill = TRUE,
              quote = "\""
            )
          })
        } else if (file_ext %in% c("xlsx", "xls")) {
          if (requireNamespace("readxl", quietly = TRUE)) {
            local_values$raw_data <- readxl::read_excel(
              input$file$datapath,
              col_names = has_header,
              .name_repair = "minimal"
            )
            local_values$raw_data <- as.data.frame(local_values$raw_data, stringsAsFactors = FALSE)
          } else {
            showNotification("readxl package required for Excel files", type = "error")
            return()
          }
        }
        
        # Check if data was read successfully
        if (is.null(local_values$raw_data) || nrow(local_values$raw_data) == 0) {
          showNotification("File appears to be empty or could not be read", type = "error")
          return()
        }
        
        cat("✅ File read successfully:", nrow(local_values$raw_data), "rows x", 
            ncol(local_values$raw_data), "columns\n")
        
        # Validate the data
        local_values$validation_result <- validate_expression_data(local_values$raw_data)
        
        if (!local_values$validation_result$valid) {
          showNotification(
            paste("❌ Validation failed:", paste(local_values$validation_result$errors, collapse = "; ")),
            type = "error",
            duration = 10
          )
        } else if (length(local_values$validation_result$warnings) > 0) {
          showNotification(
            paste("⚠️ Warnings:", paste(local_values$validation_result$warnings, collapse = "; ")),
            type = "warning",
            duration = 8
          )
        } else {
          showNotification("✅ File loaded and validated successfully!", type = "message", duration = 3)
        }
        
      }, error = function(e) {
        cat("❌ Error reading file:", e$message, "\n")
        showNotification(paste("❌ File reading error:", e$message), type = "error", duration = 15)
        local_values$raw_data <- NULL
      })
    })
    
    # File preview
    output$show_preview <- reactive({
      !is.null(local_values$raw_data) && !local_values$processing
    })
    outputOptions(output, "show_preview", suspendWhenHidden = FALSE)
    
    output$file_preview <- DT::renderDataTable({
      req(local_values$raw_data)
      
      # Show first 10 rows and columns for preview
      preview_data <- local_values$raw_data[
        1:min(10, nrow(local_values$raw_data)), 
        1:min(10, ncol(local_values$raw_data))
      ]
      
      DT::datatable(
        preview_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        )
      )
    })
    
    # Process file button
    observeEvent(input$process_file, {
      req(local_values$raw_data)
      
      if (!local_values$validation_result$valid) {
        showNotification("Cannot process invalid data. Please fix errors first.", type = "error")
        return()
      }
      
      cat("🔄 Starting v5-based data processing...\n")
      local_values$processing <- TRUE
      
      # Update progress
      update_progress <- function(value) {
        if (requireNamespace("shinyWidgets", quietly = TRUE)) {
          shinyWidgets::updateProgressBar(session, "processing_progress", value = value)
        }
      }
      
      update_progress(20)
      
      # Process data with v5 method
      processing_result <- monitor_operation("data_processing", function() {
        process_expression_data_v5(local_values$raw_data)
      })
      
      update_progress(80)
      
      if (processing_result$result$success) {
        processed_matrix <- processing_result$result$data
        
        # Store duplicate stats from processing
        local_values$duplicate_stats <- processing_result$result$duplicate_stats
        
        # Store final processed data (no gene conversion here)
        local_values$processed_data <- processed_matrix
        values$expression_data <- processed_matrix
        
        update_progress(100)
        
        # Create success message
        success_message <- paste("✅ Data processed successfully!", 
                                formatC(nrow(processed_matrix), format="d", big.mark=","),
                                "genes retained")
        
        # Add duplicate handling info
        if (!is.null(local_values$duplicate_stats) && local_values$duplicate_stats$duplicates_removed > 0) {
          duplicate_info <- paste0("(", local_values$duplicate_stats$duplicates_removed, 
                                  " duplicate gene IDs aggregated)")
          success_message <- paste(success_message, duplicate_info)
        }
        
        showNotification(
          success_message,
          type = "message",
          duration = 8
        )
        
        cat("✅ v5-based data processing completed and stored in main values\n")
      } else {
        showNotification(
          paste("❌ Processing failed:", processing_result$result$error),
          type = "error",
          duration = 10
        )
      }
      
      local_values$processing <- FALSE
    })
    
    # Processing status
    output$show_processing <- reactive({
      local_values$processing
    })
    outputOptions(output, "show_processing", suspendWhenHidden = FALSE)
    
    output$processing_status <- renderUI({
      if (local_values$processing) {
        div(
          p("Processing uploaded data with v5 duplicate handling..."),
          p("Gene conversion (if needed) will be handled in analysis modules.")
        )
      }
    })
    
    # Results summary
    output$show_results <- reactive({
      !is.null(local_values$processed_data) && !local_values$processing
    })
    outputOptions(output, "show_results", suspendWhenHidden = FALSE)
    
    output$processing_summary <- renderUI({
      req(local_values$processed_data)
      
      div(
        h6("📊 Processing Results"),
        p(paste("✅ Genes retained:", 
               formatC(nrow(local_values$processed_data), format="d", big.mark=","))),
        p(paste("✅ Samples:", ncol(local_values$processed_data))),
        p(paste("✅ Data size:", 
               round(object.size(local_values$processed_data) / (1024^2), 1), "MB")),
        
        if (!is.null(local_values$duplicate_stats) && local_values$duplicate_stats$duplicates_removed > 0) {
          div(
            h6("🔧 Duplicate Handling:"),
            p(paste("- Duplicate gene IDs found:", local_values$duplicate_stats$duplicates_removed)),
            p(paste("- Aggregated by summing counts"))
          )
        },
        
        if (!is.null(local_values$validation_result$warnings) && 
            length(local_values$validation_result$warnings) > 0) {
          div(
            h6("⚠️ Warnings:"),
            tags$ul(
              lapply(local_values$validation_result$warnings, function(w) tags$li(w))
            )
          )
        }
      )
    })
    
    # Show processed data preview
    output$show_processed_preview <- reactive({
      !is.null(local_values$processed_data) && !local_values$processing
    })
    outputOptions(output, "show_processed_preview", suspendWhenHidden = FALSE)
    
    # Processed data preview
    output$processed_preview <- DT::renderDataTable({
      req(local_values$processed_data)
      
      # Create preview data frame
      preview_data <- as.data.frame(local_values$processed_data)
      preview_data <- preview_data[1:min(15, nrow(preview_data)), 1:min(8, ncol(preview_data))]
      
      # Add Gene_ID column as first column
      preview_data <- data.frame(
        Gene_ID = rownames(preview_data),
        preview_data,
        stringsAsFactors = FALSE
      )
      
      DT::datatable(
        preview_data,
        options = list(
          scrollX = TRUE,
          pageLength = 15,
          dom = 'tip',
          columnDefs = list(
            list(width = '150px', targets = 0),
            list(className = 'dt-center', targets = 1:(ncol(preview_data)-1))
          )
        ),
        class = "display",
        rownames = FALSE
      ) %>%
      DT::formatStyle(
        'Gene_ID',
        backgroundColor = '#f8f9fa',
        fontWeight = 'bold'
      )
    })
    
    # Return reactive for processed data
    return(reactive({
      local_values$processed_data
    }))
  })
}

cat("✅ v5-based data upload module loaded\n")
cat("ℹ️ This module focuses on data loading and duplicate handling\n")
cat("ℹ️ Gene ID conversion is handled separately in analysis modules\n")