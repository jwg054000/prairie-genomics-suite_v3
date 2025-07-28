# Data Upload Module for Prairie Genomics Suite - v5 Architecture
# Based on working v5 implementation with proper separation of concerns
# 
# Author: Prairie Genomics Team
# Architecture: Clean data loading with duplicate handling, gene conversion handled separately
# Performance: Memory-optimized processing aligned with v5 approach

# Load configuration
source("config/app_config.R")

# Load required packages
required_packages <- c("dplyr", "shiny")

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
      warning = memory_mb > 500,
      critical = memory_mb > 800,
      message = if (memory_mb > 800) {
        "Critical memory usage! Consider restarting the session."
      } else if (memory_mb > 500) {
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
    
    if (n_cols < 3) {
      validation_result$valid <- FALSE
      validation_result$errors <- c(validation_result$errors, "Need at least 3 samples")
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
    
    if (estimated_memory_mb > 200) {
      validation_result$warnings <- c(validation_result$warnings, 
                                    "Large dataset detected. Processing may be slow.")
    }
    
    # Set processing info
    validation_result$info$dimensions <- c(n_rows, n_cols)
    validation_result$info$processing_method <- if (n_rows > 1000) "chunked" else "standard"
    
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
process_expression_data_v5 <- function(raw_data) {
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
dataUploadUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        h4("Upload Expression Data"),
        
        fileInput(
          ns("expression_file"),
          "Choose Expression Data File",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          width = "100%"
        ),
        
        br(),
        
        div(
          class = "alert alert-info",
          style = "font-size: 12px; padding: 8px;",
          h6("📋 Optional: pAnno/Annotation File"),
          p("Upload a separate annotation file (pAnno format) to automatically assign sample groups and metadata.", style = "margin: 0;")
        ),
        
        fileInput(
          ns("annotation_file"),
          "Choose Sample Annotation File (Optional)",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          width = "100%",
          placeholder = "pAnno file with sample metadata"
        ),
        
        radioButtons(
          ns("file_type"),
          "File Type:",
          choices = list(
            "Auto-detect" = "auto",
            "CSV (comma-separated)" = "csv", 
            "TSV (tab-separated)" = "tsv"
          ),
          selected = "auto",
          inline = TRUE
        ),
        
        checkboxInput(
          ns("has_header"),
          "File has header row",
          value = TRUE
        ),
        
        checkboxInput(
          ns("has_rownames"),
          "First column contains gene names",
          value = TRUE
        )
      ),
      
      column(
        6,
        h4("Data Processing Options"),
        
        numericInput(
          ns("min_count"),
          "Minimum read count per gene:",
          value = 10,
          min = 1,
          max = 1000,
          step = 1
        ),
        
        numericInput(
          ns("min_samples"),
          "Minimum samples with counts:",
          value = 3,
          min = 1,
          max = 50,
          step = 1
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
          ns("process_data"),
          "🚀 Process Data",
          class = "btn-primary btn-lg",
          style = "width: 100%;"
        )
      )
    ),
    
    br(),
    
    # Status messages
    conditionalPanel(
      condition = paste0("output['", ns("show_status"), "']"),
      div(id = ns("status_messages"))
    ),
    
    br(),
    
    # Data summary
    conditionalPanel(
      condition = paste0("output['", ns("show_summary"), "']"),
      fluidRow(
        column(
          3,
          valueBoxOutput(ns("n_genes"), width = 12)
        ),
        column(
          3, 
          valueBoxOutput(ns("n_samples"), width = 12)
        ),
        column(
          3,
          valueBoxOutput(ns("total_counts"), width = 12)
        ),
        column(
          3,
          valueBoxOutput(ns("median_counts"), width = 12)
        )
      )
    )
  )
}

# Data Upload Server Module (simplified, no gene conversion)
dataUploadServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    raw_data = NULL,
    processed_data = NULL,
    annotation_data = NULL,
    file_info = NULL,
    status_messages = list()
  )
  
  # File upload observer
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    tryCatch({
      # Add status message
      add_status_message("📁 Reading file...", "info")
      
      file_path <- input$expression_file$datapath
      file_ext <- tools::file_ext(input$expression_file$name)
      
      # Determine file type
      file_type <- input$file_type
      if (file_type == "auto") {
        file_type <- switch(
          tolower(file_ext),
          "csv" = "csv",
          "tsv" = "tsv", 
          "txt" = "tsv",
          "xlsx" = "excel",
          "xls" = "excel",
          "csv"  # default
        )
      }
      
      # Read data based on file type
      if (file_type == "excel") {
        if (requireNamespace("readxl", quietly = TRUE)) {
          raw_data <- readxl::read_excel(
            file_path,
            col_names = input$has_header
          )
        } else {
          add_status_message("❌ Excel support not available. Please use CSV or TSV format.", "danger")
          return()
        }
      } else if (file_type == "csv") {
        raw_data <- tryCatch({
          if (requireNamespace("readr", quietly = TRUE)) {
            readr::read_csv(
              file_path,
              col_names = input$has_header,
              show_col_types = FALSE
            )
          } else {
            read.csv(
              file_path,
              header = input$has_header,
              stringsAsFactors = FALSE
            )
          }
        }, error = function(e) {
          read.csv(
            file_path,
            header = input$has_header,
            stringsAsFactors = FALSE
          )
        })
      } else {  # tsv
        raw_data <- tryCatch({
          if (requireNamespace("readr", quietly = TRUE)) {
            readr::read_tsv(
              file_path,
              col_names = input$has_header,
              show_col_types = FALSE
            )
          } else {
            read.delim(
              file_path,
              header = input$has_header,
              stringsAsFactors = FALSE
            )
          }
        }, error = function(e) {
          read.delim(
            file_path,
            header = input$has_header,
            stringsAsFactors = FALSE
          )
        })
      }
      
      # Handle gene names and duplicates using v5 method
      if (input$has_rownames) {
        gene_names <- raw_data[[1]]
        data_only <- raw_data[, -1]
        
        # Check for duplicate gene IDs and aggregate if found
        if (any(duplicated(gene_names))) {
          n_duplicates <- sum(duplicated(gene_names))
          add_status_message(
            paste0("⚠️ Found ", n_duplicates, " duplicate gene IDs. Aggregating by summing counts..."), 
            "warning"
          )
          
          # Convert expression columns to numeric first
          data_only <- data_only %>%
            mutate_all(~ as.numeric(as.character(.)))
          
          # Combine gene names with data for aggregation
          combined_data <- cbind(gene_id = gene_names, data_only)
          
          # Aggregate duplicates by summing
          if (requireNamespace("dplyr", quietly = TRUE)) {
            aggregated <- combined_data %>%
              group_by(gene_id) %>%
              summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
            
            gene_names <- aggregated$gene_id
            raw_data <- aggregated[, -1]  # Remove gene_id column
          } else {
            # Fallback without dplyr
            unique_genes <- unique(gene_names)
            aggregated_data <- matrix(0, nrow = length(unique_genes), ncol = ncol(data_only))
            
            for (i in seq_along(unique_genes)) {
              gene <- unique_genes[i]
              matching_rows <- which(gene_names == gene)
              if (length(matching_rows) == 1) {
                aggregated_data[i, ] <- as.numeric(data_only[matching_rows, ])
              } else {
                aggregated_data[i, ] <- colSums(data_only[matching_rows, ], na.rm = TRUE)
              }
            }
            
            raw_data <- as.data.frame(aggregated_data)
            colnames(raw_data) <- colnames(data_only)
            gene_names <- unique_genes
          }
          
          add_status_message(
            paste0("✅ Aggregated duplicates. Reduced from ", length(raw_data[[1]]), 
                   " to ", length(gene_names), " unique genes"), 
            "success"
          )
        } else {
          raw_data <- data_only
        }
        
        raw_data <- as.data.frame(raw_data)
        rownames(raw_data) <- gene_names
        
      } else {
        raw_data <- as.data.frame(raw_data)
        rownames(raw_data) <- paste0("Gene_", 1:nrow(raw_data))
      }
      
      # Convert to numeric
      raw_data <- raw_data %>%
        mutate_all(~ as.numeric(as.character(.)))
      
      # Remove rows with all NA
      raw_data <- raw_data[rowSums(is.na(raw_data)) < ncol(raw_data), ]
      
      # Replace remaining NA with 0
      raw_data[is.na(raw_data)] <- 0
      
      local_values$raw_data <- raw_data
      local_values$file_info <- list(
        name = input$expression_file$name,
        size = input$expression_file$size,
        type = file_type
      )
      
      add_status_message(
        paste0("✅ Successfully loaded ", nrow(raw_data), " genes and ", 
               ncol(raw_data), " samples"), 
        "success"
      )
      
    }, error = function(e) {
      add_status_message(
        paste0("❌ Error reading file: ", e$message), 
        "danger"
      )
    })
  })
  
  # Process data when button is clicked
  observeEvent(input$process_data, {
    req(local_values$raw_data)
    
    tryCatch({
      add_status_message("🔄 Processing data with v5-based duplicate handling...", "info")
      
      # Process data with v5 method
      processing_result <- process_expression_data_v5(local_values$raw_data)
      
      if (processing_result$success) {
        processed_matrix <- processing_result$data
        
        # Filter by minimum count
        if (input$min_count > 0) {
          # Count samples with sufficient reads per gene
          sufficient_samples <- rowSums(processed_matrix >= input$min_count)
          keep_genes <- sufficient_samples >= input$min_samples
          
          filtered_genes <- sum(!keep_genes)
          processed_matrix <- processed_matrix[keep_genes, ]
          
          if (filtered_genes > 0) {
            add_status_message(
              paste0("🧹 Filtered out ", filtered_genes, " low-count genes"), 
              "info"
            )
          }
        }
        
        # Ensure all values are positive for DESeq2
        processed_matrix[processed_matrix < 0] <- 0
        processed_matrix <- round(processed_matrix)
        
        local_values$processed_data <- processed_matrix
        values$expression_data <- processed_matrix
        
        # DEBUG: Verify data was stored correctly
        cat("🔍 DEBUG: Processed matrix stored with dimensions:", nrow(processed_matrix), "x", ncol(processed_matrix), "\n")
        cat("🔍 DEBUG: Sample genes in processed data:", paste(head(rownames(processed_matrix), 3), collapse = ", "), "\n")
        
        add_status_message(
          paste0("✅ Data processing complete! Final dataset: ", 
                 nrow(processed_matrix), " genes × ", ncol(processed_matrix), " samples"),
          "success"
        )
        
      } else {
        add_status_message(
          paste0("❌ Processing failed: ", processing_result$error),
          "danger"
        )
      }
      
    }, error = function(e) {
      add_status_message(
        paste0("❌ Error processing data: ", e$message),
        "danger"
      )
    })
  })
  
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
    
    # Keep only last 5 messages
    if (length(local_values$status_messages) > 5) {
      local_values$status_messages <- tail(local_values$status_messages, 5)
    }
  }
  
  # Render status messages
  output$status_messages <- renderUI({
    if (length(local_values$status_messages) > 0) {
      message_divs <- lapply(local_values$status_messages, function(msg) {
        div(
          class = paste("alert", msg$type),
          style = "margin-bottom: 5px; padding: 8px 12px;",
          tags$small(paste0("[", msg$timestamp, "] ")),
          msg$message
        )
      })
      do.call(tagList, message_divs)
    }
  })
  
  # Control visibility of status and summary sections
  output$show_status <- reactive({
    length(local_values$status_messages) > 0
  })
  outputOptions(output, "show_status", suspendWhenHidden = FALSE)
  
  output$show_summary <- reactive({
    !is.null(local_values$processed_data)
  })
  outputOptions(output, "show_summary", suspendWhenHidden = FALSE)
  
  # Summary value boxes
  output$n_genes <- renderValueBox({
    n_genes <- if (!is.null(local_values$processed_data)) nrow(local_values$processed_data) else 0
    
    valueBox(
      value = formatC(n_genes, format = "d", big.mark = ","),
      subtitle = "Genes",
      icon = icon("dna"),
      color = "green",
      width = 12
    )
  })
  
  output$n_samples <- renderValueBox({
    n_samples <- if (!is.null(local_values$processed_data)) ncol(local_values$processed_data) else 0
    
    valueBox(
      value = n_samples,
      subtitle = "Samples",
      icon = icon("vials"),
      color = "blue", 
      width = 12
    )
  })
  
  output$total_counts <- renderValueBox({
    total_counts <- if (!is.null(local_values$processed_data)) {
      formatC(sum(local_values$processed_data), format = "d", big.mark = ",")
    } else "0"
    
    valueBox(
      value = total_counts,
      subtitle = "Total Counts",
      icon = icon("calculator"),
      color = "purple",
      width = 12
    )
  })
  
  output$median_counts <- renderValueBox({
    median_counts <- if (!is.null(local_values$processed_data)) {
      formatC(median(colSums(local_values$processed_data)), format = "d", big.mark = ",")
    } else "0"
    
    valueBox(
      value = median_counts,
      subtitle = "Median Library Size",
      icon = icon("chart-line"),
      color = "orange",
      width = 12
    )
  })
  
  # Return processed data and annotation data
  return(list(
    expression_data = reactive({ local_values$processed_data }),
    annotation_data = reactive({ local_values$annotation_data })
  ))
  })
}

cat("✅ v5-based data upload module loaded\n")
cat("ℹ️ This module focuses on data loading and duplicate handling\n")
cat("ℹ️ Gene ID conversion is handled separately in analysis modules\n")