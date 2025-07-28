# Phase 3 - Real Data Testing Server Logic
# Server-side logic for real data upload, validation, and ground truth collection
# 
# This module provides:
# - File upload handling and validation
# - Data preprocessing and standardization
# - Ground truth form generation and collection
# - System validation against expert knowledge
# - Results comparison and reporting

library(shiny)
library(readr)
library(readxl)
library(DT)

# ===========================================
# REAL DATA SERVER MODULE
# ===========================================

# Main real data testing server
real_data_testing_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive storage for uploaded data
    uploaded_data <- reactiveValues(
      expression_data = NULL,
      sample_metadata = NULL,
      ground_truth = NULL,
      validation_results = NULL,
      processing_log = list()
    )
    
    # Track upload progress
    upload_status <- reactiveValues(
      expression_uploaded = FALSE,
      metadata_uploaded = FALSE,
      validation_complete = FALSE,
      ground_truth_complete = FALSE
    )
    
    # ===========================================
    # FILE UPLOAD HANDLERS
    # ===========================================
    
    # Expression matrix upload
    observeEvent(input$expression_file, {
      req(input$expression_file)
      
      withProgress(message = "Processing expression matrix...", {
        
        tryCatch({
          # Log upload attempt
          log_entry <- list(
            timestamp = Sys.time(),
            action = "expression_upload_started",
            filename = input$expression_file$name,
            size_mb = round(input$expression_file$size / 1024 / 1024, 2)
          )
          uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
          
          # Load data based on file extension
          file_ext <- tools::file_ext(input$expression_file$name)
          expression_data <- load_expression_data(input$expression_file$datapath, file_ext)
          
          if (!is.null(expression_data)) {
            uploaded_data$expression_data <- expression_data
            upload_status$expression_uploaded <- TRUE
            
            # Log success
            log_entry <- list(
              timestamp = Sys.time(),
              action = "expression_upload_success",
              dimensions = paste(nrow(expression_data), "×", ncol(expression_data)),
              status = "success"
            )
            uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
            
            showNotification(
              paste("Expression matrix loaded successfully:", nrow(expression_data), "genes ×", ncol(expression_data), "samples"),
              type = "success",
              duration = 5
            )
            
            # Trigger validation if both files are uploaded
            if (upload_status$metadata_uploaded) {
              trigger_data_validation()
            }
          }
          
        }, error = function(e) {
          # Log error
          log_entry <- list(
            timestamp = Sys.time(),
            action = "expression_upload_error",
            error = as.character(e),
            status = "error"
          )
          uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
          
          showNotification(
            paste("Error loading expression matrix:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    })
    
    # Sample metadata upload
    observeEvent(input$metadata_file, {
      req(input$metadata_file)
      
      withProgress(message = "Processing sample metadata...", {
        
        tryCatch({
          # Log upload attempt
          log_entry <- list(
            timestamp = Sys.time(),
            action = "metadata_upload_started",
            filename = input$metadata_file$name
          )
          uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
          
          # Load metadata
          file_ext <- tools::file_ext(input$metadata_file$name)
          metadata <- load_metadata(input$metadata_file$datapath, file_ext)
          
          if (!is.null(metadata)) {
            uploaded_data$sample_metadata <- metadata
            upload_status$metadata_uploaded <- TRUE
            
            # Log success
            log_entry <- list(
              timestamp = Sys.time(),
              action = "metadata_upload_success",
              dimensions = paste(nrow(metadata), "samples ×", ncol(metadata), "columns"),
              status = "success"
            )
            uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
            
            showNotification(
              paste("Sample metadata loaded successfully:", nrow(metadata), "samples"),
              type = "success",
              duration = 5
            )
            
            # Trigger validation if both files are uploaded
            if (upload_status$expression_uploaded) {
              trigger_data_validation()
            }
          }
          
        }, error = function(e) {
          # Log error
          log_entry <- list(
            timestamp = Sys.time(),
            action = "metadata_upload_error",
            error = as.character(e),
            status = "error"
          )
          uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
          
          showNotification(
            paste("Error loading metadata:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    })
    
    # ===========================================
    # DATA VALIDATION
    # ===========================================
    
    # Trigger comprehensive data validation
    trigger_data_validation <- function() {
      withProgress(message = "Validating uploaded data...", {
        
        validation_results <- validate_uploaded_data(
          uploaded_data$expression_data,
          uploaded_data$sample_metadata
        )
        
        uploaded_data$validation_results <- validation_results
        upload_status$validation_complete <- TRUE
        
        # Log validation completion
        log_entry <- list(
          timestamp = Sys.time(),
          action = "validation_complete",
          overall_status = validation_results$overall_status$status,
          issues_found = length(validation_results$overall_status$issues)
        )
        uploaded_data$processing_log <- append(uploaded_data$processing_log, list(log_entry))
        
        # Show notification based on validation results
        if (validation_results$overall_status$status == "pass") {
          showNotification(
            "✅ Data validation passed! Your data looks ready for analysis.",
            type = "success",
            duration = 7
          )
        } else if (validation_results$overall_status$status == "warning") {
          showNotification(
            "⚠️ Data validation completed with warnings. Please review the issues found.",
            type = "warning",
            duration = 10
          )
        } else {
          showNotification(
            "❌ Data validation found critical issues. Please review and fix before proceeding.",
            type = "error",
            duration = 15
          )
        }
      })
    }
    
    # ===========================================
    # UI RENDERING
    # ===========================================
    
    # Upload status display
    output$upload_status_display <- renderUI({
      
      if (!upload_status$expression_uploaded && !upload_status$metadata_uploaded) {
        div(
          class = "alert alert-info",
          icon("info-circle"),
          "Upload your expression matrix and sample metadata to begin validation."
        )
      } else if (upload_status$expression_uploaded && !upload_status$metadata_uploaded) {
        div(
          class = "alert alert-warning",
          icon("upload"),
          "Expression matrix uploaded successfully. Please upload sample metadata to continue."
        )
      } else if (!upload_status$expression_uploaded && upload_status$metadata_uploaded) {
        div(
          class = "alert alert-warning",
          icon("upload"),
          "Sample metadata uploaded successfully. Please upload expression matrix to continue."
        )
      } else if (upload_status$expression_uploaded && upload_status$metadata_uploaded && !upload_status$validation_complete) {
        div(
          class = "alert alert-info",
          icon("spinner"),
          "Both files uploaded. Running validation..."
        )
      } else if (upload_status$validation_complete) {
        status <- uploaded_data$validation_results$overall_status$status
        
        status_config <- list(
          pass = list(class = "alert-success", icon = "check-circle", color = "#10b981"),
          warning = list(class = "alert-warning", icon = "exclamation-triangle", color = "#f59e0b"),
          fail = list(class = "alert-danger", icon = "times-circle", color = "#ef4444")
        )
        
        config <- status_config[[status]]
        
        div(
          class = paste("alert", config$class),
          icon(config$icon),
          strong("Validation Complete: "),
          paste("Status -", toupper(status)),
          if (length(uploaded_data$validation_results$overall_status$issues) > 0) {
            span(paste(" (", length(uploaded_data$validation_results$overall_status$issues), "issues found)"))
          }
        )
      }
    })
    
    # Data preview display
    output$data_preview_display <- renderUI({
      
      if (!upload_status$expression_uploaded || !upload_status$metadata_uploaded) {
        div(
          class = "alert alert-info",
          style = "text-align: center; padding: 40px;",
          h4("No Data to Preview"),
          p("Upload both expression matrix and sample metadata to see data preview.")
        )
      } else {
        
        tagList(
          # Expression data preview
          h4("📊 Expression Matrix Preview"),
          div(
            style = "margin-bottom: 30px;",
            p(paste("Dimensions:", nrow(uploaded_data$expression_data), "genes ×", 
                   ncol(uploaded_data$expression_data), "samples")),
            
            DT::dataTableOutput(session$ns("expression_preview"), height = "300px")
          ),
          
          # Metadata preview
          h4("📋 Sample Metadata Preview"),
          div(
            style = "margin-bottom: 30px;",
            p(paste("Dimensions:", nrow(uploaded_data$sample_metadata), "samples ×", 
                   ncol(uploaded_data$sample_metadata), "columns")),
            
            DT::dataTableOutput(session$ns("metadata_preview"), height = "300px")
          ),
          
          # Validation summary
          if (upload_status$validation_complete) {
            div(
              h4("✅ Validation Summary"),
              create_validation_summary_ui(uploaded_data$validation_results)
            )
          }
        )
      }
    })
    
    # Expression data preview table
    output$expression_preview <- DT::renderDataTable({
      req(uploaded_data$expression_data)
      
      # Show first 10 genes and all samples (or first 10 samples if too many)
      preview_data <- uploaded_data$expression_data[1:min(10, nrow(uploaded_data$expression_data)), 
                                                   1:min(10, ncol(uploaded_data$expression_data))]
      
      # Add gene names as first column for display
      preview_df <- data.frame(
        Gene = rownames(preview_data),
        preview_data,
        check.names = FALSE
      )
      
      DT::datatable(
        preview_df,
        options = list(
          scrollX = TRUE,
          scrollY = "250px",
          pageLength = 10,
          dom = 't'
        ),
        rownames = FALSE
      )
    })
    
    # Metadata preview table
    output$metadata_preview <- DT::renderDataTable({
      req(uploaded_data$sample_metadata)
      
      DT::datatable(
        uploaded_data$sample_metadata,
        options = list(
          scrollX = TRUE,
          scrollY = "250px", 
          pageLength = 10,
          dom = 'tp'
        ),
        rownames = FALSE
      )
    })
    
    # Ground truth form
    output$ground_truth_form <- renderUI({
      
      if (!upload_status$validation_complete) {
        div(
          class = "alert alert-info",
          style = "text-align: center; padding: 40px;",
          h4("Complete Data Upload First"),
          p("Upload and validate your data before providing ground truth information.")
        )
      } else {
        
        # Generate form based on ground truth schema
        create_ground_truth_form(uploaded_data$validation_results)
      }
    })
    
    # System validation results
    output$validation_results_display <- renderUI({
      
      if (!upload_status$ground_truth_complete) {
        div(
          class = "alert alert-info",
          style = "text-align: center; padding: 40px;",
          h4("Complete Ground Truth Collection First"),
          p("Provide your expert assessment to validate our system's performance.")
        )
      } else {
        
        # Generate validation comparison results
        create_system_validation_display(uploaded_data$ground_truth, uploaded_data$validation_results)
      }
    })
    
    # Return reactive values for integration with other modules
    return(list(
      uploaded_data = uploaded_data,
      upload_status = upload_status
    ))
  })
}

# ===========================================
# DATA LOADING FUNCTIONS
# ===========================================

# Load expression data from various formats
load_expression_data <- function(file_path, file_ext) {
  
  # Check file size and warn if very large
  file_size_mb <- file.info(file_path)$size / 1024 / 1024
  if (file_size_mb > 100) {
    warning(paste("Large file detected:", round(file_size_mb, 1), "MB. Processing may take time."))
  }
  
  switch(tolower(file_ext),
    
    "csv" = {
      tryCatch({
        # Try to read with progress indication
        data <- read_csv(file_path, show_col_types = FALSE, progress = interactive())
        
        # Validate basic structure
        if (nrow(data) == 0) {
          stop("File appears to be empty")
        }
        if (ncol(data) < 2) {
          stop("File must have at least 2 columns (genes + samples)")
        }
        
        # Convert to matrix, assuming first column is gene names
        if (is.character(data[[1]])) {
          gene_names <- data[[1]]
          
          # Check for empty gene names
          if (any(is.na(gene_names) | gene_names == "")) {
            warning("Some gene names are missing or empty")
          }
          
          # Extract numeric data
          numeric_data <- data[, -1]
          
          # Check if data is numeric
          if (!all(sapply(numeric_data, is.numeric))) {
            # Try to convert non-numeric columns
            for (i in 1:ncol(numeric_data)) {
              if (!is.numeric(numeric_data[[i]])) {
                numeric_data[[i]] <- as.numeric(as.character(numeric_data[[i]]))
              }
            }
          }
          
          data_matrix <- as.matrix(numeric_data)
          rownames(data_matrix) <- make.unique(gene_names)  # Handle duplicate gene names
          
          return(data_matrix)
        } else {
          # All columns are numeric
          data_matrix <- as.matrix(data)
          rownames(data_matrix) <- paste0("Gene_", 1:nrow(data_matrix))
          return(data_matrix)
        }
      }, error = function(e) {
        stop(paste("Error reading CSV file:", e$message))
      })
    },
    
    "tsv" = ,
    "txt" = {
      data <- read_tsv(file_path, show_col_types = FALSE)
      if (is.character(data[[1]])) {
        gene_names <- data[[1]]
        data_matrix <- as.matrix(data[, -1])
        rownames(data_matrix) <- gene_names
        return(data_matrix)
      } else {
        return(as.matrix(data))
      }
    },
    
    "xlsx" = ,
    "xls" = {
      data <- read_excel(file_path)
      if (is.character(data[[1]])) {
        gene_names <- data[[1]]
        data_matrix <- as.matrix(data[, -1])
        rownames(data_matrix) <- gene_names
        return(data_matrix)
      } else {
        return(as.matrix(data))
      }
    },
    
    "rdata" = {
      env <- new.env()
      load(file_path, envir = env)
      # Return the first matrix-like object found
      for (obj_name in ls(env)) {
        obj <- get(obj_name, envir = env)
        if (is.matrix(obj) || is.data.frame(obj)) {
          return(as.matrix(obj))
        }
      }
      stop("No matrix or data frame found in RData file")
    },
    
    "rds" = {
      data <- readRDS(file_path)
      return(as.matrix(data))
    },
    
    stop(paste("Unsupported file format:", file_ext))
  )
}

# Load sample metadata
load_metadata <- function(file_path, file_ext) {
  
  switch(tolower(file_ext),
    
    "csv" = read_csv(file_path, show_col_types = FALSE),
    "tsv" = ,
    "txt" = read_tsv(file_path, show_col_types = FALSE),
    "xlsx" = ,
    "xls" = read_excel(file_path),
    
    stop(paste("Unsupported metadata format:", file_ext))
  )
}

# ===========================================
# UI GENERATION FUNCTIONS
# ===========================================

# Create validation summary UI
create_validation_summary_ui <- function(validation_results) {
  
  if (is.null(validation_results)) return(div("No validation results available"))
  
  overall_status <- validation_results$overall_status$status
  
  status_colors <- list(
    pass = "#10b981",
    warning = "#f59e0b", 
    fail = "#ef4444"
  )
  
  color <- status_colors[[overall_status]]
  
  div(
    style = paste0("border: 2px solid ", color, "; border-radius: 8px; padding: 15px; margin: 10px 0;"),
    
    div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      div(
        style = paste0("width: 20px; height: 20px; border-radius: 50%; background: ", color, "; margin-right: 10px;")
      ),
      h5(paste("Overall Status:", toupper(overall_status)), style = "margin: 0;")
    ),
    
    if (length(validation_results$overall_status$issues) > 0) {
      div(
        h6("Issues Found:"),
        tags$ul(
          lapply(validation_results$overall_status$issues, function(issue) {
            tags$li(issue)
          })
        )
      )
    } else {
      div(
        class = "alert alert-success",
        style = "margin-top: 10px;",
        "✅ No issues detected! Your data is ready for analysis."
      )
    }
  )
}

# Create ground truth collection form
create_ground_truth_form <- function(validation_results) {
  
  # This would generate a comprehensive form based on GROUND_TRUTH_SCHEMA
  # For now, showing a simplified version
  
  div(
    div(
      class = "alert alert-info",
      style = "margin-bottom: 20px;",
      icon("lightbulb"),
      strong("Your Expertise Matters: "),
      "Help us validate our system by sharing your analysis decisions and biological knowledge."
    ),
    
    # Placeholder form - would be dynamically generated in full implementation
    div(
      style = "background: #f8fafc; padding: 20px; border-radius: 8px;",
      h4("📋 Ground Truth Collection Form"),
      p("This form will collect your expert assessment of:"),
      tags$ul(
        tags$li("Experiment type and confidence level"),
        tags$li("Analysis parameters and rationale"),
        tags$li("Expected biological outcomes"),
        tags$li("Known positive/negative controls"),
        tags$li("Previous analysis challenges and solutions")
      ),
      
      div(
        style = "margin-top: 20px; text-align: center;",
        actionButton(
          "ground_truth_form_placeholder",
          "📝 Start Ground Truth Collection",
          class = "btn btn-primary btn-lg"
        )
      )
    )
  )
}

# Create system validation display
create_system_validation_display <- function(ground_truth, validation_results) {
  
  div(
    h4("🔬 System Validation Results"),
    p("Comparison between our system's assessment and your expert knowledge"),
    
    # Placeholder for validation comparison
    div(
      class = "alert alert-info",
      "System validation comparison will be displayed here once ground truth collection is complete."
    )
  )
}

# Export functions
list(
  real_data_testing_server = real_data_testing_server,
  load_expression_data = load_expression_data,
  load_metadata = load_metadata
)