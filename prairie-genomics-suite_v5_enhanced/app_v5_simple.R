# Prairie Genomics Suite v5 Enhanced - Simplified Version
# Handles optional dependencies gracefully

# Load required packages with error handling
if (!require("pacman", quietly = TRUE)) install.packages("pacman")

# Core packages (required)
pacman::p_load(
  shiny, shinydashboard, DT, plotly, ggplot2, dplyr, readr,
  install = TRUE
)

# Increase Shiny file upload limit to 500MB (default is 5MB)
options(shiny.maxRequestSize = 500*1024^2)

# Handle Bioconductor packages separately
bioc_packages <- c("DESeq2", "sva")
for (pkg in bioc_packages) {
  tryCatch({
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat("✅", pkg, "loaded successfully\n")
    } else {
      # Try to install Bioconductor packages
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat("✅", pkg, "installed and loaded\n")
    }
  }, error = function(e) {
    cat("Note:", pkg, "not available, some features may be limited\n")
  })
}

# Handle CRAN packages
cran_packages <- c("MASS", "pheatmap", "RColorBrewer", "car", "ggrepel", "readxl", "logger")
for (pkg in cran_packages) {
  tryCatch({
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat("✅", pkg, "loaded successfully\n")
    } else {
      # Try with pacman for CRAN packages
      pacman::p_load(pkg, character.only = TRUE, install = TRUE)
      cat("✅", pkg, "installed and loaded\n")
    }
  }, error = function(e) {
    cat("Note:", pkg, "not available, some features may be limited\n")
  })
}

# Create logger fallback functions if logger not available
if (!requireNamespace("logger", quietly = TRUE)) {
  log_error <- function(msg) cat("ERROR:", msg, "\n")
  log_warn <- function(msg) cat("WARNING:", msg, "\n")
  log_info <- function(msg) cat("INFO:", msg, "\n")
}

# Check availability
bioc_available <- requireNamespace("DESeq2", quietly = TRUE)
excel_support <- requireNamespace("readxl", quietly = TRUE)
enhanced_plots <- requireNamespace("RColorBrewer", quietly = TRUE)
rgl_support <- requireNamespace("rgl", quietly = TRUE)

cat("Package status:\n")
cat("- DESeq2:", if(bioc_available) "✅" else "❌", "\n")
cat("- Excel support:", if(excel_support) "✅" else "❌", "\n")
cat("- Enhanced plots:", if(enhanced_plots) "✅" else "❌", "\n")
cat("- 3D plotting:", if(rgl_support) "✅" else "❌", "\n")

# Source modules
source("modules/enhanced_sample_annotation.R")
source("modules/enhanced_deseq2_analysis.R")
source("modules/context7_visualizations.R")
source("utils/batch_correction.R")

# Helper functions for performance and memory management
validate_file_size <- function(file_path, max_size_mb = 200) {
  tryCatch({
    size_bytes <- file.size(file_path)
    size_mb <- round(size_bytes / 1024^2, 2)
    
    result <- list(
      valid = TRUE,
      size_mb = size_mb,
      message = ""
    )
    
    if (size_mb > 500) {
      result$valid <- FALSE
      result$message <- paste("File too large:", size_mb, "MB. Maximum allowed: 500 MB")
    } else if (size_mb > max_size_mb) {
      result$message <- paste("Large file detected:", size_mb, "MB. Processing may be slow.")
    } else if (size_mb > 50) {
      result$message <- paste("Medium file detected:", size_mb, "MB. Using optimized processing.")
    }
    
    return(result)
  }, error = function(e) {
    return(list(valid = FALSE, size_mb = 0, message = paste("Error checking file size:", e$message)))
  })
}

monitor_memory <- function() {
  tryCatch({
    gc_info <- gc()
    memory_mb <- round(sum(gc_info[, 2]), 1)
    
    list(
      used_mb = memory_mb,
      warning = memory_mb > 1000,
      critical = memory_mb > 2000,
      message = if (memory_mb > 2000) {
        "Critical memory usage! Consider restarting the session."
      } else if (memory_mb > 1000) {
        "High memory usage detected. Performance may be affected."
      } else {
        paste("Memory usage:", memory_mb, "MB")
      }
    )
  }, error = function(e) {
    list(used_mb = 0, warning = FALSE, critical = FALSE, message = "Memory monitoring unavailable")
  })
}

process_large_dataset <- function(raw_data, chunk_size = 5000) {
  tryCatch({
    # Check if chunking is needed
    n_rows <- nrow(raw_data)
    n_cols <- ncol(raw_data)
    
    cat("Processing dataset:", n_rows, "x", n_cols, "\n")
    
    # Monitor memory before processing
    mem_before <- monitor_memory()
    cat("Memory before processing:", mem_before$message, "\n")
    
    # Convert first column to row names if it contains gene IDs
    if (is.character(raw_data[[1]]) || is.factor(raw_data[[1]])) {
      gene_names <- raw_data[[1]]
      expression_data <- raw_data[, -1, drop = FALSE]
    } else {
      gene_names <- paste0("Gene_", 1:n_rows)
      expression_data <- raw_data
    }
    
    # Process in chunks if dataset is large
    if (n_rows > chunk_size) {
      cat("Large dataset detected. Processing in chunks of", chunk_size, "genes.\n")
      
      # Initialize result matrix
      result_matrix <- NULL
      
      # Process in chunks
      for (start_row in seq(1, n_rows, by = chunk_size)) {
        end_row <- min(start_row + chunk_size - 1, n_rows)
        cat("Processing genes", start_row, "to", end_row, "\n")
        
        # Extract chunk
        chunk_data <- expression_data[start_row:end_row, , drop = FALSE]
        
        # Convert to numeric matrix
        chunk_matrix <- tryCatch({
          as.matrix(chunk_data)
        }, error = function(e) {
          # Try to convert character columns to numeric
          numeric_data <- lapply(chunk_data, function(x) {
            if (is.character(x) || is.factor(x)) {
              as.numeric(as.character(x))
            } else {
              as.numeric(x)
            }
          })
          as.matrix(data.frame(numeric_data))
        })
        
        # Check for numeric conversion issues
        if (!is.numeric(chunk_matrix)) {
          stop("Unable to convert data to numeric matrix")
        }
        
        # Remove all-zero genes from chunk
        keep_genes <- rowSums(chunk_matrix, na.rm = TRUE) > 0
        if (sum(keep_genes) == 0) {
          cat("Warning: All genes in chunk", start_row, "to", end_row, "have zero counts\n")
          next
        }
        
        chunk_matrix <- chunk_matrix[keep_genes, , drop = FALSE]
        chunk_gene_names <- gene_names[start_row:end_row][keep_genes]
        
        # Combine with result
        if (is.null(result_matrix)) {
          result_matrix <- chunk_matrix
          result_gene_names <- chunk_gene_names
        } else {
          result_matrix <- rbind(result_matrix, chunk_matrix)
          result_gene_names <- c(result_gene_names, chunk_gene_names)
        }
        
        # Clean up chunk data
        rm(chunk_data, chunk_matrix)
        gc()
      }
      
      # Set final rownames
      if (!is.null(result_matrix)) {
        rownames(result_matrix) <- result_gene_names
      }
      
    } else {
      # Process normally for smaller datasets
      cat("Standard processing for manageable dataset size.\n")
      
      # Convert to numeric matrix
      result_matrix <- tryCatch({
        as.matrix(expression_data)
      }, error = function(e) {
        # Try to convert character columns to numeric
        numeric_data <- lapply(expression_data, function(x) {
          if (is.character(x) || is.factor(x)) {
            as.numeric(as.character(x))
          } else {
            as.numeric(x)
          }
        })
        as.matrix(data.frame(numeric_data))
      })
      
      # Check for numeric conversion
      if (!is.numeric(result_matrix)) {
        stop("Expression data must be numeric")
      }
      
      # Remove all-zero genes
      keep_genes <- rowSums(result_matrix, na.rm = TRUE) > 0
      result_matrix <- result_matrix[keep_genes, , drop = FALSE]
      
      # Set rownames
      rownames(result_matrix) <- gene_names[keep_genes]
    }
    
    # Final validation
    if (is.null(result_matrix) || nrow(result_matrix) == 0) {
      stop("No valid genes remaining after processing")
    }
    
    # Monitor memory after processing
    mem_after <- monitor_memory()
    cat("Memory after processing:", mem_after$message, "\n")
    
    # Final cleanup
    gc()
    
    return(result_matrix)
    
  }, error = function(e) {
    cat("Error in process_large_dataset:", e$message, "\n")
    stop(paste("Data processing failed:", e$message))
  })
}

# UI
ui <- dashboardPage(
  dashboardHeader(
    title = tags$span(
      icon("dna"),
      "Prairie Genomics Suite v5",
      style = "font-size: 18px; font-weight: bold;"
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("📁 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("🧬 Sample Annotation", tabName = "annotation", icon = icon("tags")),
      menuItem("🚀 DESeq2 Analysis", tabName = "analysis", icon = icon("chart-line")),
      menuItem("🎨 Visualizations", tabName = "visualizations", icon = icon("chart-bar")),
      menuItem("📋 Export", tabName = "export", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Data Upload Tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "📁 Data Upload", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            fluidRow(
              column(
                6,
                h5("📊 Expression Data"),
                
                fileInput(
                  "expression_file",
                  "Choose Expression Data File",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx")
                ),
                
                checkboxInput(
                  "expression_has_header",
                  "File has header row",
                  value = TRUE
                ),
                
                conditionalPanel(
                  condition = "output.has_expression_file",
                  actionButton(
                    "process_expression_data",
                    "🔄 Process Data",
                    class = "btn-primary"
                  )
                )
              ),
              
              column(
                6,
                h5("📄 Test Data"),
                
                actionButton(
                  "load_simple_test",
                  "📊 Load Simple Test Data",
                  class = "btn-info",
                  style = "width: 100%; margin-bottom: 10px;"
                ),
                
                actionButton(
                  "load_complex_test",
                  "🧬 Load Multi-Group Test",
                  class = "btn-warning",
                  style = "width: 100%;"
                ),
                
                br(), br(),
                
                h6("Requirements:"),
                tags$ul(
                  tags$li("Genes as rows, samples as columns"),
                  tags$li("Raw counts (not normalized)"),
                  tags$li("Clear sample naming patterns")
                )
              )
            )
          )
        ),
        
        # Data preview
        conditionalPanel(
          condition = "output.has_processed_data",
          fluidRow(
            box(
              title = "📊 Data Preview", 
              status = "info", 
              solidHeader = TRUE, 
              width = 12,
              
              DT::dataTableOutput("data_preview")
            )
          )
        )
      ),
      
      # Sample Annotation Tab
      tabItem(
        tabName = "annotation",
        conditionalPanel(
          condition = "!output.has_processed_data",
          fluidRow(
            box(
              title = "⚠️ Data Required",
              status = "warning",
              width = 12,
              p("Please upload expression data first.")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_processed_data",
          enhancedSampleAnnotationUI("sample_annotation")
        )
      ),
      
      # DESeq2 Analysis Tab
      tabItem(
        tabName = "analysis",
        conditionalPanel(
          condition = "!output.has_annotation_data",
          fluidRow(
            box(
              title = "⚠️ Annotation Required",
              status = "warning", 
              width = 12,
              p("Please complete sample annotation first.")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_annotation_data",
          
          conditionalPanel(
            condition = "!output.deseq2_available",
            fluidRow(
              box(
                title = "⚠️ DESeq2 Not Available",
                status = "danger",
                width = 12,
                div(
                  class = "alert alert-danger",
                  h4("DESeq2 Package Required"),
                  p("Please install DESeq2 from Bioconductor:"),
                  tags$pre("BiocManager::install('DESeq2')")
                )
              )
            )
          ),
          
          conditionalPanel(
            condition = "output.deseq2_available",
            enhancedDESeq2AnalysisUI("deseq2_analysis")
          )
        )
      ),
      
      # Visualizations Tab
      tabItem(
        tabName = "visualizations",
        conditionalPanel(
          condition = "!output.has_results",
          fluidRow(
            box(
              title = "⚠️ Results Required",
              status = "warning",
              width = 12,
              p("Please complete analysis first.")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_results",
          context7VisualizationsUI("visualizations")
        )
      ),
      
      # Export Tab
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "📋 Export Results", 
            status = "primary", 
            width = 12,
            
            conditionalPanel(
              condition = "!output.has_results",
              p("No results available for export.")
            ),
            
            conditionalPanel(
              condition = "output.has_results",
              h4("Export Options"),
              
              downloadButton(
                "download_results",
                "📥 Download Results",
                class = "btn-primary",
                style = "margin: 5px;"
              ),
              
              downloadButton(
                "download_plots",
                "🎨 Download Plots",
                class = "btn-success",
                style = "margin: 5px;"
              )
            )
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    expression_data = NULL,
    annotation_data = NULL,
    deseq2_results = NULL,
    comparison_results = NULL,
    filtered_data = NULL,
    batch_data = NULL,
    raw_expression_data = NULL,
    file_info = NULL,
    processing_log = list()
  )
  
  # Enhanced data upload handling with file size validation
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    tryCatch({
      file_path <- input$expression_file$datapath
      file_name <- input$expression_file$name
      file_ext <- tools::file_ext(file_name)
      
      # Validate file size before processing
      size_check <- validate_file_size(file_path)
      values$file_info <- size_check
      
      if (!size_check$valid) {
        showNotification(size_check$message, type = "error")
        return()
      }
      
      # Show size warning if needed
      if (size_check$message != "") {
        showNotification(size_check$message, type = "warning")
      }
      
      # Show loading message for large files
      if (size_check$size_mb > 50) {
        showNotification("Processing large file... This may take a moment.", type = "message")
      }
      
      # Read file based on extension with progress indication
      cat("Reading file:", file_name, "(", size_check$size_mb, "MB )\n")
      
      if (file_ext == "csv") {
        raw_data <- readr::read_csv(file_path, col_names = input$expression_has_header, show_col_types = FALSE)
      } else if (file_ext %in% c("tsv", "txt")) {
        raw_data <- readr::read_delim(file_path, delim = "\t", col_names = input$expression_has_header, show_col_types = FALSE)
      } else if (file_ext == "xlsx" && excel_support) {
        raw_data <- readxl::read_excel(file_path, col_names = input$expression_has_header)
      } else {
        stop("Unsupported file format or readxl not available")
      }
      
      # Basic validation
      if (is.null(raw_data) || nrow(raw_data) == 0) {
        stop("File appears to be empty or unreadable")
      }
      
      if (ncol(raw_data) < 2) {
        stop("File must have at least 2 columns (genes and sample data)")
      }
      
      # Store raw data and file information
      values$raw_expression_data <- raw_data
      values$file_info <- c(size_check, list(
        rows = nrow(raw_data),
        cols = ncol(raw_data),
        file_name = file_name
      ))
      
      # Add to processing log
      log_entry <- paste("File uploaded:", file_name, "-", nrow(raw_data), "x", ncol(raw_data), 
                        "- Size:", size_check$size_mb, "MB")
      values$processing_log <- append(values$processing_log, log_entry)
      
      showNotification("File uploaded successfully! Click 'Process Data' to continue.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
      cat("Upload error:", e$message, "\n")
    })
  })
  
  # Enhanced data processing with memory management
  observeEvent(input$process_expression_data, {
    req(values$raw_expression_data)
    
    tryCatch({
      # Show processing notification
      showNotification("Processing expression data... Please wait.", type = "message")
      
      # Check memory before processing
      mem_status <- monitor_memory()
      if (mem_status$critical) {
        showNotification(mem_status$message, type = "error")
        return()
      } else if (mem_status$warning) {
        showNotification(mem_status$message, type = "warning")
      }
      
      # Process using enhanced chunked processing
      cat("Starting enhanced data processing...\n")
      expression_matrix <- process_large_dataset(values$raw_expression_data)
      
      # Final validation
      if (is.null(expression_matrix) || nrow(expression_matrix) == 0) {
        stop("No valid genes remaining after processing")
      }
      
      if (ncol(expression_matrix) < 2) {
        stop("Insufficient samples remaining after processing")
      }
      
      # Store processed data
      values$expression_data <- expression_matrix
      
      # Update processing log
      log_entry <- paste("Data processed successfully:", nrow(expression_matrix), "genes x", 
                        ncol(expression_matrix), "samples")
      values$processing_log <- append(values$processing_log, log_entry)
      
      # Final memory check
      final_mem <- monitor_memory()
      cat("Final memory status:", final_mem$message, "\n")
      
      showNotification(paste("Data processed successfully!", nrow(expression_matrix), 
                           "genes retained from", nrow(values$raw_expression_data), "original genes."), 
                      type = "message")
      
    }, error = function(e) {
      showNotification(paste("Processing error:", e$message), type = "error")
      cat("Processing error:", e$message, "\n")
      
      # Clean up on error
      gc()
    })
  })
  
  # Load test data
  observeEvent(input$load_simple_test, {
    if (file.exists("test_data.csv")) {
      test_data <- read.csv("test_data.csv", row.names = 1)
      values$expression_data <- test_data
      showNotification("Simple test data loaded!", type = "message")
    } else {
      # Create simple test data
      set.seed(42)
      test_data <- matrix(
        rnbinom(60, size = 10, mu = 100),
        nrow = 10, ncol = 6
      )
      colnames(test_data) <- c("Control_1", "Control_2", "Control_3", 
                              "Treatment_1", "Treatment_2", "Treatment_3")
      rownames(test_data) <- paste0("Gene_", 1:10)
      
      values$expression_data <- test_data
      showNotification("Simple test data generated!", type = "message")
    }
  })
  
  observeEvent(input$load_complex_test, {
    # Create complex multi-group test data
    set.seed(42)
    
    groups <- c("B1", "TransB", "aN", "DN", "SM")
    n_per_group <- 5
    n_genes <- 1000
    
    sample_names <- paste(rep(groups, each = n_per_group), 1:n_per_group, sep = "_")
    
    test_data <- matrix(
      rnbinom(n_genes * length(sample_names), size = 10, mu = 100),
      nrow = n_genes,
      ncol = length(sample_names)
    )
    
    colnames(test_data) <- sample_names
    rownames(test_data) <- paste0("Gene_", 1:n_genes)
    
    # Add some differential expression
    de_genes <- 1:100
    group_effects <- c(0, 2, -1.5, 1, -0.5)
    
    for (i in seq_along(groups)) {
      group_samples <- grep(paste0("^", groups[i]), sample_names)
      if (length(group_samples) > 0) {
        test_data[de_genes, group_samples] <- test_data[de_genes, group_samples] * 
          exp(group_effects[i])
      }
    }
    
    values$expression_data <- test_data
    showNotification("Complex multi-group test data loaded!", type = "message")
  })
  
  # Enhanced sample annotation module with error handling
  tryCatch({
    annotation_result <- callModule(enhancedSampleAnnotation, "sample_annotation", values)
    cat("✅ Sample annotation module loaded successfully\n")
    
    # Monitor annotation data changes
    observe({
      if (!is.null(values$annotation_data)) {
        cat("✅ Annotation data updated:", nrow(values$annotation_data), "samples,", 
            length(unique(values$annotation_data$Condition)), "groups\n")
      }
    })
    
  }, error = function(e) {
    cat("❌ Error loading sample annotation module:", e$message, "\n")
    showNotification("Error loading sample annotation module", type = "error")
  })
  
  # Enhanced DESeq2 analysis module (only if available)
  if (bioc_available) {
    tryCatch({
      analysis_result <- callModule(enhancedDESeq2Analysis, "deseq2_analysis", values)
      cat("✅ DESeq2 analysis module loaded successfully\n")
    }, error = function(e) {
      cat("❌ Error loading DESeq2 analysis module:", e$message, "\n")
      showNotification("Error loading DESeq2 analysis module", type = "error")
    })
  } else {
    cat("⚠️ DESeq2 not available - analysis module not loaded\n")
  }
  
  # Context7 visualizations module
  tryCatch({
    viz_result <- callModule(context7Visualizations, "visualizations", values)
    cat("✅ Visualizations module loaded successfully\n")
  }, error = function(e) {
    cat("❌ Error loading visualizations module:", e$message, "\n")
    showNotification("Error loading visualizations module", type = "error")
  })
  
  # Conditional outputs
  output$has_expression_file <- reactive({
    !is.null(input$expression_file)
  })
  outputOptions(output, "has_expression_file", suspendWhenHidden = FALSE)
  
  output$has_processed_data <- reactive({
    !is.null(values$expression_data)
  })
  outputOptions(output, "has_processed_data", suspendWhenHidden = FALSE)
  
  output$has_annotation_data <- reactive({
    !is.null(values$annotation_data)
  })
  outputOptions(output, "has_annotation_data", suspendWhenHidden = FALSE)
  
  output$has_results <- reactive({
    !is.null(values$deseq2_results) || !is.null(values$comparison_results)
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  output$deseq2_available <- reactive({
    bioc_available
  })
  outputOptions(output, "deseq2_available", suspendWhenHidden = FALSE)
  
  # Data preview
  output$data_preview <- DT::renderDataTable({
    req(values$expression_data)
    
    preview_data <- values$expression_data[1:min(10, nrow(values$expression_data)), 
                                         1:min(8, ncol(values$expression_data))]
    preview_data <- cbind(Gene = rownames(preview_data), preview_data)
    
    DT::datatable(
      preview_data,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = 'tip'
      )
    )
  })
  
  # Download handlers
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("genomics_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$deseq2_results)) {
        write.csv(values$deseq2_results, file, row.names = TRUE)
      } else if (!is.null(values$comparison_results)) {
        # Export first comparison as example
        first_result <- values$comparison_results[[1]]
        if (!is.null(first_result$results_df)) {
          write.csv(first_result$results_df, file, row.names = FALSE)
        }
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)