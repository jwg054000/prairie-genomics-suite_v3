# Simplified Data Testing App - Focused on Upload and Validation
# Minimal version for stable real data testing

library(shiny)
library(DT)
library(readr)
library(readxl)

# Source just the essential functions
source("phase3/real_data_testing.R", local = TRUE)
source("phase3/real_data_server.R", local = TRUE)

# Simple UI
ui <- fluidPage(
  titlePanel("🧬 Prairie Genomics Suite - Real Data Testing"),
  
  tags$head(
    tags$style(HTML("
      .content-wrapper, .right-side {
        background-color: #f8fafc;
      }
      .main-header .navbar {
        background-color: #3b82f6 !important;
      }
      .upload-card {
        background: white;
        border-radius: 12px;
        padding: 20px;
        margin: 15px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
    "))
  ),
  
  div(class = "container-fluid",
    fluidRow(
      column(12,
        div(
          class = "alert alert-info",
          style = "margin: 20px 0;",
          h4("🎯 Real Data Testing System"),
          p("Upload your RNA-seq data to test our scientific guardrails against your expert knowledge.")
        )
      )
    ),
    
    fluidRow(
      # Upload Section
      column(6,
        div(class = "upload-card",
          h3("📊 Data Upload"),
          
          # Expression matrix upload
          div(style = "margin-bottom: 20px;",
            h4("Expression Matrix"),
            fileInput("expression_file", 
                     "Choose CSV/TSV File",
                     accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            p("Expected format: Genes as rows, samples as columns", 
              style = "color: #64748b; font-size: 14px;")
          ),
          
          # Sample metadata upload  
          div(style = "margin-bottom: 20px;",
            h4("Sample Metadata"),
            fileInput("metadata_file",
                     "Choose CSV/TSV File", 
                     accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            p("Required columns: Sample_ID, Condition", 
              style = "color: #64748b; font-size: 14px;")
          ),
          
          # Upload status
          uiOutput("upload_status")
        )
      ),
      
      # Preview Section
      column(6,
        div(class = "upload-card",
          h3("🔍 Data Preview"),
          uiOutput("data_preview")
        )
      )
    ),
    
    # Validation Results
    fluidRow(
      column(12,
        div(class = "upload-card",
          h3("✅ Validation Results"),
          uiOutput("validation_display")
        )
      )
    )
  )
)

# Simple Server
server <- function(input, output, session) {
  
  # Reactive storage
  data_storage <- reactiveValues(
    expression_data = NULL,
    sample_metadata = NULL,
    validation_results = NULL
  )
  
  # Expression file upload
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    withProgress(message = "Loading expression data...", {
      tryCatch({
        file_ext <- tools::file_ext(input$expression_file$name)
        
        # Load data
        expr_data <- load_expression_data(input$expression_file$datapath, file_ext)
        data_storage$expression_data <- expr_data
        
        showNotification(
          paste("✅ Expression matrix loaded:", nrow(expr_data), "genes ×", ncol(expr_data), "samples"),
          type = "success",
          duration = 5
        )
        
      }, error = function(e) {
        showNotification(
          paste("❌ Error loading expression data:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
  })
  
  # Metadata file upload
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    
    withProgress(message = "Loading metadata...", {
      tryCatch({
        file_ext <- tools::file_ext(input$metadata_file$name)
        
        # Load metadata
        metadata <- load_metadata(input$metadata_file$datapath, file_ext)
        data_storage$sample_metadata <- metadata
        
        showNotification(
          paste("✅ Sample metadata loaded:", nrow(metadata), "samples"),
          type = "success", 
          duration = 5
        )
        
        # Trigger validation if both files loaded
        if (!is.null(data_storage$expression_data)) {
          validate_data()
        }
        
      }, error = function(e) {
        showNotification(
          paste("❌ Error loading metadata:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
  })
  
  # Data validation
  validate_data <- function() {
    req(data_storage$expression_data, data_storage$sample_metadata)
    
    withProgress(message = "Validating data compatibility...", {
      tryCatch({
        validation <- validate_uploaded_data(
          data_storage$expression_data,
          data_storage$sample_metadata
        )
        data_storage$validation_results <- validation
        
        showNotification(
          "✅ Data validation completed!",
          type = "success",
          duration = 5
        )
        
      }, error = function(e) {
        showNotification(
          paste("❌ Validation error:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
  }
  
  # Upload status display
  output$upload_status <- renderUI({
    expr_loaded <- !is.null(data_storage$expression_data)
    meta_loaded <- !is.null(data_storage$sample_metadata)
    
    if (!expr_loaded && !meta_loaded) {
      div(class = "alert alert-info", "📁 Ready to upload files")
    } else if (expr_loaded && !meta_loaded) {
      div(class = "alert alert-warning", "📊 Expression data loaded. Upload metadata to continue.")
    } else if (!expr_loaded && meta_loaded) {
      div(class = "alert alert-warning", "📋 Metadata loaded. Upload expression data to continue.")
    } else {
      div(class = "alert alert-success", "🎉 Both files loaded successfully!")
    }
  })
  
  # Data preview
  output$data_preview <- renderUI({
    if (is.null(data_storage$expression_data) && is.null(data_storage$sample_metadata)) {
      div(class = "alert alert-info", "Upload files to see data preview")
    } else {
      tagList(
        if (!is.null(data_storage$expression_data)) {
          div(
            h5("Expression Matrix"),
            p(paste("Dimensions:", nrow(data_storage$expression_data), "genes ×", 
                   ncol(data_storage$expression_data), "samples")),
            DT::dataTableOutput("expr_preview", height = "200px")
          )
        },
        
        if (!is.null(data_storage$sample_metadata)) {
          div(
            h5("Sample Metadata", style = "margin-top: 20px;"),
            p(paste("Samples:", nrow(data_storage$sample_metadata))),
            DT::dataTableOutput("meta_preview", height = "200px")
          )
        }
      )
    }
  })
  
  # Expression preview table
  output$expr_preview <- DT::renderDataTable({
    req(data_storage$expression_data)
    
    # Show first 5 genes and 5 samples
    preview_data <- data_storage$expression_data[1:min(5, nrow(data_storage$expression_data)), 
                                                1:min(5, ncol(data_storage$expression_data))]
    
    DT::datatable(
      cbind(Gene = rownames(preview_data), preview_data),
      options = list(dom = 't', scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  # Metadata preview table  
  output$meta_preview <- DT::renderDataTable({
    req(data_storage$sample_metadata)
    
    DT::datatable(
      data_storage$sample_metadata,
      options = list(dom = 't', scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  # Validation display
  output$validation_display <- renderUI({
    if (is.null(data_storage$validation_results)) {
      div(class = "alert alert-info", "Upload both files to see validation results")
    } else {
      validation <- data_storage$validation_results
      status <- validation$overall_status$status
      
      status_class <- switch(status,
        "pass" = "alert-success",
        "warning" = "alert-warning", 
        "fail" = "alert-danger"
      )
      
      status_icon <- switch(status,
        "pass" = "✅",
        "warning" = "⚠️",
        "fail" = "❌"
      )
      
      div(
        class = paste("alert", status_class),
        h4(paste(status_icon, "Validation Status:", toupper(status))),
        
        if (length(validation$overall_status$issues) > 0) {
          div(
            h5("Issues Found:"),
            tags$ul(
              lapply(validation$overall_status$issues, function(issue) {
                tags$li(issue)
              })
            )
          )
        } else {
          p("🎉 No issues detected! Your data is ready for analysis.")
        },
        
        div(
          style = "margin-top: 15px;",
          p(paste("Total checks:", validation$overall_status$total_checks)),
          p(paste("Passed:", validation$overall_status$passed_checks,
                 "| Warnings:", validation$overall_status$warning_checks,
                 "| Failed:", validation$overall_status$failed_checks))
        )
      )
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)