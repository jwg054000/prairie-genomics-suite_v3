# Prairie Genomics Suite - R Shiny Implementation (Simple Version)
# Minimal version for shinyapps.io deployment

# Load only essential packages
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(readr)

# Optional packages
excel_support <- FALSE
bioc_available <- FALSE
enhanced_plots <- FALSE

tryCatch({
  library(readxl)
  excel_support <- TRUE
}, error = function(e) NULL)

tryCatch({
  library(RColorBrewer)
  library(pheatmap)
  enhanced_plots <- TRUE
}, error = function(e) NULL)

# Try to load DESeq2 (may fail on shinyapps.io)
tryCatch({
  library(DESeq2)
  bioc_available <- TRUE
}, error = function(e) {
  cat("Note: DESeq2 not available - some analysis features will be limited\n")
})

# Source modules with error handling
tryCatch({
  source("modules/data_upload.R")
  source("modules/sample_annotation.R") 
  source("modules/deseq2_analysis.R")
  source("modules/visualization.R")
}, error = function(e) {
  stop(paste("Error loading modules:", e$message))
})

# Define UI (same as original)
ui <- dashboardPage(
  dashboardHeader(
    title = tags$span(
      icon("dna"), 
      "Prairie Genomics Suite",
      style = "font-size: 18px; font-weight: bold;"
    ),
    tags$li(
      class = "dropdown",
      tags$a(
        href = "#",
        tags$span("R Shiny v1.0", style = "color: #ffffff; font-size: 12px;")
      )
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("📁 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("🧬 Sample Annotation", tabName = "annotation", icon = icon("tags")),
      menuItem("🚀 DESeq2 Analysis", tabName = "analysis", icon = icon("chart-line")),
      menuItem("📊 Visualizations", tabName = "visualizations", icon = icon("chart-bar")),
      menuItem("📋 Results Export", tabName = "export", icon = icon("download")),
      
      hr(),
      
      # Analysis progress indicator
      div(
        style = "padding: 15px;",
        h4("Analysis Progress", style = "color: white; font-size: 14px;"),
        
        # Progress steps
        div(
          style = "color: #aaa; font-size: 12px;",
          div(id = "step1", "📁 Upload Data"),
          div(id = "step2", "🧬 Annotate Samples"),
          div(id = "step3", "🚀 Run Analysis"),
          div(id = "step4", "📊 View Results")
        ),
        
        br(),
        # Simple progress display without shinyWidgets
        div(
          style = "background-color: #3c8dbc; color: white; padding: 8px; border-radius: 4px; text-align: center;",
          id = "overall_progress",
          "Analysis Progress: 0%"
        )
      ),
      
      hr(),
      
      # Quick stats
      div(
        style = "padding: 15px;",
        h4("Quick Stats", style = "color: white; font-size: 14px;"),
        
        valueBoxOutput("genes_count", width = 12),
        valueBoxOutput("samples_count", width = 12),
        valueBoxOutput("groups_count", width = 12)
      )
    )
  ),
  
  dashboardBody(
    # Custom CSS
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        
        .small-box {
          margin-bottom: 10px;
        }
        
        .small-box .inner {
          padding: 8px 15px;
        }
        
        .small-box h3 {
          font-size: 20px;
          margin: 0;
        }
        
        .small-box p {
          font-size: 12px;
          margin: 0;
        }
        
        .progress {
          height: 15px;
          margin-bottom: 10px;
        }
        
        .alert-success {
          background-color: #d4edda;
          border-color: #c3e6cb;
          color: #155724;
        }
        
        .alert-warning {
          background-color: #fff3cd;
          border-color: #ffeaa7;
          color: #856404;
        }
        
        .alert-info {
          background-color: #cce5ff;
          border-color: #b8daff;
          color: #004085;
        }
      "))
    ),
    
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
            
            div(
              class = "alert alert-info",
              h4("Welcome to Prairie Genomics Suite!"),
              p("Upload your RNA-seq expression data to get started. We support CSV, TSV, and Excel formats."),
              p("Your data should have genes as rows and samples as columns.")
            ),
            
            dataUploadUI("data_upload")
          )
        ),
        
        # Data preview
        fluidRow(
          box(
            title = "📊 Data Preview", 
            status = "info", 
            solidHeader = TRUE, 
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            DT::dataTableOutput("data_preview")
          )
        )
      ),
      
      # Sample Annotation Tab
      tabItem(
        tabName = "annotation",
        fluidRow(
          box(
            title = "🧬 Sample Annotation", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_expression_data == false",
              div(
                class = "alert alert-warning",
                h4("Upload Data First"),
                p("Please upload your expression data in the Data Upload tab before proceeding with sample annotation.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_expression_data == true",
              sampleAnnotationUI("sample_annotation")
            )
          )
        )
      ),
      
      # DESeq2 Analysis Tab
      tabItem(
        tabName = "analysis",
        fluidRow(
          box(
            title = "🚀 DESeq2 Differential Expression Analysis", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_annotation_data == false",
              div(
                class = "alert alert-warning",
                h4("Complete Sample Annotation First"),
                p("Please annotate your samples in the Sample Annotation tab before running DESeq2 analysis.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_annotation_data == true",
              deseq2AnalysisUI("deseq2_analysis")
            )
          )
        )
      ),
      
      # Visualizations Tab
      tabItem(
        tabName = "visualizations",
        fluidRow(
          box(
            title = "📊 Interactive Visualizations", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_results == false",
              div(
                class = "alert alert-warning",
                h4("Run Analysis First"),
                p("Please complete your DESeq2 analysis before viewing visualizations.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_results == true",
              visualizationUI("visualization")
            )
          )
        )
      ),
      
      # Export Tab
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "📋 Results Export", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_results == false",
              div(
                class = "alert alert-warning",
                h4("No Results to Export"),
                p("Please complete your analysis before exporting results.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_results == true",
              h4("Export Options"),
              
              fluidRow(
                column(
                  6,
                  h5("📊 Analysis Results"),
                  downloadButton("download_results", "Download DESeq2 Results", 
                               class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                  downloadButton("download_significant", "Download Significant Genes Only", 
                               class = "btn-success", style = "width: 100%; margin-bottom: 10px;"),
                  downloadButton("download_filtered", "Download Filtered Expression Data", 
                               class = "btn-info", style = "width: 100%;")
                ),
                column(
                  6,
                  h5("📈 Plots"),
                  downloadButton("download_volcano", "Download Volcano Plot", 
                               class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                  downloadButton("download_heatmap", "Download Heatmap", 
                               class = "btn-success", style = "width: 100%; margin-bottom: 10px;"),
                  downloadButton("download_pca", "Download PCA Plot", 
                               class = "btn-info", style = "width: 100%;")
                )
              )
            )
          )
        )
      )
    )
  )
)

# Define Server (same as original but simplified)
server <- function(input, output, session) {
  # Reactive values to store data
  values <- reactiveValues(
    expression_data = NULL,
    annotation_data = NULL,
    deseq2_results = NULL,
    filtered_data = NULL,
    progress = 0
  )
  
  # Data Upload Module
  expression_data <- callModule(dataUpload, "data_upload", values)
  
  # Sample Annotation Module  
  annotation_data <- callModule(sampleAnnotation, "sample_annotation", values)
  
  # DESeq2 Analysis Module
  analysis_results <- callModule(deseq2Analysis, "deseq2_analysis", values)
  
  # Visualization Module
  callModule(visualization, "visualization", values)
  
  # Progress tracking
  observe({
    progress <- 0
    
    if (!is.null(values$expression_data)) {
      progress <- progress + 25
    }
    
    if (!is.null(values$annotation_data)) {
      progress <- progress + 25
    }
    
    if (!is.null(values$deseq2_results)) {
      progress <- progress + 50
    }
    
    values$progress <- progress
  })
  
  # Conditional panel outputs
  output$has_expression_data <- reactive({
    !is.null(values$expression_data)
  })
  outputOptions(output, "has_expression_data", suspendWhenHidden = FALSE)
  
  output$has_annotation_data <- reactive({
    !is.null(values$annotation_data)
  })
  outputOptions(output, "has_annotation_data", suspendWhenHidden = FALSE)
  
  output$has_results <- reactive({
    !is.null(values$deseq2_results)
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  # Quick stats value boxes
  output$genes_count <- renderValueBox({
    count <- if (!is.null(values$expression_data)) nrow(values$expression_data) else 0
    valueBox(
      value = formatC(count, format = "d", big.mark = ","),
      subtitle = "Genes",
      icon = icon("dna"),
      color = if (count > 0) "green" else "black",
      width = 12
    )
  })
  
  output$samples_count <- renderValueBox({
    count <- if (!is.null(values$expression_data)) ncol(values$expression_data) else 0
    valueBox(
      value = count,
      subtitle = "Samples", 
      icon = icon("vials"),
      color = if (count > 0) "blue" else "black",
      width = 12
    )
  })
  
  output$groups_count <- renderValueBox({
    count <- if (!is.null(values$annotation_data)) length(unique(values$annotation_data$Condition)) else 0
    valueBox(
      value = count,
      subtitle = "Groups",
      icon = icon("layer-group"),
      color = if (count > 1) "purple" else "black", 
      width = 12
    )
  })
  
  # Data preview
  output$data_preview <- DT::renderDataTable({
    if (!is.null(values$expression_data)) {
      # Show first 10 rows and columns
      preview_data <- values$expression_data[1:min(10, nrow(values$expression_data)), 
                                           1:min(10, ncol(values$expression_data))]
      preview_data <- cbind(Gene = rownames(preview_data), preview_data)
      
      DT::datatable(
        preview_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        )
      )
    }
  })
  
  # Download handlers
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("deseq2_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$deseq2_results)) {
        write.csv(values$deseq2_results, file, row.names = TRUE)
      }
    }
  )
  
  output$download_significant <- downloadHandler(
    filename = function() {
      paste0("significant_genes_", Sys.Date(), ".csv") 
    },
    content = function(file) {
      if (!is.null(values$deseq2_results)) {
        significant <- values$deseq2_results[values$deseq2_results$padj < 0.05 & 
                                           abs(values$deseq2_results$log2FoldChange) > 1, ]
        write.csv(significant, file, row.names = TRUE)
      }
    }
  )
  
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("filtered_expression_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$filtered_data)) {
        write.csv(values$filtered_data, file, row.names = TRUE)
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)