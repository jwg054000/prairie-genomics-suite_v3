# Prairie Genomics Suite - R Shiny Implementation
# Enhanced interactive genomics analysis platform
# 
# Author: Prairie Genomics Team
# Based on existing Python implementation with R Shiny for better performance

# Load required packages for shinyapps.io deployment
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(readr)

# Optional packages with graceful handling
tryCatch(library(shinyWidgets), error = function(e) NULL)
tryCatch(library(readxl), error = function(e) NULL)
tryCatch(library(RColorBrewer), error = function(e) NULL)
tryCatch(library(pheatmap), error = function(e) NULL)
tryCatch(library(ggrepel), error = function(e) NULL)

# Handle DESeq2 separately (may not be available on shinyapps.io)
bioc_available <- FALSE
tryCatch({
  library(DESeq2)
  bioc_available <- TRUE
}, error = function(e) {
  cat("DESeq2 not available - some features will be limited\n")
})

# Handle biomaRt separately (may not be available on shinyapps.io)
biomart_available <- FALSE
tryCatch({
  library(biomaRt)
  biomart_available <- TRUE
}, error = function(e) {
  cat("biomaRt not available - gene symbol conversion will be skipped\n")
})

excel_support <- requireNamespace("readxl", quietly = TRUE)
enhanced_plots <- requireNamespace("RColorBrewer", quietly = TRUE)
ggrepel_available <- requireNamespace("ggrepel", quietly = TRUE)

# Increase Shiny file upload limit to 500MB (default is 5MB)
options(shiny.maxRequestSize = 500*1024^2)

# No additional package loading needed - already handled above

# Initialize comprehensive logging system FIRST
tryCatch({
  source("logging_system.R")
  log_info("APP_INIT", "Prairie Genomics Suite starting up")
}, error = function(e) {
  cat("Warning: Logging system failed to load:", e$message, "\n")
})

# Source helper functions (from root directory)
tryCatch({
  source("data_upload.R")
  source("sample_annotation.R") 
  source("deseq2_analysis.R")
  source("pathway_analysis.R")
  source("visualization.R")
  
  # Load code visibility modules
  source("code_visibility/code_logger.R")
  source("code_visibility/code_generator.R")
  source("code_visibility/analysis_wrappers.R")
  source("code_visibility/code_display_ui.R")
  source("code_visibility/code_display_server.R")
  
  cat("✅ All modules loaded successfully\n")
  cat("✅ Code visibility system loaded\n")
}, error = function(e) {
  stop(paste("Error loading modules:", e$message))
})

# Load Phase 1 Modern Components
tryCatch({
  source("phase1/components/modern_ui_components.R")
  ui_components <- source("phase1/components/modern_ui_components.R")$value
  cat("✅ Phase 1 modern UI components loaded\n")
}, error = function(e) {
  cat("⚠️ Phase 1 components not available:", e$message, "\n")
  # Create fallback functions if Phase 1 not available
  ui_components <- list(
    modern_card = function(...) wellPanel(...),
    stats_card = function(title, value, ...) valueBox(value, title, ...),
    modern_progress = function(id, label, percentage, ...) {
      div(
        h5(label),
        div(class = "progress", 
            div(class = "progress-bar", style = paste0("width: ", percentage, "%")))
      )
    }
  )
})

# Enhanced data processing functions for large datasets
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
        
        cat("Processing chunk:", start_row, "to", end_row, "\n")
        
        # Extract chunk
        chunk_data <- expression_data[start_row:end_row, , drop = FALSE]
        chunk_genes <- gene_names[start_row:end_row]
        
        # Convert to matrix and ensure numeric
        chunk_matrix <- tryCatch({
          matrix(as.numeric(as.matrix(chunk_data)), 
                 nrow = nrow(chunk_data), 
                 ncol = ncol(chunk_data))
        }, error = function(e) {
          stop(paste("Non-numeric data found in chunk", start_row, "to", end_row))
        })
        
        # Set names
        rownames(chunk_matrix) <- chunk_genes
        colnames(chunk_matrix) <- colnames(expression_data)
        
        # Remove genes with all zeros in this chunk
        chunk_sums <- rowSums(chunk_matrix, na.rm = TRUE)
        keep_genes <- chunk_sums > 0
        
        if (sum(keep_genes) > 0) {
          chunk_matrix <- chunk_matrix[keep_genes, , drop = FALSE]
          
          # Combine with result
          if (is.null(result_matrix)) {
            result_matrix <- chunk_matrix
          } else {
            result_matrix <- rbind(result_matrix, chunk_matrix)
          }
        }
        
        # Garbage collection after each chunk
        gc()
      }
      
    } else {
      cat("Small dataset. Processing all at once.\n")
      
      # For small datasets, process normally
      expression_matrix <- tryCatch({
        matrix(as.numeric(as.matrix(expression_data)), 
               nrow = nrow(expression_data), 
               ncol = ncol(expression_data))
      }, error = function(e) {
        stop("Expression data must be numeric")
      })
      
      # Set names
      rownames(expression_matrix) <- gene_names
      colnames(expression_matrix) <- colnames(expression_data)
      
      # Remove genes with all zeros
      keep_genes <- rowSums(expression_matrix, na.rm = TRUE) > 0
      result_matrix <- expression_matrix[keep_genes, , drop = FALSE]
    }
    
    # Final memory cleanup
    gc()
    
    # Final validation
    if (is.null(result_matrix) || nrow(result_matrix) == 0) {
      stop("No valid genes remaining after processing")
    }
    
    if (ncol(result_matrix) < 2) {
      stop("Insufficient samples remaining after processing")
    }
    
    # Monitor memory after processing
    mem_after <- monitor_memory()
    cat("Memory after processing:", mem_after$message, "\n")
    
    cat("Processing completed successfully:", nrow(result_matrix), "genes retained\n")
    
    return(result_matrix)
    
  }, error = function(e) {
    cat("Error in process_large_dataset:", e$message, "\n")
    stop(paste("Data processing failed:", e$message))
  })
}

# Define UI
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
        tags$span("Phase 1 Enhanced", style = "color: #ffffff; font-size: 12px;")
      )
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("📁 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("🧬 Sample Annotation", tabName = "annotation", icon = icon("tags")),
      menuItem("🚀 DESeq2 Analysis", tabName = "analysis", icon = icon("chart-line")),
      menuItem("🧬 Pathway Analysis", tabName = "pathways", icon = icon("project-diagram")),
      menuItem("📊 Visualizations", tabName = "visualizations", icon = icon("chart-bar")),
      menuItem("📜 Code View", tabName = "code_view", icon = icon("code")),
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
        progressBar(
          id = "overall_progress",
          value = 0,
          total = 100,
          status = "primary",
          display_pct = TRUE,
          striped = TRUE
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
    # Add Phase 1 modern styling assets
    tags$head(
      # Include modern CSS components
      includeCSS("www/css/modern_components.css"),
      # Include modern JavaScript interactions
      includeScript("www/js/modern_interactions.js"),
      # Custom CSS including Phase 1 integration
      tags$style(HTML("
        /* Phase 1 Enhanced Styling */
        .content-wrapper, .right-side {
          background-color: var(--gray-50);
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
            
            # Debug/Clear data controls
            fluidRow(
              column(
                12,
                div(
                  style = "margin-bottom: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                  h6("🔧 Debug Controls", style = "margin-top: 0;"),
                  actionButton(
                    "clear_all_data",
                    "🗑️ Clear All Data",
                    class = "btn-warning btn-sm",
                    style = "margin-right: 10px;"
                  ),
                  actionButton(
                    "debug_current_data",
                    "🔍 Debug Current Data",
                    class = "btn-info btn-sm"
                  ),
                  br(),
                  tags$small("Use these controls if you see unexpected gene counts or cached data"),
                  br(),
                  div(
                    class = "alert alert-warning",
                    style = "margin-top: 10px; font-size: 12px;",
                    h6("🚨 15,357 Gene Issue Fix:"),
                    p("If you see 15,357 genes instead of your uploaded data:"),
                    tags$ol(
                      tags$li("Click 'Clear All Data' button above"),
                      tags$li("Restart your R session completely"),
                      tags$li("Re-upload your data file"),
                      tags$li("The app now blocks this problematic cached dataset")
                    )
                  )
                )
              )
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
              deseq2AnalysisUI("deseq2_analysis"),
              
              br(),
              
              # Inline code display for DESeq2 analysis
              div(
                style = "margin-top: 20px;",
                inlineCodeUI("deseq2_code", "View DESeq2 Analysis Code", "btn-outline-primary btn-sm")
              )
            )
          )
        )
      ),
      
      # Pathway Analysis Tab
      tabItem(
        tabName = "pathways",
        fluidRow(
          box(
            title = "🧬 Pathway Analysis (GO, KEGG, GSEA)", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_results == false",
              div(
                class = "alert alert-warning",
                h4("Run DESeq2 Analysis First"),
                p("Please complete your DESeq2 analysis before running pathway analysis.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_results == true",
              
              # Analysis Controls
              fluidRow(
                column(
                  4,
                  wellPanel(
                    h5("🎯 Analysis Settings"),
                    
                    # Cache status indicator
                    div(
                      id = "cache_status",
                      style = "margin-bottom: 10px; font-size: 12px;",
                      uiOutput("cache_status_text")
                    ),
                    
                    selectInput(
                      "pathway_analysis_type",
                      "Analysis Type:",
                      choices = list(
                        "Gene Ontology (GO)" = "GO",
                        "KEGG Pathways" = "KEGG", 
                        "Gene Set Enrichment (GSEA)" = "GSEA",
                        "MSigDB Collections" = "MSigDB",
                        "Reactome Pathways" = "Reactome"
                      ),
                      selected = "GO"
                    ),
                    
                    selectInput(
                      "pathway_species",
                      "Species:",
                      choices = list(
                        "Auto-detect" = "auto",
                        "Human (Homo sapiens)" = "human",
                        "Mouse (Mus musculus)" = "mouse",
                        "Rat (Rattus norvegicus)" = "rat",
                        "Zebrafish (Danio rerio)" = "zebrafish",
                        "Fruit fly (Drosophila)" = "fly",
                        "Worm (C. elegans)" = "worm",
                        "Yeast (S. cerevisiae)" = "yeast",
                        "Arabidopsis (A. thaliana)" = "arabidopsis"
                      ),
                      selected = "auto"
                    ),
                    
                    # Species availability indicator
                    div(
                      id = "species_availability",
                      style = "font-size: 11px; color: #6c757d; margin-top: 5px;",
                      htmlOutput("species_availability_text")
                    ),
                    
                    conditionalPanel(
                      condition = "input.pathway_analysis_type == 'GO'",
                      selectInput(
                        "go_ontology",
                        "GO Ontology:",
                        choices = list(
                          "Biological Process" = "BP",
                          "Molecular Function" = "MF", 
                          "Cellular Component" = "CC"
                        ),
                        selected = "BP"
                      )
                    ),
                    
                    conditionalPanel(
                      condition = "input.pathway_analysis_type == 'MSigDB' || input.pathway_analysis_type == 'GSEA'",
                      selectInput(
                        "msigdb_collection",
                        "Gene Set Collection:",
                        choices = list(
                          "Hallmark" = "H",
                          "Curated (C2)" = "C2",
                          "Ontology (C5)" = "C5",
                          "Immunologic (C7)" = "C7"
                        ),
                        selected = "H"
                      )
                    ),
                    
                    br(),
                    
                    h5("📊 Filter Settings"),
                    
                    numericInput(
                      "pathway_padj_cutoff",
                      "Adjusted p-value cutoff:",
                      value = 0.05,
                      min = 0.001,
                      max = 0.1,
                      step = 0.001
                    ),
                    
                    numericInput(
                      "pathway_fc_cutoff", 
                      "Log2 fold change cutoff:",
                      value = 1.0,
                      min = 0.1,
                      max = 5.0,
                      step = 0.1
                    ),
                    
                    br(),
                    
                    # Gene filtering preview
                    div(
                      id = "gene_filter_preview",
                      h5("📊 Gene Filter Preview"),
                      div(
                        style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                        htmlOutput("pathway_gene_preview")
                      )
                    ),
                    
                    actionButton(
                      "run_pathway_analysis",
                      "🚀 Run Pathway Analysis",
                      class = "btn-primary btn-lg",
                      style = "width: 100%;"
                    )
                  )
                ),
                
                column(
                  8,
                  # Results Display
                  conditionalPanel(
                    condition = "output.show_pathway_results == true",
                    
                    tabsetPanel(
                      id = "pathway_results_tabs",
                      
                      tabPanel(
                        "📊 Results Table",
                        br(),
                        DT::dataTableOutput("pathway_results_table")
                      ),
                      
                      tabPanel(
                        "📈 Visualizations", 
                        br(),
                        fluidRow(
                          column(
                            3,
                            # Dynamic plot type selection based on analysis type
                            conditionalPanel(
                              condition = "output.show_pathway_results == true",
                              uiOutput("dynamic_plot_type_selector")
                            ),
                            
                            numericInput(
                              "pathway_plot_top_n",
                              "Show top N pathways:",
                              value = 20,
                              min = 5,
                              max = 50,
                              step = 5
                            ),
                            
                            br(),
                            
                            downloadButton(
                              "download_pathway_plot",
                              "Download Plot",
                              class = "btn-secondary",
                              style = "width: 100%;"
                            )
                          ),
                          
                          column(
                            9,
                            plotOutput("pathway_plot", height = "600px")
                          )
                        )
                      ),
                      
                      tabPanel(
                        "📋 Summary",
                        br(),
                        verbatimTextOutput("pathway_analysis_summary")
                      ),
                      
                      tabPanel(
                        "📜 Analysis Code",
                        br(),
                        # Inline code display for pathway analysis
                        inlineCodeUI("pathway_code", "Pathway Analysis Code", "btn-outline-success btn-sm"),
                        
                        br(),
                        
                        div(
                          class = "alert alert-info",
                          style = "margin-top: 15px;",
                          h5("📋 Pathway Analysis Code"),
                          p("This tab shows the exact R code used for your pathway analysis, including gene set preparation, statistical testing, and result formatting."),
                          p("You can copy this code to reproduce the analysis independently or customize it for your specific needs."),
                          br(),
                          div(
                            class = "alert alert-success",
                            style = "margin-top: 10px; border-left: 4px solid #28a745;",
                            h6("🔧 Key Analysis Parameters (Enhanced for Accuracy):"),
                            tags$ul(
                              tags$li(HTML("<strong>pvalueCutoff = 0.01</strong> - Only pathways with p-value < 1% are considered significant")),
                              tags$li(HTML("<strong>qvalueCutoff = 0.05</strong> - Only pathways with FDR < 5% after multiple testing correction")),
                              tags$li(HTML("<strong>minGSSize = 10</strong> - Only pathways with at least 10 genes for meaningful results")),
                              tags$li(HTML("<strong>Tiered Gene Selection</strong> - Prioritizes highly significant genes (padj<0.05, |FC|>1) first"))
                            ),
                            p(HTML("<em>These stringent parameters ensure you get only the most biologically meaningful pathways, avoiding the thousands of false positives from lenient settings.</em>"))
                          )
                        )
                      )
                    )
                  ),
                  
                  # Status/Progress Display
                  conditionalPanel(
                    condition = "output.show_pathway_results == false",
                    wellPanel(
                      h4("Ready for Pathway Analysis"),
                      p("Configure your analysis settings and click 'Run Pathway Analysis' to begin."),
                      p("This analysis will identify enriched biological pathways in your significant genes."),
                      
                      h5("Analysis Types:"),
                      tags$ul(
                        tags$li("🧬 GO: Gene Ontology biological processes, molecular functions, cellular components"),
                        tags$li("🛤️  KEGG: Metabolic and signaling pathways"), 
                        tags$li("📊 GSEA: Gene Set Enrichment Analysis using ranked gene lists"),
                        tags$li("📚 MSigDB: Curated gene sets from molecular signatures database")
                      )
                    )
                  )
                )
              )
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
      
      # Code View Tab
      tabItem(
        tabName = "code_view",
        fluidRow(
          box(
            title = "📜 Code View & Export", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            div(
              class = "alert alert-info",
              h4("Scientific Reproducibility & Code Transparency"),
              p("View and download the complete R code for your analysis. This ensures full reproducibility and allows you to run the analysis independently."),
              p("All analysis steps are automatically logged with proper package loading, parameter settings, and session information.")
            ),
            
            # Main code display interface
            codeDisplayUI("main_code_display", title = "Complete Analysis Code", height = "500px"),
            
            br(),
            
            # Code export panel
            fluidRow(
              column(
                6,
                codeExportUI("code_export")
              ),
              column(
                6,
                # Code validation panel
                codeValidationUI("code_validation"),
                
                br(),
                
                # Method documentation
                methodDocumentationUI("method_docs")
              )
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

# Define Server
server <- function(input, output, session) {
  # Reactive values to store data
  values <- reactiveValues(
    expression_data = NULL,
    annotation_data = NULL,
    deseq2_results = NULL,
    deseq2_data = NULL,  # NEW: Unified DESeq2 data structure
    filtered_data = NULL,
    pathway_results = NULL,  # Current displayed results
    pathway_results_cache = list(),  # CACHE: Store results for each analysis type
    progress = 0,
    code_session_id = NULL
  )
  
  # DEBUG: Monitor all changes to expression_data
  observeEvent(values$expression_data, {
    if (!is.null(values$expression_data)) {
      cat("🔍 DEBUG: values$expression_data changed\n")
      cat("   - New dimensions:", nrow(values$expression_data), "x", ncol(values$expression_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(values$expression_data), 3), collapse = ", "), "\n")
      
      # CRITICAL: Check for the mysterious 15,357 dataset
      if (nrow(values$expression_data) == 15357) {
        cat("🚨 CRITICAL ALERT: 15,357 gene dataset detected!\n")
        cat("   - This is the known problematic dataset that overrides user uploads\n")
        cat("   - First 10 gene IDs:", paste(head(rownames(values$expression_data), 10), collapse = ", "), "\n")
        cat("   - Data class:", class(values$expression_data), "\n")
        cat("   - Attributes:", paste(names(attributes(values$expression_data)), collapse = ", "), "\n")
      }
      
      cat("   - Stack trace:\n")
      traceback()
    } else {
      cat("🔍 DEBUG: values$expression_data set to NULL\n")
    }
  }, ignoreNULL = FALSE)
  
  # Clear all data function
  clear_all_data <- function() {
    cat("🗑️ Clearing all session data\n")
    values$expression_data <- NULL
    values$annotation_data <- NULL
    values$deseq2_results <- NULL
    values$deseq2_data <- NULL  # Clear unified structure
    values$filtered_data <- NULL
    values$pathway_results <- NULL
    values$pathway_results_cache <- list()
    values$progress <- 0
    
    # CRITICAL: Clear gene conversion cache to eliminate 15,357 dataset issue
    tryCatch({
      if (exists("clear_gene_cache")) {
        clear_gene_cache()
        cat("✅ Gene conversion cache cleared\n")
      }
    }, error = function(e) {
      cat("⚠️ Could not clear gene conversion cache:", e$message, "\n")
    })
    
    # Force garbage collection to free memory
    gc()
    cat("🧹 Memory cleaned up\n")
  }
  
  # Initialize code logging session
  observe({
    if (is.null(values$code_session_id)) {
      values$code_session_id <- init_code_logger(user_name = "prairie_user")
      cat("✅ Code logging session initialized:", values$code_session_id, "\n")
    }
  })
  
  # Debug button handlers
  observeEvent(input$clear_all_data, {
    clear_all_data()
    showNotification("🗑️ All session data cleared", type = "message", duration = 3)
  })
  
  observeEvent(input$debug_current_data, {
    if (!is.null(values$expression_data)) {
      cat("🔍 DEBUG: Current expression data status:\n")
      cat("   - Dimensions:", nrow(values$expression_data), "x", ncol(values$expression_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(values$expression_data), 5), collapse = ", "), "\n")
      cat("   - Gene symbols available:", !is.null(attr(values$expression_data, "gene_symbols")), "\n")
      if (!is.null(attr(values$expression_data, "gene_symbols"))) {
        cat("   - Conversion rate:", attr(values$expression_data, "conversion_rate") %||% "unknown", "%\n")
        cat("   - Species:", attr(values$expression_data, "conversion_species") %||% "unknown", "\n")
      }
      
      showNotification(
        paste0("📊 Current data: ", nrow(values$expression_data), " genes × ", ncol(values$expression_data), " samples. Check console for details."),
        type = "message",
        duration = 5
      )
    } else {
      cat("🔍 DEBUG: No expression data loaded\n")
      showNotification("ℹ️ No expression data currently loaded", type = "message", duration = 3)
    }
  })
  
  # PERSISTENT RESULTS: Load cached results when analysis type changes
  observeEvent(input$pathway_analysis_type, {
    if (!is.null(input$pathway_analysis_type)) {
      cache_key <- paste0(input$pathway_analysis_type, "_", input$pathway_species %||% "auto")
      
      if (cache_key %in% names(values$pathway_results_cache)) {
        cat("📂 Loading cached results for", input$pathway_analysis_type, "\n")
        values$pathway_results <- values$pathway_results_cache[[cache_key]]
        values$pathway_results_updated <- Sys.time()
        
        showNotification(
          paste0("📂 Loaded cached ", input$pathway_analysis_type, " results"),
          type = "message",
          duration = 3
        )
      } else {
        cat("💡 No cached results for", input$pathway_analysis_type, "- ready for new analysis\n")
        # Clear current results when switching to uncached analysis type
        if (!is.null(values$pathway_results)) {
          values$pathway_results <- NULL
        }
      }
    }
  }, ignoreInit = TRUE)
  
  # Data Upload Module
  upload_results <- callModule(dataUpload, "data_upload", values)
  
  # Extract results from data upload module
  observe({
    if (!is.null(upload_results$expression_data())) {
      # DEBUG: Check what data is being stored in reactive values
      new_data <- upload_results$expression_data()
      cat("🔍 DEBUG: App.R - Storing expression data in reactive values\n")
      cat("   - Data dimensions:", nrow(new_data), "x", ncol(new_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(new_data), 3), collapse = ", "), "\n")
      
      values$expression_data <- new_data
      
      # Log data upload step
      if (!is.null(values$code_session_id) && !is.null(values$expression_data)) {
        file_info <- list(
          name = "expression_data.csv",
          size = paste0(object.size(values$expression_data), " bytes")
        )
        processing_info <- list(
          steps = c("remove_zero_genes", "ensure_integer"),
          final_dimensions = paste(nrow(values$expression_data), "genes x", ncol(values$expression_data), "samples")
        )
        
        log_data_upload(file_info, processing_info, values$code_session_id)
      }
    }
    if (!is.null(upload_results$annotation_data())) {
      values$panno_annotation <- upload_results$annotation_data()
    }
  })
  
  # Sample Annotation Module  
  annotation_data <- callModule(sampleAnnotation, "sample_annotation", values)
  
  # DESeq2 Analysis Module
  analysis_results <- callModule(deseq2Analysis, "deseq2_analysis", values)
  
  # Fixed Pathway Analysis Server Logic
  observeEvent(input$run_pathway_analysis, {
    # Log user action with comprehensive details
    log_user_action("PATHWAY_ANALYSIS_CLICKED", list(
      analysis_type = input$pathway_analysis_type,
      species = input$pathway_species, 
      deseq2_available = !is.null(values$deseq2_results),
      deseq2_genes = if (!is.null(values$deseq2_results)) nrow(values$deseq2_results) else 0
    ))
    
    cat("🚀 Pathway analysis button clicked\n")
    
    # Validate prerequisites with logging
    if (is.null(values$deseq2_results)) {
      log_error("PATHWAY_ANALYSIS", "No DESeq2 results available for pathway analysis", list(
        user_attempted_analysis = input$pathway_analysis_type,
        session_state = "missing_deseq2_results"
      ))
      showNotification("❌ No DESeq2 results available", type = "error")
      return()
    }
    
    log_info("PATHWAY_ANALYSIS", "Starting pathway analysis", list(
      analysis_type = input$pathway_analysis_type,
      species = input$pathway_species,
      input_genes = nrow(values$deseq2_results)
    ))
    
    cat("✅ DESeq2 results available, proceeding with pathway analysis\n")
    
    # Show progress notification
    showNotification("🔄 Running pathway analysis...", type = "message", duration = NULL, id = "pathway_progress")
    
    # Run pathway analysis
    tryCatch({
      start_time <- Sys.time()
      
      pathway_result <- run_pathway_analysis(
        deseq2_results = values$deseq2_results,
        analysis_type = input$pathway_analysis_type,
        species = input$pathway_species %||% "auto",
        ontology = if (input$pathway_analysis_type == "GO") input$go_ontology else "BP",
        padj_cutoff = input$pathway_padj_cutoff,
        fc_cutoff = input$pathway_fc_cutoff,
        gene_set_collection = if (input$pathway_analysis_type %in% c("GSEA", "MSigDB")) input$msigdb_collection else "H",
        deseq2_data = values$deseq2_data  # Pass unified structure for optimization
      )
      
      # Log successful analysis
      end_time <- Sys.time()
      duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      log_analysis_performance(
        analysis_type = input$pathway_analysis_type,
        duration_seconds = duration,
        input_genes = rownames(values$deseq2_results),
        output_pathways = if (pathway_result$success && !is.null(pathway_result$data)) nrow(pathway_result$data) else 0
      )
      
      # Debug output
      cat("📊 Pathway analysis completed. Success:", pathway_result$success, "\n")
      if (pathway_result$success && !is.null(pathway_result$data)) {
        cat("📈 Found", nrow(pathway_result$data), "enriched pathways\n")
      }
      
      # Store results - THIS IS CRITICAL
      values$pathway_results <- pathway_result
      
      # PERSISTENT RESULTS: Cache results for this analysis type
      cache_key <- paste0(input$pathway_analysis_type, "_", input$pathway_species %||% "auto")
      values$pathway_results_cache[[cache_key]] <- pathway_result
      cat("💾 Cached results for", cache_key, "\n")
      
      # Clean up memory after storing results
      gc()
      
      # Force reactive update by creating a trigger
      values$pathway_results_updated <- Sys.time()
      
      # Remove progress notification
      removeNotification("pathway_progress")
      
      # Show completion notification
      if (pathway_result$success) {
        n_pathways <- pathway_result$n_terms %||% pathway_result$n_pathways %||% pathway_result$n_gene_sets %||% 0
        if (!is.null(pathway_result$data)) {
          n_pathways <- nrow(pathway_result$data)
        }
        
        showNotification(
          paste0("✅ ", pathway_result$analysis_type, " analysis completed! Found ", n_pathways, " enriched pathways"),
          type = "message",
          duration = 5
        )
        
        cat("✅ Pathway results stored. UI should now update.\n")
      } else {
        showNotification(
          paste0("❌ Pathway analysis failed: ", pathway_result$message %||% "Unknown error"),
          type = "error",
          duration = 8
        )
      }
      
    }, error = function(e) {
      # Log the error with full context
      log_error("PATHWAY_ANALYSIS", "Pathway analysis failed", list(
        error_message = e$message,
        analysis_type = input$pathway_analysis_type,
        species = input$pathway_species,
        input_parameters = list(
          padj_cutoff = input$pathway_padj_cutoff,
          fc_cutoff = input$pathway_fc_cutoff,
          ontology = if (input$pathway_analysis_type == "GO") input$go_ontology else "BP"
        ),
        deseq2_results_available = !is.null(values$deseq2_results),
        deseq2_genes = if (!is.null(values$deseq2_results)) nrow(values$deseq2_results) else 0
      ))
      
      removeNotification("pathway_progress")
      cat("❌ Pathway analysis error:", e$message, "\n")
      
      # Store error result
      values$pathway_results <- list(
        success = FALSE,
        error = e$message,
        message = paste("Analysis failed:", e$message)
      )
      
      showNotification(
        paste0("❌ Pathway analysis error: ", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Species availability indicator
  output$species_availability_text <- renderText({
    selected_species <- input$pathway_species %||% "auto"
    
    if (selected_species == "auto") {
      "📋 Will detect species automatically from gene IDs"
    } else {
      # Check if selected species is supported
      if (exists("get_supported_species") && is.function(get_supported_species)) {
        supported_species <- get_supported_species()
        
        if (selected_species %in% names(supported_species)) {
          species_info <- supported_species[[selected_species]]
          availability_icon <- if (species_info$available) "✅" else "⚠️"
          availability_text <- if (species_info$available) "Fully supported" else "Limited support (core packages missing)"
          
          paste0(availability_icon, " ", species_info$name, " - ", availability_text)
        } else {
          "❓ Species information not available"
        }
      } else {
        "❓ Cannot check species availability"
      }
    }
  })
  
  # Real-time gene filtering preview
  output$pathway_gene_preview <- renderText({
    if (!is.null(values$deseq2_results)) {
      preview <- preview_gene_filtering(
        values$deseq2_results, 
        padj_cutoff = input$pathway_padj_cutoff %||% 0.05,
        fc_cutoff = input$pathway_fc_cutoff %||% 1.0
      )
      
      status_color <- if (preview$recommended) "success" else if (preview$significant_genes == 0) "danger" else "warning"
      status_icon <- if (preview$recommended) "✅" else if (preview$significant_genes == 0) "❌" else "⚠️"
      
      paste0(
        "<strong>", status_icon, " ", preview$significant_genes, "</strong> genes pass filters<br>",
        "<small>",
        "📈 Upregulated: ", preview$upregulated, " | ",
        "📉 Downregulated: ", preview$downregulated, "<br>",
        "📊 Filter rate: ", preview$filter_rate, "% of ", formatC(preview$total_genes, format="d", big.mark=","), " total genes",
        if (!preview$recommended && preview$significant_genes > 0) {
          "<br><span style='color: orange;'>💡 Consider adjusting thresholds for optimal results</span>"
        } else if (preview$significant_genes == 0) {
          "<br><span style='color: red;'>⚠️ No genes pass current filters - try relaxing thresholds</span>"
        } else {
          "<br><span style='color: green;'>✨ Good gene count for pathway analysis</span>"
        },
        "</small>"
      )
    } else {
      "<em>Run DESeq2 analysis first to see gene filtering preview</em>"
    }
  })
  
  # ENHANCED Pathway results table with gene details and clear explanations
  output$pathway_results_table <- DT::renderDataTable({
    req(values$pathway_results)
    
    if (values$pathway_results$success) {
      pathway_data <- values$pathway_results$data
      
      # Create enhanced display data with user-friendly columns
      display_data <- data.frame(
        "Pathway/Gene Set" = pathway_data$Description,
        stringsAsFactors = FALSE
      )
      
      # Add gene information with clear formatting
      if ("geneID" %in% colnames(pathway_data)) {
        # Format gene lists nicely
        gene_lists <- sapply(pathway_data$geneID, function(genes) {
          if (is.na(genes) || genes == "") return("No genes")
          gene_vector <- unlist(strsplit(as.character(genes), "/"))
          if (length(gene_vector) > 8) {
            paste0(paste(gene_vector[1:8], collapse = ", "), ", ... (", length(gene_vector), " total)")
          } else {
            paste(gene_vector, collapse = ", ")
          }
        })
        display_data$"Enriched Genes" <- gene_lists
      }
      
      # Add user-friendly ratio explanations
      if ("GeneRatio" %in% colnames(pathway_data)) {
        # Parse GeneRatio (e.g., "3/400" -> "3 out of 400 input genes")
        gene_ratios <- sapply(pathway_data$GeneRatio, function(ratio) {
          if (is.na(ratio)) return("N/A")
          parts <- unlist(strsplit(as.character(ratio), "/"))
          if (length(parts) == 2) {
            paste0(parts[1], " out of ", parts[2], " input genes")
          } else {
            as.character(ratio)
          }
        })
        display_data$"Your Genes in Pathway" <- gene_ratios
      }
      
      if ("BgRatio" %in% colnames(pathway_data)) {
        # Parse BgRatio (e.g., "45/18000" -> "45 out of 18,000 total database genes")
        bg_ratios <- sapply(pathway_data$BgRatio, function(ratio) {
          if (is.na(ratio)) return("N/A")
          parts <- unlist(strsplit(as.character(ratio), "/"))
          if (length(parts) == 2) {
            total_formatted <- format(as.numeric(parts[2]), big.mark = ",")
            paste0(parts[1], " out of ", total_formatted, " database genes")
          } else {
            as.character(ratio)
          }
        })
        display_data$"Total Genes in Pathway" <- bg_ratios
      }
      
      # Add statistical significance with clear formatting
      if ("pvalue" %in% colnames(pathway_data)) {
        p_values <- ifelse(pathway_data$pvalue < 0.001, 
                          "< 0.001", 
                          formatC(pathway_data$pvalue, format = "f", digits = 4))
        display_data$"P-value" <- p_values
      }
      
      if ("p.adjust" %in% colnames(pathway_data)) {
        adj_p_values <- ifelse(pathway_data$p.adjust < 0.001, 
                              "< 0.001", 
                              formatC(pathway_data$p.adjust, format = "f", digits = 4))
        display_data$"Adjusted P-value (FDR)" <- adj_p_values
      }
      
      # Add enrichment strength if available
      if ("Count" %in% colnames(pathway_data)) {
        display_data$"Gene Count" <- pathway_data$Count
      }
      
      if ("pvalue" %in% colnames(pathway_data)) {
        # Add significance stars for easy interpretation
        stars <- sapply(pathway_data$pvalue, function(p) {
          if (is.na(p)) return("")
          if (p < 0.001) return("***")
          if (p < 0.01) return("**") 
          if (p < 0.05) return("*")
          return("")
        })
        display_data$"Significance" <- stars
      }
      
      # Create enhanced datatable with explanatory caption
      dt <- DT::datatable(
        display_data,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'tip',
          columnDefs = list(
            list(width = '200px', targets = c(0)),  # Pathway column
            list(width = '300px', targets = which(colnames(display_data) == "Enriched Genes") - 1),  # Gene list column
            list(className = 'dt-center', targets = c(2, 3, 4, 5, 6))  # Center align stats
          )
        ),
        rownames = FALSE,
        caption = HTML(paste(
          "<div style='margin: 10px 0; padding: 10px; background-color: #f8f9fa; border-left: 4px solid #007bff;'>",
          "<h5 style='margin-top: 0; color: #007bff;'>📊 How to Read This Table:</h5>",
          "<p><strong>Your Genes in Pathway:</strong> Number of your significant genes found in this pathway out of total input genes</p>",
          "<p><strong>Total Genes in Pathway:</strong> Total genes known to be in this pathway in the database</p>", 
          "<p><strong>Adjusted P-value (FDR):</strong> Statistical significance after multiple testing correction (< 0.05 = significant)</p>",
          "<p><strong>Significance:</strong> *** p < 0.001, ** p < 0.01, * p < 0.05</p>",
          "</div>"
        ))
      )
      
      return(dt)
    }
  })
  
  # Dynamic plot type selector based on analysis type
  output$dynamic_plot_type_selector <- renderUI({
    req(values$pathway_results)
    
    if (values$pathway_results$success) {
      analysis_type <- values$pathway_results$analysis_type
      
      if (analysis_type == "GSEA") {
        # GSEA-specific plot options
        tagList(
          selectInput(
            "pathway_plot_type",
            "Plot Type:",
            choices = list(
              "GSEA Enrichment Plot (NES)" = "gsea",
              "Dot Plot (NES vs Significance)" = "dotplot",
              "Bar Plot" = "barplot"
            ),
            selected = "gsea"
          ),
          div(
            class = "alert alert-info",
            style = "margin-top: 10px; padding: 8px; font-size: 11px;",
            HTML("<strong>💡 GSEA Plot Guide:</strong><br/>
                  <em>Enrichment Plot</em>: Shows pathway scores (NES) with up/down regulation<br/>
                  <em>Dot Plot</em>: NES vs significance with size indicating strength")
          )
        )
      } else {
        # ORA analysis plot options (GO, KEGG, MSigDB, Reactome)
        selectInput(
          "pathway_plot_type", 
          "Plot Type:",
          choices = list(
            "Dot Plot (recommended)" = "dotplot",
            "Bar Plot" = "barplot",
            "Network Plot" = "network"
          ),
          selected = "dotplot"
        )
      }
    } else {
      # Default when no results
      selectInput(
        "pathway_plot_type",
        "Plot Type:",
        choices = list("Dot Plot" = "dotplot"),
        selected = "dotplot"
      )
    }
  })
  
  # Pathway analysis plot - SIMPLIFIED FOR SPEED
  output$pathway_plot <- renderPlot({
    req(values$pathway_results)
    
    if (values$pathway_results$success) {
      # ULTRA-FAST plotting with aggressive simplification
      plot_obj <- tryCatch({
        # Force minimal data for plotting (max 10 pathways)
        limited_results <- values$pathway_results
        if (!is.null(limited_results$data) && nrow(limited_results$data) > 10) {
          limited_results$data <- limited_results$data[1:10, ]
          cat("🔧 LIMITED to 10 pathways for fast rendering\n")
        }
        
        # Use user-selected plot type from dropdown
        selected_plot_type <- input$pathway_plot_type %||% "dotplot"
        selected_top_n <- min(input$pathway_plot_top_n %||% 10, 20)  # Limit for performance
        
        cat("🎨 Creating", selected_plot_type, "plot with", selected_top_n, "pathways\n")
        
        create_fast_pathway_plot(
          limited_results, 
          plot_type = selected_plot_type,
          top_n = selected_top_n
        )
      }, error = function(e) {
        cat("⚠️ Plot creation failed, using emergency fallback:", e$message, "\n")
        create_emergency_plot(values$pathway_results, 10)
      })
      
      if (!is.null(plot_obj)) {
        print(plot_obj)
      } else {
        # Emergency fallback plot
        plot.new()
        text(0.5, 0.5, paste("Pathway Analysis Complete\n", 
                             nrow(values$pathway_results$data), "pathways found\n",
                             "Visualization in progress..."), 
             cex = 1.2, col = "darkblue", font = 2)
      }
    }
  })
  
  # Enhanced pathway analysis summary with integration
  output$pathway_analysis_summary <- renderText({
    req(values$pathway_results)
    
    if (values$pathway_results$success) {
      # Create integrated summary if DESeq2 results are available
      if (!is.null(values$deseq2_results)) {
        integrated_summary <- create_integrated_summary(values$deseq2_results, values$pathway_results)
        if (!is.null(integrated_summary)) {
          return(integrated_summary)
        }
      }
      
      # Fallback to basic summary
      result <- values$pathway_results
      summary_text <- paste0(
        "Analysis Type: ", result$analysis_type, "\n",
        "Species: ", result$species, "\n",
        if (!is.null(result$ontology)) paste0("Ontology: ", result$ontology, "\n") else "",
        if (!is.null(result$gene_set_collection)) paste0("Gene Set Collection: ", result$gene_set_collection, "\n") else "",
        "Total Pathways Found: ", nrow(result$data), "\n",
        "Significant Pathways (p < 0.05): ", sum(result$data$pvalue < 0.05, na.rm = TRUE), "\n",
        "Analysis Status: Success"
      )
    } else {
      summary_text <- paste0(
        "Analysis Status: Failed\n",
        "Error: ", values$pathway_results$error, "\n",
        "Message: ", values$pathway_results$message
      )
    }
    
    return(summary_text)
  })
  
  # Enhanced show_pathway_results reactive with debugging
  output$show_pathway_results <- reactive({
    result <- !is.null(values$pathway_results) && 
              ("success" %in% names(values$pathway_results)) && 
              values$pathway_results$success
    
    # Debug output
    cat("🔍 show_pathway_results reactive triggered. Result:", result, "\n")
    if (!is.null(values$pathway_results)) {
      cat("   - pathway_results exists\n")
      cat("   - success field:", values$pathway_results$success %||% "MISSING", "\n")
    } else {
      cat("   - pathway_results is NULL\n")
    }
    
    return(result)
  })
  outputOptions(output, "show_pathway_results", suspendWhenHidden = FALSE)
  
  # Cache status indicator
  output$cache_status_text <- renderUI({
    if (!is.null(input$pathway_analysis_type)) {
      cache_key <- paste0(input$pathway_analysis_type, "_", input$pathway_species %||% "auto")
      
      if (cache_key %in% names(values$pathway_results_cache)) {
        cached_result <- values$pathway_results_cache[[cache_key]]
        if (cached_result$success) {
          n_pathways <- nrow(cached_result$data)
          HTML(paste0(
            "<span style='color: #28a745;'>",
            "📂 Cached: ", n_pathways, " pathways from previous analysis",
            "</span>"
          ))
        } else {
          HTML("<span style='color: #ffc107;'>⚠️ Previous analysis failed - try again</span>")
        }
      } else {
        HTML("<span style='color: #6c757d;'>💡 No cached results - run analysis to save results</span>")
      }
    } else {
      HTML("")
    }
  })
  
  # Download pathway plot
  output$download_pathway_plot <- downloadHandler(
    filename = function() {
      paste0(input$pathway_analysis_type, "_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      if (!is.null(values$pathway_results) && values$pathway_results$success) {
        plot_obj <- create_pathway_plots(
          values$pathway_results, 
          plot_type = input$pathway_plot_type %||% "dotplot",
          top_n = input$pathway_plot_top_n %||% 20
        )
        
        if (!is.null(plot_obj)) {
          ggsave(file, plot = plot_obj, width = 10, height = 8, dpi = 300)
        }
      }
    }
  )
  
  # Visualization Module
  callModule(visualization, "visualization", values)
  
  # Code Visibility Modules
  session_id_reactive <- reactive({ values$code_session_id })
  
  # Main code display server
  codeDisplayServer("main_code_display", session_id_reactive)
  
  # Code export server
  codeExportServer("code_export", session_id_reactive)
  
  # Code validation server
  codeValidationServer("code_validation", session_id_reactive)
  
  # Method documentation server
  methodDocumentationServer("method_docs")
  
  # Inline code display servers
  deseq2_code_reactive <- reactive({
    if (!is.null(values$code_session_id)) {
      get_analysis_code_snippet(values$code_session_id, category = "deseq2")
    } else {
      "# No DESeq2 analysis has been run yet"
    }
  })
  
  pathway_code_reactive <- reactive({
    if (!is.null(values$code_session_id)) {
      get_analysis_code_snippet(values$code_session_id, category = "pathway")
    } else {
      "# No pathway analysis has been run yet"
    }
  })
  
  inlineCodeServer("deseq2_code", 
                   code_snippet = deseq2_code_reactive,
                   code_title = reactive("DESeq2 Differential Expression Analysis"),
                   code_description = reactive("Complete R code for running DESeq2 analysis with your specific parameters"))
  
  inlineCodeServer("pathway_code", 
                   code_snippet = pathway_code_reactive,
                   code_title = reactive("Pathway Enrichment Analysis"),
                   code_description = reactive("R code for pathway analysis including gene set preparation and enrichment testing"))
  
  # Enhanced progress tracking with detailed steps
  observe({
    progress <- 0
    progress_details <- list()
    
    # Data upload progress (25%)
    if (!is.null(values$expression_data)) {
      progress <- progress + 25
      progress_details$data_upload <- list(
        status = "completed",
        message = paste(nrow(values$expression_data), "genes,", ncol(values$expression_data), "samples"),
        timestamp = Sys.time()
      )
    } else {
      progress_details$data_upload <- list(status = "pending", message = "Upload expression data")
    }
    
    # Sample annotation progress (25%)
    if (!is.null(values$annotation_data)) {
      progress <- progress + 25
      progress_details$sample_annotation <- list(
        status = "completed",
        message = paste(length(unique(values$annotation_data$Condition)), "groups identified"),
        timestamp = Sys.time()
      )
    } else {
      progress_details$sample_annotation <- list(status = "pending", message = "Annotate samples")
    }
    
    # DESeq2 analysis progress (35%)
    if (!is.null(values$deseq2_results)) {
      progress <- progress + 35
      significant_genes <- sum(values$deseq2_results$padj < 0.05 & abs(values$deseq2_results$log2FoldChange) > 1, na.rm = TRUE)
      progress_details$deseq2_analysis <- list(
        status = "completed",
        message = paste(significant_genes, "significant genes found"),
        timestamp = Sys.time()
      )
    } else {
      progress_details$deseq2_analysis <- list(status = "pending", message = "Run differential analysis")
    }
    
    # Pathway analysis progress (15%)
    if (!is.null(values$pathway_results) && values$pathway_results$success) {
      progress <- progress + 15
      n_pathways <- nrow(values$pathway_results$data)
      progress_details$pathway_analysis <- list(
        status = "completed",
        message = paste(n_pathways, "pathways analyzed"),
        timestamp = Sys.time()
      )
    } else {
      progress_details$pathway_analysis <- list(status = "pending", message = "Pathway analysis (optional)")
    }
    
    # Update progress bar
    updateProgressBar(session, "overall_progress", value = progress)
    
    # Update progress details in sidebar
    update_progress_sidebar(progress_details)
    
    values$progress <- progress
    values$progress_details <- progress_details
  })
  
  # Function to update progress sidebar
  update_progress_sidebar <- function(details) {
    # Update step indicators with status
    for (step_name in names(details)) {
      step_detail <- details[[step_name]]
      
      # Update step color and icon based on status
      step_color <- switch(
        step_detail$status,
        "completed" = "color: #28a745; font-weight: bold;",
        "in_progress" = "color: #ffc107; font-weight: bold;",
        "pending" = "color: #6c757d;",
        "color: #6c757d;"
      )
      
      step_icon <- switch(
        step_detail$status,
        "completed" = "✅",
        "in_progress" = "🔄",
        "pending" = "⏸️",
        "⏸️"
      )
      
      # Note: In a real implementation, you'd use updateUI or similar
      # This is a simplified version for demonstration
    }
  }
  
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
                                           1:min(6, ncol(values$expression_data))]
      
      # Add gene symbols if available
      gene_symbols <- attr(values$expression_data, "gene_symbols")
      preview_gene_ids <- rownames(preview_data)
      
      # Initialize variables for caption
      has_meaningful_symbols <- FALSE
      conversion_rate <- 0
      conversion_species <- "unknown"
      
      if (!is.null(gene_symbols)) {
        # Match gene symbols to preview gene IDs
        symbol_matches <- match(preview_gene_ids, gene_symbols$ensembl_gene_id)
        symbols <- gene_symbols$gene_symbol[symbol_matches]
        
        # Check how many genes actually have different symbols vs IDs
        valid_symbols <- !is.na(symbols) & symbols != "" & symbols != preview_gene_ids
        has_meaningful_symbols <- sum(valid_symbols) > 0
        
        if (has_meaningful_symbols) {
          # Some genes have meaningful symbols - show both columns
          # Use symbol where available, ID as fallback
          display_symbols <- ifelse(valid_symbols, symbols, preview_gene_ids)
          
          preview_data <- cbind(
            Gene_Symbol = display_symbols,
            Ensembl_ID = preview_gene_ids,
            preview_data
          )
        } else {
          # No meaningful symbols found - just show IDs
          preview_data <- cbind(
            Gene_ID = preview_gene_ids,
            preview_data
          )
        }
        
        # Get conversion info
        conversion_rate <- attr(values$expression_data, "conversion_rate") %||% 0
        conversion_species <- attr(values$expression_data, "conversion_species") %||% "unknown"
        
      } else {
        # No gene symbols available
        preview_data <- cbind(Gene_ID = preview_gene_ids, preview_data)
      }
      
      DT::datatable(
        preview_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = if (!is.null(gene_symbols)) {
          if (has_meaningful_symbols) {
            paste0("Data preview with gene symbols (", conversion_rate, "% conversion rate for ", conversion_species, " genes)")
          } else {
            paste0("Data preview showing Ensembl IDs only (", conversion_rate, "% genes had symbols, but symbols matched IDs for ", conversion_species, ")")
          }
        } else {
          "Data preview (enable gene conversion in Data Upload for symbols)"
        }
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