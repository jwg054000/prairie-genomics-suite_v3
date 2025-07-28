# Prairie Genomics Suite v5 Enhanced - Main Application
# Advanced genomics analysis platform with multi-group support and Context7 visualizations
# Based on Emory DESeq2 methodology with batch correction and comprehensive statistics

# Load required packages for enhanced version
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  shiny, shinydashboard, DT, plotly, ggplot2, dplyr, readr, 
  DESeq2, sva, MASS, pheatmap, RColorBrewer, car, ggrepel,
  readxl, logger
)

# Increase Shiny file upload limit to 500MB (default is 5MB)
options(shiny.maxRequestSize = 500*1024^2)

# Set up logging
log_appender(appender_file("prairie_genomics_v5.log"))
log_threshold(INFO)

# Source enhanced modules
source("modules/enhanced_sample_annotation.R")
source("modules/enhanced_deseq2_analysis.R") 
source("modules/context7_visualizations.R")
source("utils/batch_correction.R")

# Optional packages with graceful handling
bioc_available <- requireNamespace("DESeq2", quietly = TRUE)
excel_support <- requireNamespace("readxl", quietly = TRUE)
enhanced_plots <- requireNamespace("RColorBrewer", quietly = TRUE)

# Define Enhanced UI
ui <- dashboardPage(
  
  # Enhanced Header
  dashboardHeader(
    title = tags$span(
      icon("dna"),
      "Prairie Genomics Suite v5",
      style = "font-size: 18px; font-weight: bold;"
    ),
    tags$li(
      class = "dropdown",
      tags$a(
        href = "#",
        tags$span(
          "Enhanced R Shiny v5.0 🚀", 
          style = "color: #ffffff; font-size: 12px;"
        )
      )
    )
  ),
  
  # Enhanced Sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      
      # Main workflow tabs
      menuItem("📁 Enhanced Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("🧬 Multi-Group Annotation", tabName = "annotation", icon = icon("tags")),
      menuItem("🚀 Advanced DESeq2", tabName = "analysis", icon = icon("chart-line")),
      menuItem("🎨 Context7 Visualizations", tabName = "visualizations", icon = icon("chart-bar")),
      menuItem("📋 Results Export", tabName = "export", icon = icon("download")),
      
      hr(),
      
      # System information
      menuSubItem("💻 System Status", tabName = "system", icon = icon("info-circle")),
      menuSubItem("📖 Help & Documentation", tabName = "help", icon = icon("question-circle")),
      
      hr(),
      
      # Enhanced progress tracking
      div(
        style = "padding: 15px;",
        h4("🔄 Analysis Pipeline", style = "color: white; font-size: 14px;"),
        
        # Progress steps with status indicators
        div(
          style = "color: #aaa; font-size: 12px;",
          div(id = "step1", "📁 Data Upload", style = "margin: 5px 0;"),
          div(id = "step2", "🧬 Sample Annotation", style = "margin: 5px 0;"),
          div(id = "step3", "🔬 Batch Correction", style = "margin: 5px 0;"),
          div(id = "step4", "🚀 Statistical Analysis", style = "margin: 5px 0;"),
          div(id = "step5", "🎨 Visualization", style = "margin: 5px 0;")
        ),
        
        br(),
        
        # Overall progress
        div(
          style = "background-color: #3c8dbc; color: white; padding: 8px; border-radius: 4px; text-align: center;",
          id = "overall_progress_display",
          "Analysis Progress: 0%"
        )
      ),
      
      hr(),
      
      # Enhanced quick stats
      div(
        style = "padding: 15px;",
        h4("📊 Quick Stats", style = "color: white; font-size: 14px;"),
        
        valueBoxOutput("enhanced_genes_count", width = 12),
        valueBoxOutput("enhanced_samples_count", width = 12),
        valueBoxOutput("enhanced_groups_count", width = 12),
        valueBoxOutput("enhanced_comparisons_count", width = 12)
      )
    )
  ),
  
  # Enhanced Body
  dashboardBody(
    
    # Enhanced custom CSS
    tags$head(
      tags$style(HTML("
        /* Enhanced styling for v5 */
        .content-wrapper, .right-side {
          background-color: #f8f9fa;
        }
        
        .small-box {
          margin-bottom: 8px;
          border-radius: 8px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        .small-box .inner {
          padding: 8px 12px;
        }
        
        .small-box h3 {
          font-size: 18px;
          margin: 0;
          font-weight: bold;
        }
        
        .small-box p {
          font-size: 11px;
          margin: 0;
          font-weight: 500;
        }
        
        /* Enhanced alert styling */
        .alert {
          border-radius: 8px;
          border: none;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        .alert-success {
          background: linear-gradient(135deg, #d4edda 0%, #c3e6cb 100%);
          color: #155724;
        }
        
        .alert-warning {
          background: linear-gradient(135deg, #fff3cd 0%, #ffeaa7 100%);
          color: #856404;
        }
        
        .alert-info {
          background: linear-gradient(135deg, #d1ecf1 0%, #b8daff 100%);
          color: #0c5460;
        }
        
        .alert-danger {
          background: linear-gradient(135deg, #f8d7da 0%, #f5c6cb 100%);
          color: #721c24;
        }
        
        /* Enhanced button styling */
        .btn-primary {
          background: linear-gradient(135deg, #007bff 0%, #0056b3 100%);
          border: none;
          border-radius: 8px;
          box-shadow: 0 2px 4px rgba(0,123,255,0.3);
          transition: all 0.3s ease;
        }
        
        .btn-primary:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 8px rgba(0,123,255,0.4);
        }
        
        /* Enhanced box styling */
        .box {
          border-radius: 12px;
          box-shadow: 0 4px 6px rgba(0,0,0,0.1);
          border: none;
        }
        
        .box-header {
          border-radius: 12px 12px 0 0;
        }
        
        /* Progress indicator styling */
        #overall_progress_display {
          background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
          border-radius: 8px;
          font-size: 12px;
          font-weight: bold;
          box-shadow: 0 2px 4px rgba(40,167,69,0.3);
        }
        
        /* Enhanced table styling */
        .dataTables_wrapper {
          border-radius: 8px;
          overflow: hidden;
        }
      "))
    ),
    
    # Main content tabs
    tabItems(
      
      # Enhanced Data Upload Tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "📁 Enhanced Data Upload System", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            div(
              class = "alert alert-info",
              h4("🚀 Welcome to Prairie Genomics Suite v5!"),
              p("This enhanced version features multi-group analysis, batch effect correction, and Context7-inspired visualizations."),
              tags$ul(
                tags$li("✅ Multi-group experimental designs (GroupA, GroupB, GroupC, etc.)"),
                tags$li("✅ Automatic batch effect detection and correction"),
                tags$li("✅ All pairwise comparisons with advanced statistics"),
                tags$li("✅ Publication-quality Context7-enhanced visualizations"),
                tags$li("✅ Enhanced accessibility and performance optimization")
              )
            ),
            
            # Enhanced file upload interface
            fluidRow(
              column(
                6,
                h5("📊 Expression Data Upload"),
                
                fileInput(
                  "expression_file",
                  "Choose Expression Data File",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                  placeholder = "RNA-seq count matrix (genes × samples)"
                ),
                
                checkboxInput(
                  "expression_has_header",
                  "File has header row",
                  value = TRUE
                ),
                
                conditionalPanel(
                  condition = "output.has_expression_file",
                  
                  h6("📋 Data Validation"),
                  verbatimTextOutput("expression_file_info"),
                  
                  actionButton(
                    "process_expression_data",
                    "🔄 Process Expression Data",
                    class = "btn-primary",
                    style = "margin-top: 10px;"
                  )
                )
              ),
              
              column(
                6,
                h5("📄 Test Data & Examples"),
                
                p("New to the platform? Try our enhanced test datasets:"),
                
                fluidRow(
                  column(
                    6,
                    actionButton(
                      "load_simple_test",
                      "📊 Simple Test Data",
                      class = "btn-info btn-sm",
                      style = "width: 100%; margin-bottom: 5px;"
                    ),
                    p("2 groups, 6 samples", style = "font-size: 11px; color: #666;")
                  ),
                  column(
                    6,
                    actionButton(
                      "load_complex_test",
                      "🧬 Multi-Group Test",
                      class = "btn-warning btn-sm", 
                      style = "width: 100%; margin-bottom: 5px;"
                    ),
                    p("5 groups, 25 samples (multi-group)", style = "font-size: 11px; color: #666;")
                  )
                ),
                
                br(),
                h6("💡 Format Requirements:"),
                tags$ul(
                  tags$li("Genes as rows, samples as columns"),
                  tags$li("Raw counts (not normalized values)"),
                  tags$li("Gene IDs in first column"),
                  tags$li("Clear sample naming patterns")
                )
              )
            )
          )
        ),
        
        # Enhanced data preview with statistics
        fluidRow(
          box(
            title = "📊 Enhanced Data Preview & Quality Control", 
            status = "info", 
            solidHeader = TRUE, 
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            conditionalPanel(
              condition = "output.has_processed_data",
              
              fluidRow(
                column(
                  8,
                  h5("🔍 Data Table Preview"),
                  DT::dataTableOutput("enhanced_data_preview")
                ),
                column(
                  4,
                  h5("📈 Quality Metrics"),
                  
                  wellPanel(
                    h6("📊 Dataset Summary"),
                    uiOutput("dataset_summary"),
                    
                    br(),
                    h6("🎯 Sample Name Patterns"),
                    uiOutput("sample_patterns"),
                    
                    br(),
                    h6("⚠️ Quality Warnings"),
                    uiOutput("quality_warnings")
                  )
                )
              )
            )
          )
        )
      ),
      
      # Enhanced Sample Annotation Tab
      tabItem(
        tabName = "annotation",
        conditionalPanel(
          condition = "!output.has_processed_data",
          fluidRow(
            box(
              title = "⚠️ Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              
              div(
                class = "alert alert-warning",
                h4("📁 Upload Expression Data First"),
                p("Please upload and process your expression data in the Data Upload tab before proceeding with sample annotation.")
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_processed_data",
          enhancedSampleAnnotationUI("sample_annotation")
        )
      ),
      
      # Enhanced DESeq2 Analysis Tab
      tabItem(
        tabName = "analysis",
        conditionalPanel(
          condition = "!output.has_annotation_data",
          fluidRow(
            box(
              title = "⚠️ Annotation Required",
              status = "warning", 
              solidHeader = TRUE,
              width = 12,
              
              div(
                class = "alert alert-warning",
                h4("🧬 Complete Sample Annotation First"),
                p("Please annotate your samples in the Multi-Group Annotation tab before running enhanced DESeq2 analysis.")
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_annotation_data",
          enhancedDESeq2AnalysisUI("deseq2_analysis")
        )
      ),
      
      # Context7 Visualizations Tab
      tabItem(
        tabName = "visualizations",
        conditionalPanel(
          condition = "!output.has_results",
          fluidRow(
            box(
              title = "⚠️ Analysis Results Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              
              div(
                class = "alert alert-warning",
                h4("🚀 Run Analysis First"),
                p("Please complete your enhanced DESeq2 analysis before accessing Context7 visualizations.")
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.has_results",
          context7VisualizationsUI("visualizations")
        )
      ),
      
      # Enhanced Export Tab
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "📋 Enhanced Results Export System", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "!output.has_results",
              div(
                class = "alert alert-warning",
                h4("📊 No Results Available"),
                p("Please complete your analysis to access export options.")
              )
            ),
            
            conditionalPanel(
              condition = "output.has_results",
              
              h4("🎯 Export Options"),
              
              tabsetPanel(
                
                # Statistical results export
                tabPanel(
                  "📊 Statistical Results",
                  br(),
                  
                  fluidRow(
                    column(
                      6,
                      h5("🔢 Analysis Results"),
                      
                      downloadButton(
                        "download_all_comparisons",
                        "📥 All Pairwise Comparisons", 
                        class = "btn-primary",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_significant_only",
                        "⭐ Significant Genes Only",
                        class = "btn-success", 
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_filtered_expression",
                        "📈 Filtered Expression Matrix",
                        class = "btn-info",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_sample_annotation",
                        "🏷️ Sample Annotation",
                        class = "btn-secondary",
                        style = "width: 100%;"
                      )
                    ),
                    
                    column(
                      6,
                      h5("📊 Summary Statistics"),
                      
                      downloadButton(
                        "download_comparison_summary",
                        "📋 Comparison Summary",
                        class = "btn-primary",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_overlap_analysis",
                        "🔄 Gene Overlap Analysis",
                        class = "btn-warning",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_batch_diagnostics",
                        "🔬 Batch Effect Report",
                        class = "btn-info",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_analysis_log",
                        "📝 Analysis Log",
                        class = "btn-secondary",
                        style = "width: 100%;"
                      )
                    )
                  )
                ),
                
                # Visualization export
                tabPanel(
                  "🎨 Visualizations",
                  br(),
                  
                  p("High-resolution, publication-quality plots from Context7 visualizations."),
                  
                  fluidRow(
                    column(
                      6,
                      h5("🌋 Individual Plots"),
                      
                      p("Export individual plots for each comparison:"),
                      
                      uiOutput("individual_plot_downloads")
                    ),
                    
                    column(
                      6,
                      h5("📊 Combined Figures"),
                      
                      downloadButton(
                        "download_figure_panel",
                        "🖼️ Multi-Panel Figure",
                        class = "btn-primary",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_supplementary_figures",
                        "📚 Supplementary Figures",
                        class = "btn-success",
                        style = "width: 100%; margin-bottom: 10px;"
                      ),
                      
                      downloadButton(
                        "download_interactive_html",
                        "🌐 Interactive HTML Report",
                        class = "btn-info",
                        style = "width: 100%;"
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ),
      
      # System Status Tab
      tabItem(
        tabName = "system",
        fluidRow(
          box(
            title = "💻 System Status & Information",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            h4("🔧 Package Status"),
            verbatimTextOutput("package_status"),
            
            br(),
            h4("📊 Memory Usage"),
            verbatimTextOutput("memory_status"),
            
            br(),
            h4("⏱️ Session Information"),
            verbatimTextOutput("session_info")
          )
        )
      ),
      
      # Help & Documentation Tab
      tabItem(
        tabName = "help",
        fluidRow(
          box(
            title = "📖 Help & Documentation",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            h3("🚀 Prairie Genomics Suite v5 Enhanced"),
            
            p("Welcome to the most advanced version of Prairie Genomics Suite, featuring:"),
            
            h4("✨ New Features in v5:"),
            tags$ul(
              tags$li("🧬 Multi-group experimental design support (GroupA, GroupB, etc.)"),
              tags$li("🔬 Automatic batch effect detection and correction using ComBat-seq"),
              tags$li("🔄 Comprehensive pairwise comparison analysis"),
              tags$li("🎨 Context7-enhanced visualizations with accessibility features"),
              tags$li("📊 Advanced PCA plots with 3D support and confidence ellipses"),
              tags$li("🌋 Publication-quality volcano plots with gene highlighting"),
              tags$li("🔥 Interactive heatmaps with smart gene selection"),
              tags$li("📈 Enhanced performance for large datasets")
            ),
            
            h4("📋 Quick Start Guide:"),
            tags$ol(
              tags$li("📁 Upload your RNA-seq count matrix (genes × samples)"),
              tags$li("🧬 Use smart pattern detection or manual annotation for sample groups"),
              tags$li("🔬 Review batch effect diagnostics and apply correction if needed"),
              tags$li("🚀 Run enhanced DESeq2 analysis with multiple comparisons"),
              tags$li("🎨 Explore Context7-enhanced interactive visualizations"),
              tags$li("📥 Export publication-ready results and figures")
            ),
            
            h4("💡 Tips for Best Results:"),
            tags$ul(
              tags$li("Use clear, consistent sample naming patterns (e.g., Group_Sample_Rep)"),
              tags$li("Include at least 3 biological replicates per group"),
              tags$li("Provide raw counts, not normalized values"),
              tags$li("Review batch correction recommendations"),
              tags$li("Use gene highlighting in volcano plots for genes of interest")
            ),
            
            h4("🆘 Troubleshooting:"),
            tags$ul(
              tags$li("If pattern detection fails, try manual annotation"),
              tags$li("Check the System Status tab for package availability"),
              tags$li("Large datasets may require performance optimization"),
              tags$li("Contact support if analysis fails repeatedly")
            )
          )
        )
      )
    )
  )
)

# Define Enhanced Server
server <- function(input, output, session) {
  
  # Enhanced reactive values
  values <- reactiveValues(
    expression_data = NULL,
    annotation_data = NULL,
    deseq2_results = NULL,
    comparison_results = NULL,
    filtered_data = NULL,
    batch_data = NULL,
    progress = 0,
    analysis_log = list()
  )
  
  # Enhanced data upload handling
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    tryCatch({
      file_path <- input$expression_file$datapath
      file_ext <- tools::file_ext(input$expression_file$name)
      
      # Read file based on extension
      if (file_ext == "csv") {
        raw_data <- readr::read_csv(file_path, col_names = input$expression_has_header, show_col_types = FALSE)
      } else if (file_ext %in% c("tsv", "txt")) {
        raw_data <- readr::read_delim(file_path, delim = "\t", col_names = input$expression_has_header, show_col_types = FALSE)
      } else if (file_ext == "xlsx" && excel_support) {
        raw_data <- readxl::read_excel(file_path, col_names = input$expression_has_header)
      } else {
        stop("Unsupported file format or readxl not available")
      }
      
      # Store raw data for processing
      values$raw_expression_data <- raw_data
      
      add_to_log("Expression file uploaded successfully", "success")
      
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
      add_to_log(paste("File upload error:", e$message), "error")
    })
  })
  
  # Process expression data
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
      
      cat("Starting data processing...\n")
      
      # Process the uploaded data with chunked processing
      processed_data <- process_expression_data(values$raw_expression_data)
      
      if (!is.null(processed_data)) {
        values$expression_data <- processed_data
        values$progress <- 20
        
        add_to_log("Expression data processed successfully", "success")
        showNotification(paste("Data processed successfully!", nrow(processed_data), 
                              "genes retained from", nrow(values$raw_expression_data), "original genes."), 
                        type = "message")
        
        # Final memory check
        final_mem <- monitor_memory()
        cat("Final memory status:", final_mem$message, "\n")
      }
      
    }, error = function(e) {
      showNotification(paste("Error processing data:", e$message), type = "error")
      add_to_log(paste("Data processing error:", e$message), "error")
      cat("Processing error:", e$message, "\n")
    })
  })
  
  # Load test data
  observeEvent(input$load_simple_test, {
    # Load simple test data (existing test_data.csv)
    test_data <- read.csv("test_data.csv", row.names = 1)
    values$expression_data <- test_data
    values$progress <- 20
    
    add_to_log("Simple test data loaded", "info")
    showNotification("Simple test data loaded successfully!", type = "message")
  })
  
  observeEvent(input$load_complex_test, {
    # Create complex multi-group test data (Emory-style)
    complex_test_data <- create_complex_test_data()
    values$expression_data <- complex_test_data
    values$progress <- 20
    
    add_to_log("Complex multi-group test data loaded", "info")
    showNotification("Complex multi-group test data loaded successfully!", type = "message")
  })
  
  # Enhanced sample annotation module
  annotation_result <- callModule(enhancedSampleAnnotation, "sample_annotation", values)
  
  # Enhanced DESeq2 analysis module  
  analysis_result <- callModule(enhancedDESeq2Analysis, "deseq2_analysis", values)
  
  # Context7 visualizations module
  callModule(context7Visualizations, "visualizations", values)
  
  # Enhanced progress tracking
  observe({
    progress <- 0
    
    if (!is.null(values$expression_data)) {
      progress <- progress + 20
    }
    
    if (!is.null(values$annotation_data)) {
      progress <- progress + 20
    }
    
    if (!is.null(values$deseq2_results) || !is.null(values$comparison_results)) {
      progress <- progress + 40
    }
    
    if (progress > 60) {
      progress <- progress + 20  # Visualization step
    }
    
    values$progress <- progress
    
    # Update progress display
    session$sendCustomMessage("updateProgress", list(
      progress = progress,
      text = paste("Analysis Progress:", progress, "%")
    ))
  })
  
  # Enhanced conditional outputs
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
  
  # Enhanced value boxes
  output$enhanced_genes_count <- renderValueBox({
    count <- if (!is.null(values$expression_data)) nrow(values$expression_data) else 0
    valueBox(
      value = formatC(count, format = "d", big.mark = ","),
      subtitle = "Genes",
      icon = icon("dna", class = "fa-lg"),
      color = if (count > 0) "green" else "light-blue",
      width = 12
    )
  })
  
  output$enhanced_samples_count <- renderValueBox({
    count <- if (!is.null(values$expression_data)) ncol(values$expression_data) else 0
    valueBox(
      value = count,
      subtitle = "Samples",
      icon = icon("vials", class = "fa-lg"),
      color = if (count > 0) "blue" else "light-blue",
      width = 12
    )
  })
  
  output$enhanced_groups_count <- renderValueBox({
    count <- if (!is.null(values$annotation_data)) {
      length(unique(values$annotation_data$Condition))
    } else 0
    valueBox(
      value = count,
      subtitle = "Groups",
      icon = icon("layer-group", class = "fa-lg"),
      color = if (count > 1) "purple" else "light-blue",
      width = 12
    )
  })
  
  output$enhanced_comparisons_count <- renderValueBox({
    count <- if (!is.null(values$comparison_results)) {
      length(values$comparison_results)
    } else if (!is.null(values$deseq2_results)) {
      1
    } else 0
    valueBox(
      value = count,
      subtitle = "Comparisons",
      icon = icon("chart-line", class = "fa-lg"),
      color = if (count > 0) "yellow" else "light-blue",
      width = 12
    )
  })
  
  # Enhanced data preview
  output$enhanced_data_preview <- DT::renderDataTable({
    req(values$expression_data)
    
    # Show first 15 rows and 10 columns with gene names
    preview_data <- values$expression_data[1:min(15, nrow(values$expression_data)), 
                                         1:min(10, ncol(values$expression_data))]
    preview_data <- cbind(Gene = rownames(preview_data), preview_data)
    
    DT::datatable(
      preview_data,
      options = list(
        scrollX = TRUE,
        scrollY = "400px",
        pageLength = 15,
        dom = 'tip',
        columnDefs = list(
          list(className = 'dt-head-center', targets = "_all")
        )
      ),
      style = "bootstrap4",
      class = "table-striped table-hover"
    ) %>%
      DT::formatRound(columns = 2:ncol(preview_data), digits = 1)
  })
  
  # System status outputs
  output$package_status <- renderText({
    status <- c(
      paste("DESeq2 available:", if (bioc_available) "✅ Yes" else "❌ No"),
      paste("Excel support:", if (excel_support) "✅ Yes" else "❌ No"),
      paste("Enhanced plots:", if (enhanced_plots) "✅ Yes" else "❌ No"),
      paste("R version:", R.version.string),
      paste("Shiny version:", packageVersion("shiny"))
    )
    paste(status, collapse = "\n")
  })
  
  output$memory_status <- renderText({
    gc_info <- gc()
    paste("Memory usage:", round(sum(gc_info[, 2]), 1), "MB")
  })
  
  output$session_info <- renderText({
    paste("Session started:", Sys.time())
  })
  
  # Helper functions
  
  # Memory monitoring function
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
  
  # Robust chunked data processing function
  process_expression_data <- function(raw_data, chunk_size = 5000) {
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
      cat("Error in process_expression_data:", e$message, "\n")
      stop(paste("Data processing failed:", e$message))
    })
  }
  
  create_complex_test_data <- function() {
    # Create multi-group test data
    set.seed(42)
    
    n_genes <- 1000
    groups <- c("GroupA", "GroupB", "GroupC", "GroupD", "GroupE")
    n_per_group <- 5
    
    # Generate sample names
    sample_names <- paste(rep(groups, each = n_per_group), 1:n_per_group, sep = "_")
    
    # Generate expression data with group-specific patterns
    expression_data <- matrix(
      rnbinom(n_genes * length(sample_names), size = 10, mu = 100),
      nrow = n_genes,
      ncol = length(sample_names)
    )
    
    colnames(expression_data) <- sample_names
    rownames(expression_data) <- paste0("Gene_", 1:n_genes)
    
    # Add some differentially expressed genes
    de_genes <- 1:100
    group_effects <- c(0, 2, -1.5, 1, -0.5)  # Effect sizes for each group
    
    for (i in seq_along(groups)) {
      group_samples <- grep(paste0("^", groups[i]), sample_names)
      if (length(group_samples) > 0) {
        expression_data[de_genes, group_samples] <- expression_data[de_genes, group_samples] * 
          exp(group_effects[i])
      }
    }
    
    return(expression_data)
  }
  
  add_to_log <- function(message, type = "info") {
    log_entry <- list(
      timestamp = Sys.time(),
      message = message,
      type = type
    )
    values$analysis_log <- append(values$analysis_log, list(log_entry))
    
    # Log to file as well
    if (type == "error") {
      log_error(message)
    } else if (type == "warning") {
      log_warn(message)
    } else {
      log_info(message)
    }
  }
}

# Custom JavaScript for enhanced interactivity
js_code <- "
Shiny.addCustomMessageHandler('updateProgress', function(data) {
  var element = document.getElementById('overall_progress_display');
  if (element) {
    element.innerHTML = data.text;
    element.style.background = 'linear-gradient(135deg, #28a745 0%, #20c997 100%)';
  }
});
"

# Run the enhanced application
shinyApp(ui = ui, server = server)