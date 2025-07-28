# Prairie Genomics Suite - Optimized R Shiny Application
# Modular, memory-efficient genomics analysis platform
# 
# Author: Prairie Genomics Team - Optimized Version 2.0
# Features: Modular architecture, memory management, async processing, enhanced error handling

# =============================================================================
# INITIALIZATION AND CONFIGURATION
# =============================================================================

cat("🚀 Initializing Prairie Genomics Suite - Optimized Version\n")
cat("=" , rep("=", 60), "\n")

# Load core configuration first
source("config/app_config.R")

# Increase Shiny file upload limit
options(shiny.maxRequestSize = get_config("memory", "max_file_size_mb") * 1024^2)

# Load essential packages with error handling
load_essential_packages <- function() {
  essential_packages <- get_config("packages", "essential")
  
  for (pkg in essential_packages) {
    tryCatch({
      library(pkg, character.only = TRUE)
      cat("✅", pkg, "loaded\n")
    }, error = function(e) {
      cat("❌ Failed to load essential package:", pkg, "\n")
      stop(paste("Essential package", pkg, "is required but not available"))
    })
  }
}

load_essential_packages()

# Load optional packages with graceful handling  
load_optional_packages <- function() {
  optional_packages <- get_config("packages", "optional")
  
  for (pkg in optional_packages) {
    tryCatch({
      library(pkg, character.only = TRUE)
      cat("✅", pkg, "loaded\n")
    }, error = function(e) {
      cat("⚠️", pkg, "not available - some features may be limited\n")
    })
  }
}

load_optional_packages()

# =============================================================================
# LOAD OPTIMIZED MODULES
# =============================================================================

cat("\n📦 Loading optimized modules...\n")

# Load utilities first
tryCatch({
  source("utils/memory_manager.R")
  cat("✅ Memory management utilities loaded\n")
}, error = function(e) {
  cat("❌ Failed to load memory manager:", e$message, "\n")
  stop("Memory manager is required for optimal performance")
})

# Load data processing modules
tryCatch({
  source("modules/data_processing/file_upload.R")
  cat("✅ Data processing modules loaded\n")
}, error = function(e) {
  cat("❌ Failed to load data processing modules:", e$message, "\n")
  stop("Data processing modules are required")
})

# Load analysis modules
tryCatch({
  source("modules/analysis/deseq2_analysis.R")
  cat("✅ Analysis modules loaded\n")
}, error = function(e) {
  cat("❌ Failed to load analysis modules:", e$message, "\n")
  stop("Analysis modules are required")
})

# Load UI modules
tryCatch({
  source("modules/ui/sample_annotation.R")
  cat("✅ UI modules loaded\n")
}, error = function(e) {
  cat("❌ Failed to load UI modules:", e$message, "\n")
  stop("UI modules are required")
})

# Initialize global memory manager
init_memory_manager()

cat("✅ All modules loaded successfully - Prairie Genomics Suite ready!\n\n")

# =============================================================================
# USER INTERFACE
# =============================================================================

ui <- dashboardPage(
  
  # Header
  dashboardHeader(
    title = tags$span(
      icon("dna"), 
      get_config("app", "name"),
      style = "font-size: 18px; font-weight: bold;"
    ),
    tags$li(
      class = "dropdown",
      tags$a(
        href = "#",
        tags$span(
          paste("v", get_config("app", "version")), 
          style = "color: #ffffff; font-size: 12px;"
        )
      )
    )
  ),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      
      # Main navigation
      menuItem("📁 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("🧬 Sample Annotation", tabName = "annotation", icon = icon("tags")),
      menuItem("🚀 DESeq2 Analysis", tabName = "analysis", icon = icon("chart-line")),
      menuItem("📊 Visualizations", tabName = "visualizations", icon = icon("chart-bar")),
      menuItem("📋 Export Results", tabName = "export", icon = icon("download")),
      
      hr(),
      
      # System monitoring
      menuItem("💾 System Monitor", tabName = "monitor", icon = icon("tachometer-alt")),
      
      hr(),
      
      # Progress tracking
      div(
        style = "padding: 15px;",
        h4("Analysis Progress", style = "color: white; font-size: 14px;"),
        
        div(
          style = "color: #aaa; font-size: 12px;",
          div(id = "step1", "📁 Upload Data"),
          div(id = "step2", "🧬 Annotate Samples"), 
          div(id = "step3", "🚀 Run Analysis"),
          div(id = "step4", "📊 View Results")
        ),
        
        br(),
        
        # Progress bar with fallback
        if (requireNamespace("shinyWidgets", quietly = TRUE)) {
          shinyWidgets::progressBar(
            id = "overall_progress",
            value = 0,
            total = 100,
            status = "primary",
            display_pct = TRUE,
            striped = TRUE
          )
        } else {
          div(
            class = "progress",
            div(
              class = "progress-bar progress-bar-striped progress-bar-animated",
              role = "progressbar",
              style = "width: 0%",
              id = "overall_progress",
              "0%"
            )
          )
        }
      ),
      
      hr(),
      
      # Quick stats with memory info
      div(
        style = "padding: 15px;",
        h4("System Stats", style = "color: white; font-size: 14px;"),
        
        valueBoxOutput("genes_count", width = 12),
        valueBoxOutput("samples_count", width = 12),
        valueBoxOutput("memory_usage", width = 12)
      )
    )
  ),
  
  # Main content
  dashboardBody(
    
    # Custom CSS for optimized interface
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
          font-size: 16px;
          margin: 0;
        }
        
        .small-box p {
          font-size: 11px;
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
        
        .alert-danger {
          background-color: #f8d7da;
          border-color: #f5c6cb;
          color: #721c24;
        }
        
        .memory-critical {
          background-color: #f8d7da !important;
          border: 2px solid #dc3545 !important;
        }
        
        .memory-warning {
          background-color: #fff3cd !important;
          border: 2px solid #ffc107 !important;
        }
        
        .memory-ok {
          background-color: #d4edda !important;
          border: 2px solid #28a745 !important;
        }
      "))
    ),
    
    tabItems(
      
      # Data Upload Tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "📁 Optimized Data Upload", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            div(
              class = "alert alert-info",
              h4("Welcome to Prairie Genomics Suite - Optimized!"),
              p("Upload your RNA-seq expression data with enhanced memory management and validation."),
              p("Supported formats: CSV, TSV, Excel. Genes as rows, samples as columns."),
              p(paste("Maximum file size:", get_config("memory", "max_file_size_mb"), "MB"))
            ),
            
            dataUploadUI("data_upload"),
            
            # Development test data button
            conditionalPanel(
              condition = "typeof input.test_data_ui !== 'undefined'",
              
              br(),
              hr(),
              
              div(
                class = "alert alert-secondary",
                h6("🧪 Development Mode"),
                p("Generate test data for development and testing purposes."),
                uiOutput("test_data_ui")
              )
            )
          )
        ),
        
        # Enhanced data preview with memory monitoring
        fluidRow(
          box(
            title = "📊 Smart Data Preview", 
            status = "info", 
            solidHeader = TRUE, 
            width = 8,
            collapsible = TRUE,
            collapsed = TRUE,
            
            DT::dataTableOutput("data_preview")
          ),
          
          box(
            title = "📈 Processing Stats",
            status = "success",
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            collapsed = TRUE,
            
            uiOutput("processing_stats")
          )
        )
      ),
      
      # Sample Annotation Tab
      tabItem(
        tabName = "annotation",
        fluidRow(
          box(
            title = "🧬 Enhanced Sample Annotation", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_expression_data == false",
              div(
                class = "alert alert-warning",
                h4("Upload Data First"),
                p("Please upload your expression data in the Data Upload tab before proceeding.")
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
            title = "🚀 Optimized DESeq2 Analysis", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            conditionalPanel(
              condition = "output.has_annotation_data == false",
              div(
                class = "alert alert-warning",
                h4("Complete Sample Annotation First"),
                p("Please annotate your samples before running DESeq2 analysis.")
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
            title = "📊 Memory-Efficient Visualizations", 
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
              div(
                class = "alert alert-info",
                h4("Optimized Visualization Module"),
                p("Memory-efficient plots with timeout protection and smart rendering."),
                p("This module is being implemented with the optimized architecture.")
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
            title = "📋 Enhanced Results Export", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            div(
              class = "alert alert-info",
              h4("Export Module"),
              p("Enhanced export capabilities with memory-efficient processing."),
              p("This module is being implemented with the optimized architecture.")
            )
          )
        )
      ),
      
      # System Monitor Tab
      tabItem(
        tabName = "monitor",
        fluidRow(
          box(
            title = "💾 Real-time System Monitor", 
            status = "primary", 
            solidHeader = TRUE, 
            width = 12,
            
            fluidRow(
              column(
                6,
                h5("📊 Memory Usage"),
                uiOutput("memory_monitor"),
                
                br(),
                
                h5("🧹 Memory Management"),
                actionButton("force_gc", "Force Garbage Collection", 
                           class = "btn-warning", style = "margin-right: 10px;"),
                actionButton("memory_report", "Generate Memory Report", 
                           class = "btn-info")
              ),
              
              column(
                6,
                h5("📈 Performance Metrics"),
                uiOutput("performance_metrics"),
                
                br(),
                
                h5("⚙️ System Information"),
                verbatimTextOutput("system_info")
              )
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER LOGIC
# =============================================================================

server <- function(input, output, session) {
  
  # Initialize reactive values with memory monitoring
  values <- reactiveValues(
    expression_data = NULL,
    annotation_data = NULL,
    deseq2_results = NULL,
    progress = 0,
    memory_alerts = 0
  )
  
  # =============================================================================
  # MEMORY MONITORING
  # =============================================================================
  
  # Real-time memory monitoring
  observe({
    invalidateLater(get_config("ui", "progress_update_interval"), session)
    
    memory_status <- get_memory_status()
    
    # Update memory usage display
    output$memory_usage <- renderValueBox({
      color <- if (memory_status$critical) {
        "red"
      } else if (memory_status$warning) {
        "yellow"  
      } else {
        "green"
      }
      
      valueBox(
        value = paste0(memory_status$used_mb, "MB"),
        subtitle = "Memory",
        icon = icon("memory"),
        color = color,
        width = 12
      )
    })
    
    # Memory alerts
    if (memory_status$critical && values$memory_alerts < 3) {
      showNotification(
        "🚨 Critical memory usage detected! Consider restarting or reducing data size.",
        type = "error",
        duration = 10
      )
      values$memory_alerts <- values$memory_alerts + 1
    } else if (memory_status$warning && values$memory_alerts == 0) {
      showNotification(
        "⚠️ High memory usage detected. Monitoring performance.",
        type = "warning",
        duration = 5
      )
      values$memory_alerts <- 1
    }
  })
  
  # =============================================================================
  # DATA UPLOAD MODULE
  # =============================================================================
  
  # Data upload server module
  dataUploadServer("data_upload", values)
  
  # Sample annotation server module
  sampleAnnotationServer("sample_annotation", values)
  
  # DESeq2 analysis server module
  deseq2AnalysisServer("deseq2_analysis", values)
  
  # Optional test data for development (only when explicitly requested)
  observeEvent(input$test_data, {
    # Only load test data when button is clicked AND no real data exists
    if (is_development_mode() && !is.null(input$test_data) && input$test_data > 0) {
      
      # Check if user already has uploaded data
      if (!is.null(values$expression_data) && nrow(values$expression_data) > 0) {
        showNotification("⚠️ Real data already loaded. Clear data first to use test data.", type = "warning")
        return()
      }
      
      cat("📊 Generating test data for demonstration...\n")
      
      # Create test expression data with realistic gene IDs
      test_matrix <- matrix(
        rpois(1000 * 12, lambda = 100),
        nrow = 1000,
        ncol = 12,
        dimnames = list(
          # Use realistic Ensembl gene IDs for testing
          c(paste0("ENSG0000014", sprintf("%04d", 1:500)),   # Human genes
            paste0("ENSMUSG0000014", sprintf("%04d", 1:500))), # Mouse genes
          paste0("Sample_", 1:12)
        )
      )
      
      values$expression_data <- test_matrix
      
      # Create test annotation
      values$annotation_data <- data.frame(
        Sample = paste0("Sample_", 1:12),
        Condition = rep(c("Control", "Treatment"), each = 6),
        stringsAsFactors = FALSE
      )
      
      showNotification("✅ Test data loaded successfully (mixed human/mouse genes for testing)", type = "message")
    }
  })
  
  # =============================================================================
  # SYSTEM MONITORING OUTPUTS
  # =============================================================================
  
  output$memory_monitor <- renderUI({
    memory_status <- get_memory_status()
    
    status_class <- if (memory_status$critical) {
      "memory-critical"
    } else if (memory_status$warning) {
      "memory-warning"
    } else {
      "memory-ok"
    }
    
    div(
      class = paste("alert", status_class),
      h6("Current Memory Status"),
      p(paste("Usage:", memory_status$used_mb, "MB")),
      p(paste("Status:", if (memory_status$critical) "CRITICAL" else if (memory_status$warning) "WARNING" else "OK")),
      p(paste("Recommendation:", memory_status$recommendation))
    )
  })
  
  output$performance_metrics <- renderUI({
    memory_report <- get_memory_report()
    
    div(
      h6("Performance Summary"),
      p(paste("Uptime:", memory_report$uptime_minutes, "minutes")),
      p(paste("GC Count:", memory_report$gc_count)),
      p(paste("Efficiency:", memory_report$memory_efficiency, "%")),
      p(paste("Trend:", memory_report$current_status$trend))
    )
  })
  
  output$system_info <- renderText({
    paste(
      "R Version:", R.version.string,
      "\nPlatform:", R.version$platform,
      "\nEnvironment:", if (is_cloud_deployment()) "Cloud" else "Local",
      "\nShiny Version:", packageVersion("shiny"),
      "\nMemory Limit:", get_config("memory", "memory_critical_mb"), "MB"
    )
  })
  
  # =============================================================================
  # MANUAL MEMORY MANAGEMENT
  # =============================================================================
  
  observeEvent(input$force_gc, {
    result <- smart_gc(force = TRUE)
    
    if (result$success) {
      showNotification(
        paste("🗑️ Garbage collection completed. Freed:", round(result$freed_mb, 1), "MB"),
        type = "message"
      )
    } else {
      showNotification(result$reason, type = "warning")
    }
  })
  
  observeEvent(input$memory_report, {
    report <- get_memory_report()
    
    showModal(modalDialog(
      title = "💾 Memory Report",
      
      h5("Current Status"),
      p(paste("Memory Usage:", report$current_status$used_mb, "MB")),
      p(paste("Status:", if (report$current_status$critical) "CRITICAL" else if (report$current_status$warning) "WARNING" else "NORMAL")),
      
      h5("Performance Metrics"),
      p(paste("Session Uptime:", report$uptime_minutes, "minutes")),
      p(paste("Garbage Collections:", report$gc_count)),
      p(paste("Memory Efficiency:", report$memory_efficiency, "%")),
      
      h5("Recommendations"),
      tags$ul(
        lapply(report$recommendations, function(rec) tags$li(rec))
      ),
      
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  # =============================================================================
  # CONDITIONAL PANEL OUTPUTS
  # =============================================================================
  
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
  
  # =============================================================================
  # QUICK STATS
  # =============================================================================
  
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
  
  # =============================================================================
  # DATA PREVIEW
  # =============================================================================
  
  output$data_preview <- DT::renderDataTable({
    if (!is.null(values$expression_data)) {
      max_rows <- get_config("ui", "max_preview_rows")
      max_cols <- get_config("ui", "max_preview_cols")
      
      preview_data <- values$expression_data[
        1:min(max_rows, nrow(values$expression_data)), 
        1:min(max_cols, ncol(values$expression_data))
      ]
      preview_data <- cbind(Gene = rownames(preview_data), preview_data)
      
      DT::datatable(
        preview_data,
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          dom = 'tip'
        )
      )
    }
  })
  
  output$processing_stats <- renderUI({
    if (!is.null(values$expression_data)) {
      div(
        h6("📊 Dataset Information"),
        p(paste("Genes:", formatC(nrow(values$expression_data), format="d", big.mark=","))),
        p(paste("Samples:", ncol(values$expression_data))),
        p(paste("Size:", round(object.size(values$expression_data) / (1024^2), 1), "MB")),
        
        br(),
        
        h6("💾 Memory Impact"),
        p(paste("Current Usage:", get_memory_status()$used_mb, "MB")),
        p(paste("Status:", get_memory_status()$recommendation))
      )
    }
  })
  
  # =============================================================================
  # PROGRESS TRACKING  
  # =============================================================================
  
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
    
    # Update progress bar with fallback
    if (requireNamespace("shinyWidgets", quietly = TRUE)) {
      shinyWidgets::updateProgressBar(session, "overall_progress", value = progress)
    }
    
    values$progress <- progress
  })
  
  # Add test data and clear data buttons for demonstration
  output$test_data_ui <- renderUI({
    if (is_development_mode()) {
      div(
        actionButton("test_data", "Load Test Data", class = "btn-secondary"),
        br(), br(),
        actionButton("clear_data", "Clear All Data", class = "btn-warning")
      )
    }
  })
  
  # Clear data functionality
  observeEvent(input$clear_data, {
    if (is_development_mode()) {
      values$expression_data <- NULL
      values$annotation_data <- NULL
      showNotification("🗑️ All data cleared successfully", type = "message")
      cat("🗑️ Data cleared by user request\n")
    }
  })
}

# =============================================================================
# RUN APPLICATION
# =============================================================================

cat("🚀 Starting Prairie Genomics Suite - Optimized Version\n")
cat("💾 Memory management: ENABLED\n")
cat("⚡ Performance monitoring: ENABLED\n") 
cat("🔧 Development mode:", is_development_mode(), "\n")
cat("🌥️ Cloud deployment:", is_cloud_deployment(), "\n\n")

# Run the optimized application
shinyApp(ui = ui, server = server)