# Prairie Genomics Suite - Application Configuration
# Centralized configuration management for better maintainability
# 
# Author: Prairie Genomics Team - Optimized Version
# Date: 2025

# Global Application Configuration
APP_CONFIG <- list(
  
  # Application metadata
  app = list(
    name = "Prairie Genomics Suite",
    version = "2.0.0-optimized",
    description = "Enhanced interactive genomics analysis platform"
  ),
  
  # Memory management settings
  memory = list(
    chunk_size = 5000,                    # Default chunk size for large datasets
    memory_warning_mb = 1000,             # Warning threshold in MB
    memory_critical_mb = 2000,            # Critical threshold in MB
    max_file_size_mb = 500,               # Maximum upload file size
    gc_frequency = 10                     # Garbage collection frequency (operations)
  ),
  
  # Analysis parameters
  analysis = list(
    # DESeq2 defaults
    deseq2 = list(
      default_padj = 0.05,
      default_fc = 1.0,
      min_samples = 3,
      max_genes_for_analysis = 50000
    ),
    
    # Pathway analysis defaults
    pathway = list(
      default_species = "human",
      max_pathways_plot = 50,
      cache_duration_hours = 24,
      supported_analyses = c("GO", "KEGG", "GSEA", "MSigDB", "Reactome")
    )
  ),
  
  # UI/UX settings
  ui = list(
    max_preview_rows = 100,
    max_preview_cols = 20,
    progress_update_interval = 1000,      # milliseconds
    notification_duration = 5000,         # milliseconds
    plot_height = "600px",
    plot_width = "800px"
  ),
  
  # Performance settings
  performance = list(
    enable_async = TRUE,
    enable_caching = TRUE,
    cache_size_mb = 1000,
    plot_timeout_seconds = 30,
    analysis_timeout_minutes = 60
  ),
  
  # Package management
  packages = list(
    essential = c(
      "shiny", "shinydashboard", "DT", "ggplot2", 
      "dplyr", "readr", "plotly"
    ),
    
    bioconductor = c(
      "DESeq2", "biomaRt", "clusterProfiler", 
      "org.Hs.eg.db", "org.Mm.eg.db"
    ),
    
    optional = c(
      "shinyWidgets", "readxl", "RColorBrewer", 
      "pheatmap", "ggrepel", "enrichplot", "pathview"
    )
  ),
  
  # Cloud deployment settings
  deployment = list(
    cloud_memory_limit_mb = 1024,
    cloud_timeout_minutes = 30,
    fallback_mode = TRUE,                 # Enable simplified algorithms for cloud
    disable_heavy_features = TRUE         # Disable memory-intensive features
  ),
  
  # Logging configuration
  logging = list(
    level = "INFO",                       # DEBUG, INFO, WARN, ERROR
    log_to_file = FALSE,
    max_log_size_mb = 100,
    enable_performance_logging = TRUE
  )
)

# Environment detection functions
is_cloud_deployment <- function() {
  return(
    Sys.getenv("SHINY_SERVER") != "" || 
    Sys.getenv("SHINYAPPS_IO") != "" ||
    Sys.getenv("CONNECT_SERVER") != ""
  )
}

is_development_mode <- function() {
  return(Sys.getenv("R_ENV") == "development" || interactive())
}

# Configuration getter functions
get_config <- function(section = NULL, key = NULL) {
  if (is.null(section)) {
    return(APP_CONFIG)
  }
  
  if (is.null(key)) {
    return(APP_CONFIG[[section]])
  }
  
  return(APP_CONFIG[[section]][[key]])
}

# Dynamic configuration based on environment
configure_for_environment <- function() {
  if (is_cloud_deployment()) {
    # Adjust settings for cloud deployment
    APP_CONFIG$memory$chunk_size <<- 2000
    APP_CONFIG$memory$memory_warning_mb <<- 500
    APP_CONFIG$memory$memory_critical_mb <<- 800
    APP_CONFIG$performance$cache_size_mb <<- 200
    APP_CONFIG$ui$max_preview_rows <<- 50
    
    cat("🌥️ Cloud deployment detected - optimized configuration applied\n")
  }
  
  if (is_development_mode()) {
    APP_CONFIG$logging$level <<- "DEBUG"
    APP_CONFIG$logging$enable_performance_logging <<- TRUE
    
    cat("🔧 Development mode detected - debug configuration applied\n")
  }
}

# Initialize configuration
configure_for_environment()

# Utility functions and operators
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

# Fallback progressBar function when shinyWidgets is not available
progressBar <- function(id, value = 0, total = 100, status = "primary", 
                       display_pct = TRUE, striped = TRUE, ...) {
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::progressBar(id = id, value = value, total = total, 
                             status = status, display_pct = display_pct, 
                             striped = striped, ...)
  } else {
    # Fallback to basic HTML progress bar
    div(
      class = "progress",
      div(
        class = paste("progress-bar", if(striped) "progress-bar-striped"),
        role = "progressbar",
        style = paste0("width: ", round(100 * value / total), "%"),
        id = id,
        if (display_pct) paste0(round(100 * value / total), "%")
      )
    )
  }
}

# Export key configuration values for easy access
CHUNK_SIZE <- get_config("memory", "chunk_size")
MEMORY_WARNING_MB <- get_config("memory", "memory_warning_mb")
MEMORY_CRITICAL_MB <- get_config("memory", "memory_critical_mb")
DEFAULT_PADJ <- get_config("analysis", "deseq2")$default_padj
DEFAULT_FC <- get_config("analysis", "deseq2")$default_fc

cat("✅ Prairie Genomics Suite configuration loaded\n")
cat("   Version:", get_config("app", "version"), "\n")
cat("   Environment:", if(is_cloud_deployment()) "Cloud" else "Local", "\n")
cat("   Memory chunk size:", CHUNK_SIZE, "genes\n")