# Prairie Genomics Suite - Advanced Logging System
# Comprehensive error logging for easier debugging and maintenance
# 
# Features:
# - Multiple log levels (ERROR, WARNING, INFO, DEBUG)
# - Session tracking with unique IDs
# - User action logging
# - System information capture
# - File rotation to prevent huge logs
# - Stack trace capture for errors
# - Performance monitoring

# Load required packages
if (!require("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
  library(jsonlite)
}

# Global logging configuration
LOG_CONFIG <- list(
  enabled = TRUE,
  log_dir = "logs",
  max_file_size_mb = 10,
  max_files = 5,
  levels = c("ERROR" = 1, "WARNING" = 2, "INFO" = 3, "DEBUG" = 4),
  current_level = 3,  # INFO level by default
  session_id = paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1))
)

# Initialize logging system
initialize_logging <- function() {
  if (!LOG_CONFIG$enabled) return(invisible(NULL))
  
  # Create logs directory if it doesn't exist
  if (!dir.exists(LOG_CONFIG$log_dir)) {
    dir.create(LOG_CONFIG$log_dir, recursive = TRUE)
    cat("📁 Created logs directory:", LOG_CONFIG$log_dir, "\n")
  }
  
  # Log system startup
  log_info("SYSTEM_START", "Prairie Genomics Suite logging system initialized", list(
    session_id = LOG_CONFIG$session_id,
    r_version = R.version.string,
    platform = Sys.info()[["sysname"]],
    user = Sys.info()[["user"]],
    working_dir = getwd(),
    timestamp = Sys.time()
  ))
  
  cat("🔧 Logging system initialized - Session ID:", LOG_CONFIG$session_id, "\n")
}

# Core logging function
write_log <- function(level, category, message, details = NULL, stack_trace = NULL) {
  if (!LOG_CONFIG$enabled) return(invisible(NULL))
  
  # Check if we should log this level
  if (LOG_CONFIG$levels[[level]] > LOG_CONFIG$current_level) {
    return(invisible(NULL))
  }
  
  # Create log entry
  log_entry <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    session_id = LOG_CONFIG$session_id,
    level = level,
    category = category,
    message = message,
    details = details,
    r_version = R.version.string,
    stack_trace = stack_trace
  )
  
  # Add memory usage for performance monitoring
  if (level %in% c("ERROR", "WARNING")) {
    log_entry$memory_mb <- round(as.numeric(object.size(ls(envir = .GlobalEnv))) / 1024^2, 2)
    log_entry$system_memory <- gc()[1, 2]  # Memory usage from gc()
  }
  
  # Determine log file
  log_file <- file.path(LOG_CONFIG$log_dir, paste0("prairie_genomics_", Sys.Date(), ".log"))
  
  # Format for human readability
  formatted_entry <- sprintf(
    "[%s] %s [%s] %s: %s%s%s",
    log_entry$timestamp,
    log_entry$session_id,
    level,
    category,
    message,
    if (!is.null(details)) paste0("\n    Details: ", jsonlite::toJSON(details, auto_unbox = TRUE, pretty = TRUE)) else "",
    if (!is.null(stack_trace)) paste0("\n    Stack Trace:\n", paste("      ", stack_trace, collapse = "\n")) else ""
  )
  
  # Write to file
  tryCatch({
    cat(formatted_entry, "\n\n", file = log_file, append = TRUE)
    
    # Check file size and rotate if needed
    rotate_logs_if_needed()
    
  }, error = function(e) {
    # Fallback - write to console if file writing fails
    cat("LOG WRITE ERROR:", e$message, "\n")
    cat(formatted_entry, "\n")
  })
  
  # Also output to console for immediate feedback
  if (level %in% c("ERROR", "WARNING")) {
    cat("📝", formatted_entry, "\n")
  }
}

# Convenience functions for different log levels
log_error <- function(category, message, details = NULL) {
  # Capture stack trace for errors
  stack_trace <- capture.output(traceback())
  if (length(stack_trace) == 0) {
    stack_trace <- capture.output(sys.calls())
  }
  write_log("ERROR", category, message, details, stack_trace)
}

log_warning <- function(category, message, details = NULL) {
  write_log("WARNING", category, message, details)
}

log_info <- function(category, message, details = NULL) {
  write_log("INFO", category, message, details)
}

log_debug <- function(category, message, details = NULL) {
  write_log("DEBUG", category, message, details)
}

# User action logging
log_user_action <- function(action, inputs = NULL, results = NULL) {
  details <- list(
    action = action,
    inputs = inputs,
    results_summary = if (!is.null(results)) {
      if (is.list(results)) {
        list(
          success = results$success %||% "unknown",
          data_rows = if (!is.null(results$data)) nrow(results$data) else NULL,
          analysis_type = results$analysis_type %||% NULL
        )
      } else {
        "non-list result"
      }
    } else NULL
  )
  
  log_info("USER_ACTION", paste("User performed:", action), details)
}

# Error wrapper for critical functions
with_error_logging <- function(category, description, expr) {
  tryCatch({
    log_debug("FUNCTION_START", paste("Starting:", description))
    start_time <- Sys.time()
    
    result <- expr
    
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    log_info("FUNCTION_SUCCESS", paste("Completed:", description), list(
      duration_seconds = round(duration, 2),
      result_type = class(result)[1]
    ))
    
    result
    
  }, error = function(e) {
    log_error(category, paste("Function failed:", description), list(
      error_message = e$message,
      error_class = class(e)[1],
      function_context = description
    ))
    
    # Re-throw the error
    stop(e)
  })
}

# File rotation
rotate_logs_if_needed <- function() {
  if (!LOG_CONFIG$enabled) return(invisible(NULL))
  
  log_files <- list.files(LOG_CONFIG$log_dir, pattern = "prairie_genomics_.*\\.log$", full.names = TRUE)
  
  for (log_file in log_files) {
    file_size_mb <- file.size(log_file) / 1024^2
    
    if (file_size_mb > LOG_CONFIG$max_file_size_mb) {
      # Rotate the file
      timestamp <- format(Sys.time(), "%H%M%S")
      rotated_name <- gsub("\\.log$", paste0("_", timestamp, ".log"), log_file)
      file.rename(log_file, rotated_name)
      
      log_info("LOG_ROTATION", "Log file rotated", list(
        original_file = basename(log_file),
        rotated_file = basename(rotated_name),
        size_mb = round(file_size_mb, 2)
      ))
    }
  }
  
  # Clean up old files if too many
  all_log_files <- list.files(LOG_CONFIG$log_dir, pattern = "prairie_genomics_.*\\.log$", full.names = TRUE)
  if (length(all_log_files) > LOG_CONFIG$max_files) {
    # Sort by modification time and remove oldest
    file_info <- file.info(all_log_files)
    oldest_files <- head(order(file_info$mtime), length(all_log_files) - LOG_CONFIG$max_files)
    
    for (old_file in all_log_files[oldest_files]) {
      file.remove(old_file)
      log_info("LOG_CLEANUP", "Old log file removed", list(file = basename(old_file)))
    }
  }
}

# System information logging
log_system_info <- function() {
  system_info <- list(
    r_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()[["sysname"]],
    os_version = Sys.info()[["release"]],
    user = Sys.info()[["user"]],
    working_directory = getwd(),
    memory_limit = if (Sys.info()[["sysname"]] == "Windows") memory.limit() else "Unix system - no limit function",
    loaded_packages = (.packages()),
    sys_time = Sys.time(),
    timezone = Sys.timezone()
  )
  
  log_info("SYSTEM_INFO", "System information captured", system_info)
}

# Analysis performance logging
log_analysis_performance <- function(analysis_type, duration_seconds, input_genes, output_pathways = NULL) {
  performance_info <- list(
    analysis_type = analysis_type,
    duration_seconds = round(duration_seconds, 2),
    input_gene_count = length(input_genes),
    output_pathway_count = output_pathways,
    genes_per_second = round(length(input_genes) / duration_seconds, 2),
    memory_used_mb = round(as.numeric(object.size(ls(envir = .GlobalEnv))) / 1024^2, 2)
  )
  
  log_info("PERFORMANCE", paste("Analysis performance:", analysis_type), performance_info)
}

# Export configuration for easy access
get_log_config <- function() {
  return(LOG_CONFIG)
}

# Function to read recent logs (for debugging)
get_recent_logs <- function(lines = 50) {
  if (!LOG_CONFIG$enabled) return("Logging is disabled")
  
  log_file <- file.path(LOG_CONFIG$log_dir, paste0("prairie_genomics_", Sys.Date(), ".log"))
  
  if (!file.exists(log_file)) {
    return("No log file found for today")
  }
  
  # Read last N lines
  tryCatch({
    all_lines <- readLines(log_file)
    recent_lines <- tail(all_lines, lines)
    return(paste(recent_lines, collapse = "\n"))
  }, error = function(e) {
    return(paste("Error reading log file:", e$message))
  })
}

# Initialize when this file is sourced
initialize_logging()
log_system_info()

cat("✅ Prairie Genomics Suite Logging System Ready\n")
cat("📁 Logs directory:", LOG_CONFIG$log_dir, "\n")
cat("🆔 Session ID:", LOG_CONFIG$session_id, "\n")
cat("📊 Log level:", names(LOG_CONFIG$levels)[LOG_CONFIG$current_level], "\n")