# Debug Logs Viewer - Simple interface to view error logs
# This file provides functions to easily view and analyze logged errors

# Simple function to view recent logs
view_recent_logs <- function(lines = 100) {
  cat("📋 RECENT LOGS (Last", lines, "lines)\n")
  cat("=====================================\n\n")
  
  recent_logs <- get_recent_logs(lines)
  cat(recent_logs)
  cat("\n\n📊 End of logs\n")
}

# Function to search logs for specific patterns
search_logs <- function(pattern, context_lines = 3) {
  if (!dir.exists(LOG_CONFIG$log_dir)) {
    cat("❌ No logs directory found\n")
    return(invisible(NULL))
  }
  
  log_files <- list.files(LOG_CONFIG$log_dir, pattern = "prairie_genomics_.*\\.log$", full.names = TRUE)
  
  if (length(log_files) == 0) {
    cat("❌ No log files found\n")
    return(invisible(NULL))
  }
  
  cat("🔍 SEARCHING LOGS FOR:", pattern, "\n")
  cat("=====================================\n\n")
  
  found_matches <- FALSE
  
  for (log_file in log_files) {
    tryCatch({
      # Use grep to search with context
      cmd <- paste0("grep -i -C ", context_lines, " '", pattern, "' '", log_file, "'")
      result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
      
      if (length(result) > 0) {
        cat("📁 File:", basename(log_file), "\n")
        cat(paste(result, collapse = "\n"), "\n\n")
        found_matches <- TRUE
      }
    }, error = function(e) {
      # Fallback for Windows or if grep fails
      lines <- readLines(log_file)
      matches <- grep(pattern, lines, ignore.case = TRUE)
      
      if (length(matches) > 0) {
        cat("📁 File:", basename(log_file), "\n")
        for (match_line in matches) {
          start_line <- max(1, match_line - context_lines)
          end_line <- min(length(lines), match_line + context_lines)
          context <- lines[start_line:end_line]
          cat(paste(context, collapse = "\n"), "\n\n")
        }
        found_matches <- TRUE
      }
    })
  }
  
  if (!found_matches) {
    cat("❌ No matches found for pattern:", pattern, "\n")
  }
}

# Function to show error summary
show_error_summary <- function(days = 1) {
  if (!dir.exists(LOG_CONFIG$log_dir)) {
    cat("❌ No logs directory found\n")  
    return(invisible(NULL))
  }
  
  # Get recent log files
  cutoff_date <- Sys.Date() - days + 1
  log_files <- list.files(LOG_CONFIG$log_dir, pattern = "prairie_genomics_.*\\.log$", full.names = TRUE)
  
  cat("📊 ERROR SUMMARY (Last", days, "day(s))\n")
  cat("=====================================\n\n")
  
  total_errors <- 0
  total_warnings <- 0
  error_categories <- list()
  
  for (log_file in log_files) {
    # Extract date from filename  
    file_date_str <- regmatches(basename(log_file), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", basename(log_file)))
    
    if (length(file_date_str) > 0) {
      file_date <- as.Date(file_date_str)
      if (file_date >= cutoff_date) {
        tryCatch({
          lines <- readLines(log_file)
          
          # Count errors and warnings
          errors <- grep("\\[ERROR\\]", lines, ignore.case = TRUE)
          warnings <- grep("\\[WARNING\\]", lines, ignore.case = TRUE)
          
          total_errors <- total_errors + length(errors)
          total_warnings <- total_warnings + length(warnings)
          
          # Extract error categories
          for (error_line in lines[errors]) {
            # Extract category from log format
            category_match <- regmatches(error_line, regexpr("\\[ERROR\\] [A-Z_]+", error_line))
            if (length(category_match) > 0) {
              category <- gsub("\\[ERROR\\] ", "", category_match)
              error_categories[[category]] <- (error_categories[[category]] %||% 0) + 1
            }
          }
          
        }, error = function(e) {
          cat("⚠️ Error reading file:", basename(log_file), "\n")
        })
      }
    }
  }
  
  cat("📈 SUMMARY STATISTICS:\n")
  cat("   - Total Errors:", total_errors, "\n")
  cat("   - Total Warnings:", total_warnings, "\n")
  
  if (length(error_categories) > 0) {
    cat("\n🔍 ERROR CATEGORIES:\n")
    sorted_categories <- sort(unlist(error_categories), decreasing = TRUE)
    for (i in 1:min(10, length(sorted_categories))) {
      cat("   -", names(sorted_categories)[i], ":", sorted_categories[i], "\n")
    }
  }
  
  cat("\n💡 TIP: Use search_logs('ERROR') to see detailed error messages\n")
}

# Function to clear old logs (maintenance)
clean_old_logs <- function(days_to_keep = 7) {
  if (!dir.exists(LOG_CONFIG$log_dir)) {
    cat("❌ No logs directory found\n")
    return(invisible(NULL))
  }
  
  cutoff_date <- Sys.Date() - days_to_keep
  log_files <- list.files(LOG_CONFIG$log_dir, pattern = "prairie_genomics_.*\\.log$", full.names = TRUE)
  
  removed_count <- 0
  
  for (log_file in log_files) {
    file_date_str <- regmatches(basename(log_file), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", basename(log_file)))
    
    if (length(file_date_str) > 0) {
      file_date <- as.Date(file_date_str)
      if (file_date < cutoff_date) {
        file.remove(log_file)
        removed_count <- removed_count + 1
        cat("🗑️ Removed old log file:", basename(log_file), "\n")
      }
    }
  }
  
  cat("✅ Cleanup completed:", removed_count, "old log files removed\n")
}

# Quick debugging commands
debug_recent_errors <- function() {
  search_logs("ERROR", context_lines = 2)
}

debug_pathway_issues <- function() {
  search_logs("PATHWAY_ANALYSIS", context_lines = 3)
}

debug_user_actions <- function() {
  search_logs("USER_ACTION", context_lines = 1)
}

# Export main functions
cat("🔧 Debug Log Viewer Functions Available:\n")
cat("   - view_recent_logs(lines = 100)\n")
cat("   - search_logs('pattern', context_lines = 3)\n")
cat("   - show_error_summary(days = 1)\n")
cat("   - debug_recent_errors()\n")
cat("   - debug_pathway_issues()\n")
cat("   - debug_user_actions()\n")
cat("   - clean_old_logs(days_to_keep = 7)\n")
cat("\n💡 Example: search_logs('KEGG') to find KEGG-related issues\n")