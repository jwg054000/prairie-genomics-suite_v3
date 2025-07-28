# Memory Management Utilities - Optimized Version
# Advanced memory monitoring and management for large genomics datasets
# 
# Author: Prairie Genomics Team - Optimized Version
# Features: Smart garbage collection, memory alerts, optimization recommendations

# Load configuration
source("config/app_config.R")

# Memory Manager R6 Class
MemoryManager <- R6::R6Class("MemoryManager",
  public = list(
    
    # Initialize memory manager
    initialize = function() {
      private$start_time <- Sys.time()
      private$gc_count <- 0
      private$memory_history <- data.frame(
        timestamp = numeric(0),
        used_mb = numeric(0),
        operation = character(0)
      )
      
      cat("🧠 Memory Manager initialized\n")
    },
    
    # Get current memory status
    get_status = function() {
      tryCatch({
        gc_info <- gc()
        memory_mb <- round(sum(gc_info[, 2]), 1)
        
        status <- list(
          used_mb = memory_mb,
          warning = memory_mb > MEMORY_WARNING_MB,
          critical = memory_mb > MEMORY_CRITICAL_MB,
          recommendation = private$get_recommendation(memory_mb),
          trend = private$get_trend(),
          last_gc = private$gc_count
        )
        
        # Update history
        private$update_history(memory_mb, "status_check")
        
        return(status)
        
      }, error = function(e) {
        return(list(
          used_mb = 0,
          warning = FALSE,
          critical = FALSE,
          recommendation = "Memory monitoring unavailable",
          error = e$message
        ))
      })
    },
    
    # Smart garbage collection
    smart_gc = function(force = FALSE) {
      current_status <- self$get_status()
      
      # Decide whether to run GC
      should_gc <- force || 
                   current_status$warning || 
                   (private$gc_count %% get_config("memory", "gc_frequency") == 0)
      
      if (should_gc) {
        before_mb <- current_status$used_mb
        
        # Run garbage collection
        gc_result <- gc()
        private$gc_count <- private$gc_count + 1
        
        after_status <- self$get_status()
        freed_mb <- before_mb - after_status$used_mb
        
        private$update_history(after_status$used_mb, "garbage_collection")
        
        cat("🗑️ Garbage collection completed. Freed:", round(freed_mb, 1), "MB\n")
        
        return(list(
          success = TRUE,
          freed_mb = freed_mb,
          current_mb = after_status$used_mb,
          gc_count = private$gc_count
        ))
      } else {
        return(list(
          success = FALSE,
          reason = "Garbage collection not needed",
          current_mb = current_status$used_mb
        ))
      }
    },
    
    # Monitor operation and suggest optimizations
    monitor_operation = function(operation_name, func) {
      before_status <- self$get_status()
      before_time <- Sys.time()
      
      cat("📊 Monitoring operation:", operation_name, "\n")
      cat("   Memory before:", before_status$used_mb, "MB\n")
      
      # Execute operation
      result <- tryCatch({
        func()
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })
      
      after_status <- self$get_status()
      duration <- as.numeric(difftime(Sys.time(), before_time, units = "secs"))
      memory_delta <- after_status$used_mb - before_status$used_mb
      
      # Log operation
      private$update_history(after_status$used_mb, operation_name)
      
      # Provide feedback
      cat("   Memory after:", after_status$used_mb, "MB (Δ", 
          ifelse(memory_delta > 0, "+", ""), round(memory_delta, 1), "MB)\n")
      cat("   Duration:", round(duration, 2), "seconds\n")
      
      # Suggest optimizations if needed
      if (memory_delta > 100) {
        cat("💡 Suggestion: Operation used significant memory. Consider chunking.\n")
      }
      
      if (after_status$critical) {
        cat("⚠️ CRITICAL: Memory usage is very high. Immediate action needed.\n")
        self$smart_gc(force = TRUE)
      }
      
      return(list(
        result = result,
        memory_delta = memory_delta,
        duration = duration,
        recommendation = after_status$recommendation
      ))
    },
    
    # Get memory usage report
    get_report = function() {
      current_status <- self$get_status()
      uptime <- difftime(Sys.time(), private$start_time, units = "mins")
      
      report <- list(
        current_status = current_status,
        uptime_minutes = round(as.numeric(uptime), 1),
        gc_count = private$gc_count,
        memory_efficiency = private$calculate_efficiency(),
        recommendations = private$get_recommendations()
      )
      
      return(report)
    },
    
    # Clear memory history (for long-running sessions)
    clear_history = function() {
      private$memory_history <- data.frame(
        timestamp = numeric(0),
        used_mb = numeric(0),
        operation = character(0)
      )
      cat("🧹 Memory history cleared\n")
    }
  ),
  
  private = list(
    start_time = NULL,
    gc_count = 0,
    memory_history = NULL,
    
    # Update memory history
    update_history = function(memory_mb, operation) {
      new_entry <- data.frame(
        timestamp = as.numeric(Sys.time()),
        used_mb = memory_mb,
        operation = operation
      )
      
      private$memory_history <- rbind(private$memory_history, new_entry)
      
      # Keep only recent history (last 100 entries)
      if (nrow(private$memory_history) > 100) {
        private$memory_history <- tail(private$memory_history, 100)
      }
    },
    
    # Get memory trend
    get_trend = function() {
      if (nrow(private$memory_history) < 3) {
        return("insufficient_data")
      }
      
      recent_entries <- tail(private$memory_history, 5)
      trend_slope <- lm(used_mb ~ timestamp, data = recent_entries)$coefficients[2]
      
      if (abs(trend_slope) < 1) {
        return("stable")
      } else if (trend_slope > 0) {
        return("increasing")
      } else {
        return("decreasing")
      }
    },
    
    # Get memory recommendation
    get_recommendation = function(memory_mb) {
      if (memory_mb > MEMORY_CRITICAL_MB) {
        return("CRITICAL: Restart session or reduce data size immediately")
      } else if (memory_mb > MEMORY_WARNING_MB) {
        return("WARNING: Consider running garbage collection or processing in chunks")
      } else if (memory_mb > MEMORY_WARNING_MB * 0.7) {
        return("CAUTION: Monitor memory usage closely")
      } else {
        return("OK: Memory usage is within normal limits")
      }
    },
    
    # Calculate memory efficiency
    calculate_efficiency = function() {
      if (nrow(private$memory_history) < 2) {
        return(100)
      }
      
      # Simple efficiency metric based on memory stability
      recent_history <- tail(private$memory_history, 10)
      memory_variance <- var(recent_history$used_mb)
      
      # Lower variance = higher efficiency
      efficiency <- max(0, 100 - sqrt(memory_variance))
      return(round(efficiency, 1))
    },
    
    # Get comprehensive recommendations
    get_recommendations = function() {
      current_status <- self$get_status()
      recommendations <- character(0)
      
      if (current_status$critical) {
        recommendations <- c(recommendations, "Immediate restart recommended")
      }
      
      if (current_status$warning) {
        recommendations <- c(recommendations, "Run garbage collection")
        recommendations <- c(recommendations, "Process data in smaller chunks")
      }
      
      if (private$gc_count > 20) {
        recommendations <- c(recommendations, "Consider restarting session for optimal performance")
      }
      
      trend <- private$get_trend()
      if (trend == "increasing") {
        recommendations <- c(recommendations, "Memory usage trending upward - monitor closely")
      }
      
      if (length(recommendations) == 0) {
        recommendations <- "No specific recommendations - memory usage is optimal"
      }
      
      return(recommendations)
    }
  )
)

# Global memory manager instance
memory_manager <- NULL

# Initialize global memory manager
init_memory_manager <- function() {
  if (is.null(memory_manager)) {
    memory_manager <<- MemoryManager$new()
  }
  return(memory_manager)
}

# Convenience functions for global memory manager
get_memory_status <- function() {
  if (is.null(memory_manager)) {
    init_memory_manager()
  }
  return(memory_manager$get_status())
}

smart_gc <- function(force = FALSE) {
  if (is.null(memory_manager)) {
    init_memory_manager()
  }
  return(memory_manager$smart_gc(force))
}

monitor_operation <- function(operation_name, func) {
  if (is.null(memory_manager)) {
    init_memory_manager()
  }
  return(memory_manager$monitor_operation(operation_name, func))
}

get_memory_report <- function() {
  if (is.null(memory_manager)) {
    init_memory_manager()
  }
  return(memory_manager$get_report())
}

# Memory-aware data processing wrapper
process_with_memory_management <- function(data, processing_func, operation_name = "data_processing") {
  monitor_operation(operation_name, function() {
    processing_func(data)
  })
}

# Auto-optimization based on memory status
auto_optimize_memory <- function() {
  status <- get_memory_status()
  
  if (status$critical) {
    cat("🚨 Critical memory usage detected - running emergency cleanup\n")
    smart_gc(force = TRUE)
    
    # Additional cleanup measures
    if (exists("temp_data")) rm(temp_data, envir = .GlobalEnv)
    if (exists("cache_data")) rm(cache_data, envir = .GlobalEnv)
    
    return("emergency_cleanup_performed")
  } else if (status$warning) {
    cat("⚠️ High memory usage - running optimization\n")
    smart_gc()
    return("optimization_performed")
  } else {
    return("no_optimization_needed")
  }
}

cat("✅ Memory management utilities loaded\n")