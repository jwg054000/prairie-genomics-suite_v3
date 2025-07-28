# Phase 2 - Real-time Updates & Notifications
# Live progress tracking, notifications, and real-time UI updates
# 
# This module provides:
# - Real-time progress tracking for all operations
# - Toast notifications for completion and errors
# - Live status monitoring and updates
# - WebSocket-style real-time communication
# - Performance metrics display

library(shiny)
library(jsonlite)
library(promises)

# ===========================================
# REAL-TIME PROGRESS TRACKING SYSTEM
# ===========================================

# Create real-time progress tracker
create_realtime_progress_tracker <- function(session, tracker_id = "main_tracker") {
  
  # Progress state management
  progress_state <- reactiveValues(
    active_processes = list(),
    completed_processes = list(),
    failed_processes = list(),
    global_status = "idle"
  )
  
  # Progress tracker methods
  tracker <- list(
    
    # Start tracking a new process
    start_process = function(process_id, process_name, estimated_duration = NULL) {
      
      process_info <- list(
        id = process_id,
        name = process_name,
        status = "running",
        progress = 0,
        message = "Starting...",
        start_time = Sys.time(),
        estimated_duration = estimated_duration,
        steps = list(),
        current_step = 1
      )
      
      progress_state$active_processes[[process_id]] <- process_info
      progress_state$global_status <- "active"
      
      # Send real-time update to UI
      session$sendCustomMessage(
        type = "updateGlobalProgress",
        message = list(
          status = "active",
          active_count = length(progress_state$active_processes),
          total_count = length(progress_state$active_processes) + 
                       length(progress_state$completed_processes) +
                       length(progress_state$failed_processes)
        )
      )
      
      # Send process start notification
      session$sendCustomMessage(
        type = "showRealtimeNotification",
        message = list(
          type = "info",
          title = "Process Started",
          message = paste("Starting", process_name),
          process_id = process_id,
          duration = 3000
        )
      )
      
      cat(sprintf("[%s] Started process: %s (%s)\n", 
                 format(Sys.time(), "%H:%M:%S"), process_name, process_id))
    },
    
    # Update process progress
    update_progress = function(process_id, progress, message = NULL, step_info = NULL) {
      
      if (!process_id %in% names(progress_state$active_processes)) {
        warning(paste("Process", process_id, "not found in active processes"))
        return(FALSE)
      }
      
      # Update process information
      process <- progress_state$active_processes[[process_id]]
      process$progress <- max(0, min(100, progress))
      
      if (!is.null(message)) {
        process$message <- message
      }
      
      if (!is.null(step_info)) {
        process$steps <- append(process$steps, list(step_info))
        process$current_step <- length(process$steps)
      }
      
      # Calculate elapsed time and ETA
      elapsed_time <- as.numeric(difftime(Sys.time(), process$start_time, units = "secs"))
      
      if (progress > 5) {  # Only calculate ETA after 5% progress
        estimated_total_time <- elapsed_time / (progress / 100)
        eta_seconds <- estimated_total_time - elapsed_time
        process$eta <- if (eta_seconds > 0) eta_seconds else 0
      }
      
      process$elapsed_time <- elapsed_time
      progress_state$active_processes[[process_id]] <- process
      
      # Send real-time update to UI
      session$sendCustomMessage(
        type = "updateProcessProgress",
        message = list(
          process_id = process_id,
          progress = process$progress,
          message = process$message,
          elapsed_time = round(elapsed_time, 1),
          eta = if (!is.null(process$eta)) round(process$eta, 1) else NULL,
          current_step = process$current_step,
          total_steps = length(process$steps)
        )
      )
      
      cat(sprintf("[%s] %s: %d%% - %s\n", 
                 format(Sys.time(), "%H:%M:%S"), process$name, 
                 process$progress, process$message))
      
      return(TRUE)
    },
    
    # Complete a process successfully
    complete_process = function(process_id, final_message = "Complete", results = NULL) {
      
      if (!process_id %in% names(progress_state$active_processes)) {
        warning(paste("Process", process_id, "not found in active processes"))
        return(FALSE)
      }
      
      # Move from active to completed
      process <- progress_state$active_processes[[process_id]]
      process$status <- "completed"
      process$progress <- 100
      process$message <- final_message
      process$end_time <- Sys.time()
      process$total_duration <- as.numeric(difftime(process$end_time, process$start_time, units = "secs"))
      
      if (!is.null(results)) {
        process$results <- results
      }
      
      progress_state$completed_processes[[process_id]] <- process
      progress_state$active_processes[[process_id]] <- NULL
      
      # Update global status
      if (length(progress_state$active_processes) == 0) {
        progress_state$global_status <- "idle"
      }
      
      # Send completion notification
      session$sendCustomMessage(
        type = "showRealtimeNotification",
        message = list(
          type = "success",
          title = "Process Complete",
          message = paste(process$name, "completed successfully"),
          process_id = process_id,
          duration = 5000
        )
      )
      
      # Send final progress update
      session$sendCustomMessage(
        type = "completeProcess",
        message = list(
          process_id = process_id,
          total_duration = round(process$total_duration, 1),
          final_message = final_message
        )
      )
      
      cat(sprintf("[%s] Completed process: %s (%.1fs)\n", 
                 format(Sys.time(), "%H:%M:%S"), process$name, process$total_duration))
      
      return(TRUE)
    },
    
    # Fail a process with error
    fail_process = function(process_id, error_message, error_details = NULL) {
      
      if (!process_id %in% names(progress_state$active_processes)) {
        warning(paste("Process", process_id, "not found in active processes"))
        return(FALSE)
      }
      
      # Move from active to failed
      process <- progress_state$active_processes[[process_id]]
      process$status <- "failed"
      process$message <- error_message
      process$end_time <- Sys.time()
      process$total_duration <- as.numeric(difftime(process$end_time, process$start_time, units = "secs"))
      
      if (!is.null(error_details)) {
        process$error_details <- error_details
      }
      
      progress_state$failed_processes[[process_id]] <- process
      progress_state$active_processes[[process_id]] <- NULL
      
      # Update global status
      if (length(progress_state$active_processes) == 0) {
        progress_state$global_status <- "idle"
      }
      
      # Send error notification
      session$sendCustomMessage(
        type = "showRealtimeNotification",
        message = list(
          type = "error",
          title = "Process Failed",
          message = paste(process$name, "failed:", error_message),
          process_id = process_id,
          duration = 10000
        )
      )
      
      # Send failure update
      session$sendCustomMessage(
        type = "failProcess",
        message = list(
          process_id = process_id,
          error_message = error_message,
          total_duration = round(process$total_duration, 1)
        )
      )
      
      cat(sprintf("[%s] Failed process: %s - %s\n", 
                 format(Sys.time(), "%H:%M:%S"), process$name, error_message))
      
      return(TRUE)
    },
    
    # Get current status
    get_status = function() {
      return(list(
        global_status = progress_state$global_status,
        active_processes = progress_state$active_processes,
        completed_processes = progress_state$completed_processes,
        failed_processes = progress_state$failed_processes,
        total_active = length(progress_state$active_processes),
        total_completed = length(progress_state$completed_processes),
        total_failed = length(progress_state$failed_processes)
      ))
    },
    
    # Clear completed and failed processes
    clear_history = function() {
      progress_state$completed_processes <- list()
      progress_state$failed_processes <- list()
      
      session$sendCustomMessage(
        type = "clearProcessHistory",
        message = list()
      )
    }
  )
  
  return(tracker)
}

# ===========================================
# REAL-TIME NOTIFICATION SYSTEM
# ===========================================

# Create notification manager
create_notification_manager <- function(session) {
  
  # Notification queue
  notification_queue <- reactiveValues(
    pending = list(),
    displayed = list(),
    max_displayed = 5
  )
  
  notification_manager <- list(
    
    # Show success notification
    success = function(title, message, duration = 5000, action = NULL) {
      notification_manager$show("success", title, message, duration, action)
    },
    
    # Show error notification
    error = function(title, message, duration = 10000, action = NULL) {
      notification_manager$show("error", title, message, duration, action)
    },
    
    # Show warning notification
    warning = function(title, message, duration = 7000, action = NULL) {
      notification_manager$show("warning", title, message, duration, action)
    },
    
    # Show info notification
    info = function(title, message, duration = 4000, action = NULL) {
      notification_manager$show("info", title, message, duration, action)
    },
    
    # Show custom notification
    show = function(type, title, message, duration = 5000, action = NULL) {
      
      notification_id <- paste0("notif_", as.integer(as.numeric(Sys.time()) * 1000))
      
      notification <- list(
        id = notification_id,
        type = type,
        title = title,
        message = message,
        duration = duration,
        action = action,
        timestamp = Sys.time(),
        shown = FALSE
      )
      
      # Add to queue
      notification_queue$pending[[notification_id]] <- notification
      
      # Process queue
      notification_manager$process_queue()
    },
    
    # Process notification queue
    process_queue = function() {
      
      # Check if we can show more notifications
      displayed_count <- length(notification_queue$displayed)
      
      if (displayed_count < notification_queue$max_displayed && 
          length(notification_queue$pending) > 0) {
        
        # Get next notification
        next_notif_id <- names(notification_queue$pending)[1]
        notification <- notification_queue$pending[[next_notif_id]]
        
        # Send to UI
        session$sendCustomMessage(
          type = "showRealtimeNotification",
          message = list(
            id = notification$id,
            type = notification$type,
            title = notification$title,
            message = notification$message,
            duration = notification$duration,
            action = notification$action
          )
        )
        
        # Move from pending to displayed
        notification_queue$displayed[[next_notif_id]] <- notification
        notification_queue$pending[[next_notif_id]] <- NULL
        
        # Schedule removal
        later::later(function() {
          notification_manager$remove(next_notif_id)
        }, delay = notification$duration / 1000)
      }
    },
    
    # Remove notification
    remove = function(notification_id) {
      if (notification_id %in% names(notification_queue$displayed)) {
        notification_queue$displayed[[notification_id]] <- NULL
        
        # Process queue for next notification
        notification_manager$process_queue()
      }
    },
    
    # Clear all notifications
    clear_all = function() {
      notification_queue$pending <- list()
      notification_queue$displayed <- list()
      
      session$sendCustomMessage(
        type = "clearAllNotifications",
        message = list()
      )
    }
  )
  
  return(notification_manager)
}

# ===========================================
# SYSTEM STATUS MONITOR
# ===========================================

# Create system status monitor
create_system_monitor <- function(session, update_interval = 5) {
  
  # System metrics
  system_metrics <- reactiveValues(
    memory_usage = 0,
    memory_limit = 1000,  # MB
    active_sessions = 1,
    server_load = 0,
    last_update = Sys.time()
  )
  
  # Start monitoring loop
  observe({
    invalidateLater(update_interval * 1000, session)
    
    # Update system metrics
    tryCatch({
      # Memory usage
      gc_info <- gc()
      memory_used <- sum(gc_info[, 2])  # Used memory in MB
      system_metrics$memory_usage <- memory_used
      
      # Server load (simplified)
      system_metrics$server_load <- sample(0:100, 1)  # Mock data
      
      # Update timestamp
      system_metrics$last_update <- Sys.time()
      
      # Send to UI
      session$sendCustomMessage(
        type = "updateSystemMetrics",
        message = list(
          memory_usage = round(memory_used, 1),
          memory_percent = round(100 * memory_used / system_metrics$memory_limit, 1),
          server_load = system_metrics$server_load,
          active_sessions = system_metrics$active_sessions,
          last_update = format(system_metrics$last_update, "%H:%M:%S")
        )
      )
      
      # Alert on high memory usage
      memory_percent <- 100 * memory_used / system_metrics$memory_limit
      if (memory_percent > 80) {
        session$sendCustomMessage(
          type = "showRealtimeNotification",
          message = list(
            type = "warning",
            title = "High Memory Usage",
            message = sprintf("Memory usage at %.1f%%. Consider restarting if performance degrades.", memory_percent),
            duration = 8000
          )
        )
      }
      
    }, error = function(e) {
      cat("System monitoring error:", e$message, "\n")
    })
  })
  
  return(system_metrics)
}

# ===========================================
# REAL-TIME UI COMPONENTS
# ===========================================

# Create real-time status dashboard
create_realtime_dashboard <- function(progress_tracker, notification_manager, system_monitor) {
  
  tagList(
    # Real-time status header
    div(
      id = "realtime-status-header",
      class = "realtime-status-header",
      style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
               color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
      
      fluidRow(
        column(3,
          div(
            class = "status-metric",
            h4("System Status", style = "margin: 0; font-size: 16px;"),
            div(id = "global-status", class = "status-indicator", "●", 
                style = "color: #28a745; font-size: 20px; margin-top: 5px;")
          )
        ),
        column(3,
          div(
            class = "status-metric",
            h4("Active Processes", style = "margin: 0; font-size: 16px;"),
            div(id = "active-count", class = "metric-value", "0",
                style = "font-size: 24px; font-weight: bold; margin-top: 5px;")
          )
        ),
        column(3,
          div(
            class = "status-metric", 
            h4("Memory Usage", style = "margin: 0; font-size: 16px;"),
            div(id = "memory-usage", class = "metric-value", "0%",
                style = "font-size: 24px; font-weight: bold; margin-top: 5px;")
          )
        ),
        column(3,
          div(
            class = "status-metric",
            h4("Last Update", style = "margin: 0; font-size: 16px;"),
            div(id = "last-update", class = "metric-value", "--:--:--",
                style = "font-size: 18px; margin-top: 5px;")
          )
        )
      )
    ),
    
    # Active processes panel
    div(
      id = "active-processes-panel",
      class = "active-processes-panel",
      style = "display: none; background: white; border: 1px solid #ddd; 
               border-radius: 8px; padding: 20px; margin-bottom: 20px;",
      
      h4("Active Processes", style = "margin-top: 0; color: #333;"),
      div(id = "active-processes-list")
    ),
    
    # Notification container
    div(
      id = "realtime-notifications",
      class = "realtime-notifications",
      style = "position: fixed; top: 20px; right: 20px; z-index: 9999; width: 350px;"
    )
  )
}

# ===========================================
# EXPORT FUNCTIONS
# ===========================================

list(
  create_realtime_progress_tracker = create_realtime_progress_tracker,
  create_notification_manager = create_notification_manager,
  create_system_monitor = create_system_monitor,
  create_realtime_dashboard = create_realtime_dashboard
)