#!/usr/bin/env Rscript

# Robust launch script for Phase 3 Guardrails Demo
# with better error handling and stability

cat("🚀 Starting Prairie Genomics Suite - Phase 3 Guardrails Demo\n")
cat("📊 Real Data Testing System\n")
cat("=" , rep("=", 50), "\n", sep="")

# Load required libraries with error handling
tryCatch({
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(ggplot2)
  library(plotly)
  cat("✅ All libraries loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading libraries:", e$message, "\n")
  quit(status = 1)
})

# Check if demo file exists
if (!file.exists("phase3_guardrails_demo.R")) {
  cat("❌ Error: phase3_guardrails_demo.R not found in current directory\n")
  cat("Current directory:", getwd(), "\n")
  quit(status = 1)
}

# Set options for better stability
options(
  shiny.port = 3838,
  shiny.host = "127.0.0.1",
  shiny.maxRequestSize = 200*1024^2,  # 200MB max upload
  shiny.reactlog = FALSE,
  shiny.trace = FALSE,
  shiny.error = function() {
    cat("❌ Shiny error occurred\n")
    traceback()
  }
)

cat("🌐 App will be available at: http://127.0.0.1:3838\n")
cat("💡 Navigate to '📊 Real Data Testing' tab to upload your data\n")
cat("⏹️  Press Ctrl+C to stop the app\n")
cat("=" , rep("=", 50), "\n", sep="")

# Launch the app with error handling
tryCatch({
  shiny::runApp(
    "phase3_guardrails_demo.R",
    host = "127.0.0.1",
    port = 3838,
    launch.browser = TRUE
  )
}, error = function(e) {
  cat("❌ Error running app:", e$message, "\n")
  traceback()
}, interrupt = function(e) {
  cat("\n👋 App stopped by user\n")
})

cat("📊 Prairie Genomics Suite session ended\n")