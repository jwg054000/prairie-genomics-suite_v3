# Quick launcher for Phase 2 Demo App
# This ensures all dependencies are loaded properly

cat("🚀 Launching Prairie Genomics Suite Phase 2 Demo...\n")
cat("==================================================\n\n")

# Check for required packages
required_packages <- c("shiny", "shinydashboard", "DT", "promises", "future", "shinyjs")

# Install missing packages if needed
missing <- setdiff(required_packages, installed.packages()[, "Package"])
if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cran.rstudio.com/")
}

# Load all required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(promises)
library(future)
library(shinyjs)

# Set up async processing
plan(multisession, workers = 2)

cat("✅ All dependencies loaded successfully\n")
cat("✅ Async processing configured\n\n")

cat("📋 PHASE 2 DEMO FEATURES:\n")
cat("   • ⚡ Async Analysis - Non-blocking DESeq2 processing\n")
cat("   • 📊 Real-time Updates - Live notifications and progress\n")
cat("   • 🔐 Authentication - Firebase integration demo\n")
cat("   • 🎨 Advanced UI - Modern components and interactions\n\n")

cat("🌐 Opening demo app in your browser...\n\n")

# Run the demo app
shiny::runApp('phase2_demo_app.R', launch.browser = TRUE)