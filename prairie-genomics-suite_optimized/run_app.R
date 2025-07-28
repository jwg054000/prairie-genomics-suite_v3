# Quick Start Script for Prairie Genomics Suite - Optimized
# Run this script to start the application with minimal setup

cat("🚀 Starting Prairie Genomics Suite - Optimized Version\n")
cat("=" , rep("=", 50), "\n")

# Check if we're in the right directory
if (!file.exists("app.R")) {
  stop("❌ Please run this script from the prairie-genomics-suite_optimized directory")
}

# Set development environment
Sys.setenv(R_ENV = "development")

# Install essential packages if missing
check_and_install <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if (length(missing_packages) > 0) {
    cat("📦 Installing missing essential packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, repos = "https://cloud.r-project.org/")
  }
}

# Check essential packages
essential_packages <- c("shiny", "shinydashboard", "DT", "ggplot2", "dplyr", "readr", "R6")
check_and_install(essential_packages)

# Check performance packages for better file reading
performance_packages <- c("data.table")  # For robust file reading
check_and_install(performance_packages)

# Try to load packages
for (pkg in essential_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("❌ Essential package", pkg, "is not available\n")
    cat("   Please run: install.packages('", pkg, "')\n", sep = "")
    stop("Essential packages missing")
  }
}

cat("✅ Essential packages verified\n")

# Set options for better performance
options(
  shiny.maxRequestSize = 500 * 1024^2,  # 500MB file upload limit
  repos = c(CRAN = "https://cloud.r-project.org/")
)

cat("🎯 Configuration:\n")
cat("   - File upload limit: 500MB\n")
cat("   - Development mode: ENABLED\n")
cat("   - Memory monitoring: ENABLED\n")

# Launch the application
cat("\n🚀 Launching Prairie Genomics Suite...\n")
cat("   Open your browser to the URL shown below\n")
cat("   Test files available:\n")
cat("   - test_data.csv: Basic dataset\n")
cat("   - test_data_with_patterns.csv: Dataset with clear Control/Treatment patterns\n\n")

# Run the app
shiny::runApp("app.R", launch.browser = TRUE)