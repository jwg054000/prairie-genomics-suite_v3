# Package Installation Script - Prairie Genomics Suite Optimized
# Installs all required packages with smart dependency management
# 
# Author: Prairie Genomics Team - Optimized Version

cat("🧬 Installing Prairie Genomics Suite - Optimized Dependencies\n")
cat("=" , rep("=", 70), "\n")

# Load configuration
source("config/app_config.R")

# Set CRAN mirror
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  cat("🔧 CRAN mirror configured\n")
}

# Function to check if package is installed
is_package_installed <- function(package) {
  return(package %in% rownames(installed.packages()))
}

# Function to install package with error handling
install_package_safely <- function(package, source = "CRAN") {
  if (is_package_installed(package)) {
    cat("✅", package, "already installed\n")
    return(TRUE)
  }
  
  cat("📦 Installing", package, "from", source, "...\n")
  
  success <- tryCatch({
    if (source == "CRAN") {
      install.packages(package, dependencies = TRUE, quiet = TRUE)
    } else if (source == "Bioconductor") {
      if (!is_package_installed("BiocManager")) {
        install.packages("BiocManager", quiet = TRUE)
      }
      BiocManager::install(package, quiet = TRUE)
    }
    
    # Check if installation was successful
    if (is_package_installed(package)) {
      cat("✅", package, "installed successfully\n")
      return(TRUE)
    } else {
      cat("❌", package, "installation failed\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("❌", package, "installation error:", e$message, "\n")
    return(FALSE)
  })
  
  return(success)
}

# Install essential packages (required for basic functionality)
install_essential_packages <- function() {
  cat("\n📋 Installing Essential Packages\n")
  cat("-", rep("-", 40), "\n")
  
  essential_packages <- get_config("packages", "essential")
  
  failed_packages <- character(0)
  
  for (package in essential_packages) {
    if (!install_package_safely(package, "CRAN")) {
      failed_packages <- c(failed_packages, package)
    }
  }
  
  # Install R6 for memory manager
  if (!install_package_safely("R6", "CRAN")) {
    failed_packages <- c(failed_packages, "R6")
  }
  
  if (length(failed_packages) > 0) {
    cat("\n❌ Essential package installation failed for:", paste(failed_packages, collapse = ", "), "\n")
    cat("The application may not function properly without these packages.\n")
    return(FALSE)
  } else {
    cat("\n✅ All essential packages installed successfully!\n")
    return(TRUE)
  }
}

# Install Bioconductor packages
install_bioconductor_packages <- function() {
  cat("\n🧬 Installing Bioconductor Packages\n")
  cat("-", rep("-", 40), "\n")
  
  bioc_packages <- get_config("packages", "bioconductor")
  
  failed_packages <- character(0)
  
  for (package in bioc_packages) {
    if (!install_package_safely(package, "Bioconductor")) {
      failed_packages <- c(failed_packages, package)
    }
  }
  
  if (length(failed_packages) > 0) {
    cat("\n⚠️ Some Bioconductor packages failed to install:", paste(failed_packages, collapse = ", "), "\n")
    cat("The application will use fallback methods for these features.\n")
  } else {
    cat("\n✅ All Bioconductor packages installed successfully!\n")
  }
  
  return(failed_packages)
}

# Install optional packages
install_optional_packages <- function() {
  cat("\n🎨 Installing Optional Packages\n")
  cat("-", rep("-", 40), "\n")
  
  optional_packages <- get_config("packages", "optional")
  
  failed_packages <- character(0)
  
  for (package in optional_packages) {
    if (!install_package_safely(package, "CRAN")) {
      failed_packages <- c(failed_packages, package)
    }
  }
  
  if (length(failed_packages) > 0) {
    cat("\n⚠️ Some optional packages failed to install:", paste(failed_packages, collapse = ", "), "\n")
    cat("Some enhanced features may not be available.\n")
  } else {
    cat("\n✅ All optional packages installed successfully!\n")
  }
  
  return(failed_packages)
}

# Performance packages for optimized version
install_performance_packages <- function() {
  cat("\n⚡ Installing Performance Enhancement Packages\n")
  cat("-", rep("-", 40), "\n")
  
  performance_packages <- c(
    "future",      # Async processing
    "promises",    # Promise-based async operations
    "memoise",     # Function memoization
    "data.table",  # Fast data operations and robust file reading
    "dtplyr",      # dplyr backend for data.table
    "readr"        # Fast and robust file reading
  )
  
  failed_packages <- character(0)
  
  for (package in performance_packages) {
    if (!install_package_safely(package, "CRAN")) {
      failed_packages <- c(failed_packages, package)
    }
  }
  
  if (length(failed_packages) > 0) {
    cat("\n⚠️ Some performance packages failed to install:", paste(failed_packages, collapse = ", "), "\n")
    cat("Performance optimizations may be limited.\n")
  } else {
    cat("\n✅ All performance packages installed successfully!\n")
  }
  
  return(failed_packages)
}

# Check system requirements
check_system_requirements <- function() {
  cat("\n🔍 Checking System Requirements\n")
  cat("-", rep("-", 40), "\n")
  
  # Check R version
  r_version <- R.version.string
  cat("R Version:", r_version, "\n")
  
  if (getRversion() < "4.0.0") {
    cat("⚠️ Warning: R version 4.0.0 or higher is recommended\n")
  } else {
    cat("✅ R version is compatible\n")
  }
  
  # Check available memory
  mem_info <- tryCatch({
    gc_info <- gc()
    available_mb <- round(sum(gc_info[, 2]), 1)
    cat("Available Memory:", available_mb, "MB\n")
    
    if (available_mb < 500) {
      cat("⚠️ Warning: Low available memory. Consider closing other applications.\n")
    } else {
      cat("✅ Sufficient memory available\n")
    }
    
  }, error = function(e) {
    cat("❓ Unable to determine memory usage\n")
  })
  
  # Check platform
  platform <- R.version$platform
  cat("Platform:", platform, "\n")
  
  if (grepl("win", platform, ignore.case = TRUE)) {
    cat("💻 Windows platform detected\n")
  } else if (grepl("darwin", platform, ignore.case = TRUE)) {
    cat("🍎 macOS platform detected\n")
  } else if (grepl("linux", platform, ignore.case = TRUE)) {
    cat("🐧 Linux platform detected\n")
  }
  
  cat("✅ System check completed\n")
}

# Main installation function
main_installation <- function() {
  start_time <- Sys.time()
  
  cat("🚀 Starting installation process...\n\n")
  
  # Check system requirements
  check_system_requirements()
  
  # Install packages in order of importance
  essential_success <- install_essential_packages()
  
  if (!essential_success) {
    cat("\n❌ Installation failed! Essential packages could not be installed.\n")
    cat("Please check your internet connection and R installation.\n")
    return(FALSE)
  }
  
  # Continue with optional installations
  bioc_failed <- install_bioconductor_packages()
  optional_failed <- install_optional_packages()
  performance_failed <- install_performance_packages()
  
  # Installation summary
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 1)
  
  cat("\n" , rep("=", 70), "\n")
  cat("📋 INSTALLATION SUMMARY\n")
  cat(rep("=", 70), "\n")
  
  cat("✅ Essential packages: INSTALLED\n")
  
  if (length(bioc_failed) == 0) {
    cat("✅ Bioconductor packages: INSTALLED\n")
  } else {
    cat("⚠️ Bioconductor packages: PARTIAL (", length(bioc_failed), "failed )\n")
  }
  
  if (length(optional_failed) == 0) {
    cat("✅ Optional packages: INSTALLED\n")
  } else {
    cat("⚠️ Optional packages: PARTIAL (", length(optional_failed), "failed )\n")
  }
  
  if (length(performance_failed) == 0) {
    cat("✅ Performance packages: INSTALLED\n")
  } else {
    cat("⚠️ Performance packages: PARTIAL (", length(performance_failed), "failed )\n")
  }
  
  cat("⏱️ Installation time:", duration, "minutes\n")
  
  cat("\n🎉 Prairie Genomics Suite - Optimized is ready to use!\n")
  cat("Run the application with: shiny::runApp('app.R')\n\n")
  
  # Save installation log
  all_failed <- c(bioc_failed, optional_failed, performance_failed)
  if (length(all_failed) > 0) {
    cat("📝 Failed packages saved to 'installation_failures.txt'\n")
    writeLines(all_failed, "installation_failures.txt")
  }
  
  return(TRUE)
}

# Run installation
tryCatch({
  success <- main_installation()
  
  if (success) {
    cat("🏁 Installation completed successfully!\n")
  } else {
    cat("❌ Installation completed with errors.\n")
  }
  
}, error = function(e) {
  cat("💥 Installation script failed with error:", e$message, "\n")
  cat("Please check your R installation and internet connection.\n")
})

cat("\n📚 For troubleshooting, visit: https://github.com/your-repo/prairie-genomics-suite\n")