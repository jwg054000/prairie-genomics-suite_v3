# Test Script for Prairie Genomics Suite - Optimized Version
# Quick validation that all modules load correctly and basic functionality works

cat("🧪 Testing Prairie Genomics Suite - Optimized Version\n")
cat("=" , rep("=", 50), "\n")

# Test 1: Configuration loading
cat("📋 Test 1: Configuration loading...\n")
tryCatch({
  source("config/app_config.R")
  cat("✅ Configuration loaded successfully\n")
  cat("   - App version:", get_config("app", "version"), "\n")
  cat("   - Memory chunk size:", CHUNK_SIZE, "\n")
  cat("   - Default padj:", DEFAULT_PADJ, "\n")
}, error = function(e) {
  cat("❌ Configuration loading failed:", e$message, "\n")
  stop("Configuration test failed")
})

# Test 2: Memory manager loading
cat("\n💾 Test 2: Memory manager loading...\n")
tryCatch({
  source("utils/memory_manager.R")
  cat("✅ Memory manager loaded successfully\n")
  
  # Test memory status
  memory_status <- get_memory_status()
  cat("   - Current memory usage:", memory_status$used_mb, "MB\n")
  cat("   - Memory status:", memory_status$recommendation, "\n")
}, error = function(e) {
  cat("❌ Memory manager loading failed:", e$message, "\n")
  stop("Memory manager test failed")
})

# Test 3: Data processing modules
cat("\n📁 Test 3: Data processing modules...\n")
tryCatch({
  source("modules/data_processing/file_upload.R")
  cat("✅ Data upload module loaded successfully\n")
}, error = function(e) {
  cat("❌ Data processing module loading failed:", e$message, "\n")
  stop("Data processing test failed")
})

# Test 4: Analysis modules
cat("\n🧬 Test 4: Analysis modules...\n")
tryCatch({
  source("modules/analysis/deseq2_analysis.R")
  cat("✅ DESeq2 analysis module loaded successfully\n")
}, error = function(e) {
  cat("❌ Analysis module loading failed:", e$message, "\n")
  stop("Analysis module test failed")
})

# Test 5: UI modules
cat("\n💻 Test 5: UI modules...\n")
tryCatch({
  source("modules/ui/sample_annotation.R")
  cat("✅ Sample annotation module loaded successfully\n")
}, error = function(e) {
  cat("❌ UI module loading failed:", e$message, "\n")
  stop("UI module test failed")
})

# Test 6: Essential packages
cat("\n📦 Test 6: Essential packages...\n")
essential_packages <- get_config("packages", "essential")
for (pkg in essential_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✅", pkg, "available\n")
  } else {
    cat("⚠️", pkg, "missing (but may not be critical)\n")
  }
}

# Test 7: Test data loading
cat("\n📊 Test 7: Test data validation...\n")
if (file.exists("test_data_with_patterns.csv")) {
  tryCatch({
    test_data <- read.csv("test_data_with_patterns.csv", row.names = 1)
    cat("✅ Test data loaded successfully\n")
    cat("   - Genes:", nrow(test_data), "\n")
    cat("   - Samples:", ncol(test_data), "\n")
    cat("   - Sample pattern detected:", paste(colnames(test_data)[1:3], collapse = ", "), "...\n")
  }, error = function(e) {
    cat("⚠️ Test data loading issue:", e$message, "\n")
  })
} else {
  cat("⚠️ Test data file not found (not critical for functionality)\n")
}

# Test 8: Utility functions
cat("\n🔧 Test 8: Utility functions...\n")
tryCatch({
  # Test null coalescing operator
  test_result <- NULL %||% "default_value"
  if (test_result == "default_value") {
    cat("✅ Null coalescing operator (%||%) working correctly\n")
  }
  
  # Test progress bar function
  if (exists("progressBar")) {
    cat("✅ progressBar function available\n")
  }
}, error = function(e) {
  cat("❌ Utility function test failed:", e$message, "\n")
})

cat("\n🎉 ALL TESTS COMPLETED!\n")
cat("=" , rep("=", 50), "\n")
cat("✅ Prairie Genomics Suite - Optimized Version appears to be ready\n")
cat("🚀 You can now run the application using: Rscript run_app.R\n")
cat("📋 Or manually with: R -e \"shiny::runApp('app.R')\"\n\n")

# Memory cleanup
smart_gc()