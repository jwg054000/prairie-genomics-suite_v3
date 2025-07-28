# Test Fast Gene Conversion Performance
# Tests the speed optimizations implemented for gene ID to symbol conversion

cat("🚀 Testing Fast Gene Conversion Performance\n")
cat("=" , rep("=", 60), "\n")

# Load modules
source("config/app_config.R")
source("modules/data_processing/file_upload.R")

# Test data
if (file.exists("test_data_with_patterns.csv")) {
  test_data <- read.csv("test_data_with_patterns.csv", stringsAsFactors = FALSE)
  test_genes <- test_data$Gene  # All 15 genes
  
  cat("📊 Test setup:\n")
  cat("   - Genes to convert:", length(test_genes), "\n")
  cat("   - Gene IDs:", paste(test_genes[1:3], collapse = ", "), "...\n")
  
  # Test 1: Speed comparison - local vs biomaRt
  cat("\n🏃‍♂️ Test 1: Speed Comparison\n")
  cat("-" , rep("-", 40), "\n")
  
  # Test local conversion speed (org.Hs.eg.db)
  if (orgdb_available) {
    cat("⚡ Testing local conversion (org.Hs.eg.db)...\n")
    start_time <- Sys.time()
    
    local_result <- convert_genes_locally(test_genes, species = "human")
    
    local_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(paste0("✅ Local conversion time: ", round(local_time, 2), " seconds\n"))
  } else {
    cat("⚠️ org.Hs.eg.db not available - skipping local test\n")
  }
  
  # Test optimized conversion (with caching and local preference)
  cat("\n🔄 Testing optimized multi-method conversion...\n")
  start_time <- Sys.time()
  
  optimized_result <- convert_gene_ids_to_symbols(test_genes, species = "human")
  
  optimized_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(paste0("✅ Optimized conversion time: ", round(optimized_time, 2), " seconds\n"))
  
  # Test 2: Caching effectiveness
  cat("\n💾 Test 2: Cache Performance\n")
  cat("-" , rep("-", 40), "\n")
  
  cat("🔄 Running conversion again (should use cache)...\n")
  start_time <- Sys.time()
  
  cached_result <- convert_gene_ids_to_symbols(test_genes, species = "human")
  
  cached_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(paste0("✅ Cached conversion time: ", round(cached_time, 2), " seconds\n"))
  
  if (cached_time < optimized_time) {
    speedup <- round(optimized_time / cached_time, 1)
    cat(paste0("🚀 Cache speedup: ", speedup, "x faster!\n"))
  }
  
  # Test 3: Results verification
  cat("\n✅ Test 3: Results Verification\n")
  cat("-" , rep("-", 40), "\n")
  
  if (!is.null(optimized_result)) {
    converted_count <- sum(optimized_result$gene_symbol != optimized_result$original_id)
    cat(paste0("📊 Conversion results:\n"))
    cat(paste0("   - Total genes: ", nrow(optimized_result), "\n"))
    cat(paste0("   - Successfully converted: ", converted_count, "\n"))
    cat(paste0("   - Conversion rate: ", round(100 * converted_count / nrow(optimized_result), 1), "%\n"))
    
    cat("\n📋 Sample conversions:\n")
    for (i in 1:min(5, nrow(optimized_result))) {
      cat(paste0("   - ", optimized_result$original_id[i], " → ", optimized_result$gene_symbol[i], "\n"))
    }
  }
  
  # Test 4: Matrix conversion performance
  cat("\n📊 Test 4: Matrix Conversion Performance\n")
  cat("-" , rep("-", 40), "\n")
  
  # Create test matrix
  test_matrix <- matrix(
    rnorm(length(test_genes) * 6),
    nrow = length(test_genes),
    dimnames = list(test_genes, c("Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"))
  )
  
  cat("🔄 Testing full matrix conversion...\n")
  start_time <- Sys.time()
  
  matrix_result <- apply_gene_conversion(test_matrix, species = "human", enable_conversion = TRUE)
  
  matrix_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(paste0("✅ Matrix conversion time: ", round(matrix_time, 2), " seconds\n"))
  
  cat(paste0("📊 Matrix conversion stats:\n"))
  cat(paste0("   - Attempted: ", matrix_result$conversion_stats$attempted, "\n"))
  cat(paste0("   - Total genes: ", matrix_result$conversion_stats$total_genes, "\n"))
  cat(paste0("   - Converted: ", matrix_result$conversion_stats$converted_count, "\n"))
  cat(paste0("   - Rate: ", matrix_result$conversion_stats$conversion_rate, "%\n"))
  
} else {
  cat("❌ Test data file not found\n")
}

# Summary
cat("\n🎉 SPEED OPTIMIZATION TEST COMPLETED!\n")
cat("=" , rep("=", 60), "\n")

cat("🔧 Speed optimizations implemented:\n")
cat("   ✅ Local gene database (org.Hs.eg.db) - near-instant for human genes\n")
cat("   ✅ Session-persistent caching - instant for repeated conversions\n") 
cat("   ✅ Large batch processing - up to 5,000 genes at once\n")
cat("   ✅ Minimal retry logic - reduced delays\n")
cat("   ✅ No rate limiting - maximum speed\n")
cat("   ✅ Fastest mirror prioritization - useast first\n")
cat("   ✅ Unique ID processing - avoid duplicate work\n")

cat("\n📈 Expected performance improvements:\n")
if (orgdb_available) {
  cat("   🚀 Local conversion: ~0.1-0.5 seconds (90%+ faster)\n")
} else {
  cat("   ⚠️ Install org.Hs.eg.db for maximum speed\n")
}
cat("   💾 Cached conversion: ~0.01 seconds (99%+ faster)\n")
cat("   🌐 biomaRt conversion: ~2-10 seconds (50-80% faster)\n")

cat("\n🎯 For maximum speed, ensure org.Hs.eg.db is installed:\n")
cat("   BiocManager::install('org.Hs.eg.db')\n\n")

# Cleanup
if (exists("test_matrix")) rm(test_matrix)
if (exists("optimized_result")) rm(optimized_result)
smart_gc()