# Debug script to find the exact location of "undefined columns selected" error

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)

# Create test enrichGO result
test_genes <- c("11545", "74778", "17869", "16590", "12043", "13982", "16193", "11593", "20926", "75560")

cat("Creating test GO enrichment result...\n")
go_result <- enrichGO(
  gene = test_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

if (!is.null(go_result) && nrow(go_result@result) > 0) {
  go_df <- as.data.frame(go_result@result)
  cat("✅ GO enrichment successful - ", nrow(go_df), "pathways found\n")
  cat("Columns available: ", paste(colnames(go_df), collapse=", "), "\n\n")
  
  # Test different filtering approaches
  cat("Testing filtering approaches...\n")
  
  # Test 1: Basic filtering
  cat("\n1. Testing basic filter: go_df$p.adjust <= 0.05\n")
  tryCatch({
    test1 <- go_df[go_df$p.adjust <= 0.05, ]
    cat("   ✅ Success - ", nrow(test1), " rows\n")
  }, error = function(e) {
    cat("   ❌ Error: ", e$message, "\n")
  })
  
  # Test 2: Combined filter with &
  cat("\n2. Testing combined filter: go_df$p.adjust <= 0.05 & go_df$Count >= 3\n")
  tryCatch({
    test2 <- go_df[go_df$p.adjust <= 0.05 & go_df$Count >= 3, ]
    cat("   ✅ Success - ", nrow(test2), " rows\n")
  }, error = function(e) {
    cat("   ❌ Error: ", e$message, "\n")
  })
  
  # Test 3: With trailing comma (old style)
  cat("\n3. Testing with trailing comma: go_df[condition, ]\n")
  tryCatch({
    test3 <- go_df[
      go_df$p.adjust <= 0.05 & 
      go_df$Count >= 3
    , ]
    cat("   ✅ Success - ", nrow(test3), " rows\n")
  }, error = function(e) {
    cat("   ❌ Error: ", e$message, "\n")
  })
  
  # Test 4: With comma before closing bracket
  cat("\n4. Testing with comma before bracket: go_df[condition, ]\n")
  tryCatch({
    test4 <- go_df[
      go_df$p.adjust <= 0.05 & 
      go_df$Count >= 3,
      ]
    cat("   ✅ Success - ", nrow(test4), " rows\n")
  }, error = function(e) {
    cat("   ❌ Error: ", e$message, "\n")
  })
  
  # Test 5: Check if columns exist
  cat("\n5. Checking column existence:\n")
  cat("   - 'p.adjust' exists: ", "p.adjust" %in% colnames(go_df), "\n")
  cat("   - 'Count' exists: ", "Count" %in% colnames(go_df), "\n")
  cat("   - 'padj' exists: ", "padj" %in% colnames(go_df), "\n")
  
} else {
  cat("❌ No GO results obtained\n")
}

cat("\nDebug complete!\n")