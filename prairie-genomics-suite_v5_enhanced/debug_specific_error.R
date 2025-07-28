# Debug specific undefined columns error
library(clusterProfiler)
library(org.Mm.eg.db)

# Simulate what's happening in the pathway analysis
cat("Creating test GO result...\n")
test_genes <- c("11545", "74778", "17869", "16590", "12043", "13982", "16193", "11593", "20926", "75560")

go_result <- enrichGO(
  gene = test_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 300
)

cat("GO enrichment completed\n")
go_df <- as.data.frame(go_result@result)
cat("Converted to data frame - ", nrow(go_df), "rows\n")
cat("Columns: ", paste(colnames(go_df), collapse=", "), "\n\n")

# Simulate the filtering process
cat("Testing filtering steps...\n")

padj_cutoff <- 0.05

# Step 1: Check individual column operations
cat("\n1. Testing column existence and operations:\n")
cat("   - Testing go_df$p.adjust <= padj_cutoff...\n")
tryCatch({
  test1 <- go_df$p.adjust <= padj_cutoff
  cat("     ✅ Success\n")
}, error = function(e) {
  cat("     ❌ Error:", e$message, "\n")
})

cat("   - Testing go_df$Count >= 3...\n")
tryCatch({
  test2 <- go_df$Count >= 3
  cat("     ✅ Success\n")
}, error = function(e) {
  cat("     ❌ Error:", e$message, "\n")
})

# Step 2: Test combined condition
cat("\n2. Testing combined condition:\n")
tryCatch({
  combined <- go_df$p.adjust <= padj_cutoff & go_df$Count >= 3
  cat("   ✅ Combined condition successful\n")
  cat("   - TRUE values:", sum(combined, na.rm=TRUE), "\n")
}, error = function(e) {
  cat("   ❌ Error in combined condition:", e$message, "\n")
})

# Step 3: Test subsetting
cat("\n3. Testing data frame subsetting:\n")
tryCatch({
  standard_filter <- go_df[
    go_df$p.adjust <= padj_cutoff & 
    go_df$Count >= 3,
    ]
  cat("   ✅ Subsetting successful - ", nrow(standard_filter), "rows\n")
}, error = function(e) {
  cat("   ❌ Error in subsetting:", e$message, "\n")
  cat("   Trying alternative syntax...\n")
  
  # Try alternative
  tryCatch({
    indices <- which(go_df$p.adjust <= padj_cutoff & go_df$Count >= 3)
    standard_filter <- go_df[indices, ]
    cat("   ✅ Alternative worked - ", nrow(standard_filter), "rows\n")
  }, error = function(e2) {
    cat("   ❌ Alternative also failed:", e2$message, "\n")
  })
})

# Step 4: Test metadata addition
cat("\n4. Testing metadata addition:\n")
if (nrow(go_df) > 0) {
  tryCatch({
    go_df$test_column <- "test"
    cat("   ✅ Can add new columns\n")
  }, error = function(e) {
    cat("   ❌ Cannot add columns:", e$message, "\n")
  })
}

cat("\nDebug complete!\n")