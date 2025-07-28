# GENE SPECIES DIAGNOSTIC
# Determine if gene IDs are human or mouse and what type they are

cat("🔍 GENE SPECIES DIAGNOSTIC\n")
cat("=" , rep("=", 60), "\n")

# Function to analyze gene IDs and determine species
analyze_gene_ids <- function(gene_ids) {
  cat("📊 Analyzing", length(gene_ids), "gene IDs...\n")
  
  # Sample of gene IDs for display
  sample_genes <- head(gene_ids, 10)
  cat("📋 Sample gene IDs:\n", paste(sample_genes, collapse = ", "), "\n\n")
  
  # Check for different gene ID patterns
  patterns <- list(
    human_ensembl = "^ENSG\\d{11}",
    mouse_ensembl = "^ENSMUSG\\d{11}", 
    general_ensembl = "^ENS[A-Z]*G\\d{11}",
    human_symbol = "^[A-Z][A-Z0-9-]*$",
    mouse_symbol = "^[A-Z][a-z0-9-]*$",
    entrez = "^\\d+$",
    refseq = "^[NX][MR]_\\d+"
  )
  
  results <- list()
  
  for (pattern_name in names(patterns)) {
    pattern <- patterns[[pattern_name]]
    matches <- sum(grepl(pattern, gene_ids))
    percentage <- round(100 * matches / length(gene_ids), 1)
    results[[pattern_name]] <- list(matches = matches, percentage = percentage)
    
    cat(sprintf("%-15s: %d/%d genes (%.1f%%)\n", pattern_name, matches, length(gene_ids), percentage))
  }
  
  # Determine most likely species and ID type
  cat("\n🎯 ANALYSIS RESULTS:\n")
  
  if (results$human_ensembl$percentage > 80) {
    cat("✅ DETECTED: Human Ensembl gene IDs (ENSG...)\n")
    cat("✅ RECOMMENDED: Use species='human' with org.Hs.eg.db\n")
    return(list(species = "human", id_type = "ensembl", confidence = "high"))
  } else if (results$mouse_ensembl$percentage > 80) {
    cat("✅ DETECTED: Mouse Ensembl gene IDs (ENSMUSG...)\n") 
    cat("✅ RECOMMENDED: Use species='mouse' with org.Mm.eg.db\n")
    return(list(species = "mouse", id_type = "ensembl", confidence = "high"))
  } else if (results$general_ensembl$percentage > 50) {
    cat("⚠️ DETECTED: Mixed or other species Ensembl IDs\n")
    cat("⚠️ RECOMMENDED: Check specific species in your data\n")
    return(list(species = "unknown", id_type = "ensembl", confidence = "medium"))
  } else if (results$human_symbol$percentage > 70) {
    cat("✅ DETECTED: Human gene symbols (uppercase)\n")
    cat("ℹ️ RECOMMENDED: No conversion needed\n")
    return(list(species = "human", id_type = "symbol", confidence = "high"))
  } else if (results$mouse_symbol$percentage > 70) {
    cat("✅ DETECTED: Mouse gene symbols (mixed case)\n")
    cat("ℹ️ RECOMMENDED: No conversion needed\n") 
    return(list(species = "mouse", id_type = "symbol", confidence = "high"))
  } else {
    cat("❓ UNKNOWN: Mixed or unrecognized gene ID format\n")
    cat("⚠️ RECOMMENDED: Manual inspection needed\n")
    return(list(species = "unknown", id_type = "unknown", confidence = "low"))
  }
}

# Test with some common examples
cat("🧪 TESTING WITH EXAMPLE GENE IDs:\n\n")

# Test 1: Human Ensembl
cat("TEST 1: Human Ensembl IDs\n")
human_ensembl <- c("ENSG00000139618", "ENSG00000012048", "ENSG00000111640", "ENSG00000075624")
result1 <- analyze_gene_ids(human_ensembl)
cat("\n")

# Test 2: Mouse Ensembl  
cat("TEST 2: Mouse Ensembl IDs\n")
mouse_ensembl <- c("ENSMUSG00000051951", "ENSMUSG00000020122", "ENSMUSG00000057666", "ENSMUSG00000030636")
result2 <- analyze_gene_ids(mouse_ensembl)
cat("\n")

# Test 3: Human symbols
cat("TEST 3: Human Gene Symbols\n")
human_symbols <- c("BRCA1", "BRCA2", "TP53", "GAPDH", "ACTB")
result3 <- analyze_gene_ids(human_symbols)
cat("\n")

# Test 4: Mouse symbols
cat("TEST 4: Mouse Gene Symbols\n") 
mouse_symbols <- c("Brca1", "Brca2", "Trp53", "Gapdh", "Actb")
result4 <- analyze_gene_ids(mouse_symbols)
cat("\n")

cat("🔧 TO DIAGNOSE YOUR ACTUAL DATA:\n")
cat("1. Look at your gene IDs in the first column of your data file\n")
cat("2. Compare them to the patterns above\n")
cat("3. If they start with ENSMUSG, you need mouse conversion\n")
cat("4. If they start with ENSG, you need human conversion\n")
cat("5. Make sure to select the correct species in the UI dropdown\n\n")

cat("✅ Diagnostic completed - check your gene ID format above\n")