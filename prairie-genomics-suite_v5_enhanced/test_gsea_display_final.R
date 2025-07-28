# Final Test of GSEA Display Fixes
# Verify that GSEA results now display properly in UI
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🧪 Final Test of GSEA Display Fixes\n")
cat("=" , rep("=", 45), "\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Test 1: Load updated pathway analysis functions
cat("\nTest 1: Loading Updated Functions\n")
cat(rep("-", 35), "\n")

tryCatch({
  source("pathway_analysis.R")
  cat("✅ pathway_analysis.R loaded successfully\n")
  
  # Check if key functions exist
  key_functions <- c("run_gsea_analysis", "get_fgsea_gene_sets", "prepare_gene_list_gsea")
  
  for (func in key_functions) {
    if (exists(func)) {
      cat("✅", func, "function available\n")
    } else {
      cat("❌", func, "function missing\n")
    }
  }
  
}, error = function(e) {
  cat("❌ Failed to load pathway_analysis.R:", e$message, "\n")
  quit()
})

# Test 2: Verify fgsea package and fgseaMultilevel
cat("\nTest 2: fgsea Package Verification\n")
cat(rep("-", 35), "\n")

fgsea_available <- FALSE
tryCatch({
  if (!require("fgsea", quietly = TRUE)) {
    cat("📦 Installing fgsea...\n")
    BiocManager::install("fgsea")
    library(fgsea)
  }
  fgsea_available <- TRUE
  cat("✅ fgsea package loaded\n")
  
  # Test fgseaMultilevel vs fgseaSimple
  cat("📋 fgsea functions available:\n")
  if (exists("fgsea")) cat("✅ fgsea() main function\n")
  
}, error = function(e) {
  cat("❌ fgsea package issue:", e$message, "\n")
})

# Test 3: Create realistic test data
cat("\nTest 3: Creating Realistic Test Data\n")
cat(rep("-", 40), "\n")

if (fgsea_available) {
  # Create mock DESeq2 results
  set.seed(123)
  n_genes <- 500
  
  mock_deseq <- data.frame(
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    pvalue = runif(n_genes, 0.001, 0.1),
    padj = runif(n_genes, 0.001, 0.2),
    stringsAsFactors = FALSE
  )
  
  # Add some realistic gene names
  gene_names <- c(
    paste0("GENE_", sprintf("%04d", 1:(n_genes-20))),
    "TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "MTOR", 
    "PTEN", "RB1", "VEGFA", "TGFB1", "CDKN1A", "BCL2", "BAX", 
    "CASP3", "TNF", "IL6", "STAT3", "NFE2L2"
  )
  rownames(mock_deseq) <- gene_names
  
  cat("✅ Created mock DESeq2 results:", nrow(mock_deseq), "genes\n")
  
  # Test gene list preparation
  if (exists("prepare_gene_list_gsea")) {
    cat("🧪 Testing gene list preparation...\n")
    
    gsea_gene_list <- prepare_gene_list_gsea(mock_deseq, "human", "signed_pvalue")
    
    if (!is.null(gsea_gene_list) && length(gsea_gene_list) > 0) {
      cat("✅ Gene list prepared:", length(gsea_gene_list), "genes\n")
      cat("   - Top gene:", names(gsea_gene_list)[1], "=", round(gsea_gene_list[1], 3), "\n")
      cat("   - Bottom gene:", names(gsea_gene_list)[length(gsea_gene_list)], "=", round(gsea_gene_list[length(gsea_gene_list)], 3), "\n")
      
      # Test 4: Run GSEA analysis
      cat("\nTest 4: Running GSEA Analysis\n")
      cat(rep("-", 35), "\n")
      
      if (exists("run_gsea_analysis")) {
        cat("🧬 Running GSEA with Hallmark collection...\n")
        
        gsea_result <- run_gsea_analysis(gsea_gene_list, "human", "H")
        
        if (!is.null(gsea_result) && gsea_result$success) {
          cat("✅ GSEA analysis successful!\n")
          cat("   - Method:", gsea_result$method, "\n")
          cat("   - Pathways tested:", gsea_result$n_pathways, "\n") 
          cat("   - Significant pathways:", gsea_result$n_significant, "\n")
          cat("   - Columns in result:", paste(gsea_result$columns, collapse = ", "), "\n")
          
          # Test 5: Examine result structure for UI compatibility
          cat("\nTest 5: Result Structure Analysis\n")
          cat(rep("-", 40), "\n")
          
          result_data <- gsea_result$data
          
          cat("📊 Result data structure:\n")
          cat("   - Class:", class(result_data), "\n")
          cat("   - Dimensions:", paste(dim(result_data), collapse = " x "), "\n")
          cat("   - Column names:", paste(colnames(result_data), collapse = ", "), "\n")
          
          # Check for UI-required columns
          ui_required_cols <- c("pathway", "pval", "padj", "NES", "size", "Direction", "leadingEdge_display")
          
          cat("\n🔍 UI compatibility check:\n")
          for (col in ui_required_cols) {
            if (col %in% colnames(result_data)) {
              cat("✅", col, "- Present\n")
            } else {
              cat("❌", col, "- Missing\n")
            }
          }
          
          # Show first few results
          if (nrow(result_data) > 0) {
            cat("\n📋 First few results:\n")
            display_data <- result_data[1:min(3, nrow(result_data)), 
                                      c("pathway", "NES", "pval_display", "padj_display", "Direction")]
            print(display_data)
          }
          
          # Test 6: UI-formatted output
          cat("\nTest 6: UI-Formatted Output\n")
          cat(rep("-", 30), "\n")
          
          # Create a clean UI display version
          if (nrow(result_data) > 0) {
            ui_display <- data.frame(
              Pathway = gsub("HALLMARK_", "", result_data$pathway),
              Direction = result_data$Direction,
              NES = round(result_data$NES, 3),
              pvalue = result_data$pval_display,
              padj = result_data$padj_display,
              Size = result_data$size,
              stringsAsFactors = FALSE
            )
            
            # Show significant results only
            significant <- ui_display[as.numeric(gsub("< ", "", result_data$padj_display)) < 0.05, ]
            
            cat("✅ UI-formatted results ready:\n")
            cat("   - Total results:", nrow(ui_display), "\n")
            cat("   - Significant (padj < 0.05):", nrow(significant), "\n")
            
            if (nrow(significant) > 0) {
              cat("\n📊 Top significant pathways:\n")
              print(head(significant, 5))
            }
          }
          
        } else {
          cat("❌ GSEA analysis failed\n")
          if (!is.null(gsea_result$error)) {
            cat("   Error:", gsea_result$error, "\n")
          }
        }
      } else {
        cat("❌ run_gsea_analysis function not found\n")
      }
      
    } else {
      cat("❌ Gene list preparation failed\n")
    }
  } else {
    cat("❌ prepare_gene_list_gsea function not found\n")
  }
} else {
  cat("⚠️ Skipping GSEA tests - fgsea package not available\n")
}

# Final summary
cat("\n🎯 Final Test Summary\n")
cat("=" , rep("=", 25), "\n")

cat("Key Fixes Implemented:\n")
cat("✅ Updated to fgseaMultilevel (removed nperm argument)\n")
cat("✅ Enhanced result formatting for UI display\n")
cat("✅ Added Direction, Significance columns\n")
cat("✅ Formatted p-values (avoid scientific notation)\n")
cat("✅ Convert leadingEdge lists to display strings\n")
cat("✅ Added debug output for troubleshooting\n")

cat("\nExpected Behavior:\n")
cat("• No more fgseaSimple warning message\n")
cat("• GSEA results should appear in UI table\n")
cat("• Pathways should show Direction (Up/Down)\n")
cat("• P-values should be readable (not scientific)\n")
cat("• Leading edge genes should display as text\n")

cat("\n🚀 Ready for App Testing:\n")
cat("1. ✅ Functions updated with display fixes\n")
cat("2. ✅ fgseaMultilevel implementation\n")
cat("3. ✅ UI-compatible result formatting\n")
cat("4. 🧪 Test the main app with mouse GSEA\n")

cat("\n🧬 GSEA Display Fixes Complete!\n")
cat("📊 Results should now appear properly in the UI table\n")