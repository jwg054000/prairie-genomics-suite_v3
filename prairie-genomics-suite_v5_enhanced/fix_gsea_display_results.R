# Fix GSEA Results Display in UI
# GSEA is running successfully but results aren't showing in the UI table
# Author: Prairie Genomics Suite Development Team
# Date: January 24, 2025

cat("🔧 Fixing GSEA Results Display Issues\n")
cat("=" , rep("=", 50), "\n")

# Issue Analysis: GSEA runs successfully but results don't show in UI
cat("📊 Issue Analysis:\n")
cat("✅ GSEA analysis completes successfully (48 pathways, 25 significant)\n")
cat("✅ Console shows correct statistics\n")
cat("❌ Results table empty in UI\n")
cat("❌ Visualizations not displaying\n")

cat("\n🔍 Possible Causes:\n")
cat("• fgsea result format incompatible with UI expectations\n")
cat("• Column name mismatches between fgsea output and display code\n")
cat("• fgseaSimple vs fgseaMultilevel format differences\n")
cat("• Data frame conversion issues\n")
cat("• UI table filtering problems\n")

# Issue 1: Update to fgseaMultilevel
cat("\n⚡ Issue 1: Update to fgseaMultilevel\n")
cat(rep("-", 40), "\n")

cat("Warning message indicates using deprecated fgseaSimple\n")
cat("Solution: Remove nperm argument to use recommended fgseaMultilevel\n")

improved_fgsea_function <- '
# IMPROVED: GSEA analysis using fgseaMultilevel (recommended)
run_gsea_analysis <- function(gene_list, species, gene_set_collection = "H") {
  tryCatch({
    cat("🔬 Running GSEA analysis with fgseaMultilevel\\n")
    
    # Ensure fgsea is available
    if (!require("fgsea", quietly = TRUE)) {
      cat("📦 Installing fgsea package...\\n")
      if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
      }
      BiocManager::install("fgsea")
      library(fgsea)
    }
    
    # Get gene sets in proper format for fgsea
    pathways <- get_fgsea_gene_sets(species, gene_set_collection)
    
    if (is.null(pathways) || length(pathways) == 0) {
      return(list(
        success = FALSE,
        error = "No gene sets retrieved",
        message = paste("Could not get", gene_set_collection, "gene sets for", species)
      ))
    }
    
    # Filter gene sets by size (15-500 genes as per reference)
    pathway_sizes <- sapply(pathways, length)
    valid_pathways <- pathways[pathway_sizes >= 15 & pathway_sizes <= 500]
    
    cat("📋 Gene set filtering:\\n")
    cat("   - Original pathways:", length(pathways), "\\n")
    cat("   - After size filtering (15-500):", length(valid_pathways), "\\n")
    
    if (length(valid_pathways) == 0) {
      return(list(
        success = FALSE,
        error = "No valid gene sets after filtering",
        message = "All gene sets were outside the 15-500 gene size range"
      ))
    }
    
    # Run fgseaMultilevel (RECOMMENDED - removed nperm argument)
    cat("🧬 Running fgseaMultilevel with", length(valid_pathways), "pathways...\\n")
    
    fgsea_results <- fgsea(
      pathways = valid_pathways,
      stats = gene_list,
      minSize = 15,
      maxSize = 500
      # Removed nperm argument to use fgseaMultilevel
    )
    
    if (is.null(fgsea_results) || nrow(fgsea_results) == 0) {
      return(list(
        success = FALSE,
        error = "No significant pathways found",
        message = "Try different gene set collection or check gene ranking"
      ))
    }
    
    # Convert to data.frame and ensure proper column names
    fgsea_results <- as.data.frame(fgsea_results)
    
    # Sort by adjusted p-value
    fgsea_results <- fgsea_results[order(fgsea_results$padj), ]
    
    # Add interpretation columns for UI display
    fgsea_results$Direction <- ifelse(fgsea_results$NES > 0, "Upregulated", "Downregulated") 
    fgsea_results$Significance <- ifelse(fgsea_results$padj < 0.05, "Significant", "Not Significant")
    fgsea_results$AbsNES <- abs(fgsea_results$NES)
    
    # Format p-values for display (avoid scientific notation)
    fgsea_results$pval_display <- ifelse(fgsea_results$pval < 0.001, 
                                        "< 0.001", 
                                        sprintf("%.3f", fgsea_results$pval))
    fgsea_results$padj_display <- ifelse(fgsea_results$padj < 0.001, 
                                        "< 0.001", 
                                        sprintf("%.3f", fgsea_results$padj))
    
    # Convert leadingEdge to character for display
    if ("leadingEdge" %in% colnames(fgsea_results)) {
      fgsea_results$leadingEdge_display <- sapply(fgsea_results$leadingEdge, function(x) {
        if (length(x) > 5) {
          paste(c(x[1:5], "..."), collapse = ", ")
        } else {
          paste(x, collapse = ", ")
        }
      })
    }
    
    cat("✅ GSEA completed successfully:\\n")
    cat("   - Total pathways tested:", nrow(fgsea_results), "\\n")
    cat("   - Significant pathways (padj < 0.05):", sum(fgsea_results$padj < 0.05), "\\n")
    cat("   - Upregulated pathways:", sum(fgsea_results$NES > 0 & fgsea_results$padj < 0.05), "\\n")
    cat("   - Downregulated pathways:", sum(fgsea_results$NES < 0 & fgsea_results$padj < 0.05), "\\n")
    
    # Debug: Print column names and first few rows
    cat("\\n🔍 Debug - Result structure:\\n")
    cat("   - Columns:", paste(colnames(fgsea_results), collapse = ", "), "\\n")
    cat("   - First pathway:", fgsea_results$pathway[1], "\\n")
    cat("   - First NES:", round(fgsea_results$NES[1], 3), "\\n")
    cat("   - First padj:", sprintf("%.3e", fgsea_results$padj[1]), "\\n")
    
    return(list(
      success = TRUE,
      data = fgsea_results,
      analysis_type = "GSEA",
      method = "fgseaMultilevel", 
      species = species,
      gene_set_collection = gene_set_collection,
      n_pathways = nrow(fgsea_results),
      n_significant = sum(fgsea_results$padj < 0.05),
      columns = colnames(fgsea_results)  # For debugging
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message,
      message = paste("GSEA analysis failed:", e$message)
    ))
  })
}
'

cat("✅ Created improved GSEA function with fgseaMultilevel\n")

# Issue 2: Create UI-compatible result formatting
cat("\n📋 Issue 2: UI-Compatible Result Formatting\n")
cat(rep("-", 45), "\n")

ui_formatting_function <- '
# Format GSEA results for UI display
format_gsea_results_for_ui <- function(gsea_results) {
  if (is.null(gsea_results) || !gsea_results$success || nrow(gsea_results$data) == 0) {
    return(data.frame(
      Pathway = character(0),
      Direction = character(0),
      NES = numeric(0),
      pvalue = character(0),
      padj = character(0),
      size = numeric(0),
      leadingEdge = character(0)
    ))
  }
  
  df <- gsea_results$data
  
  # Create clean display dataframe
  display_df <- data.frame(
    Pathway = df$pathway,
    Direction = ifelse(df$NES > 0, "↑ Upregulated", "↓ Downregulated"),
    NES = round(df$NES, 3),
    pvalue = ifelse(df$pval < 0.001, "< 0.001", sprintf("%.3f", df$pval)),
    padj = ifelse(df$padj < 0.001, "< 0.001", sprintf("%.3f", df$padj)),
    size = df$size,
    leadingEdge = sapply(df$leadingEdge, function(x) {
      if (length(x) > 3) {
        paste(c(x[1:3], paste("...", length(x) - 3, "more")), collapse = ", ")
      } else {
        paste(x, collapse = ", ")
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Sort by significance
  display_df <- display_df[order(display_df$padj), ]
  
  return(display_df)
}

# Test result formatting with mock data
test_gsea_formatting <- function() {
  cat("🧪 Testing GSEA result formatting...\\n")
  
  # Create mock fgsea results
  mock_results <- data.frame(
    pathway = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_GLYCOLYSIS"),
    pval = c(0.001, 0.0001, 0.02),
    padj = c(0.01, 0.002, 0.05),
    NES = c(2.1, -1.8, 1.4),
    size = c(200, 150, 180),
    leadingEdge = list(c("COX1", "COX2", "ATP5A"), c("MYC", "MAX", "MXD1"), c("HK1", "HK2", "PFKM")),
    stringsAsFactors = FALSE
  )
  
  # Format for UI
  formatted <- format_gsea_results_for_ui(list(success = TRUE, data = mock_results))
  
  cat("✅ Formatted", nrow(formatted), "results for UI display\\n")
  print(formatted)
  
  return(formatted)
}
'

cat("✅ Created UI formatting functions\n")

# Issue 3: Debug current result structure
cat("\n🔍 Issue 3: Debug Current Results\n")
cat(rep("-", 35), "\n")

debug_function <- '
# Debug GSEA results to identify display issues
debug_gsea_results <- function(gsea_result) {
  cat("🔍 Debugging GSEA Results Structure\\n")
  cat(rep("=", 40), "\\n")
  
  if (is.null(gsea_result)) {
    cat("❌ gsea_result is NULL\\n")
    return()
  }
  
  if (!is.list(gsea_result)) {
    cat("❌ gsea_result is not a list\\n")
    return()
  }
  
  cat("📋 gsea_result structure:\\n")
  cat("   - success:", gsea_result$success, "\\n")
  cat("   - error:", gsea_result$error %||% "None", "\\n")
  cat("   - analysis_type:", gsea_result$analysis_type %||% "Unknown", "\\n")
  
  if ("data" %in% names(gsea_result)) {
    data <- gsea_result$data
    cat("\\n📊 Data structure:\\n")
    cat("   - Class:", class(data), "\\n")
    cat("   - Dimensions:", paste(dim(data), collapse = " x "), "\\n")
    cat("   - Columns:", paste(colnames(data), collapse = ", "), "\\n")
    
    if (nrow(data) > 0) {
      cat("\\n📝 First few rows:\\n")
      print(head(data, 3))
      
      # Check for required columns
      required_cols <- c("pathway", "pval", "padj", "NES", "size")
      missing_cols <- required_cols[!required_cols %in% colnames(data)]
      
      if (length(missing_cols) > 0) {
        cat("\\n❌ Missing required columns:", paste(missing_cols, collapse = ", "), "\\n")
      } else {
        cat("\\n✅ All required columns present\\n")
      }
      
      # Check for significant results
      if ("padj" %in% colnames(data)) {
        sig_count <- sum(data$padj < 0.05, na.rm = TRUE)
        cat("   - Significant results (padj < 0.05):", sig_count, "\\n")
      }
    } else {
      cat("\\n⚠️ Data frame is empty\\n")
    }
  } else {
    cat("\\n❌ No data element found in gsea_result\\n")
  }
}
'

cat("✅ Created debugging function\n")

# Write all functions to file
complete_fix <- paste(
  "# Complete GSEA Display Fix",
  "# Generated:", Sys.time(),
  "",
  "# Fix 1: Updated GSEA function with fgseaMultilevel",
  improved_fgsea_function,
  "",
  "# Fix 2: UI formatting functions", 
  ui_formatting_function,
  "",
  "# Fix 3: Debug functions",
  debug_function,
  "",
  "# Test the formatting",
  "cat('🧪 Testing GSEA result formatting...\\n')",
  "test_result <- test_gsea_formatting()",
  "",
  "cat('✅ GSEA display fixes loaded successfully\\n')",
  sep = "\n"
)

writeLines(complete_fix, "gsea_display_fixes.R")
cat("✅ Complete fix saved to gsea_display_fixes.R\n")

# Create a quick test script
cat("\n🧪 Creating Quick Test Script\n")
cat(rep("-", 35), "\n")

quick_test <- '
# Quick test of GSEA display fix
cat("🧪 Quick GSEA Display Test\\n")

# Load the pathway analysis functions
if (file.exists("pathway_analysis.R")) {
  source("pathway_analysis.R")
  cat("✅ pathway_analysis.R loaded\\n")
} else {
  cat("❌ pathway_analysis.R not found\\n")
  quit()
}

# Test with mock data
mock_gene_list <- c(2.5, 1.8, -1.2, 3.1, -2.0, 1.5, -1.8, 2.2)
names(mock_gene_list) <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "MTOR")
mock_gene_list <- sort(mock_gene_list, decreasing = TRUE)

cat("🧬 Testing GSEA with mock gene list...\\n")
cat("   - Genes:", length(mock_gene_list), "\\n")
cat("   - Range:", round(min(mock_gene_list), 2), "to", round(max(mock_gene_list), 2), "\\n")

# Test GSEA function
if (exists("run_gsea_analysis")) {
  result <- run_gsea_analysis(mock_gene_list, "human", "H")
  
  if (!is.null(result) && result$success) {
    cat("✅ GSEA function works!\\n")
    cat("   - Pathways found:", result$n_pathways, "\\n")
    cat("   - Significant:", result$n_significant, "\\n")
    cat("   - Columns:", paste(result$columns, collapse = ", "), "\\n")
  } else {
    cat("❌ GSEA function failed\\n") 
    if (!is.null(result$error)) {
      cat("   Error:", result$error, "\\n")
    }
  }
} else {
  cat("❌ run_gsea_analysis function not found\\n")
}
'

writeLines(quick_test, "test_gsea_display_quick.R")
cat("✅ Quick test script saved\n")

# Summary and recommendations
cat("\n🎯 Fix Summary\n")
cat("=" , rep("=", 20), "\n")

cat("Issues Identified:\n")
cat("• Using deprecated fgseaSimple instead of fgseaMultilevel\n")
cat("• Results format may not match UI expectations\n")
cat("• Possible column name mismatches\n")
cat("• leadingEdge list format issues\n")

cat("\nSolutions Implemented:\n")
cat("✅ Updated to fgseaMultilevel (removed nperm argument)\n")
cat("✅ Enhanced result formatting for UI display\n")
cat("✅ Added debug functions to identify issues\n")
cat("✅ Improved column formatting and p-value display\n")
cat("✅ leadingEdge list converted to display strings\n")

cat("\n🚀 Next Steps:\n")
cat("1. Update pathway_analysis.R with new GSEA function\n")
cat("2. Run test_gsea_display_quick.R to verify fixes\n")
cat("3. Check UI table rendering with formatted results\n")
cat("4. Test with real data in the app\n")

cat("\n💡 Key Changes:\n")
cat("• fgsea() call: Removed nperm argument → uses fgseaMultilevel\n")
cat("• Added Direction, Significance columns for UI\n")
cat("• P-value formatting: Scientific notation → readable format\n")
cat("• leadingEdge: List → comma-separated string\n")
cat("• Enhanced debug output for troubleshooting\n")

cat("\n🧬 GSEA Display Issues Should Be Fixed!\n")
cat("📊 Results should now appear properly in the UI table\n")