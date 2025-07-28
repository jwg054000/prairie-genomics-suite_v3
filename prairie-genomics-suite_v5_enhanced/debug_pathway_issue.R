# Pathway Analysis Issue Diagnostic Tool
# Help identify and resolve current pathway analysis problems
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🔍 PATHWAY ANALYSIS DIAGNOSTIC TOOL\n")
cat("===================================\n\n")

# Function to capture detailed error information
diagnose_pathway_issue <- function(deseq2_results = NULL, analysis_type = "GO", 
                                  species = "human", debug_level = "verbose") {
  
  cat("🩺 COMPREHENSIVE PATHWAY ANALYSIS DIAGNOSIS\n")
  cat("===========================================\n\n")
  
  # Check system setup
  cat("1. SYSTEM ENVIRONMENT CHECK\n")
  cat("---------------------------\n")
  cat("R Version:", R.version.string, "\n")
  cat("Working Directory:", getwd(), "\n")
  cat("Available Memory:", round(memory.limit()/1024, 2), "GB\n")
  
  # Check package availability
  cat("\n2. PACKAGE AVAILABILITY\n")
  cat("-----------------------\n")
  
  packages_to_check <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db", 
                        "org.Mm.eg.db", "msigdbr", "fgsea", "pathview")
  
  for (pkg in packages_to_check) {
    status <- if (requireNamespace(pkg, quietly = TRUE)) "✅ Available" else "❌ Missing"
    cat(sprintf("%-20s: %s\n", pkg, status))
  }
  
  # Check pathway analysis functions
  cat("\n3. FUNCTION AVAILABILITY CHECK\n")
  cat("------------------------------\n")
  
  functions_to_check <- c("run_pathway_analysis", "prepare_gene_list_ora", 
                         "run_go_analysis", "run_kegg_analysis", "run_gsea_analysis")
  
  for (func in functions_to_check) {
    status <- if (exists(func)) "✅ Loaded" else "❌ Missing"
    cat(sprintf("%-25s: %s\n", func, status))
  }
  
  # If DESeq2 results provided, test the workflow
  if (!is.null(deseq2_results)) {
    cat("\n4. DATA VALIDATION\n")
    cat("------------------\n")
    
    # Check data structure
    cat("Data class:", class(deseq2_results), "\n")
    cat("Dimensions:", nrow(deseq2_results), "x", ncol(deseq2_results), "\n")
    cat("Column names:", paste(colnames(deseq2_results), collapse = ", "), "\n")
    
    # Check for required columns
    required_cols <- c("log2FoldChange", "padj", "pvalue")
    missing_cols <- required_cols[!required_cols %in% colnames(deseq2_results)]
    
    if (length(missing_cols) > 0) {
      cat("❌ Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    } else {
      cat("✅ All required columns present\n")
    }
    
    # Check data quality
    cat("\nData Quality:\n")
    cat("- Genes with padj values:", sum(!is.na(deseq2_results$padj)), "\n")
    cat("- Significant genes (padj < 0.05):", sum(deseq2_results$padj < 0.05, na.rm = TRUE), "\n")
    cat("- High FC genes (|FC| > 1):", sum(abs(deseq2_results$log2FoldChange) > 1, na.rm = TRUE), "\n")
    
    # Test gene list preparation
    cat("\n5. GENE LIST PREPARATION TEST\n")
    cat("-----------------------------\n")
    
    tryCatch({
      cat("Testing prepare_gene_list_ora...\n")
      gene_list <- prepare_gene_list_ora(deseq2_results, 0.05, 1.0, species)
      
      if (!is.null(gene_list) && length(gene_list) > 0) {
        cat("✅ Gene list preparation successful\n")
        cat("- Genes prepared:", length(gene_list), "\n")
        cat("- Sample genes:", paste(head(gene_list, 3), collapse = ", "), "\n")
        
        # Test pathway analysis
        cat("\n6. PATHWAY ANALYSIS TEST\n")
        cat("------------------------\n")
        
        if (analysis_type == "GO") {
          cat("Testing GO analysis...\n")
          result <- tryCatch({
            run_go_analysis(gene_list, species, "BP")
          }, error = function(e) {
            list(success = FALSE, error = e$message)
          })
          
          if (result$success) {
            cat("✅ GO analysis successful\n")
            cat("- Terms found:", nrow(result$data), "\n")
            cat("- Execution time:", round(result$execution_time, 2), "seconds\n")
          } else {
            cat("❌ GO analysis failed:", result$error, "\n")
          }
          
        } else if (analysis_type == "KEGG") {
          cat("Testing KEGG analysis...\n")
          result <- tryCatch({
            run_kegg_analysis(gene_list, species)
          }, error = function(e) {
            list(success = FALSE, error = e$message)
          })
          
          if (result$success) {
            cat("✅ KEGG analysis successful\n")
            cat("- Pathways found:", nrow(result$data), "\n")
          } else {
            cat("❌ KEGG analysis failed:", result$error, "\n")
          }
        }
        
      } else {
        cat("❌ Gene list preparation failed - no genes returned\n")
      }
      
    }, error = function(e) {
      cat("❌ Gene list preparation error:", e$message, "\n")
    })
  }
  
  # Test with sample data if no data provided
  if (is.null(deseq2_results)) {
    cat("\n4. SAMPLE DATA TEST\n")
    cat("-------------------\n")
    cat("Creating sample DESeq2 results for testing...\n")
    
    set.seed(123)
    sample_results <- data.frame(
      baseMean = runif(100, 10, 1000),
      log2FoldChange = rnorm(100, 0, 2),
      lfcSE = runif(100, 0.1, 0.5),
      stat = rnorm(100, 0, 3),
      pvalue = runif(100, 0, 1),
      padj = runif(100, 0, 1)
    )
    
    # Make some genes significant
    sample_results$padj[1:20] <- runif(20, 0.001, 0.049)
    sample_results$log2FoldChange[1:20] <- rnorm(20, 0, 3)
    rownames(sample_results) <- paste0("ENSG", sprintf("%011d", 1:100))
    
    # Recursive call with sample data
    cat("Running diagnosis with sample data...\n")
    diagnose_pathway_issue(sample_results, analysis_type, species, "brief")
  }
  
  cat("\n💡 TROUBLESHOOTING SUGGESTIONS\n")
  cat("==============================\n")
  cat("If you're experiencing issues, check:\n")
  cat("1. ✅ Package installation: install.packages(c('clusterProfiler', 'org.Hs.eg.db'))\n")
  cat("2. ✅ Internet connection: KEGG analysis requires online access\n")
  cat("3. ✅ Data format: Ensure DESeq2 results have required columns\n")
  cat("4. ✅ Gene IDs: Check if gene IDs are in correct format (Ensembl, Symbol, etc.)\n")
  cat("5. ✅ Significance filters: Ensure some genes pass padj and FC thresholds\n")
  cat("6. ✅ Species annotation: Verify correct species databases are installed\n")
  
  cat("\n📞 For specific errors, please provide:\n")
  cat("- The exact error message\n")
  cat("- Analysis type you're trying to run (GO, KEGG, GSEA)\n")
  cat("- Species you're analyzing\n")
  cat("- Size of your dataset\n")
  cat("- Any specific steps where it fails\n\n")
}

# Quick diagnosis function for immediate use
quick_pathway_diagnosis <- function() {
  cat("🚀 QUICK PATHWAY DIAGNOSIS\n")
  cat("=========================\n\n")
  
  # Check if pathway_analysis.R is loaded
  if (!exists("run_pathway_analysis")) {
    cat("❌ pathway_analysis.R not loaded\n")
    cat("💡 Run: source('pathway_analysis.R')\n\n")
  } else {
    cat("✅ pathway_analysis.R loaded\n\n")
  }
  
  # Test basic functionality
  cat("Testing basic pathway functions...\n")
  
  # Test clusterProfiler
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    cat("✅ clusterProfiler available\n")
  } else {
    cat("❌ clusterProfiler missing - install with: BiocManager::install('clusterProfiler')\n")
  }
  
  # Test organism databases
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cat("✅ Human annotations available\n")
  } else {
    cat("❌ Human annotations missing - install with: BiocManager::install('org.Hs.eg.db')\n")
  }
  
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    cat("✅ Mouse annotations available\n")
  } else {
    cat("⚠️ Mouse annotations missing - install with: BiocManager::install('org.Mm.eg.db')\n")
  }
  
  cat("\n🔧 To run full diagnosis with your data:\n")
  cat("diagnose_pathway_issue(your_deseq2_results, 'GO', 'human')\n\n")
}

# Error capture function
capture_pathway_error <- function(expression) {
  tryCatch({
    eval(expression)
  }, error = function(e) {
    cat("❌ ERROR CAPTURED:\n")
    cat("Message:", e$message, "\n")
    cat("Call:", deparse(e$call), "\n")
    
    # Specific error handling
    if (grepl("clusterProfiler", e$message)) {
      cat("💡 This appears to be a clusterProfiler issue\n")
      cat("Try: BiocManager::install('clusterProfiler')\n")
    } else if (grepl("org\\.", e$message)) {
      cat("💡 This appears to be an organism database issue\n")
      cat("Try: BiocManager::install('org.Hs.eg.db')\n")
    } else if (grepl("network|internet", e$message)) {
      cat("💡 This appears to be a network connectivity issue\n")
      cat("Check your internet connection for KEGG analysis\n")
    } else if (grepl("timeout", e$message)) {
      cat("💡 This appears to be a timeout issue\n")
      cat("Try reducing gene list size or using different analysis type\n")
    }
    
    return(list(error = TRUE, message = e$message))
  })
}

# Load and run quick diagnosis
cat("🔧 Pathway Analysis Diagnostic Tool Loaded\n")
cat("==========================================\n\n")

cat("Available functions:\n")
cat("• quick_pathway_diagnosis() - Quick system check\n")
cat("• diagnose_pathway_issue(data, type, species) - Full diagnosis\n")
cat("• capture_pathway_error(expression) - Error analysis\n\n")

cat("💡 To start, run: quick_pathway_diagnosis()\n\n")

# Auto-run quick diagnosis
quick_pathway_diagnosis()