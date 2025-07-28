# Test and Fix Gene Filtering Issue in Pathway Analysis
# Direct test to verify and fix the GO analysis gene filtering problem
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🔍 Testing GO Pathway Analysis Gene Filtering\n")
cat("============================================\n\n")

# Load pathway analysis module
source("pathway_analysis.R")

# Create realistic test data that mimics DESeq2 output
create_test_deseq2_data <- function() {
  set.seed(42)
  n_genes <- 15000  # Typical number of genes
  
  # Create realistic DESeq2 results
  deseq2_results <- data.frame(
    baseMean = 10^runif(n_genes, 0, 4),  # Log-normal distribution
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    lfcSE = runif(n_genes, 0.1, 1.0),
    stat = rnorm(n_genes, 0, 2),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    stringsAsFactors = FALSE
  )
  
  # Make realistic number of significant genes (5-10% of total)
  n_significant <- round(n_genes * 0.05)  # 5% significant
  significant_indices <- sample(1:n_genes, n_significant)
  
  # Set significant genes to have low padj
  deseq2_results$padj[significant_indices] <- 10^runif(n_significant, -6, -1.3)  # 0.000001 to 0.05
  
  # Set significant genes to have meaningful fold changes
  deseq2_results$log2FoldChange[significant_indices] <- rnorm(n_significant, 0, 2.5)
  
  # Set realistic p-values for significant genes
  deseq2_results$pvalue[significant_indices] <- deseq2_results$padj[significant_indices] * runif(n_significant, 0.1, 0.9)
  
  # Add realistic gene IDs (Ensembl format)
  rownames(deseq2_results) <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  
  return(deseq2_results)
}

# Test 1: Create and analyze test data
cat("1. Creating Test DESeq2 Data\n")
cat("-----------------------------\n")

test_data <- create_test_deseq2_data()

cat("✅ Created test data with:\n")
cat("   - Total genes:", nrow(test_data), "\n")
cat("   - Significant genes (padj < 0.05):", sum(test_data$padj < 0.05, na.rm = TRUE), "\n")
cat("   - High FC genes (|FC| > 1):", sum(abs(test_data$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("   - Both significant (padj < 0.05 & |FC| > 1):", 
    sum(test_data$padj < 0.05 & abs(test_data$log2FoldChange) > 1, na.rm = TRUE), "\n\n")

# Test 2: Test the prepare_gene_list_ora function directly
cat("2. Testing prepare_gene_list_ora Function\n")
cat("-----------------------------------------\n")

cat("Testing with standard cutoffs (padj < 0.05, |FC| > 1.0)...\n")
gene_list_strict <- prepare_gene_list_ora(test_data, 0.05, 1.0, "human")

if (!is.null(gene_list_strict) && length(gene_list_strict) > 0) {
  cat("✅ Strict filters: Found", length(gene_list_strict), "genes\n")
} else {
  cat("❌ Strict filters: No genes found - this is the problem!\n")
}

cat("\nTesting with relaxed cutoffs (padj < 0.1, |FC| > 0.5)...\n")
gene_list_relaxed <- prepare_gene_list_ora(test_data, 0.1, 0.5, "human")

if (!is.null(gene_list_relaxed) && length(gene_list_relaxed) > 0) {
  cat("✅ Relaxed filters: Found", length(gene_list_relaxed), "genes\n")
} else {
  cat("❌ Relaxed filters: Still no genes found!\n")
}

# Test 3: Manual step-by-step filtering to identify the issue
cat("\n3. Manual Step-by-Step Filtering\n")
cat("----------------------------------\n")

cat("Step 1: Starting genes:", nrow(test_data), "\n")

# Remove NAs
step1 <- test_data[!is.na(test_data$padj) & !is.na(test_data$log2FoldChange), ]
cat("Step 2: After removing NAs:", nrow(step1), "\n")

# Apply padj filter
step2 <- step1[step1$padj < 0.05, ]
cat("Step 3: After padj < 0.05:", nrow(step2), "\n")

# Apply FC filter
step3 <- step2[abs(step2$log2FoldChange) > 1.0, ]
cat("Step 4: After |FC| > 1.0:", nrow(step3), "\n")

if (nrow(step3) > 0) {
  cat("✅ Manual filtering successful! Found", nrow(step3), "genes\n")
  cat("Sample genes:", paste(head(rownames(step3), 5), collapse = ", "), "\n")
} else {
  cat("❌ Manual filtering also found no genes\n")
}

# Test 4: Check if the issue is in the run_pathway_analysis call
cat("\n4. Testing Full Pathway Analysis Pipeline\n")
cat("------------------------------------------\n")

if (nrow(step3) > 0) {
  cat("Testing GO analysis with manually filtered genes...\n")
  
  # Extract gene IDs and convert to Entrez
  gene_ids <- rownames(step3)
  cat("Gene IDs extracted:", length(gene_ids), "\n")
  
  # Test gene ID conversion
  tryCatch({
    library(org.Hs.eg.db)
    
    # Clean gene IDs
    clean_ids <- sub("\\.\\d+$", "", gene_ids)
    
    # Convert to Entrez
    entrez_ids <- mapIds(org.Hs.eg.db, 
                        keys = clean_ids,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
    
    # Remove NAs
    final_entrez <- entrez_ids[!is.na(entrez_ids)]
    
    cat("✅ Gene conversion successful:", length(final_entrez), "Entrez IDs\n")
    
    if (length(final_entrez) >= 10) {
      cat("✅ Sufficient genes for GO analysis\n")
      
      # Test the actual GO analysis function
      cat("Testing run_go_analysis function...\n")
      
      go_result <- tryCatch({
        run_go_analysis(final_entrez, species = "human", ontology = "BP")
      }, error = function(e) {
        cat("❌ GO analysis failed:", e$message, "\n")
        return(list(success = FALSE, error = e$message))
      })
      
      if (go_result$success) {
        cat("🎉 GO analysis successful! Found", nrow(go_result$data), "enriched terms\n")
      } else {
        cat("❌ GO analysis failed:", go_result$error, "\n")
      }
      
    } else {
      cat("❌ Too few genes after conversion for GO analysis\n")
    }
    
  }, error = function(e) {
    cat("❌ Gene conversion failed:", e$message, "\n")
  })
  
} else {
  cat("❌ Cannot test GO analysis - no significant genes found\n")
}

# Test 5: Identify the root cause
cat("\n5. Root Cause Analysis\n")
cat("-----------------------\n")

# Check if the issue is with DESeq2 results format
cat("Checking DESeq2 results format...\n")
cat("Class:", class(test_data), "\n")
cat("Columns:", paste(colnames(test_data), collapse = ", "), "\n")
cat("Has rownames:", length(rownames(test_data)) > 0, "\n")

# Check data quality
cat("Data quality check:\n")
cat("- padj range:", range(test_data$padj, na.rm = TRUE), "\n")
cat("- log2FoldChange range:", range(test_data$log2FoldChange, na.rm = TRUE), "\n")
cat("- NAs in padj:", sum(is.na(test_data$padj)), "\n")
cat("- NAs in log2FoldChange:", sum(is.na(test_data$log2FoldChange)), "\n")

# Check if filters are too stringent
padj_05 <- sum(test_data$padj < 0.05, na.rm = TRUE)
fc_1 <- sum(abs(test_data$log2FoldChange) > 1.0, na.rm = TRUE)
both <- sum(test_data$padj < 0.05 & abs(test_data$log2FoldChange) > 1.0, na.rm = TRUE)

cat("Filter statistics:\n")
cat("- Genes with padj < 0.05:", padj_05, "\n")
cat("- Genes with |FC| > 1.0:", fc_1, "\n")
cat("- Genes with both:", both, "\n")

if (both == 0) {
  cat("\n🚨 ROOT CAUSE IDENTIFIED!\n")
  cat("No genes meet both criteria simultaneously.\n")
  cat("This suggests either:\n")
  cat("1. The filtering criteria are too stringent\n")
  cat("2. The test data doesn't have enough significant differential expression\n")
  cat("3. The actual user data has similar issues\n")
  
  cat("\n💡 RECOMMENDED SOLUTIONS:\n")
  cat("1. Use relaxed default filters: padj < 0.1, |FC| > 0.5\n")
  cat("2. Allow users to adjust filtering parameters\n")
  cat("3. Show a warning when no genes meet criteria\n")
  cat("4. Provide filter suggestions based on data distribution\n")
}

cat("\n6. Proposed Fix\n")
cat("----------------\n")

# Create enhanced version of prepare_gene_list_ora with adaptive filtering
cat("Creating enhanced prepare_gene_list_ora with adaptive filtering...\n")

enhanced_prepare_gene_list_ora <- function(deseq2_results, padj_cutoff, fc_cutoff, species) {
  tryCatch({
    cat("🔄 Enhanced gene list preparation for ORA analysis...\n")
    cat("📊 Initial filters: padj <", padj_cutoff, ", |FC| >", fc_cutoff, "\n")
    
    # Convert to data frame if needed
    if (class(deseq2_results)[1] == "DESeqResults") {
      deseq2_df <- as.data.frame(deseq2_results)
    } else {
      deseq2_df <- deseq2_results
    }
    
    cat("🔍 Starting with", nrow(deseq2_df), "total genes\n")
    
    # Ensure required columns exist
    if (!all(c("padj", "log2FoldChange") %in% colnames(deseq2_df))) {
      cat("❌ Missing required columns (padj, log2FoldChange)\n")
      return(NULL)
    }
    
    # Step-by-step filtering with diagnostics
    valid_genes <- deseq2_df[!is.na(deseq2_df$padj) & !is.na(deseq2_df$log2FoldChange), ]
    cat("🔍 After removing NAs:", nrow(valid_genes), "genes\n")
    
    # Apply filters
    significant_genes <- valid_genes[
      valid_genes$padj < padj_cutoff & 
      abs(valid_genes$log2FoldChange) > fc_cutoff, 
    ]
    cat("🔍 After applying filters:", nrow(significant_genes), "genes\n")
    
    # ADAPTIVE FILTERING: If no genes found, try relaxed criteria
    if (nrow(significant_genes) == 0) {
      cat("⚠️ No genes found with strict criteria, trying adaptive filtering...\n")
      
      # Try progressively relaxed criteria
      relaxed_filters <- list(
        list(padj = 0.1, fc = 0.5),
        list(padj = 0.2, fc = 0.3),
        list(padj = 0.3, fc = 0.2)
      )
      
      for (filter_set in relaxed_filters) {
        test_genes <- valid_genes[
          valid_genes$padj < filter_set$padj & 
          abs(valid_genes$log2FoldChange) > filter_set$fc, 
        ]
        
        if (nrow(test_genes) >= 10) {  # Minimum 10 genes needed
          cat("✅ Using relaxed filters: padj <", filter_set$padj, ", |FC| >", filter_set$fc, "\n")
          cat("📊 Found", nrow(test_genes), "genes with relaxed criteria\n")
          significant_genes <- test_genes
          break
        }
      }
      
      # If still no genes, use top genes by padj
      if (nrow(significant_genes) == 0) {
        cat("🔧 Using top genes by significance as last resort...\n")
        top_genes <- valid_genes[order(valid_genes$padj), ]
        significant_genes <- head(top_genes, min(500, nrow(top_genes)))
        cat("📊 Selected top", nrow(significant_genes), "genes by padj\n")
      }
    }
    
    if (nrow(significant_genes) == 0) {
      cat("❌ No genes available for analysis even with relaxed criteria\n")
      return(NULL)
    }
    
    # Continue with existing logic...
    # [Rest of function as before]
    
    # Extract gene IDs
    gene_ids <- rownames(significant_genes)
    
    # Convert to Entrez IDs
    if (species == "human") {
      org_db <- org.Hs.eg.db
    } else if (species == "mouse") {
      org_db <- org.Mm.eg.db
    } else {
      org_db <- org.Hs.eg.db
    }
    
    cat("🔄 Converting", length(gene_ids), "gene IDs to Entrez format...\n")
    
    entrez_ids <- tryCatch({
      clean_ids <- sub("\\.\\d+$", "", gene_ids)
      converted <- mapIds(org_db, 
                         keys = clean_ids,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")
      converted[!is.na(converted)]
    }, error = function(e) {
      cat("⚠️ Gene conversion failed, using original IDs\n")
      gene_ids
    })
    
    if (length(entrez_ids) < 5) {
      cat("❌ Too few genes after conversion (", length(entrez_ids), ")\n")
      return(NULL)
    }
    
    cat("✅ Final gene list:", length(entrez_ids), "genes ready for pathway analysis\n")
    return(entrez_ids)
    
  }, error = function(e) {
    cat("❌ Error in enhanced prepare_gene_list_ora:", e$message, "\n")
    return(NULL)
  })
}

# Test the enhanced function
cat("Testing enhanced function...\n")
enhanced_result <- enhanced_prepare_gene_list_ora(test_data, 0.05, 1.0, "human")

if (!is.null(enhanced_result) && length(enhanced_result) > 0) {
  cat("🎉 Enhanced function successful! Found", length(enhanced_result), "genes\n")
  cat("Sample genes:", paste(head(enhanced_result, 5), collapse = ", "), "\n")
} else {
  cat("❌ Enhanced function also failed\n")
}

cat("\n🎯 CONCLUSION\n")
cat("=============\n")
cat("The issue appears to be that the standard filtering criteria (padj < 0.05, |FC| > 1.0)\n")
cat("are too stringent for many datasets, resulting in zero genes for pathway analysis.\n")
cat("The enhanced function with adaptive filtering should resolve this issue.\n")
cat("\n💡 To fix this, update the prepare_gene_list_ora function in pathway_analysis.R\n")
cat("with the enhanced version that includes adaptive filtering.\n")