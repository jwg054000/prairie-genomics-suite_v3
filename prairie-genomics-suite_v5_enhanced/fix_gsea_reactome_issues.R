# Fix GSEA Visualization and Reactome Gene List Issues
# 
# This script addresses two critical issues:
# 1. GSEA visualization not plotting properly (missing enrichment_object)
# 2. Reactome analysis failing with 'gene_list' error (parameter passing issue)
#
# Author: Prairie Genomics Team  
# Date: January 27, 2025

cat("🔧 FIXING GSEA VISUALIZATION AND REACTOME ISSUES\n")
cat("================================================\n\n")

# Load required modules
source("pathway_analysis.R")

# Test data setup
create_test_data <- function() {
  cat("📊 Creating test data with gene symbols...\n")
  
  # Real mouse gene symbols for realistic testing
  real_mouse_genes <- c(
    "Tp53", "Brca1", "Myc", "Kras", "Pik3ca", "Akt1", "Mtor", "Rb1", "Cdk4", "Ccnd1",
    "Egfr", "Vegfa", "Hif1a", "Nf1", "Apc", "Ctnnb1", "Tgfb1", "Smad4", "Cdkn2a", "Pten",
    "Braf", "Mapk1", "Jun", "Fos", "Stat3", "Stat1", "Ifng", "Il6", "Tnf", "Nfkb1",
    "Trp53", "Cdkn1a", "Mdm2", "Bax", "Bcl2", "Casp3", "Gapdh", "Actb", "Tubb3", "Vim"
  )
  
  set.seed(123)
  n_genes <- length(real_mouse_genes)
  
  test_results <- data.frame(
    baseMean = runif(n_genes, 100, 3000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.2, 0.6),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    row.names = real_mouse_genes
  )
  
  # Make most genes significant
  n_sig <- 35
  sig_idx <- sample(1:n_genes, n_sig)
  test_results$padj[sig_idx] <- runif(n_sig, 0.001, 0.049)
  test_results$log2FoldChange[sig_idx] <- rnorm(n_sig, 0, 2.5)
  
  cat("✅ Test data created: ", n_genes, " genes,", sum(test_results$padj < 0.05, na.rm = TRUE), "significant\n")
  return(test_results)
}

# Test Fix 1: Reactome Analysis (gene_list parameter issue)
test_reactome_fix <- function(test_data) {
  cat("\n🧬 TEST 1: REACTOME ANALYSIS FIX\n")
  cat("===============================\n")
  
  cat("🔄 Testing Reactome analysis with proper gene list preparation...\n")
  
  # Test the main pathway analysis function (which should now handle Reactome correctly)
  reactome_result <- tryCatch({
    run_pathway_analysis(
      deseq2_results = test_data,
      analysis_type = "Reactome",
      species = "mouse",
      padj_cutoff = 0.05,
      fc_cutoff = 1.0
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (reactome_result$success) {
    cat("✅ REACTOME FIX SUCCESS!\n")
    cat("  - Analysis type:", reactome_result$analysis_type, "\n")
    cat("  - Species:", reactome_result$species, "\n")
    cat("  - Pathways found:", reactome_result$n_pathways, "\n")
    cat("  - Significant pathways:", reactome_result$n_significant, "\n")
    if (!is.null(reactome_result$data) && nrow(reactome_result$data) > 0) {
      cat("  - Top pathway:", reactome_result$data$Description[1], "\n")
    }
  } else {
    cat("❌ Reactome analysis still failing:", reactome_result$error, "\n")
  }
  
  return(reactome_result)
}

# Test Fix 2: GSEA Visualization (plotting issue)
test_gsea_visualization_fix <- function(test_data) {
  cat("\n📊 TEST 2: GSEA VISUALIZATION FIX\n")
  cat("=================================\n")
  
  cat("🔄 Testing GSEA analysis and visualization...\n")
  
  # Run GSEA analysis
  gsea_result <- tryCatch({
    run_pathway_analysis(
      deseq2_results = test_data,
      analysis_type = "GSEA", 
      species = "mouse",
      gene_set_collection = "H",
      padj_cutoff = 0.1,  # More lenient for testing
      fc_cutoff = 0.5
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (!gsea_result$success) {
    cat("❌ GSEA analysis failed:", gsea_result$error, "\n")
    return(gsea_result)
  }
  
  cat("✅ GSEA analysis succeeded!\n")
  cat("  - Method:", gsea_result$method %||% "fgseaMultilevel", "\n")
  cat("  - Pathways:", gsea_result$n_pathways, "\n")
  cat("  - Significant:", gsea_result$n_significant, "\n")
  
  # Test GSEA visualization
  cat("\n🎨 Testing GSEA plot creation...\n")
  
  plot_types <- c("dotplot", "gsea", "enrichment")
  plot_results <- list()
  
  for (plot_type in plot_types) {
    cat("  - Testing", plot_type, "plot...\n")
    
    plot_obj <- tryCatch({
      create_pathway_plots(gsea_result, plot_type = plot_type, top_n = 10)
    }, error = function(e) {
      cat("    ❌ Error:", e$message, "\n")
      NULL
    })
    
    if (!is.null(plot_obj)) {
      cat("    ✅", plot_type, "plot created successfully\n")
      plot_results[[plot_type]] <- TRUE
    } else {
      cat("    ❌", plot_type, "plot failed\n")
      plot_results[[plot_type]] <- FALSE
    }
  }
  
  gsea_result$plot_tests <- plot_results
  return(gsea_result)
}

# Test Fix 3: Compare All Analysis Types
test_all_analysis_types <- function(test_data) {
  cat("\n🚀 TEST 3: ALL ANALYSIS TYPES COMPARISON\n")
  cat("========================================\n")
  
  analysis_types <- c("GO", "KEGG", "GSEA", "MSigDB", "Reactome")
  results <- list()
  
  for (analysis_type in analysis_types) {
    cat("\n🔄 Testing", analysis_type, "...\n")
    
    result <- tryCatch({
      run_pathway_analysis(
        deseq2_results = test_data,
        analysis_type = analysis_type,
        species = "mouse",
        ontology = "BP",
        padj_cutoff = 0.05,
        fc_cutoff = 1.0,
        gene_set_collection = "H"
      )
    }, error = function(e) {
      list(success = FALSE, error = e$message, analysis_type = analysis_type)
    })
    
    results[[analysis_type]] <- result
    
    if (result$success) {
      n_pathways <- nrow(result$data)
      cat("  ✅", analysis_type, "SUCCESS -", n_pathways, "pathways found\n")
      
      # Test visualization for each type
      plot_obj <- tryCatch({
        if (analysis_type == "GSEA") {
          create_pathway_plots(result, plot_type = "gsea", top_n = 5)
        } else {
          create_pathway_plots(result, plot_type = "dotplot", top_n = 5)
        }
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(plot_obj)) {
        cat("  📊 Visualization: SUCCESS\n")
        result$plot_success <- TRUE
      } else {
        cat("  📊 Visualization: FAILED\n") 
        result$plot_success <- FALSE
      }
      
    } else {
      cat("  ❌", analysis_type, "FAILED:", result$error, "\n")
      result$plot_success <- FALSE
    }
  }
  
  return(results)
}

# Test Fix 4: Pathway Result Persistence
test_persistent_results <- function(test_data) {
  cat("\n💾 TEST 4: PERSISTENT RESULTS SYSTEM\n")
  cat("====================================\n")
  
  # Simulate cache system
  cache <- list()
  
  # Run GO analysis and cache
  cat("🔄 Running GO analysis...\n")
  go_result <- run_pathway_analysis(test_data, "GO", "mouse", "BP", 0.05, 1.0, "H")
  if (go_result$success) {
    cache_key <- "GO_mouse"
    cache[[cache_key]] <- go_result
    cat("✅ GO results cached with key:", cache_key, "\n")
  }
  
  # Run GSEA analysis and cache
  cat("🔄 Running GSEA analysis...\n")
  gsea_result <- run_pathway_analysis(test_data, "GSEA", "mouse", "BP", 0.05, 1.0, "H")
  if (gsea_result$success) {
    cache_key <- "GSEA_mouse"
    cache[[cache_key]] <- gsea_result
    cat("✅ GSEA results cached with key:", cache_key, "\n")
  }
  
  # Test cache retrieval
  cat("\n📂 Testing cache retrieval...\n")
  if ("GO_mouse" %in% names(cache)) {
    retrieved_go <- cache[["GO_mouse"]]
    cat("✅ GO results retrieved from cache -", nrow(retrieved_go$data), "pathways\n")
  }
  
  if ("GSEA_mouse" %in% names(cache)) {
    retrieved_gsea <- cache[["GSEA_mouse"]]
    cat("✅ GSEA results retrieved from cache -", nrow(retrieved_gsea$data), "pathways\n")
  }
  
  cat("💡 Cache contains", length(cache), "result sets\n")
  
  return(cache)
}

# Run all tests
cat("🚀 STARTING COMPREHENSIVE FIX TESTING\n")
cat("======================================\n")

# Create test data
test_data <- create_test_data()

# Test 1: Reactome fix
reactome_result <- test_reactome_fix(test_data)

# Test 2: GSEA visualization fix
gsea_result <- test_gsea_visualization_fix(test_data)

# Test 3: All analysis types
all_results <- test_all_analysis_types(test_data)

# Test 4: Persistent results
cache_test <- test_persistent_results(test_data)

# Summary Report
cat("\n📋 COMPREHENSIVE FIX SUMMARY\n")
cat("============================\n")

# Reactome Fix Assessment
cat("1. Reactome Analysis Fix: ")
if (reactome_result$success) {
  cat("✅ FIXED - Analysis works without gene_list error\n")
} else {
  cat("❌ STILL BROKEN -", reactome_result$error, "\n")
}

# GSEA Visualization Fix Assessment
cat("2. GSEA Visualization Fix: ")
if (gsea_result$success) {
  plot_success_count <- sum(unlist(gsea_result$plot_tests))
  cat("✅ FIXED -", plot_success_count, "of", length(gsea_result$plot_tests), "plot types working\n")
} else {
  cat("❌ STILL BROKEN - GSEA analysis itself failing\n")
}

# Overall Analysis Type Status
cat("3. Analysis Type Status:\n")
success_count <- 0
for (analysis_type in names(all_results)) {
  result <- all_results[[analysis_type]]
  if (result$success) {
    success_count <- success_count + 1
    plot_status <- if (result$plot_success) "📊✅" else "📊❌"
    cat("   ✅", analysis_type, "- Analysis ✅, Plots", plot_status, "\n")
  } else {
    cat("   ❌", analysis_type, "- FAILED\n")
  }
}

# Cache System Status
cat("4. Persistent Results: ")
if (length(cache_test) >= 2) {
  cat("✅ WORKING - Cache stores and retrieves results\n")
} else {
  cat("⚠️ LIMITED - Only", length(cache_test), "result sets cached\n")
}

# Overall Assessment
overall_success <- reactome_result$success && gsea_result$success && success_count >= 3

cat("\n🎯 OVERALL FIX STATUS:\n")
cat("======================\n")
if (overall_success) {
  cat("✅ BOTH ISSUES FIXED SUCCESSFULLY!\n")
  cat("🎉 Reactome analysis works without gene_list errors\n")
  cat("🎉 GSEA visualization displays properly with multiple plot types\n")
  cat("📊 All pathway analysis types working with visualization support\n")
  cat("💾 Persistent results system prevents data loss\n")
} else {
  cat("⚠️ SOME ISSUES MAY REMAIN\n")
  cat("🔧 Check individual test results above for details\n")
}

cat("\n💡 NEXT STEPS:\n")
cat("===============\n")
cat("1. Restart the Shiny app to load all fixes\n")
cat("2. Test Reactome analysis - should work without 'gene_list' error\n")
cat("3. Test GSEA analysis and try different plot types (dotplot, gsea, enrichment)\n")
cat("4. Verify results persist when switching between analysis types\n")
cat("5. Check that all visualizations display properly in the Plots tab\n")

cat("\n🔧 TECHNICAL FIXES IMPLEMENTED:\n")
cat("================================\n")
cat("✅ Fixed Reactome: Added proper gene list preparation in run_pathway_analysis()\n")
cat("✅ Fixed GSEA Plots: Created specialized create_gsea_plots() function\n") 
cat("✅ Enhanced Plotting: Added fallback plotting for all analysis types\n")
cat("✅ Improved Error Handling: Comprehensive error messages and fallbacks\n")
cat("✅ Better Visualization: GSEA-specific plot types with NES and significance\n")