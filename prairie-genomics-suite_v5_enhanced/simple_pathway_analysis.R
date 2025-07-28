# 🧬 SIMPLIFIED PATHWAY ANALYSIS - PHASE 4B
# Expert-Validated Pathway Analysis for All Comparisons
# Using only base R and available Bioconductor packages

library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)

cat("🧬 PHASE 4B: SIMPLIFIED PATHWAY ANALYSIS PIPELINE\n")
cat("=================================================\n")
cat("Expert validation: 'they are all completely accurate!' ✅\n\n")

# Expert-validated parameters
PADJ_THRESHOLD <- 0.05
FC_THRESHOLD <- 1.5

# Define comparison files
comparison_files <- c(
  "MC9_vs_M1245_detailed_results.csv",
  "MC9_vs_M242_detailed_results.csv", 
  "MLM_vs_M1245_detailed_results.csv",
  "MLM_vs_M242_detailed_results.csv",
  "M1245_vs_M242_detailed_results.csv"
)

comparison_names <- c(
  "MC9_vs_M1245",
  "MC9_vs_M242",
  "MLM_vs_M1245", 
  "MLM_vs_M242",
  "M1245_vs_M242"
)

# Create output directory
output_dir <- "pathway_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("📁 Created output directory:", output_dir, "\n\n")
}

# Initialize results summary
all_results_summary <- data.frame(
  Comparison = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  GO_BP_Pathways = integer(),
  GO_MF_Pathways = integer(),
  KEGG_Pathways = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

# Process each comparison
for (i in seq_along(comparison_files)) {
  
  file_path <- comparison_files[i]
  comp_name <- comparison_names[i]
  
  cat("🔬 Processing:", comp_name, "\n")
  cat("─────────────────────────────────\n")
  
  if (!file.exists(file_path)) {
    cat("❌ File not found:", file_path, "\n\n")
    next
  }
  
  tryCatch({
    
    # Load data
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("📊 Loaded", nrow(data), "genes from", file_path, "\n")
    
    # Filter significant genes using expert-validated thresholds
    significant_data <- data[!is.na(data$padj) & 
                            !is.na(data$external_gene_name) &
                            data$padj < PADJ_THRESHOLD & 
                            abs(data$log2FoldChange) >= log2(FC_THRESHOLD), ]
    
    if (nrow(significant_data) == 0) {
      cat("⚠️ No significant genes found for", comp_name, "\n\n")
      next
    }
    
    cat("✅ Found", nrow(significant_data), "significant genes\n")
    
    # Extract gene symbols
    gene_symbols <- unique(significant_data$external_gene_name)
    gene_symbols <- gene_symbols[gene_symbols != "" & !is.na(gene_symbols)]
    
    cat("🧬 Analyzing", length(gene_symbols), "unique gene symbols\n")
    
    # Initialize pathway counts
    go_bp_count <- 0
    go_mf_count <- 0
    kegg_count <- 0
    
    # GO Biological Process Analysis
    cat("🎯 Running GO Biological Process analysis...\n")
    go_bp_result <- NULL
    tryCatch({
      go_bp_result <- enrichGO(
        gene = gene_symbols,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
      
      if (!is.null(go_bp_result) && nrow(go_bp_result@result) > 0) {
        go_bp_count <- nrow(go_bp_result@result)
        cat("   ✅ Found", go_bp_count, "GO BP pathways\n")
        
        # Save top results
        top_go_bp <- head(go_bp_result@result, 10)
        write.csv(top_go_bp, file.path(output_dir, paste0(comp_name, "_GO_BP_top10.csv")), row.names = FALSE)
      } else {
        cat("   ⚠️ No significant GO BP pathways found\n")
      }
    }, error = function(e) {
      cat("   ❌ GO BP analysis error:", e$message, "\n")
    })
    
    # GO Molecular Function Analysis
    cat("🔬 Running GO Molecular Function analysis...\n")
    go_mf_result <- NULL
    tryCatch({
      go_mf_result <- enrichGO(
        gene = gene_symbols,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
      
      if (!is.null(go_mf_result) && nrow(go_mf_result@result) > 0) {
        go_mf_count <- nrow(go_mf_result@result)
        cat("   ✅ Found", go_mf_count, "GO MF pathways\n")
        
        # Save top results
        top_go_mf <- head(go_mf_result@result, 10)
        write.csv(top_go_mf, file.path(output_dir, paste0(comp_name, "_GO_MF_top10.csv")), row.names = FALSE)
      } else {
        cat("   ⚠️ No significant GO MF pathways found\n")
      }
    }, error = function(e) {
      cat("   ❌ GO MF analysis error:", e$message, "\n")
    })
    
    # KEGG Pathway Analysis
    cat("🛤️ Running KEGG pathway analysis...\n")
    kegg_result <- NULL
    tryCatch({
      # Convert to Entrez IDs for KEGG
      entrez_conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      
      if (nrow(entrez_conversion) > 0) {
        kegg_result <- enrichKEGG(
          gene = entrez_conversion$ENTREZID,
          organism = "mmu",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          minGSSize = 10,
          maxGSSize = 500
        )
        
        if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
          # Convert back to gene symbols
          kegg_result <- setReadable(kegg_result, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
          kegg_count <- nrow(kegg_result@result)
          cat("   ✅ Found", kegg_count, "KEGG pathways\n")
          
          # Save top results
          top_kegg <- head(kegg_result@result, 10)
          write.csv(top_kegg, file.path(output_dir, paste0(comp_name, "_KEGG_top10.csv")), row.names = FALSE)
        } else {
          cat("   ⚠️ No significant KEGG pathways found\n")
        }
      } else {
        cat("   ⚠️ No genes converted to Entrez IDs for KEGG\n")
      }
    }, error = function(e) {
      cat("   ❌ KEGG analysis error:", e$message, "\n")
    })
    
    # Create visualizations if results exist
    if (!is.null(go_bp_result) && nrow(go_bp_result@result) > 0) {
      tryCatch({
        png(file.path(output_dir, paste0(comp_name, "_GO_BP_dotplot.png")), width = 1200, height = 800, res = 150)
        print(dotplot(go_bp_result, showCategory = 15, title = paste("GO BP -", comp_name)))
        dev.off()
        cat("   📊 Saved GO BP dotplot\n")
      }, error = function(e) {
        cat("   ⚠️ Could not create GO BP plot:", e$message, "\n")
      })
    }
    
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      tryCatch({
        png(file.path(output_dir, paste0(comp_name, "_KEGG_dotplot.png")), width = 1200, height = 800, res = 150)
        print(dotplot(kegg_result, showCategory = 15, title = paste("KEGG Pathways -", comp_name)))
        dev.off()
        cat("   📊 Saved KEGG dotplot\n")
      }, error = function(e) {
        cat("   ⚠️ Could not create KEGG plot:", e$message, "\n")
      })
    }
    
    # Add to summary
    all_results_summary <- rbind(all_results_summary, data.frame(
      Comparison = comp_name,
      Total_Genes = nrow(data),
      Significant_Genes = nrow(significant_data),
      GO_BP_Pathways = go_bp_count,
      GO_MF_Pathways = go_mf_count,
      KEGG_Pathways = kegg_count,
      Status = "COMPLETE",
      stringsAsFactors = FALSE
    ))
    
    cat("✅ Completed", comp_name, "\n")
    cat("   - GO BP:", go_bp_count, "| GO MF:", go_mf_count, "| KEGG:", kegg_count, "\n\n")
    
  }, error = function(e) {
    cat("❌ Error processing", comp_name, ":", e$message, "\n\n")
    
    # Add error to summary
    all_results_summary <- rbind(all_results_summary, data.frame(
      Comparison = comp_name,
      Total_Genes = 0,
      Significant_Genes = 0,
      GO_BP_Pathways = 0,
      GO_MF_Pathways = 0,
      KEGG_Pathways = 0,
      Status = "ERROR",
      stringsAsFactors = FALSE
    ))
  })
}

# Save comprehensive summary
write.csv(all_results_summary, file.path(output_dir, "PATHWAY_ANALYSIS_SUMMARY.csv"), row.names = FALSE)

# Create summary report
sink(file.path(output_dir, "PATHWAY_ANALYSIS_REPORT.txt"))

cat("🧬 PHASE 4B PATHWAY ANALYSIS SUMMARY\n")
cat("====================================\n\n")
cat("Analysis Date:", format(Sys.Date(), "%B %d, %Y"), "\n")
cat("Expert Validation: Joshua Garton - 'they are all completely accurate!' ✅\n")
cat("Parameters: p.adj < 0.05, FC >= 1.5x\n\n")

cat("RESULTS OVERVIEW:\n")
cat("=================\n")
for (i in 1:nrow(all_results_summary)) {
  row <- all_results_summary[i, ]
  cat(sprintf("%-15s | %4d genes | %4d sig | %3d GO-BP | %3d GO-MF | %3d KEGG | %s\n",
              row$Comparison, row$Total_Genes, row$Significant_Genes, 
              row$GO_BP_Pathways, row$GO_MF_Pathways, row$KEGG_Pathways, row$Status))
}

cat("\nTOTAL PATHWAYS IDENTIFIED:\n")
cat("==========================\n")
cat("GO Biological Process:", sum(all_results_summary$GO_BP_Pathways), "\n")
cat("GO Molecular Function:", sum(all_results_summary$GO_MF_Pathways), "\n")
cat("KEGG Pathways:", sum(all_results_summary$KEGG_Pathways), "\n")
cat("Total Pathways:", sum(all_results_summary$GO_BP_Pathways + all_results_summary$GO_MF_Pathways + all_results_summary$KEGG_Pathways), "\n\n")

cat("FILES GENERATED:\n")
cat("================\n")
cat("- PATHWAY_ANALYSIS_SUMMARY.csv - Complete results table\n")
cat("- *_GO_BP_top10.csv - Top 10 GO Biological Process pathways per comparison\n")
cat("- *_GO_MF_top10.csv - Top 10 GO Molecular Function pathways per comparison\n")
cat("- *_KEGG_top10.csv - Top 10 KEGG pathways per comparison\n")
cat("- *_GO_BP_dotplot.png - GO BP visualization plots\n")
cat("- *_KEGG_dotplot.png - KEGG visualization plots\n\n")

cat("🎯 PHASE 4B COMPLETE!\n")
cat("Expert-validated pathway analysis applied to all comparisons\n")
cat("Results ready for biological interpretation and Phase 4C\n")

sink()

cat("🎉 PHASE 4B PATHWAY ANALYSIS COMPLETE!\n")
cat("======================================\n")
cat("Total pathways identified across all comparisons:\n")
cat("- GO Biological Process:", sum(all_results_summary$GO_BP_Pathways), "\n")
cat("- GO Molecular Function:", sum(all_results_summary$GO_MF_Pathways), "\n") 
cat("- KEGG Pathways:", sum(all_results_summary$KEGG_Pathways), "\n")
cat("- Grand Total:", sum(all_results_summary$GO_BP_Pathways + all_results_summary$GO_MF_Pathways + all_results_summary$KEGG_Pathways), "pathways\n\n")

cat("📁 All results saved in:", output_dir, "\n")
cat("📊 Summary report: PATHWAY_ANALYSIS_REPORT.txt\n")
cat("📈 Individual CSV files and plots generated for each comparison\n\n")

cat("🎯 Ready for expert review and Phase 4C!\n")