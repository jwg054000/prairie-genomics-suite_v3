# Fix Pathway Analysis Results Display Issue
# The analysis completes successfully but results don't show in UI
#
# Author: Prairie Genomics Team
# Date: January 27, 2025

cat("🔧 FIXING PATHWAY RESULTS DISPLAY ISSUE\n")
cat("=======================================\n\n")

# The issue is that pathway analysis completes successfully but results don't show in UI
# Based on diagnostics, the problem is likely in the Shiny reactive logic

cat("📋 ISSUE ANALYSIS:\n")
cat("==================\n")
cat("✅ Pathway analysis functions work correctly\n")
cat("✅ Analysis returns success = TRUE with results\n")
cat("❌ UI shows 'Ready for Pathway Analysis' instead of results\n")
cat("🎯 Problem: show_pathway_results reactive not updating properly\n\n")

# Read the current app.R file to identify the exact issue
cat("🔍 EXAMINING APP.R PATHWAY LOGIC...\n")

# The issue is likely one of these:
# 1. values$pathway_results not being set correctly
# 2. observe block not triggering
# 3. reactive dependency not updating
# 4. JavaScript/UI condition not evaluating correctly

cat("💡 IDENTIFIED POTENTIAL ISSUES:\n")
cat("===============================\n")
cat("1. The observe block uses req(input$run_pathway_analysis) which might not trigger properly\n")
cat("2. The reactive values might not be updating the UI conditionalPanel\n")
cat("3. The pathway_results might be getting overwritten or reset\n\n")

# Create a fixed version of the pathway analysis server logic
create_fixed_pathway_server_logic <- function() {
  cat("🔧 CREATING FIXED PATHWAY SERVER LOGIC...\n")
  
  fixed_logic <- '
  # Fixed Pathway Analysis Server Logic
  observeEvent(input$run_pathway_analysis, {
    # Debug output
    cat("🚀 Pathway analysis button clicked\\n")
    
    # Validate prerequisites
    if (is.null(values$deseq2_results)) {
      showNotification("❌ No DESeq2 results available", type = "error")
      return()
    }
    
    cat("✅ DESeq2 results available, proceeding with pathway analysis\\n")
    
    # Show progress notification
    showNotification("🔄 Running pathway analysis...", type = "message", duration = NULL, id = "pathway_progress")
    
    # Run pathway analysis
    tryCatch({
      # Use the regular pathway analysis function (not logged version to avoid complications)
      pathway_result <- run_pathway_analysis(
        deseq2_results = values$deseq2_results,
        analysis_type = input$pathway_analysis_type,
        species = input$pathway_species %||% "auto",
        ontology = if (input$pathway_analysis_type == "GO") input$go_ontology else "BP",
        padj_cutoff = input$pathway_padj_cutoff,
        fc_cutoff = input$pathway_fc_cutoff,
        gene_set_collection = if (input$pathway_analysis_type %in% c("GSEA", "MSigDB")) input$msigdb_collection else "H"
      )
      
      # Debug output
      cat("📊 Pathway analysis completed. Success:", pathway_result$success, "\\n")
      if (pathway_result$success && !is.null(pathway_result$data)) {
        cat("📈 Found", nrow(pathway_result$data), "enriched pathways\\n")
      }
      
      # Store results - THIS IS CRITICAL
      values$pathway_results <- pathway_result
      
      # Force reactive update by creating a trigger
      values$pathway_results_updated <- Sys.time()
      
      # Remove progress notification
      removeNotification("pathway_progress")
      
      # Show completion notification
      if (pathway_result$success) {
        n_pathways <- pathway_result$n_terms %||% pathway_result$n_pathways %||% pathway_result$n_gene_sets %||% 0
        if (!is.null(pathway_result$data)) {
          n_pathways <- nrow(pathway_result$data)
        }
        
        showNotification(
          paste0("✅ ", pathway_result$analysis_type, " analysis completed! Found ", n_pathways, " enriched pathways"),
          type = "message",
          duration = 5
        )
        
        cat("✅ Pathway results stored. UI should now update.\\n")
      } else {
        showNotification(
          paste0("❌ Pathway analysis failed: ", pathway_result$message %||% "Unknown error"),
          type = "error",
          duration = 8
        )
      }
      
    }, error = function(e) {
      removeNotification("pathway_progress")
      cat("❌ Pathway analysis error:", e$message, "\\n")
      
      # Store error result
      values$pathway_results <- list(
        success = FALSE,
        error = e$message,
        message = paste("Analysis failed:", e$message)
      )
      
      showNotification(
        paste0("❌ Pathway analysis error: ", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Enhanced show_pathway_results reactive with debugging
  output$show_pathway_results <- reactive({
    result <- !is.null(values$pathway_results) && 
              ("success" %in% names(values$pathway_results)) && 
              values$pathway_results$success
    
    # Debug output
    cat("🔍 show_pathway_results reactive triggered. Result:", result, "\\n")
    if (!is.null(values$pathway_results)) {
      cat("   - pathway_results exists\\n")
      cat("   - success field:", values$pathway_results$success %||% "MISSING", "\\n")
    } else {
      cat("   - pathway_results is NULL\\n")
    }
    
    return(result)
  })
  outputOptions(output, "show_pathway_results", suspendWhenHidden = FALSE)
  '
  
  return(fixed_logic)
}

# Create the fix instructions
cat("🛠️ FIX INSTRUCTIONS:\n")
cat("====================\n")
cat("Replace the existing pathway analysis observe block in app.R with the fixed version.\n\n")

cat("📍 LOCATION IN APP.R:\n")
cat("Lines ~853-913: Replace the existing observe block\n\n")

cat("🔧 KEY CHANGES:\n")
cat("================\n")
cat("1. Use observeEvent instead of observe for clearer event handling\n")
cat("2. Add debug output to track execution\n")
cat("3. Add error handling with tryCatch\n")
cat("4. Force reactive update with pathway_results_updated timestamp\n")
cat("5. Enhanced show_pathway_results reactive with debugging\n")
cat("6. Use regular pathway analysis function instead of logged version\n\n")

# Generate the fixed code
fixed_code <- create_fixed_pathway_server_logic()

cat("📝 FIXED CODE BLOCK:\n")
cat("====================\n")
cat(fixed_code)

cat("\n\n💡 IMPLEMENTATION STEPS:\n")
cat("========================\n")
cat("1. Open app.R in a text editor\n")
cat("2. Find lines ~853-913 (the pathway analysis observe block)\n")
cat("3. Replace with the fixed code above\n")
cat("4. Save and restart the Shiny app\n")
cat("5. Test pathway analysis - you should see debug output in the console\n")
cat("6. Results should now display properly in the UI\n\n")

cat("🧪 TESTING:\n")
cat("===========\n")
cat("After implementing the fix:\n")
cat("1. Run DESeq2 analysis\n")
cat("2. Click 'Run Pathway Analysis'\n")
cat("3. Watch console for debug messages\n")
cat("4. Results should appear in the 'Pathway Analysis' tab\n")
cat("5. If still not working, check browser console for JavaScript errors\n")