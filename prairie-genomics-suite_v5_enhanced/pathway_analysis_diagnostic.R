# Diagnostic version of pathway analysis to find the exact error location
# This wraps every potentially problematic operation in detailed error tracking

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Create a diagnostic wrapper that tracks exactly where errors occur
diagnostic_wrapper <- function(operation_name, operation_func) {
  cat("\n🔍 DIAGNOSTIC: Attempting", operation_name, "...\n")
  result <- tryCatch({
    res <- operation_func()
    cat("   ✅ SUCCESS:", operation_name, "completed\n")
    res
  }, error = function(e) {
    cat("   ❌ ERROR in", operation_name, ":\n")
    cat("      Message:", e$message, "\n")
    cat("      Call:", deparse(e$call), "\n")
    
    # Try to get more context
    if (exists(".Last.value")) {
      cat("      Last value type:", class(.Last.value), "\n")
    }
    
    # Return NULL to continue diagnosis
    NULL
  }, warning = function(w) {
    cat("   ⚠️  WARNING in", operation_name, ":", w$message, "\n")
    invokeRestart("muffleWarning")
  })
  
  return(result)
}

# Diagnostic GO analysis function
run_go_analysis_diagnostic <- function(gene_list, species = "mouse", ontology = "BP") {
  cat("\n========================================\n")
  cat("🏥 DIAGNOSTIC GO ANALYSIS\n")
  cat("========================================\n")
  
  # Step 1: Check input
  diagnostic_wrapper("input validation", function() {
    cat("   Gene list length:", length(gene_list), "\n")
    cat("   Gene list sample:", paste(head(gene_list, 3), collapse=", "), "\n")
    cat("   Species:", species, "\n")
    cat("   Ontology:", ontology, "\n")
    TRUE
  })
  
  # Step 2: Determine organism database
  org_db <- diagnostic_wrapper("organism database selection", function() {
    if (species == "mouse") {
      org.Mm.eg.db
    } else if (species == "human") {
      org.Hs.eg.db
    } else {
      stop("Unsupported species")
    }
  })
  
  if (is.null(org_db)) return(NULL)
  
  # Step 3: Run enrichGO
  go_results <- diagnostic_wrapper("enrichGO execution", function() {
    enrichGO(
      gene = gene_list,
      OrgDb = org_db,
      ont = ontology,
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05,
      minGSSize = 10,
      maxGSSize = 300,
      readable = FALSE
    )
  })
  
  if (is.null(go_results)) return(NULL)
  
  # Step 4: Check result structure
  diagnostic_wrapper("result structure check", function() {
    cat("   Result class:", class(go_results), "\n")
    cat("   Has @result slot:", .hasSlot(go_results, "result"), "\n")
    if (.hasSlot(go_results, "result")) {
      cat("   @result class:", class(go_results@result), "\n")
      cat("   @result dimensions:", dim(go_results@result), "\n")
    }
    TRUE
  })
  
  # Step 5: Extract data frame
  go_df <- diagnostic_wrapper("data frame extraction", function() {
    as.data.frame(go_results@result)
  })
  
  if (is.null(go_df)) return(NULL)
  
  # Step 6: Check data frame structure
  diagnostic_wrapper("data frame structure check", function() {
    cat("   Dimensions:", nrow(go_df), "x", ncol(go_df), "\n")
    cat("   Column names:", paste(colnames(go_df), collapse=", "), "\n")
    cat("   Column classes:\n")
    for (col in colnames(go_df)) {
      cat("     ", col, ":", class(go_df[[col]]), "\n")
    }
    TRUE
  })
  
  # Step 7: Test individual column access
  diagnostic_wrapper("column access test", function() {
    for (col in c("p.adjust", "Count", "pvalue")) {
      if (col %in% colnames(go_df)) {
        cat("   Testing go_df$", col, "... ", sep="")
        test <- go_df[[col]]
        cat("OK (", length(test), " values)\n", sep="")
      } else {
        cat("   Column '", col, "' NOT FOUND\n", sep="")
      }
    }
    TRUE
  })
  
  # Step 8: Test logical operations
  diagnostic_wrapper("logical operations test", function() {
    if ("p.adjust" %in% colnames(go_df)) {
      cat("   Testing go_df$p.adjust <= 0.05...\n")
      test1 <- go_df$p.adjust <= 0.05
      cat("     Result length:", length(test1), "\n")
      cat("     TRUE count:", sum(test1, na.rm=TRUE), "\n")
    }
    
    if ("Count" %in% colnames(go_df)) {
      cat("   Testing go_df$Count >= 3...\n")
      test2 <- go_df$Count >= 3
      cat("     Result length:", length(test2), "\n")
      cat("     TRUE count:", sum(test2, na.rm=TRUE), "\n")
    }
    TRUE
  })
  
  # Step 9: Test subsetting with different methods
  cat("\n🔬 TESTING DIFFERENT SUBSETTING METHODS:\n")
  
  # Method 1: Direct logical indexing
  subset1 <- diagnostic_wrapper("Method 1: Direct logical indexing", function() {
    go_df[go_df$p.adjust <= 0.05, ]
  })
  
  # Method 2: which() indexing
  subset2 <- diagnostic_wrapper("Method 2: which() indexing", function() {
    indices <- which(go_df$p.adjust <= 0.05)
    go_df[indices, ]
  })
  
  # Method 3: subset() function
  subset3 <- diagnostic_wrapper("Method 3: subset() function", function() {
    subset(go_df, p.adjust <= 0.05)
  })
  
  # Method 4: Combined conditions
  subset4 <- diagnostic_wrapper("Method 4: Combined conditions", function() {
    go_df[go_df$p.adjust <= 0.05 & go_df$Count >= 3, ]
  })
  
  # Method 5: Step-by-step filtering
  subset5 <- diagnostic_wrapper("Method 5: Step-by-step filtering", function() {
    # Create condition vectors separately
    cond1 <- go_df$p.adjust <= 0.05
    cond2 <- go_df$Count >= 3
    combined <- cond1 & cond2
    
    cat("     Condition 1 TRUE:", sum(cond1, na.rm=TRUE), "\n")
    cat("     Condition 2 TRUE:", sum(cond2, na.rm=TRUE), "\n")
    cat("     Combined TRUE:", sum(combined, na.rm=TRUE), "\n")
    
    # Now subset
    go_df[combined, ]
  })
  
  # Method 6: dplyr-style filtering (if available)
  subset6 <- diagnostic_wrapper("Method 6: Base R filter", function() {
    # Use base R to avoid dplyr dependency
    keep_rows <- rep(TRUE, nrow(go_df))
    keep_rows <- keep_rows & (go_df$p.adjust <= 0.05)
    keep_rows <- keep_rows & (go_df$Count >= 3)
    go_df[keep_rows, ]
  })
  
  cat("\n📊 DIAGNOSTIC SUMMARY:\n")
  cat("Original rows:", nrow(go_df), "\n")
  if (!is.null(subset1)) cat("Method 1 rows:", nrow(subset1), "\n")
  if (!is.null(subset2)) cat("Method 2 rows:", nrow(subset2), "\n")
  if (!is.null(subset3)) cat("Method 3 rows:", nrow(subset3), "\n")
  if (!is.null(subset4)) cat("Method 4 rows:", nrow(subset4), "\n")
  if (!is.null(subset5)) cat("Method 5 rows:", nrow(subset5), "\n")
  if (!is.null(subset6)) cat("Method 6 rows:", nrow(subset6), "\n")
  
  return(list(
    go_results = go_results,
    go_df = go_df,
    subsets = list(
      method1 = subset1,
      method2 = subset2,
      method3 = subset3,
      method4 = subset4,
      method5 = subset5,
      method6 = subset6
    )
  ))
}

# Test the diagnostic function
cat("🧪 RUNNING DIAGNOSTIC TEST...\n")
test_genes <- c("11545", "74778", "17869", "16590", "12043", "13982", "16193", "11593", "20926", "75560")

result <- run_go_analysis_diagnostic(test_genes)

if (!is.null(result)) {
  cat("\n✅ Diagnostic completed successfully!\n")
  cat("We can now identify which subsetting method works and apply it to the main code.\n")
} else {
  cat("\n❌ Diagnostic failed - check error messages above.\n")
}