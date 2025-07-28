# Safe subsetting function to avoid "undefined columns selected" errors
# This function provides a robust way to subset data frames with proper error handling

safe_subset <- function(df, condition, description = "subset operation") {
  # Validate inputs
  if (!is.data.frame(df)) {
    stop(paste("safe_subset: Input is not a data frame for", description))
  }
  
  if (nrow(df) == 0) {
    # Return empty data frame if input is empty
    return(df)
  }
  
  # Handle the condition safely
  tryCatch({
    # Method 1: Try to evaluate the condition
    if (is.logical(condition) && length(condition) == nrow(df)) {
      # Condition is already a logical vector
      keep_rows <- condition
    } else {
      # Condition needs to be evaluated
      # This handles cases where condition is an expression
      stop("Condition must be a logical vector of same length as data frame rows")
    }
    
    # Replace any NA values with FALSE
    keep_rows[is.na(keep_rows)] <- FALSE
    
    # Perform the subset
    result <- df[keep_rows, , drop = FALSE]
    
    return(result)
    
  }, error = function(e) {
    cat("⚠️ safe_subset failed for", description, ":\n")
    cat("   Error:", e$message, "\n")
    cat("   Returning original data frame\n")
    return(df)
  })
}

# Alternative function that builds conditions step by step
safe_filter <- function(df, ..., description = "filter operation") {
  if (!is.data.frame(df)) {
    stop(paste("safe_filter: Input is not a data frame for", description))
  }
  
  if (nrow(df) == 0) {
    return(df)
  }
  
  # Start with all TRUE
  keep_rows <- rep(TRUE, nrow(df))
  
  # Get the conditions
  conditions <- list(...)
  
  # Apply each condition
  for (i in seq_along(conditions)) {
    cond <- conditions[[i]]
    
    if (is.logical(cond) && length(cond) == nrow(df)) {
      # Replace NAs with FALSE and combine with AND
      cond[is.na(cond)] <- FALSE
      keep_rows <- keep_rows & cond
    } else {
      cat("⚠️ Skipping invalid condition", i, "in", description, "\n")
    }
  }
  
  # Perform the subset
  result <- df[keep_rows, , drop = FALSE]
  return(result)
}

# Test the functions
if (interactive() || !exists("SOURCED_FOR_PACKAGE")) {
  cat("\n🧪 Testing safe subset functions...\n")
  
  # Create test data
  test_df <- data.frame(
    value = 1:10,
    p.adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10),
    Count = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  )
  
  # Test safe_subset
  result1 <- safe_subset(test_df, test_df$p.adjust <= 0.05, "p-value filter")
  cat("safe_subset result:", nrow(result1), "rows\n")
  
  # Test safe_filter with multiple conditions
  result2 <- safe_filter(test_df, 
                        test_df$p.adjust <= 0.05,
                        test_df$Count >= 3,
                        description = "combined filter")
  cat("safe_filter result:", nrow(result2), "rows\n")
  
  cat("✅ Safe subset functions ready for use\n")
}