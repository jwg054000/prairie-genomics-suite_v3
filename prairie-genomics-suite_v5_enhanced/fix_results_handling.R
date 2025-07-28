# Script to identify and fix the results handling issue

# Search for where we might be returning enrichResult objects directly
cat("Searching for potential issues in pathway_analysis.R...\n\n")

# Read the file
lines <- readLines("pathway_analysis.R")

# Find functions that might return enrichResult objects
function_starts <- grep("^\\w+.*<-\\s*function", lines)
return_statements <- grep("return\\(", lines)

# Check each function
for (func_start in function_starts) {
  # Get function name
  func_name <- gsub("\\s*<-.*", "", lines[func_start])
  func_name <- trimws(func_name)
  
  # Find the end of this function (next function start or end of file)
  func_end <- min(c(function_starts[function_starts > func_start], length(lines)))
  
  # Get return statements in this function
  func_returns <- return_statements[return_statements > func_start & return_statements < func_end]
  
  if (length(func_returns) > 0) {
    cat("\nFunction:", func_name, "\n")
    cat("Returns at lines:", func_returns, "\n")
    
    # Check what's being returned
    for (ret_line in func_returns) {
      return_content <- lines[ret_line]
      if (grepl("return\\(result\\)", return_content) || 
          grepl("return\\(go_result", return_content) ||
          grepl("return\\(kegg_result", return_content) ||
          grepl("return\\(enrichment_result", return_content)) {
        cat("  ⚠️  Line", ret_line, "might return enrichResult object:", trimws(return_content), "\n")
      }
    }
  }
}

# Look for specific problem patterns
cat("\n\nSearching for problematic patterns...\n")

# Pattern 1: Direct access to $success on enrichResult
success_checks <- grep("\\$success", lines)
cat("\nChecking $success usage:\n")
for (line_num in success_checks) {
  line_content <- lines[line_num]
  # Check what variable is being accessed
  var_name <- gsub(".*\\b(\\w+)\\$success.*", "\\1", line_content)
  if (var_name %in% c("results", "result", "go_results", "kegg_results")) {
    cat("  Line", line_num, ":", trimws(line_content), "\n")
    
    # Check context
    if (line_num > 5) {
      # Look back to see where this variable comes from
      for (i in (line_num-5):(line_num-1)) {
        if (grepl(paste0(var_name, "\\s*<-"), lines[i])) {
          cat("    <- Assigned at line", i, ":", trimws(lines[i]), "\n")
        }
      }
    }
  }
}

cat("\n✅ Analysis complete. Key findings:\n")
cat("1. The error occurs when checking results$success on an enrichResult object\n")
cat("2. This happens because individual analysis functions return different formats\n")
cat("3. Need to ensure consistent return format from all pathway analysis functions\n")