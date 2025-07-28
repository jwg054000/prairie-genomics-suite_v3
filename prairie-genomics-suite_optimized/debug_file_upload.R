# Debug File Upload Issues
# Run this script to test file reading outside of Shiny for debugging

cat("🔍 Prairie Genomics Suite - File Upload Debugging\n")
cat("=" , rep("=", 50), "\n")

# Test file reading with your problematic file
test_file_path <- "MC9.raw.counts.test.csv"  # Replace with your actual file path

if (!file.exists(test_file_path)) {
  cat("❌ Test file not found:", test_file_path, "\n")
  cat("💡 Please update the test_file_path variable with the correct path\n")
  cat("💡 Or place your file in the current directory\n")
  quit()
}

cat("📁 Testing file:", test_file_path, "\n")
cat("📊 File size:", file.size(test_file_path), "bytes\n")

# Try to read first few lines
cat("\n📋 First 3 lines of file:\n")
tryCatch({
  first_lines <- readLines(test_file_path, n = 3)
  for (i in seq_along(first_lines)) {
    cat("Line", i, ":", substr(first_lines[i], 1, 100), "\n")
  }
}, error = function(e) {
  cat("❌ Could not read file lines:", e$message, "\n")
})

# Auto-detect separator
cat("\n🔍 Detecting separator...\n")
first_line <- tryCatch({
  readLines(test_file_path, n = 1)
}, error = function(e) {
  cat("❌ Could not read first line:", e$message, "\n")
  return("")
})

if (first_line != "") {
  separators <- c("," = "COMMA", "\t" = "TAB", ";" = "SEMICOLON")
  detected_seps <- character(0)
  
  for (sep in names(separators)) {
    if (grepl(sep, first_line, fixed = TRUE)) {
      detected_seps <- c(detected_seps, separators[sep])
    }
  }
  
  if (length(detected_seps) > 0) {
    cat("✅ Detected separators:", paste(detected_seps, collapse = ", "), "\n")
  } else {
    cat("⚠️ No common separators detected\n")
  }
}

# Test different reading methods
cat("\n🧪 Testing different file reading methods...\n")

# Method 1: Base R read.csv
cat("\n1️⃣ Testing base R read.csv...\n")
tryCatch({
  data1 <- read.csv(test_file_path, nrows = 5, check.names = FALSE)
  cat("✅ Base R read.csv: SUCCESS -", nrow(data1), "rows x", ncol(data1), "columns\n")
  cat("   Column names:", paste(head(colnames(data1), 3), collapse = ", "), "...\n")
}, error = function(e) {
  cat("❌ Base R read.csv: FAILED -", e$message, "\n")
})

# Method 2: data.table::fread
if (requireNamespace("data.table", quietly = TRUE)) {
  cat("\n2️⃣ Testing data.table::fread...\n")
  tryCatch({
    data2 <- data.table::fread(test_file_path, nrows = 5, data.table = FALSE)
    cat("✅ data.table::fread: SUCCESS -", nrow(data2), "rows x", ncol(data2), "columns\n")
    cat("   Column names:", paste(head(colnames(data2), 3), collapse = ", "), "...\n")
  }, error = function(e) {
    cat("❌ data.table::fread: FAILED -", e$message, "\n")
  })
} else {
  cat("\n2️⃣ data.table package not available - installing...\n")
  tryCatch({
    install.packages("data.table")
    library(data.table)
    data2 <- fread(test_file_path, nrows = 5, data.table = FALSE)
    cat("✅ data.table::fread (after install): SUCCESS -", nrow(data2), "rows x", ncol(data2), "columns\n")
  }, error = function(e) {
    cat("❌ data.table installation/usage: FAILED -", e$message, "\n")
  })
}

# Method 3: readr::read_delim
if (requireNamespace("readr", quietly = TRUE)) {
  cat("\n3️⃣ Testing readr::read_delim...\n")
  tryCatch({
    data3 <- readr::read_delim(test_file_path, n_max = 5, show_col_types = FALSE)
    data3 <- as.data.frame(data3)
    cat("✅ readr::read_delim: SUCCESS -", nrow(data3), "rows x", ncol(data3), "columns\n")
    cat("   Column names:", paste(head(colnames(data3), 3), collapse = ", "), "...\n")
  }, error = function(e) {
    cat("❌ readr::read_delim: FAILED -", e$message, "\n")
  })
} else {
  cat("\n3️⃣ readr package not available\n")
}

# Method 4: read.table with different parameters
cat("\n4️⃣ Testing read.table with various parameters...\n")
tryCatch({
  data4 <- read.table(test_file_path, nrows = 5, header = TRUE, sep = ",", 
                      stringsAsFactors = FALSE, fill = TRUE, quote = "")
  cat("✅ read.table (quote=''): SUCCESS -", nrow(data4), "rows x", ncol(data4), "columns\n")
}, error = function(e) {
  cat("❌ read.table (quote=''): FAILED -", e$message, "\n")
  
  # Try without quotes
  tryCatch({
    data4b <- read.table(test_file_path, nrows = 5, header = TRUE, sep = ",", 
                         stringsAsFactors = FALSE, fill = TRUE, quote = "\"")
    cat("✅ read.table (quote='\"'): SUCCESS -", nrow(data4b), "rows x", ncol(data4b), "columns\n")
  }, error = function(e2) {
    cat("❌ read.table (quote='\"'): FAILED -", e2$message, "\n")
  })
})

cat("\n" , rep("=", 50), "\n")
cat("🎯 DEBUGGING SUMMARY\n")
cat(rep("=", 50), "\n")
cat("✅ If any method succeeded, the file is readable\n")
cat("❌ If all methods failed, there may be a file format issue\n")
cat("💡 Check for:\n")
cat("   - Unusual characters or encoding issues\n")
cat("   - Inconsistent number of columns per row\n")  
cat("   - Mixed data types in columns\n")
cat("   - Special characters in headers\n")
cat("\n📋 Try opening the file in a text editor to inspect formatting\n")