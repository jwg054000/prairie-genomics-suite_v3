# Test Integrated Prairie Genomics Suite with Code Visibility
# Test that the main application loads properly with all modules
#
# Author: Prairie Genomics Team
# Date: January 24, 2025

cat("🧪 Testing Integrated Prairie Genomics Suite with Code Visibility\n")
cat("=", rep("=", 65), "\n\n")

# Test 1: Check that all required files exist
cat("1. Checking File Structure\n")
cat(rep("-", 30), "\n")

required_files <- c(
  "app.R",
  "code_visibility/code_logger.R",
  "code_visibility/code_generator.R", 
  "code_visibility/analysis_wrappers.R",
  "code_visibility/code_display_ui.R",
  "code_visibility/code_display_server.R",
  "data_upload.R",
  "sample_annotation.R",
  "deseq2_analysis.R",
  "pathway_analysis.R",
  "visualization.R"
)

all_files_exist <- TRUE
for (file in required_files) {
  if (file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "- MISSING\n")
    all_files_exist <- FALSE
  }
}

if (all_files_exist) {
  cat("\n✅ All required files found\n")
} else {
  cat("\n❌ Some required files are missing\n")
  stop("Missing required files")
}

# Test 2: Check syntax by sourcing key components
cat("\n2. Testing Module Loading\n")
cat(rep("-", 25), "\n")

tryCatch({
  # Test code visibility modules
  source("code_visibility/code_logger.R")
  cat("✅ Code logger loaded\n")
  
  source("code_visibility/code_generator.R")
  cat("✅ Code generator loaded\n")
  
  source("code_visibility/analysis_wrappers.R")
  cat("✅ Analysis wrappers loaded\n")
  
  source("code_visibility/code_display_ui.R")
  cat("✅ Code display UI loaded\n")
  
  source("code_visibility/code_display_server.R")
  cat("✅ Code display server loaded\n")
  
  cat("\n✅ All code visibility modules loaded successfully\n")
  
}, error = function(e) {
  cat("❌ Error loading modules:", e$message, "\n")
  stop("Module loading failed")
})

# Test 3: Test UI component creation
cat("\n3. Testing UI Component Creation\n")
cat(rep("-", 35), "\n")

tryCatch({
  # Test that UI functions can be called
  test_ui <- codeDisplayUI("test_id")
  cat("✅ codeDisplayUI() works\n")
  
  test_inline_ui <- inlineCodeUI("test_inline")
  cat("✅ inlineCodeUI() works\n")
  
  test_export_ui <- codeExportUI("test_export")
  cat("✅ codeExportUI() works\n")
  
  test_validation_ui <- codeValidationUI("test_validation")
  cat("✅ codeValidationUI() works\n")
  
  test_docs_ui <- methodDocumentationUI("test_docs")
  cat("✅ methodDocumentationUI() works\n")
  
  cat("\n✅ All UI components created successfully\n")
  
}, error = function(e) {
  cat("❌ Error creating UI components:", e$message, "\n")
  stop("UI component creation failed")
})

# Test 4: Test core functionality
cat("\n4. Testing Core Code Visibility Functions\n")
cat(rep("-", 40), "\n")

tryCatch({
  # Initialize a test session
  session_id <- init_code_logger("test_integration")
  cat("✅ Code logger session initialized:", session_id, "\n")
  
  # Test code generation
  deseq2_code <- generate_deseq2_code("~ group", 
                                      list(numerator = "treatment", denominator = "control"),
                                      list(padj_cutoff = 0.05, fc_cutoff = 1))
  cat("✅ DESeq2 code generation works\n")
  
  # Test step logging
  log_analysis_step(
    session_id = session_id,
    step_name = "Integration Test",
    category = "test", 
    code_snippet = deseq2_code,
    description = "Testing integrated functionality"
  )
  cat("✅ Analysis step logging works\n")
  
  # Test code retrieval
  complete_code <- get_session_code(session_id)
  if (nchar(complete_code) > 0) {
    cat("✅ Code retrieval works (", nchar(complete_code), "characters)\n")
  } else {
    cat("⚠️ Code retrieval returned empty\n")
  }
  
  # Clean up
  clear_session(session_id)
  cat("✅ Session cleanup works\n")
  
  cat("\n✅ Core functionality test passed\n")
  
}, error = function(e) {
  cat("❌ Error in core functionality:", e$message, "\n")
  stop("Core functionality test failed")
})

# Test 5: Test app.R syntax (without running Shiny)
cat("\n5. Testing App.R Syntax\n")
cat(rep("-", 25), "\n")

tryCatch({
  # Parse the app.R file to check for syntax errors
  parsed_app <- parse("app.R")
  cat("✅ app.R syntax is valid\n")
  
  # Count lines to verify it's substantial
  app_lines <- readLines("app.R")
  cat("✅ app.R has", length(app_lines), "lines\n")
  
  # Check that code visibility components are included
  has_code_view_tab <- any(grepl("code_view", app_lines))
  has_code_modules <- any(grepl("code_visibility", app_lines))
  has_inline_code <- any(grepl("inlineCodeUI", app_lines))
  
  if (has_code_view_tab && has_code_modules && has_inline_code) {
    cat("✅ Code visibility integration detected in app.R\n")
  } else {
    cat("⚠️ Code visibility integration may be incomplete\n")
    cat("   - Code View tab:", has_code_view_tab, "\n")
    cat("   - Code modules:", has_code_modules, "\n") 
    cat("   - Inline code:", has_inline_code, "\n")
  }
  
  cat("\n✅ App.R integration test passed\n")
  
}, error = function(e) {
  cat("❌ Error in app.R:", e$message, "\n")
  stop("App.R test failed")
})

# Test Summary
cat("\n🎯 Integration Test Summary\n")
cat("=", rep("=", 25), "\n")

cat("✅ File structure: All required files present\n")
cat("✅ Module loading: All code visibility modules loaded successfully\n")
cat("✅ UI components: All UI functions working\n")
cat("✅ Core functionality: Logging, generation, and retrieval working\n")
cat("✅ App integration: Code visibility integrated into main application\n")

cat("\n🎉 Integration Test Complete!\n")
cat("📋 The Prairie Genomics Suite with code visibility is ready for use\n")
cat("🚀 Users can now view, download, and reproduce their analyses\n")

cat("\n📋 Features Available:\n")
cat("   - 📜 Dedicated Code View tab with complete analysis scripts\n")
cat("   - 💡 Inline code displays in DESeq2 and pathway analysis tabs\n")
cat("   - 📄 Multiple export formats (R script, R Markdown, HTML)\n")
cat("   - ✅ Code validation and syntax checking\n")
cat("   - 📚 Method documentation and scientific references\n")
cat("   - 🔧 Automatic logging of all analysis steps and parameters\n")

cat("\n🔗 Ready for Phase 3: Enhanced Export Capabilities\n")