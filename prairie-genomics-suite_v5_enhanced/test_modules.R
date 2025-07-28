# Test script to check if modules load correctly

# Test loading each module individually
cat("Testing module loading...\n")

tryCatch({
  cat("Loading enhanced_sample_annotation.R...\n")
  source("modules/enhanced_sample_annotation.R")
  cat("✅ enhanced_sample_annotation.R loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading enhanced_sample_annotation.R:", e$message, "\n")
})

tryCatch({
  cat("Loading enhanced_deseq2_analysis.R...\n")
  source("modules/enhanced_deseq2_analysis.R")
  cat("✅ enhanced_deseq2_analysis.R loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading enhanced_deseq2_analysis.R:", e$message, "\n")
})

tryCatch({
  cat("Loading context7_visualizations.R...\n")
  source("modules/context7_visualizations.R")
  cat("✅ context7_visualizations.R loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading context7_visualizations.R:", e$message, "\n")
})

tryCatch({
  cat("Loading batch_correction.R...\n")
  source("utils/batch_correction.R")
  cat("✅ batch_correction.R loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading batch_correction.R:", e$message, "\n")
})

cat("Module loading test complete.\n")