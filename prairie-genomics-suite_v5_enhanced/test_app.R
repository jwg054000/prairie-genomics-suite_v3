# Quick test of R Shiny Prairie Genomics Suite
# Tests basic functionality and package availability

cat("🧬 Testing R Shiny Prairie Genomics Suite\n")
cat("=========================================\n")

# Test key packages
packages_to_test <- c(
  "shiny", "shinydashboard", "DT", 
  "ggplot2", "plotly", "dplyr"
)

all_present <- TRUE
for (pkg in packages_to_test) {
  cat(paste0("Testing ", pkg, "... "))
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✅\n")
  } else {
    cat("❌\n") 
    all_present <- FALSE
  }
}

if (!all_present) {
  cat("\n⚠️  Some packages missing. Run install.R first.\n")
  quit(status = 1)
}

cat("\n📱 Testing basic Shiny UI creation...\n")

# Test basic UI creation
library(shiny)
library(shinydashboard)

# Create minimal test UI
test_ui <- dashboardPage(
  dashboardHeader(title = "Test"),
  dashboardSidebar(),
  dashboardBody(
    h1("Prairie Genomics Suite Test"),
    p("If you see this, the basic UI works!")
  )
)

cat("✅ UI creation successful\n")

cat("\n🔬 Testing data processing functions...\n")

# Test basic data operations
library(dplyr)
test_data <- data.frame(
  gene = paste0("Gene_", 1:10),
  sample1 = rnorm(10, 100, 20),
  sample2 = rnorm(10, 120, 25),
  stringsAsFactors = FALSE
)

processed_data <- test_data %>%
  mutate(mean_expr = (sample1 + sample2) / 2) %>%
  arrange(desc(mean_expr))

cat("✅ Data processing successful\n")

cat("\n📊 Testing visualization packages...\n")

library(ggplot2)
library(plotly)

# Create test plot
test_plot <- ggplot(test_data, aes(x = sample1, y = sample2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Test Scatter Plot")

cat("✅ ggplot2 plot creation successful\n")

# Test plotly conversion
test_plotly <- ggplotly(test_plot)
cat("✅ plotly conversion successful\n")

cat("\n🎉 All tests passed! R Shiny environment is ready.\n")
cat("\n🚀 To run the full application:\n")
cat("   R -e \"shiny::runApp('app.R')\"\n")
cat("\n📖 Or in RStudio:\n")  
cat("   1. Open app.R\n")
cat("   2. Click 'Run App' button\n")

cat("\n📋 System Info:\n")
cat(paste0("R version: ", R.version.string, "\n"))
cat(paste0("Platform: ", R.version$platform, "\n"))

# Package versions
cat("\n📦 Package Versions:\n")
for (pkg in packages_to_test) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(paste0("   ", pkg, ": ", version, "\n"))
  }
}