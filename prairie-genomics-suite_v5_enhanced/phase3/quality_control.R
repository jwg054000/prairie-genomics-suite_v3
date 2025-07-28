# Phase 3 - Automated Quality Control System
# Real-time data quality assessment with educational feedback
# 
# This module provides:
# - Automated quality control checks
# - Traffic light system (green/yellow/red)
# - Educational explanations for each metric
# - Interactive quality visualizations
# - Pre-flight analysis validation

library(shiny)
library(ggplot2)
library(plotly)
library(pheatmap)
library(RColorBrewer)

# ===========================================
# QUALITY CONTROL CHECKS CATALOG
# ===========================================

# Define comprehensive QC checks with educational context
QC_CHECKS <- list(
  
  # Sample quality checks
  sample_quality = list(
    library_size = list(
      name = "Library Size Distribution",
      description = "Checks if samples have similar total read counts",
      icon = "📚",
      thresholds = list(good = 2, warning = 5, fail = 10),  # CV thresholds
      explanation = "Large differences in library size can indicate technical problems",
      what_it_means = "Similar library sizes suggest consistent sequencing depth",
      fix_if_bad = "Consider normalizing or removing low-quality samples"
    ),
    
    feature_detection = list(
      name = "Gene Detection Rate", 
      description = "Percentage of genes detected in each sample",
      icon = "🔍",
      thresholds = list(good = 70, warning = 50, fail = 30),  # % genes detected
      explanation = "Low detection rates may indicate poor RNA quality or sequencing",
      what_it_means = "More detected genes = better RNA quality and sequencing depth",
      fix_if_bad = "Check RNA integrity scores and sequencing metrics"
    ),
    
    outlier_detection = list(
      name = "Sample Outliers",
      description = "Identifies samples that don't cluster with their group",
      icon = "🎯",
      thresholds = list(good = 0, warning = 1, fail = 2),  # Number of outliers
      explanation = "Outliers may indicate mislabeled samples or technical failures",
      what_it_means = "Samples should cluster with their biological replicates",
      fix_if_bad = "Verify sample labels and consider removing true outliers"
    )
  ),
  
  # Experimental design checks
  design_quality = list(
    replicate_adequacy = list(
      name = "Replicate Count",
      description = "Sufficient biological replicates for statistical power",
      icon = "🔢",
      thresholds = list(good = 4, warning = 3, fail = 2),  # Min replicates per condition
      explanation = "More replicates = better statistical power and reliability",
      what_it_means = "Need enough replicates to distinguish signal from noise",
      fix_if_bad = "Add more biological replicates if possible"
    ),
    
    batch_effects = list(
      name = "Batch Effect Detection",
      description = "Checks if technical batches confound biological conditions",
      icon = "🏭",
      thresholds = list(good = 0.1, warning = 0.3, fail = 0.5),  # Proportion variance explained
      explanation = "Batch effects can create false discoveries or mask real biology",
      what_it_means = "Processing date/person shouldn't explain more variance than biology",
      fix_if_bad = "Include batch as covariate or balance batches across conditions"
    ),
    
    condition_balance = list(
      name = "Condition Balance",
      description = "Similar number of samples across conditions",
      icon = "⚖️",
      thresholds = list(good = 1.2, warning = 2.0, fail = 3.0),  # Max ratio between groups
      explanation = "Balanced designs have better statistical properties",
      what_it_means = "Each condition should have similar sample sizes",
      fix_if_bad = "Add samples to smaller groups or use appropriate statistical methods"
    )
  ),
  
  # Data quality checks
  data_quality = list(
    count_distribution = list(
      name = "Count Distribution",
      description = "Raw count data follows expected patterns",
      icon = "📊",
      thresholds = list(good = 0.8, warning = 0.6, fail = 0.4),  # Correlation with expected
      explanation = "RNA-seq counts should follow negative binomial distribution",
      what_it_means = "Proper count distribution is essential for DESeq2 assumptions",
      fix_if_bad = "Check for count processing errors or contamination"
    ),
    
    zero_inflation = list(
      name = "Zero Count Analysis",
      description = "Reasonable proportion of zero counts per gene",
      icon = "🕳️",
      thresholds = list(good = 30, warning = 50, fail = 70),  # % genes with >50% zeros
      explanation = "Too many zeros can indicate poor sequencing depth or dropouts",
      what_it_means = "Some zeros are normal, but excess suggests technical issues",
      fix_if_bad = "Increase sequencing depth or remove poorly expressed genes"
    ),
    
    normalization_factors = list(
      name = "Size Factor Variation",
      description = "DESeq2 normalization factors within reasonable range",
      icon = "📏",
      thresholds = list(good = 2, warning = 5, fail = 10),  # Max fold difference
      explanation = "Large size factors suggest very different library compositions",
      what_it_means = "Size factors correct for library size and composition",
      fix_if_bad = "Investigate samples with extreme size factors"
    )
  ),
  
  # Biological validation checks
  biological_validation = list(
    housekeeping_stability = list(
      name = "Housekeeping Gene Stability",
      description = "Reference genes show consistent expression",
      icon = "🏠",
      thresholds = list(good = 0.2, warning = 0.5, fail = 1.0),  # CV of housekeeping genes
      explanation = "Housekeeping genes should be stable across conditions",
      what_it_means = "Stable reference genes validate your experimental system",
      fix_if_bad = "Check experimental conditions or choose different reference genes"
    ),
    
    positive_controls = list(
      name = "Positive Control Validation",
      description = "Known markers behave as expected",
      icon = "✅",
      thresholds = list(good = 80, warning = 60, fail = 40),  # % controls behaving correctly
      explanation = "Positive controls validate your experimental manipulation",
      what_it_means = "Known responders should show expected changes",
      fix_if_bad = "Verify experimental conditions and control gene selection"
    ),
    
    pathway_coherence = list(
      name = "Pathway Coherence",
      description = "Related genes show coordinated expression changes",
      icon = "🛤️",
      thresholds = list(good = 0.7, warning = 0.5, fail = 0.3),  # Pathway correlation score
      explanation = "Genes in same pathways often change together",
      what_it_means = "Coherent pathways suggest biologically meaningful results",
      fix_if_bad = "Consider experimental artifacts or batch effects"
    )
  )
)

# ===========================================
# QUALITY CONTROL ENGINE
# ===========================================

# Run comprehensive quality control analysis
run_qc_analysis <- function(expression_data, sample_metadata, dds_object = NULL) {
  
  qc_results <- list()
  
  # Sample quality checks
  qc_results$sample_quality <- check_sample_quality(expression_data, sample_metadata)
  
  # Design quality checks
  qc_results$design_quality <- check_design_quality(sample_metadata)
  
  # Data quality checks (requires DESeq2 object)
  if (!is.null(dds_object)) {
    qc_results$data_quality <- check_data_quality(dds_object)
    qc_results$biological_validation <- check_biological_validation(dds_object, sample_metadata)
  }
  
  # Calculate overall quality score
  qc_results$overall_score <- calculate_overall_qc_score(qc_results)
  
  # Generate recommendations
  qc_results$recommendations <- generate_qc_recommendations(qc_results)
  
  return(qc_results)
}

# Check sample quality metrics
check_sample_quality <- function(expression_data, sample_metadata) {
  
  results <- list()
  
  # Library size distribution
  library_sizes <- colSums(expression_data)
  cv_library_size <- sd(library_sizes) / mean(library_sizes) * 100
  
  results$library_size <- list(
    value = cv_library_size,
    status = classify_qc_status(cv_library_size, QC_CHECKS$sample_quality$library_size$thresholds, "lower_better"),
    raw_data = library_sizes,
    interpretation = interpret_library_sizes(cv_library_size)
  )
  
  # Gene detection rate
  detection_rates <- apply(expression_data > 0, 2, mean) * 100
  min_detection <- min(detection_rates)
  
  results$feature_detection <- list(
    value = min_detection,
    status = classify_qc_status(min_detection, QC_CHECKS$sample_quality$feature_detection$thresholds, "higher_better"),
    raw_data = detection_rates,
    interpretation = interpret_detection_rates(min_detection)
  )
  
  # Outlier detection using PCA
  outlier_analysis <- detect_sample_outliers(expression_data, sample_metadata)
  
  results$outlier_detection <- list(
    value = outlier_analysis$n_outliers,
    status = classify_qc_status(outlier_analysis$n_outliers, QC_CHECKS$sample_quality$outlier_detection$thresholds, "lower_better"),
    raw_data = outlier_analysis,
    interpretation = interpret_outliers(outlier_analysis$n_outliers)
  )
  
  return(results)
}

# Check experimental design quality
check_design_quality <- function(sample_metadata) {
  
  results <- list()
  
  # Replicate adequacy
  condition_counts <- table(sample_metadata$Condition)
  min_replicates <- min(condition_counts)
  
  results$replicate_adequacy <- list(
    value = min_replicates,
    status = classify_qc_status(min_replicates, QC_CHECKS$design_quality$replicate_adequacy$thresholds, "higher_better"),
    raw_data = condition_counts,
    interpretation = interpret_replicates(min_replicates)
  )
  
  # Condition balance
  max_ratio <- max(condition_counts) / min(condition_counts)
  
  results$condition_balance <- list(
    value = max_ratio,
    status = classify_qc_status(max_ratio, QC_CHECKS$design_quality$condition_balance$thresholds, "lower_better"),
    raw_data = condition_counts,
    interpretation = interpret_balance(max_ratio)
  )
  
  # Batch effect detection (if batch info available)
  if ("Batch" %in% colnames(sample_metadata)) {
    batch_analysis <- analyze_batch_effects(sample_metadata)
    results$batch_effects <- list(
      value = batch_analysis$confounding_score,
      status = classify_qc_status(batch_analysis$confounding_score, QC_CHECKS$design_quality$batch_effects$thresholds, "lower_better"),
      raw_data = batch_analysis,
      interpretation = interpret_batch_effects(batch_analysis$confounding_score)
    )
  }
  
  return(results)
}

# ===========================================
# QC STATUS CLASSIFICATION
# ===========================================

# Classify QC metric into traffic light system
classify_qc_status <- function(value, thresholds, direction = "higher_better") {
  
  if (direction == "higher_better") {
    if (value >= thresholds$good) {
      return("good")
    } else if (value >= thresholds$warning) {
      return("warning")
    } else {
      return("fail")
    }
  } else {  # lower_better
    if (value <= thresholds$good) {
      return("good")
    } else if (value <= thresholds$warning) {
      return("warning")  
    } else {
      return("fail")
    }
  }
}

# ===========================================
# QC DASHBOARD UI
# ===========================================

# Create quality control dashboard UI
create_qc_dashboard_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # QC Header
    div(
      class = "qc-header",
      style = "background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
               color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
      
      div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        
        div(
          h2("🚦 Quality Control Dashboard", style = "margin: 0;"),
          p("Real-time validation of your RNA-seq data", style = "margin: 5px 0 0 0; opacity: 0.9;")
        ),
        
        div(
          id = ns("overall_status"),
          class = "overall-qc-status",
          style = "text-align: center;",
          uiOutput(ns("overall_status_display"))
        )
      )
    ),
    
    # QC Categories Tabs
    tabsetPanel(
      id = ns("qc_tabs"),
      type = "pills",
      
      # Sample Quality Tab
      tabPanel(
        "📚 Sample Quality",
        value = "sample_quality",
        
        div(style = "padding: 20px 0;",
          fluidRow(
            column(12,
              uiOutput(ns("sample_quality_cards"))
            )
          ),
          
          fluidRow(
            column(6,
              div(
                class = "qc-plot-container",
                style = "background: white; padding: 20px; border-radius: 12px; 
                         border: 1px solid #e5e7eb;",
                h4("📊 Library Size Distribution"),
                plotlyOutput(ns("library_size_plot"), height = "300px")
              )
            ),
            
            column(6,
              div(
                class = "qc-plot-container", 
                style = "background: white; padding: 20px; border-radius: 12px; 
                         border: 1px solid #e5e7eb;",
                h4("🔍 Gene Detection Rates"),
                plotlyOutput(ns("detection_rate_plot"), height = "300px")
              )
            )
          )
        )
      ),
      
      # Design Quality Tab
      tabPanel(
        "⚖️ Design Quality",
        value = "design_quality",
        
        div(style = "padding: 20px 0;",
          uiOutput(ns("design_quality_cards")),
          
          fluidRow(
            column(12,
              div(
                class = "qc-plot-container",
                style = "background: white; padding: 20px; border-radius: 12px; 
                         border: 1px solid #e5e7eb; margin-top: 20px;",
                h4("👥 Sample Groups Overview"),
                plotlyOutput(ns("sample_groups_plot"), height = "400px")
              )
            )
          )
        )
      ),
      
      # Data Quality Tab
      tabPanel(
        "📊 Data Quality",
        value = "data_quality",
        
        div(style = "padding: 20px 0;",
          uiOutput(ns("data_quality_cards")),
          
          conditionalPanel(
            condition = "output.data_available == true",
            ns = ns,
            
            fluidRow(
              column(6,
                div(
                  class = "qc-plot-container",
                  style = "background: white; padding: 20px; border-radius: 12px; 
                           border: 1px solid #e5e7eb;",
                  h4("📈 Count Distribution"),
                  plotlyOutput(ns("count_distribution_plot"), height = "300px")
                )
              ),
              
              column(6,
                div(
                  class = "qc-plot-container",
                  style = "background: white; padding: 20px; border-radius: 12px; 
                           border: 1px solid #e5e7eb;",
                  h4("📏 Normalization Factors"),
                  plotlyOutput(ns("size_factors_plot"), height = "300px")
                )
              )
            )
          )
        )
      ),
      
      # Recommendations Tab
      tabPanel(
        "💡 Recommendations",
        value = "recommendations",
        
        div(style = "padding: 20px 0;",
          uiOutput(ns("recommendations_display"))
        )
      )
    )
  )
}

# ===========================================
# QC DASHBOARD SERVER
# ===========================================

# Quality control dashboard server
qc_dashboard_server <- function(id, expression_data, sample_metadata, dds_object = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Run QC analysis reactively
    qc_results <- reactive({
      req(expression_data(), sample_metadata())
      
      withProgress(message = "Running quality control...", {
        run_qc_analysis(expression_data(), sample_metadata(), dds_object())
      })
    })
    
    # Overall status display
    output$overall_status_display <- renderUI({
      req(qc_results())
      
      overall <- qc_results()$overall_score
      
      status_color <- switch(overall$status,
                           good = "#10b981",
                           warning = "#f59e0b", 
                           fail = "#ef4444")
      
      status_icon <- switch(overall$status,
                          good = "✅",
                          warning = "⚠️",
                          fail = "❌")
      
      div(
        div(
          style = paste0("font-size: 48px; color: ", status_color, ";"),
          status_icon
        ),
        div(
          style = "font-weight: 600; margin-top: 5px;",
          paste("Overall:", toupper(overall$status))
        ),
        div(
          style = "font-size: 14px; opacity: 0.9; margin-top: 2px;",
          paste("Score:", overall$score, "/100")
        )
      )
    })
    
    # Sample quality cards
    output$sample_quality_cards <- renderUI({
      req(qc_results())
      
      sample_qc <- qc_results()$sample_quality
      
      div(
        class = "qc-cards-container",
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
                 gap: 20px; margin-bottom: 20px;",
        
        lapply(names(sample_qc), function(check_name) {
          create_qc_card(check_name, sample_qc[[check_name]], "sample_quality")
        })
      )
    })
    
    # Library size plot
    output$library_size_plot <- renderPlotly({
      req(qc_results())
      
      lib_sizes <- qc_results()$sample_quality$library_size$raw_data
      
      p <- ggplot(data.frame(Sample = names(lib_sizes), Size = lib_sizes), 
                  aes(x = Sample, y = Size)) +
        geom_col(fill = "#3b82f6", alpha = 0.7) +
        labs(title = "Library Sizes", x = "Sample", y = "Total Counts") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # Return QC results for other modules
    return(qc_results)
  })
}

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

# Create individual QC status card
create_qc_card <- function(check_name, check_result, category) {
  
  check_info <- QC_CHECKS[[category]][[check_name]]
  
  status_colors <- list(
    good = list(bg = "#f0fdf4", border = "#10b981", text = "#065f46"),
    warning = list(bg = "#fffbeb", border = "#f59e0b", text = "#92400e"),
    fail = list(bg = "#fef2f2", border = "#ef4444", text = "#991b1b")
  )
  
  colors <- status_colors[[check_result$status]]
  
  div(
    class = "qc-status-card",
    style = paste0(
      "background: ", colors$bg, "; ",
      "border: 2px solid ", colors$border, "; ",
      "border-radius: 12px; padding: 20px; ",
      "transition: all 0.3s ease;"
    ),
    
    # Header
    div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      
      div(
        style = "font-size: 24px; margin-right: 10px;",
        check_info$icon
      ),
      
      div(
        style = "flex: 1;",
        h5(check_info$name, style = paste0("margin: 0; color: ", colors$text, ";")),
        p(check_info$description, style = "margin: 0; font-size: 14px; opacity: 0.8;")
      ),
      
      div(
        style = paste0("font-weight: bold; font-size: 18px; color: ", colors$text, ";"),
        toupper(check_result$status)
      )
    ),
    
    # Value display
    div(
      style = "margin-bottom: 15px;",
      div(
        style = paste0("font-size: 24px; font-weight: bold; color: ", colors$text, ";"),
        format_qc_value(check_result$value, check_name)
      )
    ),
    
    # Interpretation
    div(
      style = "font-size: 14px; line-height: 1.4;",
      check_result$interpretation
    ),
    
    # Learn more button
    div(
      style = "margin-top: 15px;",
      tags$button(
        class = "btn btn-sm btn-outline-info",
        style = "font-size: 12px;",
        "Learn More"
      )
    )
  )
}

# Format QC values for display
format_qc_value <- function(value, check_name) {
  
  formatters <- list(
    library_size = function(x) paste0(round(x, 1), "% CV"),
    feature_detection = function(x) paste0(round(x, 1), "%"),
    outlier_detection = function(x) paste(x, "outliers"),
    replicate_adequacy = function(x) paste(x, "replicates"),
    condition_balance = function(x) paste0(round(x, 1), ":1 ratio")
  )
  
  formatter <- formatters[[check_name]]
  if (!is.null(formatter)) {
    return(formatter(value))
  } else {
    return(as.character(round(value, 2)))
  }
}

# Export functions
list(
  create_qc_dashboard_ui = create_qc_dashboard_ui,
  qc_dashboard_server = qc_dashboard_server,
  run_qc_analysis = run_qc_analysis,
  QC_CHECKS = QC_CHECKS
)