# Phase 3 - Smart Error Prevention System
# Real-time error interception and educational guidance
# 
# This module provides:
# - Proactive mistake detection and prevention
# - Educational error messages with solutions
# - Common pitfall interception 
# - Statistical assumption validation
# - Real-time parameter safety checks

library(shiny)
library(ggplot2)
library(plotly)
library(DT)

# ===========================================
# ERROR PREVENTION CATALOG
# ===========================================

# Define comprehensive error prevention rules with educational context
ERROR_PREVENTION_RULES <- list(
  
  # Statistical errors
  statistical_errors = list(
    
    insufficient_replicates = list(
      name = "Insufficient Biological Replicates",
      detect_function = "check_sample_size",
      severity = "critical",
      icon = "⚠️",
      trigger_condition = "min_samples_per_group < 3",
      user_message = "DESeq2 requires at least 3 biological replicates per condition for reliable statistical inference.",
      explanation = "With fewer than 3 replicates, you cannot distinguish biological variation from technical noise. This makes your results unreliable and non-reproducible.",
      consequences = list(
        "Inflated false discovery rates",
        "Unable to estimate biological variation properly", 
        "Results won't replicate in independent experiments",
        "Reviewers will reject your manuscript"
      ),
      solutions = list(
        primary = "Add more biological replicates (4-6 recommended)",
        alternative = "Treat as exploratory analysis only",
        workaround = "Use more conservative statistical thresholds"
      ),
      educational_links = list(
        "Why replicates matter" = "https://example.com/replicates",
        "Sample size calculation" = "https://example.com/power"
      ),
      prevent_action = TRUE,
      allow_override = FALSE
    ),
    
    confounded_design = list(
      name = "Confounded Experimental Design",
      detect_function = "check_design_confounding",
      severity = "critical", 
      icon = "🔄",
      trigger_condition = "batch_treatment_confounded == TRUE",
      user_message = "Your experimental design has batch effects confounded with treatment groups.",
      explanation = "When treatment groups are processed in different batches, you cannot distinguish treatment effects from batch effects. This is a fundamental design flaw that cannot be fixed computationally.",
      consequences = list(
        "Cannot distinguish treatment from batch effects",
        "Results may be entirely due to processing differences", 
        "False discoveries or missed true positives",
        "Experiment results will not be reproducible"
      ),
      solutions = list(
        primary = "Redesign experiment with balanced batches",
        alternative = "Include batch as covariate (limited effectiveness)",
        workaround = "Acknowledge limitation in interpretation"
      ),
      educational_links = list(
        "Experimental design principles" = "https://example.com/design",
        "Batch effects explained" = "https://example.com/batch"
      ),
      prevent_action = TRUE,
      allow_override = FALSE
    ),
    
    inappropriate_multiple_testing = list(
      name = "Inappropriate Multiple Testing Correction",
      detect_function = "check_multiple_testing",
      severity = "high",
      icon = "📊",
      trigger_condition = "too_many_comparisons OR wrong_correction_method",
      user_message = "Your analysis approach will lead to excessive false discoveries due to multiple testing issues.",
      explanation = "When testing thousands of genes simultaneously, some will appear significant by chance alone. Without proper correction, most of your 'discoveries' will be false positives.",
      consequences = list(
        "95% of significant results may be false positives",
        "Resources wasted on false leads",
        "Failed replication attempts",
        "Publication retraction risk"
      ),
      solutions = list(
        primary = "Use appropriate FDR correction (Benjamini-Hochberg)",
        alternative = "Focus on pathway-level analysis",
        workaround = "Use more stringent significance thresholds"
      ),
      educational_links = list(
        "Multiple testing explained" = "https://example.com/multiple",
        "FDR vs FWER" = "https://example.com/fdr"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    ),
    
    violation_of_assumptions = list(
      name = "Statistical Assumptions Violated",
      detect_function = "check_deseq2_assumptions", 
      severity = "high",
      icon = "📐",
      trigger_condition = "assumptions_violated == TRUE",
      user_message = "Your data violates key statistical assumptions required for DESeq2.",
      explanation = "DESeq2 assumes counts follow a negative binomial distribution and that dispersion can be estimated reliably. Violations lead to incorrect p-values and effect sizes.",
      consequences = list(
        "Incorrect statistical inference",
        "Unreliable p-values and fold changes",
        "Poor reproducibility across experiments",
        "Invalid biological conclusions"
      ),
      solutions = list(
        primary = "Use alternative analysis method appropriate for your data",
        alternative = "Transform data to meet assumptions",
        workaround = "Acknowledge limitations and use conservative interpretation"
      ),
      educational_links = list(
        "DESeq2 assumptions" = "https://example.com/assumptions",
        "Alternative methods" = "https://example.com/alternatives"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    )
  ),
  
  # Data quality errors
  data_quality_errors = list(
    
    outlier_samples = list(
      name = "Extreme Outlier Samples Detected",
      detect_function = "check_sample_outliers",
      severity = "medium",
      icon = "🎯",
      trigger_condition = "extreme_outliers_present == TRUE",
      user_message = "Some samples are extreme outliers that may skew your analysis results.",
      explanation = "Outlier samples with very different expression patterns can dominate statistical tests and lead to false discoveries or mask true biological signals.",
      consequences = list(
        "Skewed statistical results",
        "False positive discoveries",
        "Reduced power to detect true effects",
        "Misrepresentation of biological reality"
      ),
      solutions = list(
        primary = "Investigate outliers - check sample processing notes",
        alternative = "Remove confirmed technical outliers",
        workaround = "Use robust statistical methods"
      ),
      educational_links = list(
        "Outlier detection" = "https://example.com/outliers",
        "Sample QC" = "https://example.com/sampleqc"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    ),
    
    low_sequencing_depth = list(
      name = "Insufficient Sequencing Depth",
      detect_function = "check_sequencing_depth",
      severity = "medium", 
      icon = "📊",
      trigger_condition = "median_lib_size < 1000000",
      user_message = "Your samples have very low read counts, which limits statistical power.",
      explanation = "Low sequencing depth means many genes will have zero or very low counts, reducing your ability to detect differential expression and increasing uncertainty in fold change estimates.",
      consequences = list(
        "Reduced statistical power",
        "Increased uncertainty in results",
        "Many genes undetectable",
        "Poor reproducibility"
      ),
      solutions = list(
        primary = "Increase sequencing depth for better results",
        alternative = "Focus analysis on highly expressed genes",
        workaround = "Use appropriate statistical methods for low counts"
      ),
      educational_links = list(
        "Sequencing depth guidelines" = "https://example.com/depth",
        "Low count methods" = "https://example.com/lowcounts"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    ),
    
    contamination_detected = list(
      name = "Potential Sample Contamination",
      detect_function = "check_contamination",
      severity = "high",
      icon = "🦠",
      trigger_condition = "contamination_markers_present == TRUE",
      user_message = "Your samples may be contaminated with DNA from other organisms.",
      explanation = "Contamination can introduce spurious signals and confound your biological interpretation. It's particularly problematic in single-cell or low-input RNA-seq experiments.",
      consequences = list(
        "False biological signals",
        "Misidentification of cell types",
        "Incorrect pathway enrichment",
        "Non-reproducible results"
      ),
      solutions = list(
        primary = "Investigate contamination source and re-process samples",
        alternative = "Filter out contaminating sequences", 
        workaround = "Acknowledge contamination in interpretation"
      ),
      educational_links = list(
        "Contamination detection" = "https://example.com/contamination",
        "Quality control" = "https://example.com/qc"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    )
  ),
  
  # Interpretation errors
  interpretation_errors = list(
    
    fold_change_misinterpretation = list(
      name = "Fold Change Threshold Too Low",
      detect_function = "check_fold_change_threshold",
      severity = "medium",
      icon = "📈",
      trigger_condition = "fc_cutoff < 1.5",
      user_message = "Your fold change threshold may be too low to represent biological significance.",
      explanation = "Very small fold changes (e.g., <1.5x) may be statistically significant but biologically meaningless due to technical noise and biological variation.",
      consequences = list(
        "Many false positive 'discoveries'",
        "Results don't replicate in validation",
        "Wasted resources on unimportant changes",
        "Difficulty distinguishing signal from noise"
      ),
      solutions = list(
        primary = "Use fold change threshold ≥1.5x for most experiments",
        alternative = "Focus on highly significant changes first",
        workaround = "Validate all small changes experimentally"
      ),
      educational_links = list(
        "Effect size interpretation" = "https://example.com/effectsize",
        "Biological significance" = "https://example.com/biological"
      ),
      prevent_action = FALSE,
      allow_override = TRUE
    ),
    
    pvalue_fishing = list(
      name = "P-value Threshold Manipulation",
      detect_function = "check_pvalue_fishing",
      severity = "high",
      icon = "🎣",
      trigger_condition = "multiple_pvalue_thresholds_tested == TRUE",
      user_message = "Testing multiple p-value thresholds to get desired results is statistically invalid.",
      explanation = "Trying different significance thresholds until you get the 'right' number of significant genes invalidates your statistical inference and inflates false discovery rates.",
      consequences = list(
        "Inflated false discovery rates",
        "Non-reproducible results", 
        "Statistical malpractice",
        "Manuscript rejection"
      ),
      solutions = list(
        primary = "Pre-specify significance threshold before analysis",
        alternative = "Report results at standard threshold (0.05)",
        workaround = "Acknowledge exploratory nature of threshold selection"
      ),
      educational_links = list(
        "P-hacking explained" = "https://example.com/phacking",
        "Pre-registration" = "https://example.com/prereg"
      ),
      prevent_action = TRUE,
      allow_override = FALSE
    )
  )
)

# ===========================================
# ERROR DETECTION ENGINE
# ===========================================

# Main error detection function
detect_analysis_errors <- function(expression_data, sample_metadata, analysis_parameters) {
  
  detected_errors <- list()
  
  # Check statistical errors
  statistical_checks <- check_statistical_errors(expression_data, sample_metadata, analysis_parameters)
  if (length(statistical_checks) > 0) {
    detected_errors$statistical <- statistical_checks
  }
  
  # Check data quality errors
  quality_checks <- check_data_quality_errors(expression_data, sample_metadata)
  if (length(quality_checks) > 0) {
    detected_errors$data_quality <- quality_checks
  }
  
  # Check interpretation errors
  interpretation_checks <- check_interpretation_errors(analysis_parameters)
  if (length(interpretation_checks) > 0) {
    detected_errors$interpretation <- interpretation_checks
  }
  
  # Prioritize errors by severity
  detected_errors$priority_order <- prioritize_errors(detected_errors)
  
  return(detected_errors)
}

# Check for statistical errors
check_statistical_errors <- function(expression_data, sample_metadata, parameters) {
  
  errors <- list()
  
  # Check sample size adequacy
  condition_counts <- table(sample_metadata$Condition)
  min_samples <- min(condition_counts)
  
  if (min_samples < 3) {
    errors$insufficient_replicates <- list(
      rule = ERROR_PREVENTION_RULES$statistical_errors$insufficient_replicates,
      detected_values = list(
        min_samples_per_group = min_samples,
        condition_counts = condition_counts
      ),
      severity_level = "critical",
      timestamp = Sys.time()
    )
  }
  
  # Check for design confounding
  if ("Batch" %in% colnames(sample_metadata)) {
    confounding_analysis <- analyze_confounding(sample_metadata)
    if (confounding_analysis$is_confounded) {
      errors$confounded_design <- list(
        rule = ERROR_PREVENTION_RULES$statistical_errors$confounded_design,
        detected_values = confounding_analysis,
        severity_level = "critical",
        timestamp = Sys.time()
      )
    }
  }
  
  # Check multiple testing approach
  if (!is.null(parameters$multiple_testing_method)) {
    mt_check <- validate_multiple_testing_approach(expression_data, parameters)
    if (!is.null(mt_check) && mt_check$inappropriate) {
      errors$inappropriate_multiple_testing <- list(
        rule = ERROR_PREVENTION_RULES$statistical_errors$inappropriate_multiple_testing,
        detected_values = mt_check,
        severity_level = "high",
        timestamp = Sys.time()
      )
    }
  }
  
  return(errors)
}

# Check for data quality errors
check_data_quality_errors <- function(expression_data, sample_metadata) {
  
  errors <- list()
  
  # Check for extreme outliers
  outlier_analysis <- detect_extreme_outliers(expression_data, sample_metadata)
  if (outlier_analysis$has_extreme_outliers) {
    errors$outlier_samples <- list(
      rule = ERROR_PREVENTION_RULES$data_quality_errors$outlier_samples,
      detected_values = outlier_analysis,
      severity_level = "medium",
      timestamp = Sys.time()
    )
  }
  
  # Check sequencing depth
  library_sizes <- colSums(expression_data)
  median_lib_size <- median(library_sizes)
  
  if (median_lib_size < 1000000) {
    errors$low_sequencing_depth <- list(
      rule = ERROR_PREVENTION_RULES$data_quality_errors$low_sequencing_depth,
      detected_values = list(
        median_library_size = median_lib_size,
        library_sizes = library_sizes
      ),
      severity_level = "medium",
      timestamp = Sys.time()
    )
  }
  
  return(errors)
}

# Check for interpretation errors
check_interpretation_errors <- function(parameters) {
  
  errors <- list()
  
  # Check fold change threshold
  if (!is.null(parameters$fc_cutoff) && parameters$fc_cutoff < 1.5) {
    errors$fold_change_misinterpretation <- list(
      rule = ERROR_PREVENTION_RULES$interpretation_errors$fold_change_misinterpretation,
      detected_values = list(
        current_threshold = parameters$fc_cutoff,
        recommended_minimum = 1.5
      ),
      severity_level = "medium", 
      timestamp = Sys.time()
    )
  }
  
  return(errors)
}

# ===========================================
# ERROR PREVENTION UI COMPONENTS
# ===========================================

# Create error prevention dashboard UI
create_error_prevention_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Error prevention header
    div(
      class = "error-prevention-header",
      style = "background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%); 
               color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
      
      div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        
        div(
          h2("🛡️ Smart Error Prevention", style = "margin: 0;"),
          p("Protecting your analysis from common mistakes", style = "margin: 5px 0 0 0; opacity: 0.9;")
        ),
        
        div(
          id = ns("prevention_status"),
          class = "prevention-status",
          style = "text-align: center;",
          uiOutput(ns("prevention_status_display"))
        )
      )
    ),
    
    # Error alerts container
    div(
      id = ns("error_alerts_container"),
      class = "error-alerts-container",
      uiOutput(ns("error_alerts"))
    ),
    
    # Prevention dashboard
    conditionalPanel(
      condition = paste0("output['", ns("has_errors"), "'] == true"),
      
      tabsetPanel(
        id = ns("prevention_tabs"),
        type = "pills",
        
        # Critical Errors Tab
        tabPanel(
          "🚨 Critical Issues",
          value = "critical",
          
          div(style = "padding: 20px 0;",
            h4("Issues That Must Be Fixed", style = "color: #dc2626; margin-bottom: 20px;"),
            uiOutput(ns("critical_errors_display"))
          )
        ),
        
        # Warnings Tab  
        tabPanel(
          "⚠️ Warnings",
          value = "warnings",
          
          div(style = "padding: 20px 0;",
            h4("Issues That Should Be Addressed", style = "color: #d97706; margin-bottom: 20px;"),
            uiOutput(ns("warning_errors_display"))
          )
        ),
        
        # Suggestions Tab
        tabPanel(
          "💡 Suggestions",
          value = "suggestions",
          
          div(style = "padding: 20px 0;",
            h4("Recommended Improvements", style = "color: #0284c7; margin-bottom: 20px;"),
            uiOutput(ns("suggestion_errors_display"))
          )
        )
      )
    ),
    
    # Prevention summary (when no errors)
    conditionalPanel(
      condition = paste0("output['", ns("has_errors"), "'] == false"),
      
      div(
        class = "prevention-success",
        style = "background: #f0fdf4; border: 2px solid #10b981; border-radius: 12px; 
                 padding: 30px; text-align: center; margin-top: 20px;",
        
        div(style = "font-size: 48px; margin-bottom: 15px;", "✅"),
        h3("No Issues Detected!", style = "color: #065f46; margin-bottom: 10px;"),
        p("Your analysis setup looks good to proceed.", style = "color: #047857; font-size: 16px;")
      )
    )
  )
}

# ===========================================
# ERROR PREVENTION SERVER
# ===========================================

# Error prevention server logic
error_prevention_server <- function(id, expression_data, sample_metadata, analysis_parameters) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive error detection
    detected_errors <- reactive({
      req(expression_data(), sample_metadata())
      
      withProgress(message = "Checking for potential issues...", {
        detect_analysis_errors(
          expression_data(), 
          sample_metadata(), 
          analysis_parameters()
        )
      })
    })
    
    # Check if errors exist
    output$has_errors <- reactive({
      errors <- detected_errors()
      length(errors) > 0 && 
        (length(errors$statistical) > 0 || 
         length(errors$data_quality) > 0 || 
         length(errors$interpretation) > 0)
    })
    outputOptions(output, "has_errors", suspendWhenHidden = FALSE)
    
    # Prevention status display
    output$prevention_status_display <- renderUI({
      errors <- detected_errors()
      
      # Count errors by severity
      critical_count <- count_errors_by_severity(errors, "critical")
      warning_count <- count_errors_by_severity(errors, "high")
      suggestion_count <- count_errors_by_severity(errors, "medium")
      
      total_issues <- critical_count + warning_count + suggestion_count
      
      if (total_issues == 0) {
        div(
          div(style = "font-size: 48px; color: #10b981;", "✅"),
          div(style = "font-weight: 600; margin-top: 5px;", "ALL CLEAR"),
          div(style = "font-size: 14px; opacity: 0.9; margin-top: 2px;", "No issues detected")
        )
      } else {
        status_color <- if (critical_count > 0) "#ef4444" else if (warning_count > 0) "#f59e0b" else "#0284c7"
        status_icon <- if (critical_count > 0) "🚨" else if (warning_count > 0) "⚠️" else "💡"
        
        div(
          div(style = paste0("font-size: 48px; color: ", status_color, ";"), status_icon),
          div(style = "font-weight: 600; margin-top: 5px;", 
              paste(total_issues, "ISSUES")),
          div(style = "font-size: 14px; opacity: 0.9; margin-top: 2px;",
              paste(critical_count, "critical,", warning_count, "warnings"))
        )
      }
    })
    
    # Error alerts display
    output$error_alerts <- renderUI({
      errors <- detected_errors()
      
      if (length(errors) == 0) return(NULL)
      
      # Show only most critical errors in alerts
      critical_errors <- get_errors_by_severity(errors, "critical")
      
      if (length(critical_errors) > 0) {
        lapply(critical_errors, function(error) {
          create_error_alert_card(error, "critical")
        })
      }
    })
    
    # Critical errors display
    output$critical_errors_display <- renderUI({
      errors <- detected_errors()
      critical_errors <- get_errors_by_severity(errors, "critical")
      
      if (length(critical_errors) == 0) {
        div(
          class = "alert alert-success",
          icon("check"),
          "No critical issues detected!"
        )
      } else {
        lapply(critical_errors, function(error) {
          create_detailed_error_card(error, "critical")
        })
      }
    })
    
    # Warning errors display
    output$warning_errors_display <- renderUI({
      errors <- detected_errors()
      warning_errors <- get_errors_by_severity(errors, "high")
      
      if (length(warning_errors) == 0) {
        div(
          class = "alert alert-success",
          icon("check"),
          "No warnings detected!"
        )
      } else {
        lapply(warning_errors, function(error) {
          create_detailed_error_card(error, "warning")
        })
      }
    })
    
    # Return error status for other modules
    return(list(
      has_critical_errors = reactive({
        errors <- detected_errors()
        count_errors_by_severity(errors, "critical") > 0
      }),
      error_summary = detected_errors
    ))
  })
}

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

# Create error alert card
create_error_alert_card <- function(error, severity) {
  rule <- error$rule
  
  severity_styles <- list(
    critical = list(bg = "#fef2f2", border = "#ef4444", text = "#991b1b"),
    warning = list(bg = "#fffbeb", border = "#f59e0b", text = "#92400e"),
    suggestion = list(bg = "#f0f9ff", border = "#0284c7", text = "#1e40af")
  )
  
  style <- severity_styles[[severity]]
  
  div(
    class = "error-alert-card",
    style = paste0(
      "background: ", style$bg, "; ",
      "border: 2px solid ", style$border, "; ",
      "border-radius: 12px; padding: 20px; margin-bottom: 15px; ",
      "box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
    ),
    
    # Header
    div(
      style = "display: flex; align-items: center; margin-bottom: 15px;",
      
      div(style = "font-size: 32px; margin-right: 15px;", rule$icon),
      
      div(
        style = "flex: 1;",
        h4(rule$name, style = paste0("margin: 0; color: ", style$text, ";")),
        p(rule$user_message, style = "margin: 5px 0 0 0; font-size: 16px; font-weight: 500;")
      ),
      
      if (rule$prevent_action) {
        div(
          class = "prevention-badge",
          style = "background: #dc2626; color: white; padding: 5px 10px; 
                   border-radius: 6px; font-size: 12px; font-weight: 600;",
          "BLOCKS ANALYSIS"
        )
      }
    ),
    
    # Quick explanation
    div(
      style = "margin-bottom: 15px; padding: 15px; background: white; 
               border-radius: 8px; border-left: 4px solid #3b82f6;",
      h5("Why This Matters:", style = "margin: 0 0 10px 0; color: #374151;"),
      p(rule$explanation, style = "margin: 0; line-height: 1.6; color: #4b5563;")
    ),
    
    # Primary solution
    div(
      style = "display: flex; align-items: center; justify-content: space-between;",
      
      div(
        style = "flex: 1;",
        h5("Recommended Solution:", style = "margin: 0 0 5px 0; color: #059669;"),
        p(rule$solutions$primary, style = "margin: 0; color: #047857; font-weight: 500;")
      ),
      
      actionButton(
        paste0("learn_more_", gsub("[^A-Za-z0-9]", "_", rule$name)),
        "Learn More",
        class = "btn btn-outline-info btn-sm"
      )
    )
  )
}

# Count errors by severity
count_errors_by_severity <- function(errors, severity) {
  count <- 0
  
  for (category in names(errors)) {
    if (category == "priority_order") next
    
    for (error in errors[[category]]) {
      if (error$severity_level == severity) {
        count <- count + 1
      }
    }
  }
  
  return(count)
}

# Get errors by severity
get_errors_by_severity <- function(errors, severity) {
  filtered_errors <- list()
  
  for (category in names(errors)) {
    if (category == "priority_order") next
    
    for (error_name in names(errors[[category]])) {
      error <- errors[[category]][[error_name]]
      if (error$severity_level == severity) {
        filtered_errors[[paste(category, error_name, sep = "_")]] <- error
      }
    }
  }
  
  return(filtered_errors)
}

# Analyze design confounding
analyze_confounding <- function(sample_metadata) {
  if (!"Batch" %in% colnames(sample_metadata)) {
    return(list(is_confounded = FALSE))
  }
  
  # Check if treatment groups are completely separated by batch
  contingency <- table(sample_metadata$Condition, sample_metadata$Batch)
  
  # If any row or column has all zeros except one cell, it's confounded
  is_confounded <- any(rowSums(contingency > 0) == 1) || any(colSums(contingency > 0) == 1)
  
  return(list(
    is_confounded = is_confounded,
    contingency_table = contingency,
    confounding_score = if (is_confounded) 1.0 else 0.0
  ))
}

# Detect extreme outliers
detect_extreme_outliers <- function(expression_data, sample_metadata) {
  # Simple PCA-based outlier detection
  pca_data <- prcomp(t(log2(expression_data + 1)), scale. = TRUE)
  pc1_scores <- pca_data$x[, 1]
  pc2_scores <- pca_data$x[, 2]
  
  # Define outliers as samples >3 standard deviations from mean
  pc1_outliers <- abs(scale(pc1_scores)) > 3
  pc2_outliers <- abs(scale(pc2_scores)) > 3
  
  extreme_outliers <- pc1_outliers | pc2_outliers
  
  return(list(
    has_extreme_outliers = any(extreme_outliers),
    outlier_samples = names(which(extreme_outliers)),
    outlier_count = sum(extreme_outliers),
    pca_scores = data.frame(
      Sample = rownames(pca_data$x),
      PC1 = pc1_scores,
      PC2 = pc2_scores,
      Is_Outlier = extreme_outliers
    )
  ))
}

# Validate multiple testing approach (placeholder implementation)
validate_multiple_testing_approach <- function(expression_data, parameters) {
  # Simple validation - can be expanded later
  n_genes <- nrow(expression_data)
  method <- parameters$multiple_testing_method
  
  # Basic checks
  if (is.null(method)) {
    return(list(inappropriate = FALSE, reason = "No method specified"))
  }
  
  # Very large gene sets might need special consideration
  if (n_genes > 50000 && method == "bonferroni") {
    return(list(
      inappropriate = TRUE,
      reason = "Bonferroni correction too conservative for very large gene sets",
      recommendation = "Consider Benjamini-Hochberg (BH) correction instead"
    ))
  }
  
  return(list(inappropriate = FALSE, reason = "Multiple testing approach appears appropriate"))
}

# Prioritize errors by severity
prioritize_errors <- function(detected_errors) {
  
  priority_order <- list()
  
  # Extract all errors from different categories
  all_errors <- list()
  
  for (category in names(detected_errors)) {
    if (category == "priority_order") next  # Skip the priority_order field itself
    
    if (is.list(detected_errors[[category]])) {
      for (error_name in names(detected_errors[[category]])) {
        error_info <- detected_errors[[category]][[error_name]]
        all_errors[[paste(category, error_name, sep = "_")]] <- list(
          category = category,
          name = error_name,
          severity = error_info$severity_level,
          error_info = error_info
        )
      }
    }
  }
  
  # Sort by severity: critical > high > medium > low
  severity_order <- c("critical", "high", "medium", "low")
  
  sorted_errors <- all_errors[order(match(sapply(all_errors, function(x) x$severity), severity_order))]
  
  return(names(sorted_errors))
}

# Export functions
list(
  create_error_prevention_ui = create_error_prevention_ui,
  error_prevention_server = error_prevention_server,
  detect_analysis_errors = detect_analysis_errors,
  validate_multiple_testing_approach = validate_multiple_testing_approach,
  prioritize_errors = prioritize_errors,
  ERROR_PREVENTION_RULES = ERROR_PREVENTION_RULES
)