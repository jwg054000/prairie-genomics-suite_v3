# Phase 3 - Results Validation Engine
# Biological plausibility checking and publication readiness assessment
# 
# This module provides:
# - Positive/negative control validation
# - Biological plausibility checks
# - Pathway coherence analysis
# - Publication readiness assessment
# - Suspicious pattern detection

library(shiny)
library(ggplot2)
library(plotly)
library(DT)

# ===========================================
# VALIDATION RULES CATALOG
# ===========================================

# Define comprehensive validation rules with biological context
VALIDATION_RULES <- list(
  
  # Positive control validation
  positive_controls = list(
    
    housekeeping_genes = list(
      name = "Housekeeping Gene Stability",
      description = "Reference genes should remain stable across conditions",
      genes = c("ACTB", "GAPDH", "RPL13A", "HPRT1", "TBP", "YWHAZ"),
      expected_behavior = "stable",
      tolerance = list(
        fold_change_max = 1.5,
        pvalue_min = 0.05
      ),
      validation_function = "validate_housekeeping_stability",
      importance = "critical",
      explanation = "Housekeeping genes are used as internal controls. If they're changing significantly, it suggests experimental problems or inappropriate normalization.",
      pass_criteria = "≥80% of housekeeping genes show <1.5-fold change and p>0.05",
      fail_consequences = list(
        "Questions about data quality",
        "Potential normalization issues", 
        "May indicate technical problems",
        "Reviewers will be skeptical"
      ),
      remedies = list(
        "Check sample processing consistency",
        "Verify normalization method",
        "Consider alternative reference genes",
        "Investigate technical confounders"
      )
    ),
    
    treatment_markers = list(
      name = "Treatment Response Markers",
      description = "Known treatment-responsive genes should show expected changes",
      validation_function = "validate_treatment_markers",
      importance = "high",
      explanation = "If known positive controls don't respond as expected, it questions whether the treatment worked or was properly applied.",
      context_dependent = TRUE,  # Depends on experiment type
      treatment_signatures = list(
        drug_treatment = list(
          apoptosis = c("BAX", "BBC3", "CDKN1A", "TP53"),
          cell_cycle = c("CCNA2", "CCNB1", "CDK1", "MKI67"),
          stress_response = c("ATF3", "FOS", "JUN", "EGR1")
        ),
        immune_stimulation = list(
          interferon_response = c("ISG15", "IFIT1", "STAT1", "IRF1"),
          inflammation = c("TNF", "IL1B", "IL6", "NFKB1"),
          complement = c("C3", "CFB", "C1QA", "C1QB")
        ),
        differentiation = list(
          stem_cell_markers = c("NANOG", "OCT4", "SOX2", "KLF4"),
          differentiation_markers = c("GATA1", "RUNX1", "PAX6", "MYOD1")
        )
      )
    ),
    
    negative_controls = list(
      name = "Negative Control Validation",
      description = "Genes that shouldn't change must remain stable",
      validation_function = "validate_negative_controls",
      importance = "medium",
      explanation = "Negative controls help identify false positives and validate experimental specificity.",
      examples = list(
        "Unrelated pathway genes",
        "Tissue-specific genes not expressed in your system",
        "Genes known to be unaffected by your treatment"
      )
    )
  ),
  
  # Biological plausibility checks
  biological_plausibility = list(
    
    pathway_coherence = list(
      name = "Pathway Coherence Analysis",
      description = "Genes in the same pathway should show coordinated changes",
      validation_function = "validate_pathway_coherence",
      importance = "high",
      explanation = "Biological pathways are interconnected networks. Random changes suggest technical artifacts rather than biological responses.",
      coherence_thresholds = list(
        high_coherence = 0.7,    # >70% of pathway genes change together
        moderate_coherence = 0.5, # 50-70% coherence
        low_coherence = 0.3      # <30% coherence (suspicious)
      ),
      major_pathways = list(
        "Cell Cycle" = c("CCNA2", "CCNB1", "CDK1", "CDK2", "CDKN1A", "CDKN1B"),
        "Apoptosis" = c("BAX", "BCL2", "CASP3", "CASP8", "TP53", "MDM2"),
        "DNA Repair" = c("BRCA1", "BRCA2", "ATM", "PARP1", "RAD51", "XPC"),
        "Metabolism" = c("PFKM", "ALDOA", "ENO1", "PKM", "LDHA", "G6PD"),
        "Immune Response" = c("TNF", "IL1B", "IL6", "STAT1", "IRF1", "NFKB1")
      )
    ),
    
    effect_size_plausibility = list(
      name = "Effect Size Biological Plausibility",
      description = "Fold changes should be within biologically reasonable ranges",
      validation_function = "validate_effect_sizes",
      importance = "medium",
      explanation = "Extremely large fold changes (>100x) are often technical artifacts, while very small changes (<1.2x) may not be biologically meaningful.",
      plausible_ranges = list(
        typical_range = c(1.5, 10),     # Most real changes
        high_but_possible = c(10, 50),  # Stress responses, etc.
        suspicious_range = c(50, 1000), # Likely technical
        impossible_range = 1000         # Almost certainly artifacts
      ),
      context_adjustments = list(
        stress_response = list(multiplier = 2, reason = "Stress can cause large changes"),
        development = list(multiplier = 3, reason = "Developmental switches can be dramatic"),
        cancer = list(multiplier = 2, reason = "Oncogenes can show large effects")
      )
    ),
    
    expression_distribution = list(
      name = "Expression Pattern Analysis",
      description = "Overall expression changes should follow expected patterns",
      validation_function = "validate_expression_patterns",
      importance = "medium",
      explanation = "Real biological responses typically affect 5-15% of genes. If too many or too few genes change, it suggests problems.",
      expected_patterns = list(
        typical_response = list(
          percent_changed = c(5, 15),
          up_down_ratio = c(0.5, 2.0),  # Roughly balanced up/down
          description = "Normal biological response"
        ),
        strong_response = list(
          percent_changed = c(15, 25),
          up_down_ratio = c(0.3, 3.0),
          description = "Strong treatment effect"
        ),
        suspicious_response = list(
          percent_changed = c(25, 100),
          description = "Too many genes changing - possible artifacts"
        )
      )
    )
  ),
  
  # Technical validation
  technical_validation = list(
    
    pvalue_distribution = list(
      name = "P-value Distribution Analysis",
      description = "P-values should follow expected statistical distributions",
      validation_function = "validate_pvalue_distribution",
      importance = "high",
      explanation = "Uniform p-value distribution under null hypothesis is expected. Deviations suggest problems with statistical assumptions or multiple testing.",
      expected_patterns = list(
        uniform_null = "P-values should be roughly uniform between 0.1-1.0",
        enriched_small = "Small p-values (0.0-0.1) can be enriched if real effects exist",
        suspicious_patterns = c("Too many p-values near 0.05", "Bimodal distribution", "No small p-values")
      )
    ),
    
    fold_change_distribution = list(
      name = "Fold Change Distribution",
      description = "Fold changes should show realistic distribution patterns",
      validation_function = "validate_fc_distribution",
      importance = "medium",
      explanation = "Real biological fold changes typically follow log-normal distributions. Unusual patterns suggest technical problems.",
      warning_signs = list(
        "Too many genes with identical fold changes",
        "Artificial cutoffs in the distribution", 
        "No genes with small fold changes",
        "Extremely skewed up vs down regulation"
      )
    )
  )
)

# ===========================================
# RESULTS VALIDATION ENGINE
# ===========================================

# Main validation function
validate_analysis_results <- function(deseq_results, expression_data, sample_metadata, analysis_parameters) {
  
  validation_results <- list()
  
  # Run positive control validation
  validation_results$positive_controls <- validate_positive_controls(
    deseq_results, expression_data, sample_metadata, analysis_parameters
  )
  
  # Run biological plausibility checks
  validation_results$biological_plausibility <- validate_biological_plausibility(
    deseq_results, expression_data, analysis_parameters
  )
  
  # Run technical validation
  validation_results$technical_validation <- validate_technical_aspects(
    deseq_results, expression_data
  )
  
  # Calculate overall validation score
  validation_results$overall_assessment <- calculate_validation_score(validation_results)
  
  # Generate recommendations
  validation_results$recommendations <- generate_validation_recommendations(validation_results)
  
  # Publication readiness assessment
  validation_results$publication_readiness <- assess_publication_readiness(validation_results)
  
  return(validation_results)
}

# Validate positive controls
validate_positive_controls <- function(deseq_results, expression_data, sample_metadata, parameters) {
  
  controls <- list()
  
  # Housekeeping gene validation
  hk_validation <- validate_housekeeping_genes(deseq_results)
  controls$housekeeping_genes <- list(
    rule = VALIDATION_RULES$positive_controls$housekeeping_genes,
    result = hk_validation,
    status = if (hk_validation$pass_rate >= 0.8) "pass" else "fail",
    timestamp = Sys.time()
  )
  
  # Treatment marker validation (context-dependent)
  if (!is.null(parameters$experiment_type)) {
    treatment_validation <- validate_treatment_response_markers(
      deseq_results, parameters$experiment_type
    )
    
    controls$treatment_markers <- list(
      rule = VALIDATION_RULES$positive_controls$treatment_markers,
      result = treatment_validation,
      status = treatment_validation$overall_status,
      timestamp = Sys.time()
    )
  }
  
  return(controls)
}

# Validate housekeeping gene stability
validate_housekeeping_genes <- function(deseq_results) {
  
  hk_genes <- VALIDATION_RULES$positive_controls$housekeeping_genes$genes
  tolerance <- VALIDATION_RULES$positive_controls$housekeeping_genes$tolerance
  
  # Find housekeeping genes in results
  hk_in_results <- deseq_results[rownames(deseq_results) %in% hk_genes, ]
  
  if (nrow(hk_in_results) == 0) {
    return(list(
      status = "warning",
      message = "No housekeeping genes found in results",
      genes_tested = 0,
      pass_rate = NA
    ))
  }
  
  # Check stability criteria
  stable_genes <- 0
  gene_details <- list()
  
  for (gene in rownames(hk_in_results)) {
    gene_data <- hk_in_results[gene, ]
    
    # Convert log2FoldChange to fold change
    fold_change <- 2^abs(gene_data$log2FoldChange)
    pvalue <- gene_data$pvalue
    
    is_stable <- fold_change <= tolerance$fold_change_max && 
                 (is.na(pvalue) || pvalue >= tolerance$pvalue_min)
    
    if (is_stable) stable_genes <- stable_genes + 1
    
    gene_details[[gene]] <- list(
      fold_change = fold_change,
      pvalue = pvalue,
      is_stable = is_stable,
      status = if (is_stable) "stable" else "changed"
    )
  }
  
  pass_rate <- stable_genes / nrow(hk_in_results)
  
  return(list(
    genes_tested = nrow(hk_in_results),
    stable_genes = stable_genes,
    pass_rate = pass_rate,
    gene_details = gene_details,
    overall_status = if (pass_rate >= 0.8) "pass" else "fail",
    interpretation = interpret_housekeeping_results(pass_rate, stable_genes, nrow(hk_in_results))
  ))
}

# Validate pathway coherence
validate_pathway_coherence <- function(deseq_results) {
  
  pathways <- VALIDATION_RULES$biological_plausibility$pathway_coherence$major_pathways
  coherence_results <- list()
  
  for (pathway_name in names(pathways)) {
    pathway_genes <- pathways[[pathway_name]]
    
    # Find pathway genes in results
    pathway_in_results <- deseq_results[rownames(deseq_results) %in% pathway_genes, ]
    
    if (nrow(pathway_in_results) < 3) {
      coherence_results[[pathway_name]] <- list(
        status = "insufficient_data",
        genes_found = nrow(pathway_in_results),
        coherence_score = NA
      )
      next
    }
    
    # Calculate coherence (simplified - could be more sophisticated)
    significant_genes <- pathway_in_results[pathway_in_results$padj < 0.05 & !is.na(pathway_in_results$padj), ]
    
    if (nrow(significant_genes) == 0) {
      coherence_score <- 0
    } else {
      # Check if significant genes change in same direction
      up_genes <- sum(significant_genes$log2FoldChange > 0)
      down_genes <- sum(significant_genes$log2FoldChange < 0)
      total_sig <- nrow(significant_genes)
      
      # Coherence = proportion of genes changing in majority direction
      coherence_score <- max(up_genes, down_genes) / total_sig
    }
    
    # Classify coherence level
    thresholds <- VALIDATION_RULES$biological_plausibility$pathway_coherence$coherence_thresholds
    coherence_level <- if (coherence_score >= thresholds$high_coherence) {
      "high"
    } else if (coherence_score >= thresholds$moderate_coherence) {
      "moderate"
    } else {
      "low"
    }
    
    coherence_results[[pathway_name]] <- list(
      genes_found = nrow(pathway_in_results),
      significant_genes = nrow(significant_genes),
      coherence_score = coherence_score,
      coherence_level = coherence_level,
      up_genes = if (exists("up_genes")) up_genes else 0,
      down_genes = if (exists("down_genes")) down_genes else 0,
      status = if (coherence_level %in% c("high", "moderate")) "pass" else "warning"
    )
  }
  
  return(coherence_results)
}

# ===========================================
# VALIDATION UI COMPONENTS
# ===========================================

# Create results validation dashboard UI
create_results_validation_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Validation header
    div(
      class = "validation-header",
      style = "background: linear-gradient(135deg, #059669 0%, #047857 100%); 
               color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
      
      div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        
        div(
          h2("✅ Results Validation Engine", style = "margin: 0;"),
          p("Ensuring biological plausibility and publication readiness", style = "margin: 5px 0 0 0; opacity: 0.9;")
        ),
        
        div(
          id = ns("validation_status"),
          class = "validation-status",
          style = "text-align: center;",
          uiOutput(ns("validation_status_display"))
        )
      )
    ),
    
    # Validation dashboard
    tabsetPanel(
      id = ns("validation_tabs"),
      type = "pills",
      
      # Positive Controls Tab
      tabPanel(
        "🎯 Positive Controls",
        value = "controls",
        
        div(style = "padding: 20px 0;",
          h4("Validation of Known Controls", style = "color: #059669; margin-bottom: 20px;"),
          p("These checks verify that genes with known behavior are acting as expected.", 
            style = "color: #64748b; margin-bottom: 20px;"),
          uiOutput(ns("positive_controls_display"))
        )
      ),
      
      # Biological Plausibility Tab
      tabPanel(
        "🧬 Biological Plausibility",  
        value = "biology",
        
        div(style = "padding: 20px 0;",
          h4("Biological Coherence Analysis", style = "color: #059669; margin-bottom: 20px;"),
          p("These checks ensure your results make biological sense.", 
            style = "color: #64748b; margin-bottom: 20px;"),
          uiOutput(ns("biology_validation_display"))
        )
      ),
      
      # Technical Validation Tab
      tabPanel(
        "⚙️ Technical Validation",
        value = "technical",
        
        div(style = "padding: 20px 0;",
          h4("Statistical & Technical Checks", style = "color: #059669; margin-bottom: 20px;"),
          p("These checks validate the statistical properties of your analysis.", 
            style = "color: #64748b; margin-bottom: 20px;"),
          uiOutput(ns("technical_validation_display"))
        )
      ),
      
      # Publication Readiness Tab
      tabPanel(
        "📄 Publication Ready?",
        value = "publication",
        
        div(style = "padding: 20px 0;",
          h4("Publication Readiness Assessment", style = "color: #059669; margin-bottom: 20px;"),
          uiOutput(ns("publication_assessment_display"))
        )
      )
    )
  )
}

# ===========================================
# VALIDATION SERVER
# ===========================================

# Results validation server logic
results_validation_server <- function(id, deseq_results, expression_data, sample_metadata, analysis_parameters) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive validation
    validation_results <- reactive({
      req(deseq_results(), expression_data(), sample_metadata())
      
      withProgress(message = "Validating results...", {
        validate_analysis_results(
          deseq_results(),
          expression_data(), 
          sample_metadata(),
          analysis_parameters()
        )
      })
    })
    
    # Validation status display
    output$validation_status_display <- renderUI({
      results <- validation_results()
      
      if (is.null(results)) {
        div(
          div(style = "font-size: 48px; color: #64748b;", "⏳"),
          div(style = "font-weight: 600; margin-top: 5px;", "PENDING"),
          div(style = "font-size: 14px; margin-top: 2px;", "Waiting for results")
        )
      } else {
        overall <- results$overall_assessment
        
        status_color <- switch(overall$status,
                             excellent = "#059669",
                             good = "#10b981", 
                             warning = "#f59e0b",
                             poor = "#ef4444")
        
        status_icon <- switch(overall$status,
                            excellent = "🏆",
                            good = "✅", 
                            warning = "⚠️",
                            poor = "❌")
        
        div(
          div(style = paste0("font-size: 48px; color: ", status_color, ";"), status_icon),
          div(style = "font-weight: 600; margin-top: 5px;", 
              paste("VALIDATION:", toupper(overall$status))),
          div(style = "font-size: 14px; opacity: 0.9; margin-top: 2px;",
              paste("Score:", overall$score, "/100"))
        )
      }
    })
    
    # Positive controls display
    output$positive_controls_display <- renderUI({
      results <- validation_results()
      if (is.null(results)) return(div("Loading..."))
      
      controls <- results$positive_controls
      
      lapply(names(controls), function(control_name) {
        create_validation_result_card(controls[[control_name]], control_name)
      })
    })
    
    # Return validation results for other modules
    return(validation_results)
  })
}

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

# Create validation result card
create_validation_result_card <- function(validation_result, test_name) {
  
  status_styles <- list(
    pass = list(bg = "#f0fdf4", border = "#10b981", text = "#065f46"),
    warning = list(bg = "#fffbeb", border = "#f59e0b", text = "#92400e"),
    fail = list(bg = "#fef2f2", border = "#ef4444", text = "#991b1b")
  )
  
  style <- status_styles[[validation_result$status]]
  
  div(
    class = "validation-result-card",
    style = paste0(
      "background: ", style$bg, "; ",
      "border: 2px solid ", style$border, "; ",
      "border-radius: 12px; padding: 20px; margin-bottom: 15px;"
    ),
    
    # Header
    div(
      style = "display: flex; align-items: center; justify-content: between; margin-bottom: 15px;",
      
      div(
        style = "flex: 1;",
        h5(validation_result$rule$name, style = paste0("margin: 0; color: ", style$text, ";")),
        p(validation_result$rule$description, style = "margin: 5px 0 0 0; font-size: 14px; opacity: 0.8;")
      ),
      
      div(
        class = "validation-status-badge",
        style = paste0("background: ", style$border, "; color: white; padding: 5px 12px; ",
                      "border-radius: 6px; font-size: 12px; font-weight: 600;"),
        toupper(validation_result$status)
      )
    ),
    
    # Results summary
    if (!is.null(validation_result$result)) {
      div(
        style = "background: white; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
        renderValidationDetails(validation_result$result, test_name)
      )
    },
    
    # Explanation
    div(
      style = "font-size: 14px; line-height: 1.6; color: #374151;",
      validation_result$rule$explanation
    )
  )
}

# Render validation details (simplified)
renderValidationDetails <- function(result, test_name) {
  
  if (test_name == "housekeeping_genes") {
    div(
      h6("Results:", style = "margin: 0 0 10px 0;"),
      p(paste("Genes tested:", result$genes_tested)),
      p(paste("Stable genes:", result$stable_genes)),
      p(paste("Pass rate:", round(result$pass_rate * 100, 1), "%")),
      if (!is.null(result$interpretation)) {
        p(result$interpretation, style = "font-style: italic; color: #64748b;")
      }
    )
  } else {
    div("Detailed results available")
  }
}

# Calculate overall validation score (simplified)
calculate_validation_score <- function(validation_results) {
  
  # Simplified scoring system
  score <- 100
  status <- "excellent"
  
  # Check positive controls
  if (!is.null(validation_results$positive_controls)) {
    for (control in validation_results$positive_controls) {
      if (control$status == "fail") {
        score <- score - 30
        status <- "poor"
      } else if (control$status == "warning") {
        score <- score - 15
        if (status == "excellent") status <- "warning"
      }
    }
  }
  
  # Adjust final status based on score
  if (score >= 90) status <- "excellent"
  else if (score >= 75) status <- "good"
  else if (score >= 60) status <- "warning"
  else status <- "poor"
  
  return(list(
    score = max(0, score),
    status = status,
    timestamp = Sys.time()
  ))
}

# Interpret housekeeping results
interpret_housekeeping_results <- function(pass_rate, stable_genes, total_genes) {
  if (pass_rate >= 0.8) {
    paste("Excellent! ", stable_genes, " of ", total_genes, " housekeeping genes are stable as expected.")
  } else if (pass_rate >= 0.6) {
    paste("Warning: Only ", stable_genes, " of ", total_genes, " housekeeping genes are stable. Check for normalization issues.")
  } else {
    paste("Concern: Only ", stable_genes, " of ", total_genes, " housekeeping genes are stable. This suggests data quality problems.")
  }
}

# Export functions
list(
  create_results_validation_ui = create_results_validation_ui,
  results_validation_server = results_validation_server,
  validate_analysis_results = validate_analysis_results,
  VALIDATION_RULES = VALIDATION_RULES
)