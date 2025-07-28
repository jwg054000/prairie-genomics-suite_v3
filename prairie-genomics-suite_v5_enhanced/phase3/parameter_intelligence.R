# Phase 3 - Smart Parameter Selection System
# Intelligent parameter configuration for different experiment types
# 
# This module provides:
# - Context-aware parameter recommendations
# - Educational explanations for each setting
# - Interactive parameter previews
# - Validation of parameter combinations
# - Statistical power calculations

library(shiny)
library(ggplot2)
library(plotly)

# ===========================================
# PARAMETER INTELLIGENCE ENGINE
# ===========================================

# Define parameter rules for different experiment types
PARAMETER_RULES <- list(
  
  # Statistical significance thresholds
  significance = list(
    cell_line_treatment = list(
      padj_cutoff = 0.05,
      explanation = "Standard threshold for cell line experiments with good power",
      alternatives = c(0.01, 0.05, 0.1),
      why = "Cell lines have less biological variation than primary samples"
    ),
    
    tissue_comparison = list(
      padj_cutoff = 0.01,
      explanation = "Stricter threshold needed due to tissue heterogeneity",
      alternatives = c(0.001, 0.01, 0.05),
      why = "Tissue samples have more biological variation requiring stringent filtering"
    ),
    
    patient_samples = list(
      padj_cutoff = 0.05,
      explanation = "Balanced threshold for clinical heterogeneity",
      alternatives = c(0.01, 0.05, 0.1),
      why = "Patient samples need balance between sensitivity and specificity"
    ),
    
    time_course = list(
      padj_cutoff = 0.05,
      explanation = "Standard threshold with temporal FDR correction",
      alternatives = c(0.01, 0.05, 0.1),
      why = "Time course data benefits from moderate stringency"
    )
  ),
  
  # Effect size filtering
  effect_size = list(
    cell_line_treatment = list(
      fc_cutoff = 2.0,
      explanation = "2-fold change filters noise while retaining biology",
      alternatives = c(1.5, 2.0, 3.0),
      why = "Cell lines show strong treatment responses, higher threshold reduces noise"
    ),
    
    tissue_comparison = list(
      fc_cutoff = 1.5,
      explanation = "Lower threshold captures subtle tissue differences",
      alternatives = c(1.2, 1.5, 2.0),
      why = "Tissue differences can be subtle but biologically meaningful"
    ),
    
    patient_samples = list(
      fc_cutoff = 1.5,
      explanation = "Moderate threshold balances effect size and patient variation",
      alternatives = c(1.2, 1.5, 2.0),
      why = "Patient variation requires balanced effect size filtering"
    ),
    
    time_course = list(
      fc_cutoff = 1.5,
      explanation = "Captures gradual changes over time",
      alternatives = c(1.2, 1.5, 2.0),
      why = "Temporal changes can be gradual but consistent"
    )
  ),
  
  # Sample size requirements
  sample_size = list(
    minimum_per_group = list(
      cell_line_treatment = 3,
      tissue_comparison = 4,
      patient_samples = 10,
      time_course = 3
    ),
    
    recommended_per_group = list(
      cell_line_treatment = 4,
      tissue_comparison = 6,
      patient_samples = 20,
      time_course = 4
    ),
    
    power_analysis = TRUE
  )
)

# ===========================================
# PARAMETER RECOMMENDATION ENGINE
# ===========================================

# Generate smart parameter recommendations
recommend_parameters <- function(experiment_type, sample_data, expression_data = NULL) {
  
  # Get base recommendations for experiment type
  base_params <- list()
  
  # Extract significance parameters
  if (!is.null(experiment_type) && is.character(experiment_type) && 
      experiment_type %in% names(PARAMETER_RULES$significance)) {
    sig_rule <- PARAMETER_RULES$significance[[experiment_type]]
    base_params$padj_cutoff <- sig_rule$padj_cutoff
    base_params$padj_explanation <- sig_rule$explanation
    base_params$padj_alternatives <- sig_rule$alternatives
    base_params$padj_reasoning <- sig_rule$why
  }
  
  # Extract effect size parameters  
  if (!is.null(experiment_type) && is.character(experiment_type) && 
      experiment_type %in% names(PARAMETER_RULES$effect_size)) {
    fc_rule <- PARAMETER_RULES$effect_size[[experiment_type]]
    base_params$fc_cutoff <- fc_rule$fc_cutoff
    base_params$fc_explanation <- fc_rule$explanation
    base_params$fc_alternatives <- fc_rule$alternatives
    base_params$fc_reasoning <- fc_rule$why
  }
  
  # Calculate sample size adequacy
  if (!is.null(sample_data)) {
    sample_analysis <- analyze_sample_adequacy(experiment_type, sample_data)
    base_params <- c(base_params, sample_analysis)
  }
  
  # Add data-driven adjustments if expression data available
  if (!is.null(expression_data)) {
    data_adjustments <- calculate_data_driven_adjustments(expression_data, experiment_type)
    if (!is.null(data_adjustments)) {
      base_params <- modifyList(base_params, data_adjustments)
    }
  }
  
  return(base_params)
}

# Analyze sample size adequacy
analyze_sample_adequacy <- function(experiment_type, sample_data) {
  
  # Count samples per condition
  condition_counts <- table(sample_data$Condition)
  min_per_group <- min(condition_counts)
  n_conditions <- length(condition_counts)
  
  # Get requirements for experiment type
  min_required <- PARAMETER_RULES$sample_size$minimum_per_group[[experiment_type]]
  recommended <- PARAMETER_RULES$sample_size$recommended_per_group[[experiment_type]]
  
  # Assess adequacy
  adequacy_status <- if (min_per_group >= recommended) {
    "excellent"
  } else if (min_per_group >= min_required) {
    "adequate"
  } else {
    "insufficient"
  }
  
  # Calculate statistical power
  power_est <- calculate_statistical_power(min_per_group, n_conditions, experiment_type)
  
  return(list(
    sample_count = as.list(condition_counts),
    min_per_group = min_per_group,
    adequacy_status = adequacy_status,
    adequacy_color = switch(adequacy_status,
                           excellent = "success",
                           adequate = "warning", 
                           insufficient = "danger"),
    power_estimate = power_est,
    sample_recommendation = generate_sample_recommendation(min_per_group, recommended, adequacy_status)
  ))
}

# Calculate statistical power estimate
calculate_statistical_power <- function(n_per_group, n_conditions, experiment_type) {
  
  # Simplified power calculation based on typical effect sizes
  typical_effect_sizes <- list(
    cell_line_treatment = 1.5,  # Strong treatment effects
    tissue_comparison = 1.0,    # Moderate tissue differences
    patient_samples = 0.8,      # Smaller clinical effects
    time_course = 1.2           # Temporal changes
  )
  
  effect_size <- typical_effect_sizes[[experiment_type]] %||% 1.0
  
  # Use simplified power calculation
  # Power increases with sample size and effect size
  base_power <- min(0.95, (n_per_group * effect_size) / 10)
  
  # Adjust for multiple comparisons
  if (n_conditions > 2) {
    base_power <- base_power * 0.9  # Slight reduction
  }
  
  power_percentage <- round(base_power * 100)
  
  return(list(
    power = power_percentage,
    effect_size = effect_size,
    status = if (power_percentage >= 80) "good" else if (power_percentage >= 60) "marginal" else "low"
  ))
}

# Generate sample recommendation text
generate_sample_recommendation <- function(current, recommended, status) {
  if (status == "insufficient") {
    return(paste0("⚠️ Need at least ", recommended, " samples per group for reliable results (currently ", current, ")"))
  } else if (status == "adequate") {
    return(paste0("✓ Sample size is adequate, but ", recommended, " per group would be optimal"))
  } else {
    return("✅ Excellent sample size for robust statistical analysis")
  }
}

# ===========================================
# INTERACTIVE PARAMETER UI
# ===========================================

# Create parameter selection UI with educational components
create_parameter_ui <- function(id, experiment_type) {
  ns <- NS(id)
  
  tagList(
    h3("⚙️ Analysis Parameters", style = "margin-bottom: 20px;"),
    
    p("I've selected optimal parameters for your ", tags$strong(experiment_type), " experiment. 
      You can adjust them if needed:", style = "color: #64748b; margin-bottom: 30px;"),
    
    fluidRow(
      # Left column - Parameter controls
      column(6,
        # Significance threshold
        div(
          class = "parameter-group",
          style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          
          div(
            style = "display: flex; align-items: center; margin-bottom: 10px;",
            h4("📊 Significance Threshold", style = "margin: 0; flex: 1;"),
            actionButton(ns("help_padj"), "?", class = "btn btn-sm btn-outline-info",
                        style = "border-radius: 50%; width: 30px; height: 30px;")
          ),
          
          div(
            class = "parameter-control",
            sliderInput(
              ns("padj_cutoff"),
              "Adjusted p-value cutoff:",
              min = 0.001,
              max = 0.1,
              value = 0.05,
              step = 0.001,
              width = "100%"
            ),
            
            div(
              id = ns("padj_explanation"),
              class = "parameter-explanation",
              style = "font-size: 14px; color: #64748b; margin-top: 10px; 
                       padding: 10px; background: white; border-radius: 6px;",
              "Standard threshold for cell line experiments with good statistical power"
            )
          )
        ),
        
        # Effect size threshold
        div(
          class = "parameter-group",
          style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          
          div(
            style = "display: flex; align-items: center; margin-bottom: 10px;",
            h4("📈 Effect Size Filter", style = "margin: 0; flex: 1;"),
            actionButton(ns("help_fc"), "?", class = "btn btn-sm btn-outline-info",
                        style = "border-radius: 50%; width: 30px; height: 30px;")
          ),
          
          div(
            class = "parameter-control",
            sliderInput(
              ns("fc_cutoff"),
              "Minimum fold change:",
              min = 1.1,
              max = 5.0,
              value = 2.0,
              step = 0.1,
              width = "100%"
            ),
            
            div(
              id = ns("fc_explanation"),
              class = "parameter-explanation",
              style = "font-size: 14px; color: #64748b; margin-top: 10px; 
                       padding: 10px; background: white; border-radius: 6px;",
              "2-fold change filters noise while retaining biological significance"
            )
          )
        )
      ),
      
      # Right column - Preview and validation
      column(6,
        # Sample adequacy check
        div(
          class = "validation-panel",
          style = "background: white; padding: 20px; border-radius: 12px; 
                   border: 1px solid #e5e7eb; margin-bottom: 20px;",
          
          h4("🧪 Sample Size Analysis", style = "margin-bottom: 15px;"),
          uiOutput(ns("sample_adequacy_display"))
        ),
        
        # Parameter preview
        div(
          class = "preview-panel",
          style = "background: white; padding: 20px; border-radius: 12px; 
                   border: 1px solid #e5e7eb;",
          
          h4("🔍 Expected Results Preview", style = "margin-bottom: 15px;"),
          uiOutput(ns("results_preview"))
        )
      )
    ),
    
    # Advanced options (collapsed by default)
    div(
      class = "advanced-options",
      style = "margin-top: 30px;",
      
      actionButton(
        ns("toggle_advanced"),
        "🔧 Advanced Options",
        class = "btn btn-outline-secondary",
        style = "margin-bottom: 15px;"
      ),
      
      conditionalPanel(
        condition = paste0("input['", ns("show_advanced"), "'] == true"),
        
        div(
          style = "background: #fef3c7; padding: 20px; border-radius: 12px; 
                   border-left: 4px solid #f59e0b;",
          
          div(
            class = "alert alert-warning",
            style = "margin-bottom: 20px;",
            icon("exclamation-triangle"),
            strong("Advanced Mode: "),
            "These settings are for experienced users. Changing them may affect result quality."
          ),
          
          fluidRow(
            column(6,
              checkboxInput(
                ns("independent_filtering"),
                "Independent filtering",
                value = TRUE
              ),
              
              numericInput(
                ns("min_counts"),
                "Minimum counts per gene:",
                value = 10,
                min = 1,
                max = 100
              )
            ),
            
            column(6,
              selectInput(
                ns("shrinkage_type"),
                "Effect size shrinkage:",
                choices = list(
                  "Automatic (recommended)" = "auto",
                  "Normal shrinkage" = "normal", 
                  "No shrinkage" = "none"
                ),
                selected = "auto"
              ),
              
              checkboxInput(
                ns("batch_correction"),
                "Detect batch effects",
                value = TRUE
              )
            )
          )
        )
      )
    )
  )
}

# ===========================================
# PARAMETER SERVER LOGIC
# ===========================================

# Parameter selection server
parameter_server <- function(id, experiment_type, sample_data, expression_data = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Initialize with smart recommendations
    recommendations <- reactive({
      recommend_parameters(experiment_type, sample_data(), expression_data())
    })
    
    # Update UI with recommendations
    observe({
      rec <- recommendations()
      
      # Update significance threshold
      if (!is.null(rec$padj_cutoff)) {
        updateSliderInput(session, "padj_cutoff", value = rec$padj_cutoff)
      }
      
      # Update effect size threshold
      if (!is.null(rec$fc_cutoff)) {
        updateSliderInput(session, "fc_cutoff", value = rec$fc_cutoff)
      }
      
      # Update explanations
      output$padj_explanation <- renderText({
        rec$padj_explanation %||% "Significance threshold for differential expression"
      })
      
      output$fc_explanation <- renderText({
        rec$fc_explanation %||% "Minimum fold change to consider biologically relevant"
      })
    })
    
    # Sample adequacy display
    output$sample_adequacy_display <- renderUI({
      rec <- recommendations()
      
      if (!is.null(rec$sample_count)) {
        adequacy_color <- switch(rec$adequacy_status,
                               excellent = "#10b981",
                               adequate = "#f59e0b",
                               insufficient = "#ef4444")
        
        tagList(
          # Sample counts
          div(
            style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(100px, 1fr)); 
                     gap: 10px; margin-bottom: 15px;",
            
            lapply(names(rec$sample_count), function(condition) {
              div(
                style = "text-align: center; padding: 10px; background: #f8fafc; 
                         border-radius: 6px;",
                div(style = "font-size: 20px; font-weight: bold; color: #374151;", 
                    rec$sample_count[[condition]]),
                div(style = "font-size: 12px; color: #64748b;", condition)
              )
            })
          ),
          
          # Adequacy status
          div(
            class = "alert",
            style = paste0("background: ", adequacy_color, "15; border: 1px solid ", 
                          adequacy_color, "50; color: ", adequacy_color, "; margin-bottom: 15px;"),
            
            div(style = "font-weight: 600; margin-bottom: 5px;",
                paste("Status:", toupper(rec$adequacy_status))),
            div(style = "font-size: 14px;", rec$sample_recommendation)
          ),
          
          # Power estimate
          if (!is.null(rec$power_estimate)) {
            div(
              style = "background: #f0f9ff; padding: 10px; border-radius: 6px;",
              div(style = "font-weight: 600; margin-bottom: 5px;", "Statistical Power"),
              div(style = "font-size: 14px;",
                  paste0("Estimated power: ", rec$power_estimate$power, "%")),
              div(style = "font-size: 12px; color: #64748b; margin-top: 5px;",
                  "Power ≥80% is considered adequate for reliable results")
            )
          }
        )
      } else {
        div(
          class = "alert alert-info",
          "Upload sample data to see adequacy analysis"
        )
      }
    })
    
    # Results preview
    output$results_preview <- renderUI({
      if (!is.null(expression_data())) {
        n_genes <- nrow(expression_data())
        
        # Estimate how many genes will be significant based on parameters
        est_significant <- estimate_significant_genes(
          n_genes, 
          input$padj_cutoff, 
          input$fc_cutoff,
          experiment_type
        )
        
        tagList(
          div(
            style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px;",
            
            div(
              style = "text-align: center; padding: 15px; background: #f0f9ff; 
                       border-radius: 8px;",
              div(style = "font-size: 24px; font-weight: bold; color: #0284c7;",
                  format(est_significant$total, big.mark = ",")),
              div(style = "font-size: 14px; color: #64748b;", "Expected significant genes")
            ),
            
            div(
              style = "text-align: center; padding: 15px; background: #f0fdf4; 
                       border-radius: 8px;",
              div(style = "font-size: 24px; font-weight: bold; color: #16a34a;",
                  paste0(round(est_significant$percent, 1), "%")),
              div(style = "font-size: 14px; color: #64748b;", "Of total genes")
            )
          ),
          
          div(
            style = "margin-top: 15px; font-size: 14px; color: #64748b; 
                     text-align: center;",
            "Estimates based on typical ", experiment_type, " experiments"
          )
        )
      } else {
        div(
          class = "alert alert-info",
          style = "text-align: center;",
          "Upload expression data to see results preview"
        )
      }
    })
    
    # Return current parameter settings
    return(reactive({
      list(
        padj_cutoff = input$padj_cutoff,
        fc_cutoff = input$fc_cutoff,
        independent_filtering = input$independent_filtering %||% TRUE,
        min_counts = input$min_counts %||% 10,
        shrinkage_type = input$shrinkage_type %||% "auto",
        batch_correction = input$batch_correction %||% TRUE
      )
    }))
  })
}

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

# Estimate number of significant genes based on parameters
estimate_significant_genes <- function(total_genes, padj_cutoff, fc_cutoff, experiment_type) {
  
  # Typical discovery rates for different experiment types
  base_rates <- list(
    cell_line_treatment = 0.08,  # ~8% genes change with treatment
    tissue_comparison = 0.15,    # ~15% genes differ between tissues
    patient_samples = 0.05,      # ~5% genes differ in clinical samples
    time_course = 0.10           # ~10% genes show temporal changes
  )
  
  base_rate <- base_rates[[experiment_type]] %||% 0.08
  
  # Adjust for stringency
  stringency_factor <- (0.05 / padj_cutoff) * (fc_cutoff / 2.0)
  adjusted_rate <- base_rate / stringency_factor
  
  # Ensure reasonable bounds
  adjusted_rate <- max(0.01, min(0.25, adjusted_rate))
  
  estimated_count <- round(total_genes * adjusted_rate)
  
  return(list(
    total = estimated_count,
    percent = adjusted_rate * 100
  ))
}

# Calculate data-driven parameter adjustments
calculate_data_driven_adjustments <- function(expression_data, experiment_type) {
  
  # This is a placeholder for more sophisticated data-driven adjustments
  # For now, return NULL to indicate no adjustments
  
  tryCatch({
    # Could analyze data properties like:
    # - Distribution of expression values
    # - Number of zero counts
    # - Dynamic range of expression
    # - Sample clustering patterns
    
    # For now, just return a simple adjustment based on data size
    n_genes <- nrow(expression_data)
    n_samples <- ncol(expression_data)
    
    adjustments <- list()
    
    # Adjust significance threshold for very large datasets
    if (n_genes > 30000) {
      adjustments$padj_adjustment_note <- "Consider more stringent threshold for large gene sets"
    }
    
    # Adjust for small sample sizes
    if (n_samples < 6) {
      adjustments$sample_size_note <- "Small sample size detected - results may have lower power"
    }
    
    return(adjustments)
    
  }, error = function(e) {
    # Return NULL if any errors occur
    return(NULL)
  })
}

# Export functions
list(
  create_parameter_ui = create_parameter_ui,
  parameter_server = parameter_server,
  recommend_parameters = recommend_parameters,
  calculate_data_driven_adjustments = calculate_data_driven_adjustments,
  PARAMETER_RULES = PARAMETER_RULES
)