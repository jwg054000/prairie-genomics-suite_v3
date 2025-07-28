# Phase 3 - Guided Analysis Workflow Wizard
# Intelligent step-by-step guidance for wet-lab scientists
# 
# This module provides:
# - Auto-detection of experimental design
# - Step-by-step workflow with validation
# - Educational tooltips and explanations
# - Smart defaults based on experiment type
# - Progress tracking with checkpoints

library(shiny)
library(shinydashboard)

# ===========================================
# WORKFLOW WIZARD CORE
# ===========================================

# Define standard experiment types and their optimal settings
EXPERIMENT_PRESETS <- list(
  cell_line_treatment = list(
    name = "Cell Line Treatment",
    description = "Comparing treated vs untreated cell lines",
    icon = "🧫",
    min_samples = 3,
    recommended_samples = 4,
    typical_genes = 15000,
    key_parameters = list(
      padj_cutoff = 0.05,
      fc_cutoff = 2.0,
      outlier_detection = "auto",
      batch_correction = "check"
    ),
    validation_checks = c("replicate_correlation", "housekeeping_genes", "treatment_response")
  ),
  
  tissue_comparison = list(
    name = "Tissue Comparison", 
    description = "Comparing different tissue types or conditions",
    icon = "🫁",
    min_samples = 4,
    recommended_samples = 6,
    typical_genes = 20000,
    key_parameters = list(
      padj_cutoff = 0.01,
      fc_cutoff = 1.5,
      outlier_detection = "strict",
      batch_correction = "required"
    ),
    validation_checks = c("tissue_markers", "batch_effects", "sample_clustering")
  ),
  
  time_course = list(
    name = "Time Course",
    description = "Gene expression changes over time",
    icon = "⏱️",
    min_samples = 3,
    recommended_samples = 4,
    typical_genes = 15000,
    key_parameters = list(
      padj_cutoff = 0.05,
      fc_cutoff = 1.5,
      outlier_detection = "lenient",
      batch_correction = "check",
      special = "time_series_analysis"
    ),
    validation_checks = c("temporal_patterns", "circadian_genes", "progression_markers")
  ),
  
  patient_samples = list(
    name = "Patient/Clinical Samples",
    description = "Comparing patient groups or conditions",
    icon = "👥",
    min_samples = 10,
    recommended_samples = 20,
    typical_genes = 20000,
    key_parameters = list(
      padj_cutoff = 0.05,
      fc_cutoff = 1.5,
      outlier_detection = "adaptive",
      batch_correction = "required",
      special = "clinical_covariates"
    ),
    validation_checks = c("clinical_confounders", "batch_effects", "population_structure")
  ),
  
  custom = list(
    name = "Custom Analysis",
    description = "I'll configure my own settings",
    icon = "⚙️",
    min_samples = 3,
    recommended_samples = 4,
    typical_genes = 15000,
    key_parameters = list(
      padj_cutoff = 0.05,
      fc_cutoff = 2.0,
      outlier_detection = "auto",
      batch_correction = "check"
    ),
    validation_checks = c("basic_qc")
  )
)

# ===========================================
# EXPERIMENTAL DESIGN AUTO-DETECTION
# ===========================================

# Intelligently detect experiment type from sample names and metadata
detect_experiment_type <- function(sample_names, metadata = NULL) {
  
  # Initialize detection scores
  scores <- list(
    cell_line_treatment = 0,
    tissue_comparison = 0,
    time_course = 0,
    patient_samples = 0
  )
  
  # Pattern analysis
  patterns <- list(
    treatment = c("treat", "control", "ctrl", "vehicle", "dmso", "untreated", "mock"),
    cell_line = c("hela", "hek293", "mcf7", "a549", "jurkat", "k562", "rep", "replicate"),
    tissue = c("liver", "brain", "heart", "kidney", "muscle", "lung", "cortex", "tumor"),
    time = c("hr", "hour", "day", "week", "month", "time", "t0", "t1", "t2"),
    patient = c("patient", "subject", "donor", "sample", "case", "control", "pt", "ctrl")
  )
  
  # Convert to lowercase for matching
  names_lower <- tolower(sample_names)
  
  # Check for treatment patterns (cell line experiment)
  treatment_matches <- sum(sapply(patterns$treatment, function(p) any(grepl(p, names_lower))))
  cell_line_matches <- sum(sapply(patterns$cell_line, function(p) any(grepl(p, names_lower))))
  if (treatment_matches > 0 || cell_line_matches > 1) {
    scores$cell_line_treatment <- treatment_matches + cell_line_matches
  }
  
  # Check for tissue patterns
  tissue_matches <- sum(sapply(patterns$tissue, function(p) any(grepl(p, names_lower))))
  if (tissue_matches > 0) {
    scores$tissue_comparison <- tissue_matches * 2
  }
  
  # Check for time patterns
  time_matches <- sum(sapply(patterns$time, function(p) any(grepl(p, names_lower))))
  # Look for sequential numbering (T0, T1, T2...)
  if (any(grepl("t\\d+", names_lower)) || any(grepl("\\d+h", names_lower))) {
    time_matches <- time_matches + 3
  }
  if (time_matches > 0) {
    scores$time_course <- time_matches * 2
  }
  
  # Check for patient/clinical patterns
  patient_matches <- sum(sapply(patterns$patient, function(p) any(grepl(p, names_lower))))
  # Many unique sample names suggest patient samples
  if (length(unique(sample_names)) > 10) {
    patient_matches <- patient_matches + 2
  }
  if (patient_matches > 0) {
    scores$patient_samples <- patient_matches
  }
  
  # Use metadata if available for better detection
  if (!is.null(metadata) && ncol(metadata) > 1) {
    col_names <- tolower(colnames(metadata))
    
    if (any(grepl("treatment|drug|compound", col_names))) {
      scores$cell_line_treatment <- scores$cell_line_treatment + 3
    }
    if (any(grepl("tissue|organ", col_names))) {
      scores$tissue_comparison <- scores$tissue_comparison + 3
    }
    if (any(grepl("time|hour|day", col_names))) {
      scores$time_course <- scores$time_course + 3
    }
    if (any(grepl("patient|subject|clinical|age|sex", col_names))) {
      scores$patient_samples <- scores$patient_samples + 3
    }
  }
  
  # Determine most likely experiment type
  max_score <- max(unlist(scores))
  if (max_score == 0) {
    return("custom")
  }
  
  detected_type <- names(scores)[which.max(scores)]
  
  # Return detection result with confidence
  confidence <- max_score / sum(unlist(scores))
  
  return(list(
    type = detected_type,
    confidence = round(confidence * 100),
    scores = scores,
    recommendation = EXPERIMENT_PRESETS[[detected_type]]
  ))
}

# ===========================================
# WORKFLOW WIZARD UI
# ===========================================

# Create the workflow wizard UI
create_workflow_wizard_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Wizard header with progress
    div(
      id = ns("wizard_header"),
      class = "wizard-header",
      style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
               color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
      
      h2("🧬 Guided RNA-seq Analysis", style = "margin: 0 0 10px 0;"),
      p("Let me help you analyze your data correctly", style = "margin: 0; opacity: 0.9;"),
      
      # Progress indicator
      div(
        id = ns("progress_indicator"),
        class = "wizard-progress",
        style = "margin-top: 20px;",
        uiOutput(ns("progress_display"))
      )
    ),
    
    # Wizard content area
    div(
      id = ns("wizard_content"),
      class = "wizard-content",
      style = "background: white; padding: 30px; border-radius: 12px; 
               box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
      
      # Step content will be dynamically loaded here
      uiOutput(ns("step_content"))
    ),
    
    # Navigation buttons
    div(
      id = ns("wizard_navigation"),
      class = "wizard-navigation",
      style = "margin-top: 20px; display: flex; justify-content: space-between;",
      
      actionButton(
        ns("prev_step"),
        "← Previous",
        class = "btn btn-default btn-lg",
        style = "display: none;"
      ),
      
      div(
        style = "flex: 1; text-align: center;",
        uiOutput(ns("step_validation"))
      ),
      
      actionButton(
        ns("next_step"),
        "Next →",
        class = "btn btn-primary btn-lg"
      )
    )
  )
}

# ===========================================
# WORKFLOW WIZARD SERVER
# ===========================================

# Workflow wizard server logic
workflow_wizard_server <- function(id, expression_data, annotation_data) {
  moduleServer(id, function(input, output, session) {
    
    # Wizard state
    wizard_state <- reactiveValues(
      current_step = 1,
      total_steps = 5,
      experiment_type = NULL,
      detected_type = NULL,
      parameters = list(),
      validation_passed = list(),
      analysis_ready = FALSE
    )
    
    # Define wizard steps
    wizard_steps <- list(
      list(name = "Upload Review", icon = "📁", validator = "validate_upload"),
      list(name = "Experiment Type", icon = "🔬", validator = "validate_experiment"),
      list(name = "Sample Setup", icon = "🧪", validator = "validate_samples"),
      list(name = "Parameters", icon = "⚙️", validator = "validate_parameters"),
      list(name = "Review & Run", icon = "🚀", validator = "validate_final")
    )
    
    # Progress display
    output$progress_display <- renderUI({
      steps_ui <- lapply(1:wizard_state$total_steps, function(i) {
        step <- wizard_steps[[i]]
        status_class <- if (i < wizard_state$current_step) {
          "step-complete"
        } else if (i == wizard_state$current_step) {
          "step-active"
        } else {
          "step-pending"
        }
        
        div(
          class = paste("wizard-step", status_class),
          style = "display: inline-block; margin: 0 10px; text-align: center;",
          
          div(
            class = "step-icon",
            style = "font-size: 24px; margin-bottom: 5px;",
            step$icon
          ),
          div(
            class = "step-name",
            style = "font-size: 12px; font-weight: 600;",
            step$name
          )
        )
      })
      
      div(
        style = "display: flex; justify-content: center; align-items: center;",
        steps_ui
      )
    })
    
    # Step content renderer
    output$step_content <- renderUI({
      step_name <- wizard_steps[[wizard_state$current_step]]$name
      
      switch(step_name,
        "Upload Review" = render_upload_review_step(),
        "Experiment Type" = render_experiment_type_step(),
        "Sample Setup" = render_sample_setup_step(),
        "Parameters" = render_parameters_step(),
        "Review & Run" = render_review_run_step()
      )
    })
    
    # Step 1: Upload Review
    render_upload_review_step <- function() {
      ns <- session$ns
      
      data_summary <- if (!is.null(expression_data())) {
        samples <- colnames(expression_data())
        n_samples <- length(samples)
        n_genes <- nrow(expression_data())
        
        # Auto-detect experiment type
        detection <- detect_experiment_type(samples, annotation_data())
        wizard_state$detected_type <- detection
        
        tagList(
          h3("📊 Data Upload Summary"),
          
          div(
            class = "data-summary-cards",
            style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0;",
            
            div(
              class = "summary-card",
              style = "background: #f0f9ff; padding: 20px; border-radius: 8px; text-align: center;",
              h4(n_samples, style = "margin: 0; color: #0284c7; font-size: 36px;"),
              p("Samples", style = "margin: 5px 0 0 0; color: #64748b;")
            ),
            
            div(
              class = "summary-card", 
              style = "background: #f0fdf4; padding: 20px; border-radius: 8px; text-align: center;",
              h4(format(n_genes, big.mark = ","), style = "margin: 0; color: #16a34a; font-size: 36px;"),
              p("Genes", style = "margin: 5px 0 0 0; color: #64748b;")
            ),
            
            div(
              class = "summary-card",
              style = "background: #fef3c7; padding: 20px; border-radius: 8px; text-align: center;",
              h4(paste0(detection$confidence, "%"), style = "margin: 0; color: #d97706; font-size: 36px;"),
              p("Detection Confidence", style = "margin: 5px 0 0 0; color: #64748b;")
            )
          ),
          
          div(
            class = "alert alert-info",
            style = "margin-top: 20px;",
            icon("lightbulb"),
            strong("Auto-Detection: "),
            paste("This looks like a", detection$recommendation$name, "experiment"),
            if (detection$confidence < 80) {
              span(" (Low confidence - please verify)", style = "color: #d97706;")
            }
          ),
          
          h4("Sample Names Detected:"),
          div(
            style = "background: #f8fafc; padding: 15px; border-radius: 8px; 
                     max-height: 200px; overflow-y: auto;",
            tags$code(paste(samples, collapse = ", "))
          )
        )
      } else {
        div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          "No data uploaded yet. Please upload your expression data first."
        )
      }
      
      data_summary
    }
    
    # Step 2: Experiment Type Selection
    render_experiment_type_step <- function() {
      ns <- session$ns
      
      tagList(
        h3("🔬 Select Your Experiment Type"),
        p("Choose the option that best describes your RNA-seq experiment:", 
          style = "color: #64748b; margin-bottom: 20px;"),
        
        div(
          class = "experiment-type-grid",
          style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px;",
          
          lapply(names(EXPERIMENT_PRESETS), function(exp_type) {
            preset <- EXPERIMENT_PRESETS[[exp_type]]
            is_detected <- (!is.null(wizard_state$detected_type) && 
                           wizard_state$detected_type$type == exp_type)
            
            div(
              class = "experiment-type-card",
              style = paste0(
                "border: 2px solid ", ifelse(is_detected, "#3b82f6", "#e5e7eb"), "; ",
                "border-radius: 12px; padding: 20px; cursor: pointer; ",
                "transition: all 0.3s ease; position: relative;"
              ),
              onclick = paste0("Shiny.setInputValue('", ns("selected_experiment"), "', '", exp_type, "', {priority: 'event'})"),
              
              if (is_detected) {
                div(
                  style = "position: absolute; top: 10px; right: 10px; 
                          background: #3b82f6; color: white; padding: 4px 8px; 
                          border-radius: 4px; font-size: 12px;",
                  "Auto-detected"
                )
              },
              
              div(style = "font-size: 48px; margin-bottom: 10px;", preset$icon),
              h4(preset$name, style = "margin: 0 0 10px 0;"),
              p(preset$description, style = "color: #64748b; font-size: 14px;"),
              
              div(
                style = "margin-top: 15px; font-size: 12px; color: #94a3b8;",
                paste("Min samples:", preset$min_samples, "| Recommended:", preset$recommended_samples)
              )
            )
          })
        ),
        
        # Educational info
        div(
          class = "alert alert-info",
          style = "margin-top: 20px;",
          icon("info-circle"),
          strong("Not sure? "),
          "Hover over each option to learn more about when to use it."
        )
      )
    }
    
    # Navigation handlers
    observeEvent(input$next_step, {
      if (wizard_state$current_step < wizard_state$total_steps) {
        # Validate current step
        validator <- wizard_steps[[wizard_state$current_step]]$validator
        if (validate_step(validator)) {
          wizard_state$current_step <- wizard_state$current_step + 1
          updateActionButton(session, "prev_step", style = "display: inline-block;")
          
          if (wizard_state$current_step == wizard_state$total_steps) {
            updateActionButton(session, "next_step", label = "Run Analysis 🚀")
          }
        }
      } else {
        # Final step - run analysis
        wizard_state$analysis_ready <- TRUE
      }
    })
    
    observeEvent(input$prev_step, {
      if (wizard_state$current_step > 1) {
        wizard_state$current_step <- wizard_state$current_step - 1
        
        if (wizard_state$current_step == 1) {
          updateActionButton(session, "prev_step", style = "display: none;")
        }
        
        updateActionButton(session, "next_step", label = "Next →")
      }
    })
    
    # Return wizard state for integration
    return(wizard_state)
  })
}

# ===========================================
# STEP VALIDATION FUNCTIONS
# ===========================================

validate_step <- function(validator_name) {
  # Placeholder for step validation logic
  # Each step will have specific validation rules
  return(TRUE)
}

# Export functions
list(
  create_workflow_wizard_ui = create_workflow_wizard_ui,
  workflow_wizard_server = workflow_wizard_server,
  detect_experiment_type = detect_experiment_type,
  EXPERIMENT_PRESETS = EXPERIMENT_PRESETS
)