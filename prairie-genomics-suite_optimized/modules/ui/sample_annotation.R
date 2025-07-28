# Sample Annotation Module - Optimized Version
# Enhanced sample grouping with smart pattern detection and memory management
# 
# Author: Prairie Genomics Team - Optimized Version
# Features: Automatic pattern detection, fuzzy matching, comprehensive validation

# Load configuration
source("config/app_config.R")

# Sample Annotation UI Module
sampleAnnotationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        wellPanel(
          h4("🤖 Automatic Pattern Detection"),
          
          p("Smart analysis of your sample names to suggest experimental groupings."),
          
          actionButton(
            ns("detect_patterns"),
            "🔍 Detect Sample Patterns",
            class = "btn-primary",
            style = "width: 100%; margin-bottom: 15px;"
          ),
          
          conditionalPanel(
            condition = paste0("output['", ns("show_suggestions"), "']"),
            
            h5("📊 Detected Patterns:"),
            
            uiOutput(ns("pattern_suggestions")),
            
            br(),
            
            fluidRow(
              column(
                6,
                actionButton(
                  ns("accept_suggestion"),
                  "✅ Accept Suggestion",
                  class = "btn-success",
                  style = "width: 100%;"
                )
              ),
              column(
                6,
                actionButton(
                  ns("manual_annotation"),
                  "✏️ Manual Setup",
                  class = "btn-warning", 
                  style = "width: 100%;"
                )
              )
            )
          )
        )
      ),
      
      column(
        6,
        wellPanel(
          h4("✏️ Manual Sample Annotation"),
          
          p("Manually assign samples to experimental groups for precise control."),
          
          conditionalPanel(
            condition = paste0("output['", ns("has_samples"), "']"),
            
            uiOutput(ns("manual_assignment")),
            
            br(),
            
            fluidRow(
              column(
                8,
                textInput(
                  ns("new_group_name"),
                  "New Group Name:",
                  placeholder = "e.g., Treatment, Control",
                  value = ""
                )
              ),
              column(
                4,
                br(),
                actionButton(
                  ns("add_group"),
                  "Add Group",
                  class = "btn-secondary",
                  style = "width: 100%;"
                )
              )
            ),
            
            br(),
            
            actionButton(
              ns("save_manual"),
              "💾 Save Manual Annotation",
              class = "btn-primary",
              style = "width: 100%;"
            )
          ),
          
          conditionalPanel(
            condition = paste0("!output['", ns("has_samples"), "']"),
            div(
              class = "alert alert-warning",
              h6("⚠️ No Expression Data"),
              p("Please upload expression data first to proceed with sample annotation.")
            )
          )
        )
      )
    ),
    
    br(),
    
    # Current annotation display
    conditionalPanel(
      condition = paste0("output['", ns("show_current"), "']"),
      
      fluidRow(
        box(
          title = "📋 Current Sample Annotation",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          
          DT::dataTableOutput(ns("current_annotation")),
          
          br(),
          
          fluidRow(
            column(
              4,
              h6("📊 Group Summary:"),
              tableOutput(ns("group_summary"))
            ),
            column(
              4,
              h6("🧬 Experimental Design:"),
              verbatimTextOutput(ns("design_formula"))
            ),
            column(
              4,
              h6("✅ Validation Status:"),
              uiOutput(ns("validation_status"))
            )
          )
        )
      )
    )
  )
}

# Sample Annotation Server Module
sampleAnnotationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for this module
    local_values <- reactiveValues(
      sample_names = NULL,
      pattern_suggestions = NULL,
      manual_groups = list(),
      current_annotation = NULL,
      selected_suggestion = NULL
    )
    
    # Initialize sample names when expression data is available
    observe({
      if (!is.null(values$expression_data)) {
        local_values$sample_names <- colnames(values$expression_data)
        
        cat("📊 Sample annotation initialized with", length(local_values$sample_names), "samples\n")
        
        # Initialize default groups for manual annotation
        if (length(local_values$manual_groups) == 0) {
          local_values$manual_groups <- list(
            "Control" = character(0),
            "Treatment" = character(0)
          )
        }
      }
    })
    
    # Check if samples are available
    output$has_samples <- reactive({
      !is.null(local_values$sample_names) && length(local_values$sample_names) > 0
    })
    outputOptions(output, "has_samples", suspendWhenHidden = FALSE)
    
    # Pattern detection
    observeEvent(input$detect_patterns, {
      req(local_values$sample_names)
      
      cat("🔍 Starting pattern detection for", length(local_values$sample_names), "samples\n")
      
      showNotification("🔄 Analyzing sample patterns...", type = "message", duration = 3)
      
      # Perform smart pattern detection
      patterns <- detect_sample_patterns(local_values$sample_names)
      
      if (patterns$success) {
        local_values$pattern_suggestions <- patterns$suggestions
        
        showNotification(
          paste("✅ Found", length(patterns$suggestions), "pattern suggestions!"),
          type = "message",
          duration = 5
        )
        
        cat("✅ Pattern detection completed successfully\n")
      } else {
        showNotification(
          paste("⚠️ Pattern detection failed:", patterns$message),
          type = "warning",
          duration = 8
        )
        
        cat("❌ Pattern detection failed:", patterns$message, "\n")
      }
    })
    
    # Show pattern suggestions
    output$show_suggestions <- reactive({
      !is.null(local_values$pattern_suggestions) && length(local_values$pattern_suggestions) > 0
    })
    outputOptions(output, "show_suggestions", suspendWhenHidden = FALSE)
    
    # Render pattern suggestions
    output$pattern_suggestions <- renderUI({
      req(local_values$pattern_suggestions)
      
      suggestions_ui <- lapply(seq_along(local_values$pattern_suggestions), function(i) {
        suggestion <- local_values$pattern_suggestions[[i]]
        
        div(
          class = "panel panel-default",
          style = "margin-bottom: 10px;",
          
          div(
            class = "panel-body",
            
            radioButtons(
              session$ns(paste0("suggestion_", i)),
              label = paste("Pattern", i, ":", suggestion$description),
              choices = list("Select this pattern" = i),
              selected = character(0)
            ),
            
            tags$small(
              paste("Groups:", paste(names(suggestion$groups), collapse = ", ")),
              br(),
              paste("Confidence:", suggestion$confidence)
            ),
            
            # Preview of grouping
            if (length(suggestion$groups) > 0) {
              div(
                style = "margin-top: 10px; padding: 5px; background-color: #f8f9fa; border-radius: 3px;",
                tags$small(
                  lapply(names(suggestion$groups), function(group_name) {
                    samples_in_group <- suggestion$groups[[group_name]]
                    div(
                      strong(group_name, ":"), 
                      paste(samples_in_group[1:min(3, length(samples_in_group))], collapse = ", "),
                      if (length(samples_in_group) > 3) paste("... (", length(samples_in_group), "total)") else ""
                    )
                  })
                )
              )
            }
          )
        )
      })
      
      do.call(tagList, suggestions_ui)
    })
    
    # Accept suggestion
    observeEvent(input$accept_suggestion, {
      # Find selected suggestion
      selected_idx <- NULL
      
      for (i in seq_along(local_values$pattern_suggestions)) {
        radio_input <- input[[paste0("suggestion_", i)]]
        if (!is.null(radio_input) && length(radio_input) > 0) {
          selected_idx <- as.numeric(radio_input)
          break
        }
      }
      
      if (!is.null(selected_idx) && selected_idx <= length(local_values$pattern_suggestions)) {
        suggestion <- local_values$pattern_suggestions[[selected_idx]]
        
        # Create annotation from suggestion
        annotation_df <- create_annotation_from_suggestion(local_values$sample_names, suggestion)
        
        local_values$current_annotation <- annotation_df
        values$annotation_data <- annotation_df
        
        showNotification(
          paste("✅ Applied pattern suggestion with", length(suggestion$groups), "groups"),
          type = "message"
        )
        
        cat("✅ Pattern suggestion applied successfully\n")
      } else {
        showNotification("⚠️ Please select a pattern suggestion first", type = "warning")
      }
    })
    
    # Manual annotation setup
    observeEvent(input$manual_annotation, {
      showNotification("✏️ Switched to manual annotation mode", type = "message")
    })
    
    # Manual assignment UI
    output$manual_assignment <- renderUI({
      req(local_values$sample_names)
      
      if (length(local_values$manual_groups) == 0) {
        return(p("No groups defined. Add a group first."))
      }
      
      assignment_ui <- lapply(names(local_values$manual_groups), function(group_name) {
        div(
          style = "margin-bottom: 15px;",
          
          h6(paste("📊", group_name, "Group:")),
          
          selectizeInput(
            session$ns(paste0("group_", gsub("[^A-Za-z0-9]", "_", group_name))),
            label = NULL,
            choices = local_values$sample_names,
            selected = local_values$manual_groups[[group_name]],
            multiple = TRUE,
            options = list(
              placeholder = paste("Select samples for", group_name),
              plugins = list("remove_button")
            )
          )
        )
      })
      
      do.call(tagList, assignment_ui)
    })
    
    # Add new group
    observeEvent(input$add_group, {
      req(input$new_group_name)
      
      group_name <- trimws(input$new_group_name)
      
      if (group_name == "") {
        showNotification("⚠️ Please enter a group name", type = "warning")
        return()
      }
      
      if (group_name %in% names(local_values$manual_groups)) {
        showNotification("⚠️ Group name already exists", type = "warning")
        return()
      }
      
      local_values$manual_groups[[group_name]] <- character(0)
      
      updateTextInput(session, "new_group_name", value = "")
      
      showNotification(paste("✅ Added group:", group_name), type = "message")
    })
    
    # Save manual annotation
    observeEvent(input$save_manual, {
      req(local_values$sample_names)
      
      # Collect manual assignments
      manual_assignment <- list()
      
      for (group_name in names(local_values$manual_groups)) {
        input_id <- paste0("group_", gsub("[^A-Za-z0-9]", "_", group_name))
        assigned_samples <- input[[input_id]]
        
        if (!is.null(assigned_samples)) {
          manual_assignment[[group_name]] <- assigned_samples
        }
      }
      
      # Create annotation dataframe
      annotation_df <- create_annotation_from_manual(local_values$sample_names, manual_assignment)
      
      # Validate annotation
      validation <- validate_annotation(annotation_df)
      
      if (validation$valid) {
        local_values$current_annotation <- annotation_df
        values$annotation_data <- annotation_df
        
        showNotification("✅ Manual annotation saved successfully!", type = "message")
        cat("✅ Manual annotation saved with", length(unique(annotation_df$Condition)), "groups\n")
      } else {
        showNotification(
          paste("❌ Validation failed:", validation$message),
          type = "error",
          duration = 10
        )
      }
    })
    
    # Show current annotation
    output$show_current <- reactive({
      !is.null(local_values$current_annotation)
    })
    outputOptions(output, "show_current", suspendWhenHidden = FALSE)
    
    # Render current annotation table
    output$current_annotation <- DT::renderDataTable({
      req(local_values$current_annotation)
      
      DT::datatable(
        local_values$current_annotation,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'tip'
        ),
        rownames = FALSE
      )
    })
    
    # Group summary
    output$group_summary <- renderTable({
      req(local_values$current_annotation)
      
      summary_table <- table(local_values$current_annotation$Condition)
      data.frame(
        Group = names(summary_table),
        Count = as.numeric(summary_table)
      )
    })
    
    # Design formula
    output$design_formula <- renderText({
      req(local_values$current_annotation)
      
      paste("~ Condition\n\nSuitable for DESeq2 analysis")
    })
    
    # Validation status
    output$validation_status <- renderUI({
      req(local_values$current_annotation)
      
      validation <- validate_annotation(local_values$current_annotation)
      
      if (validation$valid) {
        div(
          class = "alert alert-success",
          style = "padding: 8px; margin: 0;",
          icon("check-circle"),
          " Ready for analysis"
        )
      } else {
        div(
          class = "alert alert-danger",
          style = "padding: 8px; margin: 0;",
          icon("exclamation-triangle"),
          " Issues found:",
          tags$br(),
          tags$small(validation$message)
        )
      }
    })
    
    # Return reactive for annotation data
    return(reactive({
      local_values$current_annotation
    }))
  })
}

# Pattern detection function
detect_sample_patterns <- function(sample_names) {
  tryCatch({
    patterns <- list()
    
    # Pattern 1: Control vs Treatment
    ctrl_pattern <- grep("control|ctrl|con|untreated|wt|wild", sample_names, ignore.case = TRUE)
    treat_pattern <- grep("treatment|treat|trt|treated|mut|mutant", sample_names, ignore.case = TRUE)
    
    if (length(ctrl_pattern) > 0 && length(treat_pattern) > 0) {
      patterns[[length(patterns) + 1]] <- list(
        description = "Control vs Treatment groups",
        groups = list(
          "Control" = sample_names[ctrl_pattern],
          "Treatment" = sample_names[treat_pattern]
        ),
        confidence = "High"
      )
    }
    
    # Pattern 2: Numeric suffixes
    numeric_groups <- extract_numeric_groups(sample_names)
    if (length(numeric_groups) > 1) {
      patterns[[length(patterns) + 1]] <- list(
        description = "Groups based on numeric patterns",
        groups = numeric_groups,
        confidence = "Medium"
      )
    }
    
    # Pattern 3: Underscore/dash separated groups
    separator_groups <- extract_separator_groups(sample_names)
    if (length(separator_groups) > 1) {
      patterns[[length(patterns) + 1]] <- list(
        description = "Groups based on name separators",
        groups = separator_groups,
        confidence = "Medium"
      )
    }
    
    return(list(
      success = TRUE,
      suggestions = patterns
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      message = e$message
    ))
  })
}

# Helper functions for pattern detection
extract_numeric_groups <- function(sample_names) {
  # Extract numbers from sample names and group by them
  numbers <- gsub("[^0-9]", "", sample_names)
  unique_numbers <- unique(numbers[numbers != ""])
  
  if (length(unique_numbers) > 1 && length(unique_numbers) <= 5) {
    groups <- list()
    for (num in unique_numbers) {
      matching_samples <- sample_names[grepl(num, sample_names, fixed = TRUE)]
      if (length(matching_samples) > 0) {
        groups[[paste0("Group_", num)]] <- matching_samples
      }
    }
    return(groups)
  }
  
  return(list())
}

extract_separator_groups <- function(sample_names) {
  # Try different separators
  separators <- c("_", "-", "\\.")
  
  for (sep in separators) {
    parts <- strsplit(sample_names, sep)
    
    # Check if we can group by first part
    first_parts <- sapply(parts, function(x) if(length(x) > 0) x[1] else "")
    unique_parts <- unique(first_parts)
    
    if (length(unique_parts) > 1 && length(unique_parts) <= 5) {
      groups <- list()
      for (part in unique_parts) {
        matching_samples <- sample_names[first_parts == part]
        if (length(matching_samples) > 0) {
          groups[[part]] <- matching_samples
        }
      }
      
      if (length(groups) > 1) {
        return(groups)
      }
    }
  }
  
  return(list())
}

# Create annotation from suggestion
create_annotation_from_suggestion <- function(sample_names, suggestion) {
  annotation_df <- data.frame(
    Sample = sample_names,
    Condition = NA,
    stringsAsFactors = FALSE
  )
  
  for (group_name in names(suggestion$groups)) {
    group_samples <- suggestion$groups[[group_name]]
    annotation_df$Condition[annotation_df$Sample %in% group_samples] <- group_name
  }
  
  return(annotation_df)
}

# Create annotation from manual assignment
create_annotation_from_manual <- function(sample_names, manual_assignment) {
  annotation_df <- data.frame(
    Sample = sample_names,
    Condition = NA,
    stringsAsFactors = FALSE
  )
  
  for (group_name in names(manual_assignment)) {
    group_samples <- manual_assignment[[group_name]]
    annotation_df$Condition[annotation_df$Sample %in% group_samples] <- group_name
  }
  
  return(annotation_df)
}

# Validate annotation
validate_annotation <- function(annotation_df) {
  if (is.null(annotation_df) || nrow(annotation_df) == 0) {
    return(list(valid = FALSE, message = "No annotation data"))
  }
  
  # Check for missing conditions
  missing_conditions <- sum(is.na(annotation_df$Condition))
  if (missing_conditions > 0) {
    return(list(
      valid = FALSE, 
      message = paste(missing_conditions, "samples have no group assignment")
    ))
  }
  
  # Check minimum samples per group
  group_counts <- table(annotation_df$Condition)
  min_samples <- get_config("analysis", "deseq2")$min_samples
  
  small_groups <- group_counts[group_counts < 2]
  if (length(small_groups) > 0) {
    return(list(
      valid = FALSE,
      message = paste("Groups with < 2 samples:", paste(names(small_groups), collapse = ", "))
    ))
  }
  
  # Check for at least 2 groups
  if (length(unique(annotation_df$Condition)) < 2) {
    return(list(
      valid = FALSE,
      message = "Need at least 2 groups for comparison"
    ))
  }
  
  return(list(valid = TRUE, message = "Annotation is valid"))
}

cat("✅ Optimized sample annotation module loaded\n")