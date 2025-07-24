# Sample Annotation Module for Prairie Genomics Suite
# Handles sample grouping and experimental design setup

# UI Function
sampleAnnotationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        h4("Automatic Pattern Detection"),
        
        div(
          class = "alert alert-info",
          p("We'll analyze your sample names to suggest groupings automatically.")
        ),
        
        actionButton(
          ns("detect_patterns"),
          "🔍 Detect Sample Patterns",
          class = "btn-primary",
          style = "width: 100%; margin-bottom: 15px;"
        ),
        
        conditionalPanel(
          condition = paste0("output['", ns("show_suggestions"), "']"),
          
          h5("Detected Patterns:"),
          
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
      ),
      
      column(
        6,
        h4("Manual Sample Annotation"),
        
        div(
          id = ns("manual_section"),
          
          p("Or manually assign samples to experimental groups:"),
          
          conditionalPanel(
            condition = paste0("output['", ns("has_samples"), "']"),
            
            uiOutput(ns("manual_assignment")),
            
            br(),
            
            fluidRow(
              column(
                6,
                textInput(
                  ns("new_group_name"),
                  "New Group Name:",
                  placeholder = "e.g., Treatment",
                  value = ""
                )
              ),
              column(
                6,
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
          title = "Current Sample Annotation",
          status = "success",
          solidHeader = TRUE,
          width = 12,
          
          DT::dataTableOutput(ns("current_annotation")),
          
          br(),
          
          fluidRow(
            column(
              4,
              h5("Group Summary:"),
              tableOutput(ns("group_summary"))
            ),
            column(
              8,
              h5("Experimental Design:"),
              verbatimTextOutput(ns("design_formula"))
            )
          )
        )
      )
    )
  )
}

# Server Function
sampleAnnotation <- function(input, output, session, values) {
  ns <- session$ns
  
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
      
      # Initialize default groups for manual annotation
      if (length(local_values$manual_groups) == 0) {
        local_values$manual_groups <- list(
          "Control" = character(0),
          "Treatment" = character(0)
        )
      }
    }
  })
  
  # Pattern detection
  observeEvent(input$detect_patterns, {
    req(local_values$sample_names)
    
    sample_names <- local_values$sample_names
    
    # Analyze sample name patterns
    suggestions <- detect_sample_patterns(sample_names)
    local_values$pattern_suggestions <- suggestions
    
    showNotification(
      paste0("Found ", length(suggestions), " potential grouping patterns"),
      type = "default"
    )
  })
  
  # Function to detect patterns in sample names
  detect_sample_patterns <- function(sample_names) {
    suggestions <- list()
    
    # Pattern 1: Common prefixes/suffixes
    # Look for patterns like Control_1, Control_2, Treatment_1, Treatment_2
    prefix_pattern <- extract_prefix_groups(sample_names)
    if (length(unique(prefix_pattern$groups)) > 1) {
      suggestions[["prefix"]] <- list(
        name = "Prefix-based Grouping",
        confidence = calculate_confidence(prefix_pattern$groups),
        groups = prefix_pattern$groups,
        description = paste0("Groups based on sample name prefixes: ",
                           paste(unique(prefix_pattern$groups), collapse = ", "))
      )
    }
    
    # Pattern 2: Suffix patterns  
    suffix_pattern <- extract_suffix_groups(sample_names)
    if (length(unique(suffix_pattern$groups)) > 1) {
      suggestions[["suffix"]] <- list(
        name = "Suffix-based Grouping",
        confidence = calculate_confidence(suffix_pattern$groups),
        groups = suffix_pattern$groups,
        description = paste0("Groups based on sample name suffixes: ",
                           paste(unique(suffix_pattern$groups), collapse = ", "))
      )
    }
    
    # Pattern 3: Numeric patterns (time points, doses, etc.)
    numeric_pattern <- extract_numeric_groups(sample_names)
    if (length(unique(numeric_pattern$groups)) > 1) {
      suggestions[["numeric"]] <- list(
        name = "Numeric-based Grouping",
        confidence = calculate_confidence(numeric_pattern$groups),
        groups = numeric_pattern$groups,
        description = paste0("Groups based on numeric patterns: ",
                           paste(unique(numeric_pattern$groups), collapse = ", "))
      )
    }
    
    # Sort by confidence
    suggestions <- suggestions[order(sapply(suggestions, function(x) x$confidence), decreasing = TRUE)]
    
    return(suggestions)
  }
  
  # Helper functions for pattern detection
  extract_prefix_groups <- function(sample_names) {
    # Extract common prefixes before numbers or separators
    prefixes <- sub("([A-Za-z]+)[_\\-\\.]*\\d*.*", "\\1", sample_names)
    
    # If all prefixes are the same, try splitting differently
    if (length(unique(prefixes)) == 1) {
      # Try splitting on first number
      prefixes <- sub("^([A-Za-z]+\\d*)[_\\-\\.].*", "\\1", sample_names)
    }
    
    return(list(groups = prefixes))
  }
  
  extract_suffix_groups <- function(sample_names) {
    # Extract common suffixes after numbers or separators
    suffixes <- sub(".*[_\\-\\.]([A-Za-z]+)\\d*$", "\\1", sample_names)
    
    return(list(groups = suffixes))
  }
  
  extract_numeric_groups <- function(sample_names) {
    # Extract numeric patterns that might represent conditions
    numbers <- stringr::str_extract(sample_names, "\\d+")
    
    # Group similar numbers
    if (!all(is.na(numbers))) {
      numeric_values <- as.numeric(numbers)
      # Create groups based on numeric ranges
      groups <- cut(numeric_values, breaks = 3, labels = c("Low", "Medium", "High"))
      groups <- as.character(groups)
      groups[is.na(groups)] <- "Unknown"
      
      return(list(groups = groups))
    }
    
    return(list(groups = rep("Group1", length(sample_names))))
  }
  
  calculate_confidence <- function(groups) {
    # Calculate confidence based on group balance and distinctness
    group_counts <- table(groups)
    
    # Penalize very unbalanced groups
    balance_score <- 1 - abs(max(group_counts) - min(group_counts)) / sum(group_counts)
    
    # Reward clear group separation
    separation_score <- length(unique(groups)) / length(groups)
    
    # Reward having multiple groups
    group_score <- min(1, length(unique(groups)) / 4)
    
    confidence <- (balance_score + separation_score + group_score) / 3
    return(round(confidence, 2))
  }
  
  # Render pattern suggestions
  output$pattern_suggestions <- renderUI({
    req(local_values$pattern_suggestions)
    
    suggestions <- local_values$pattern_suggestions
    
    if (length(suggestions) == 0) {
      return(div(
        class = "alert alert-warning",
        "No clear patterns detected. Please use manual annotation."
      ))
    }
    
    suggestion_cards <- lapply(names(suggestions), function(name) {
      suggestion <- suggestions[[name]]
      
      confidence_color <- if (suggestion$confidence > 0.7) "success" else 
                         if (suggestion$confidence > 0.4) "warning" else "danger"
      
      wellPanel(
        style = paste0("border-left: 4px solid ", 
                      switch(confidence_color, 
                             "success" = "#28a745",
                             "warning" = "#ffc107", 
                             "danger" = "#dc3545")),
        
        fluidRow(
          column(
            8,
            h5(suggestion$name),
            p(suggestion$description),
            tags$small(paste0("Confidence: ", suggestion$confidence * 100, "%"))
          ),
          column(
            4,
            radioButtons(
              paste0("select_", name),
              "Select:",
              choices = list("Use this" = name),
              selected = if (name == names(suggestions)[1]) name else character(0)
            )
          )
        ),
        
        # Preview of grouping
        tags$small("Sample preview:"),
        tags$ul(
          lapply(1:min(5, length(suggestion$groups)), function(i) {
            tags$li(paste0(local_values$sample_names[i], " → ", suggestion$groups[i]))
          })
        )
      )
    })
    
    do.call(tagList, suggestion_cards)
  })
  
  # Accept suggestion
  observeEvent(input$accept_suggestion, {
    req(local_values$pattern_suggestions)
    
    # Find selected suggestion
    selected_name <- NULL
    for (name in names(local_values$pattern_suggestions)) {
      if (!is.null(input[[paste0("select_", name)]]) && 
          input[[paste0("select_", name)]] == name) {
        selected_name <- name
        break
      }
    }
    
    if (!is.null(selected_name)) {
      suggestion <- local_values$pattern_suggestions[[selected_name]]
      
      # Create annotation data frame
      annotation_df <- data.frame(
        Sample = local_values$sample_names,
        Condition = suggestion$groups,
        stringsAsFactors = FALSE
      )
      
      local_values$current_annotation <- annotation_df
      values$annotation_data <- annotation_df
      
    showNotification(
      paste0("Applied ", suggestion$name, " with ", length(unique(suggestion$groups)), " groups"),
      type = "message"
    )
    } else {
      showNotification("Please select a suggestion first", type = "error")
    }
  })
  
  # Manual annotation setup
  observeEvent(input$manual_annotation, {
    # Initialize manual groups with default
    if (length(local_values$manual_groups) == 0) {
      local_values$manual_groups <- list("Group1" = character(0))
    }
  })
  
  # Add new group
  observeEvent(input$add_group, {
    req(input$new_group_name, input$new_group_name != "")
    
    group_name <- input$new_group_name
    if (!group_name %in% names(local_values$manual_groups)) {
      local_values$manual_groups[[group_name]] <- character(0)
      
      # Clear input
      updateTextInput(session, "new_group_name", value = "")
      
      showNotification(paste0("Added group: ", group_name), type = "message")
    } else {
      showNotification("Group already exists", type = "error")
    }
  })
  
  # Render manual assignment interface
  output$manual_assignment <- renderUI({
    req(local_values$sample_names)
    req(length(local_values$manual_groups) > 0)
    
    sample_names <- local_values$sample_names
    groups <- names(local_values$manual_groups)
    
    # Create dropdown for each sample
    sample_inputs <- lapply(1:length(sample_names), function(i) {
      sample_name <- sample_names[i]
      
      # Try to guess initial assignment based on sample name
      initial_assignment <- ""
      if (grepl("control", tolower(sample_name))) {
        initial_assignment <- "Control"
      } else if (grepl("treatment|treat", tolower(sample_name))) {
        initial_assignment <- "Treatment"
      }
      
      fluidRow(
        column(
          6,
          tags$label(sample_name, style = "font-weight: bold;")
        ),
        column(
          6,
          selectInput(
            ns(paste0("sample_group_", i)),
            NULL,
            choices = c("Select group..." = "", groups),
            selected = if (initial_assignment %in% groups) initial_assignment else "",
            width = "100%"
          )
        )
      )
    })
    
    do.call(tagList, sample_inputs)
  })
  
  # Save manual annotation
  observeEvent(input$save_manual, {
    req(local_values$sample_names)
    
    sample_names <- local_values$sample_names
    
    # Collect group assignments
    group_assignments <- character(length(sample_names))
    
    for (i in 1:length(sample_names)) {
      assignment <- input[[paste0("sample_group_", i)]]
      group_assignments[i] <- if (is.null(assignment) || assignment == "") "Unassigned" else assignment
    }
    
    # Debug information
    cat("Sample assignments:\n")
    for (i in 1:length(sample_names)) {
      cat(paste0("  ", sample_names[i], " -> ", group_assignments[i], "\n"))
    }
    
    # Check if all samples are assigned
    unassigned_samples <- sum(group_assignments == "Unassigned")
    if (unassigned_samples > 0) {
      showNotification(paste0("Please assign all samples to groups. ", unassigned_samples, " samples still unassigned."), type = "error")
      return()
    }
    
    # Create annotation data frame
    annotation_df <- data.frame(
      Sample = sample_names,
      Condition = group_assignments,
      stringsAsFactors = FALSE
    )
    
    local_values$current_annotation <- annotation_df
    values$annotation_data <- annotation_df
    
    showNotification(
      paste0("Manual annotation saved with ", length(unique(group_assignments)), " groups"),
      type = "message"
    )
  })
  
  # Control visibility
  output$show_suggestions <- reactive({
    !is.null(local_values$pattern_suggestions) && length(local_values$pattern_suggestions) > 0
  })
  outputOptions(output, "show_suggestions", suspendWhenHidden = FALSE)
  
  output$has_samples <- reactive({
    !is.null(local_values$sample_names) && length(local_values$sample_names) > 0
  })
  outputOptions(output, "has_samples", suspendWhenHidden = FALSE)
  
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
        dom = 'tip',
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      rownames = FALSE
    )
  })
  
  # Group summary
  output$group_summary <- renderTable({
    req(local_values$current_annotation)
    
    group_counts <- table(local_values$current_annotation$Condition)
    
    data.frame(
      Group = names(group_counts),
      Samples = as.integer(group_counts),
      stringsAsFactors = FALSE
    )
  }, rownames = FALSE)
  
  # Design formula
  output$design_formula <- renderText({
    req(local_values$current_annotation)
    
    n_groups <- length(unique(local_values$current_annotation$Condition))
    
    if (n_groups == 2) {
      return("Design: ~ Condition\nAnalysis: Two-group comparison (Control vs Treatment)")
    } else if (n_groups > 2) {
      return(paste0("Design: ~ Condition\nAnalysis: Multi-group comparison (", n_groups, " groups)"))
    } else {
      return("Warning: Only one group detected. Need at least 2 groups for analysis.")
    }
  })
  
  # Return annotation data
  return(reactive({ local_values$current_annotation }))
}