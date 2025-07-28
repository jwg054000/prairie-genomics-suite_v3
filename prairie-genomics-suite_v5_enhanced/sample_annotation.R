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
        
        # Status display for pAnno data
        conditionalPanel(
          condition = paste0("output['", ns("panno_status"), "']"),
          uiOutput(ns("panno_info"))
        ),
        
        # Standard pattern detection info
        conditionalPanel(
          condition = paste0("!output['", ns("panno_status"), "']"),
          div(
            class = "alert alert-info",
            p("We'll analyze your sample names to suggest groupings automatically.")
          )
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
      
      # Check for pAnno annotation data
      if (!is.null(values$panno_annotation)) {
        apply_panno_annotation()
      } else {
        # Initialize default groups for manual annotation
        if (length(local_values$manual_groups) == 0) {
          local_values$manual_groups <- list(
            "Control" = character(0),
            "Treatment" = character(0)
          )
        }
      }
    }
  })
  
  # Enhanced pAnno annotation with fuzzy matching and comprehensive validation
  apply_panno_annotation <- function() {
    if (is.null(values$panno_annotation) || is.null(local_values$sample_names)) {
      return()
    }
    
    tryCatch({
      panno_data <- values$panno_annotation
      expression_samples <- local_values$sample_names
      
      cat("Applying pAnno annotation with enhanced matching...\n")
      cat("Expression samples:", length(expression_samples), "\n")
      cat("pAnno samples:", nrow(panno_data), "\n")
      
      # Perform intelligent sample matching with fuzzy logic
      matching_result <- perform_intelligent_sample_matching(expression_samples, panno_data)
      
      if (!matching_result$success) {
        showNotification(
          paste0("⚠️ pAnno matching failed: ", matching_result$message),
          type = "warning",
          duration = 8
        )
        return()
      }
      
      # Create comprehensive annotation dataframe
      annotation_df <- matching_result$annotation_data
      
      # Validate annotation for DESeq2 requirements
      validation_result <- validate_annotation_for_deseq2(annotation_df)
      
      if (!validation_result$valid) {
        showNotification(
          paste0("❌ pAnno annotation validation failed: ", validation_result$message),
          type = "error",
          duration = 10
        )
        return()
      }
      
      # Apply the validated annotation
      local_values$current_annotation <- annotation_df
      values$annotation_data <- annotation_df
      
      # Create detailed success message
      matched_samples <- sum(!is.na(annotation_df$Condition))
      unique_groups <- unique(annotation_df$Condition[!is.na(annotation_df$Condition)])
      group_counts <- table(annotation_df$Condition[!is.na(annotation_df$Condition)])
      
      success_message <- paste0(
        "🎯 pAnno auto-applied successfully!\n",
        "✅ ", matched_samples, "/", nrow(annotation_df), " samples matched (", 
        round(100 * matched_samples / nrow(annotation_df), 1), "%)\n",
        "📊 ", length(unique_groups), " groups: ", paste(names(group_counts), "(", group_counts, ")", collapse = ", "),
        "\n🧬 Ready for DESeq2 analysis!"
      )
      
      showNotification(
        success_message,
        type = "message",
        duration = 12
      )
      
      # Skip pattern detection since we have validated pAnno data
      local_values$pattern_suggestions <- list()
      
      # Log detailed results
      cat("pAnno application successful:\n")
      cat("- Total samples:", nrow(annotation_df), "\n")
      cat("- Matched samples:", matched_samples, "\n")
      cat("- Groups created:", paste(unique_groups, collapse = ", "), "\n")
      cat("- Matching method:", matching_result$method, "\n")
      
    }, error = function(e) {
      showNotification(
        paste0("❌ Critical error in pAnno annotation: ", e$message),
        type = "error",
        duration = 10
      )
      cat("Error in apply_panno_annotation:", e$message, "\n")
    })
  }
  
  # Intelligent sample matching with multiple strategies
  perform_intelligent_sample_matching <- function(expression_samples, panno_data) {
    tryCatch({
      # Strategy 1: Exact matching
      exact_matches <- intersect(expression_samples, panno_data$Sample)
      
      if (length(exact_matches) == length(expression_samples)) {
        cat("Perfect exact matching achieved\n")
        return(create_matched_annotation(expression_samples, panno_data, "exact"))
      }
      
      # Strategy 2: Case-insensitive matching  
      panno_samples_lower <- tolower(panno_data$Sample)
      expression_samples_lower <- tolower(expression_samples)
      
      case_matches <- sum(expression_samples_lower %in% panno_samples_lower)
      if (case_matches == length(expression_samples)) {
        cat("Perfect case-insensitive matching achieved\n")
        return(create_matched_annotation_fuzzy(expression_samples, panno_data, "case_insensitive"))
      }
      
      # Strategy 3: Fuzzy matching with common variations
      fuzzy_result <- perform_fuzzy_matching(expression_samples, panno_data)
      if (fuzzy_result$success && fuzzy_result$match_rate >= 0.8) {
        cat("High-quality fuzzy matching achieved:", fuzzy_result$match_rate, "\n")
        return(fuzzy_result)
      }
      
      # Strategy 4: Partial matching as last resort
      if (length(exact_matches) >= length(expression_samples) * 0.5) {
        cat("Partial matching achieved with", length(exact_matches), "exact matches\n")
        return(create_matched_annotation(expression_samples, panno_data, "partial"))
      }
      
      return(list(
        success = FALSE,
        message = paste0("Unable to match samples. Found ", length(exact_matches), 
                        " exact matches out of ", length(expression_samples), " required."),
        method = "failed"
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        message = paste0("Matching error: ", e$message),
        method = "error"
      ))
    })
  }
  
  # Create matched annotation dataframe
  create_matched_annotation <- function(expression_samples, panno_data, method) {
    annotation_df <- data.frame(
      Sample = expression_samples,
      Condition = NA,
      stringsAsFactors = FALSE
    )
    
    # Add all columns from pAnno data
    for (col in colnames(panno_data)) {
      if (col != "Sample" && col != "Condition") {
        annotation_df[[col]] <- NA
      }
    }
    
    for (i in seq_len(nrow(annotation_df))) {
      sample_name <- expression_samples[i]
      panno_match <- panno_data[panno_data$Sample == sample_name, ]
      
      if (nrow(panno_match) > 0) {
        annotation_df$Condition[i] <- as.character(panno_match$Condition[1])
        
        # Copy additional columns
        for (col in colnames(panno_data)) {
          if (col != "Sample" && col != "Condition" && col %in% colnames(annotation_df)) {
            annotation_df[[col]][i] <- panno_match[[col]][1]
          }
        }
      }
    }
    
    matched_count <- sum(!is.na(annotation_df$Condition))
    
    return(list(
      success = TRUE,
      annotation_data = annotation_df,
      method = method,
      match_rate = matched_count / length(expression_samples),
      message = paste0("Matched ", matched_count, "/", length(expression_samples), " samples")
    ))
  }
  
  # Create matched annotation dataframe with fuzzy/case-insensitive matching
  create_matched_annotation_fuzzy <- function(expression_samples, panno_data, method) {
    annotation_df <- data.frame(
      Sample = expression_samples,
      Condition = NA,
      stringsAsFactors = FALSE
    )
    
    # Add all columns from pAnno data
    for (col in colnames(panno_data)) {
      if (col != "Sample" && col != "Condition") {
        annotation_df[[col]] <- NA
      }
    }
    
    # Create matching mappings based on method
    if (method == "case_insensitive") {
      # Case-insensitive matching
      panno_samples_lower <- tolower(panno_data$Sample)
      expression_samples_lower <- tolower(expression_samples)
      
      for (i in seq_len(length(expression_samples))) {
        exp_sample <- expression_samples[i]
        exp_sample_lower <- expression_samples_lower[i]
        
        # Find case-insensitive match
        match_idx <- which(panno_samples_lower == exp_sample_lower)
        
        if (length(match_idx) > 0) {
          panno_match <- panno_data[match_idx[1], ]
          
          annotation_df$Condition[i] <- as.character(panno_match$Condition)
          
          # Copy additional columns
          for (col in colnames(panno_data)) {
            if (col != "Sample" && col != "Condition" && col %in% colnames(annotation_df)) {
              annotation_df[[col]][i] <- panno_match[[col]]
            }
          }
        }
      }
    } else if (method == "fuzzy") {
      # Fuzzy matching using similarity scoring
      for (i in seq_len(length(expression_samples))) {
        exp_sample <- expression_samples[i]
        best_match_idx <- NULL
        best_score <- 0
        
        for (j in seq_len(nrow(panno_data))) {
          panno_sample <- panno_data$Sample[j]
          score <- calculate_sample_similarity(exp_sample, panno_sample)
          
          if (score > best_score && score >= 0.8) {  # 80% similarity threshold
            best_score <- score
            best_match_idx <- j
          }
        }
        
        if (!is.null(best_match_idx)) {
          panno_match <- panno_data[best_match_idx, ]
          
          annotation_df$Condition[i] <- as.character(panno_match$Condition)
          
          # Copy additional columns
          for (col in colnames(panno_data)) {
            if (col != "Sample" && col != "Condition" && col %in% colnames(annotation_df)) {
              annotation_df[[col]][i] <- panno_match[[col]]
            }
          }
        }
      }
    }
    
    matched_count <- sum(!is.na(annotation_df$Condition))
    
    return(list(
      success = TRUE,
      annotation_data = annotation_df,
      method = method,
      match_rate = matched_count / length(expression_samples),
      message = paste0("Fuzzy matched ", matched_count, "/", length(expression_samples), " samples using ", method, " method")
    ))
  }
  
  # Fuzzy matching implementation
  perform_fuzzy_matching <- function(expression_samples, panno_data) {
    matches <- list()
    
    for (exp_sample in expression_samples) {
      best_match <- NULL
      best_score <- 0
      
      for (panno_sample in panno_data$Sample) {
        score <- calculate_sample_similarity(exp_sample, panno_sample)
        
        if (score > best_score && score >= 0.8) {  # 80% similarity threshold
          best_score <- score
          best_match <- panno_sample
        }
      }
      
      if (!is.null(best_match)) {
        matches[[exp_sample]] <- best_match
      }
    }
    
    if (length(matches) >= length(expression_samples) * 0.8) {
      # Create annotation using fuzzy matches
      annotation_df <- data.frame(
        Sample = expression_samples,
        Condition = NA,
        stringsAsFactors = FALSE
      )
      
      for (col in colnames(panno_data)) {
        if (col != "Sample" && col != "Condition") {
          annotation_df[[col]] <- NA
        }
      }
      
      for (exp_sample in names(matches)) {
        panno_sample <- matches[[exp_sample]]
        panno_match <- panno_data[panno_data$Sample == panno_sample, ]
        
        if (nrow(panno_match) > 0) {
          idx <- which(annotation_df$Sample == exp_sample)
          annotation_df$Condition[idx] <- as.character(panno_match$Condition[1])
          
          for (col in colnames(panno_data)) {
            if (col != "Sample" && col != "Condition" && col %in% colnames(annotation_df)) {
              annotation_df[[col]][idx] <- panno_match[[col]][1]
            }
          }
        }
      }
      
      return(list(
        success = TRUE,
        annotation_data = annotation_df,
        method = "fuzzy",
        match_rate = length(matches) / length(expression_samples),
        message = paste0("Fuzzy matched ", length(matches), "/", length(expression_samples), " samples")
      ))
    }
    
    return(list(
      success = FALSE,
      method = "fuzzy_failed",
      match_rate = length(matches) / length(expression_samples),
      message = "Insufficient fuzzy matches"
    ))
  }
  
  # Calculate similarity between two sample names
  calculate_sample_similarity <- function(sample1, sample2) {
    # Exact match
    if (sample1 == sample2) return(1.0)
    
    # Case insensitive
    if (tolower(sample1) == tolower(sample2)) return(0.95)
    
    # Remove common separators and compare
    clean1 <- gsub("[._-]", "", tolower(sample1))
    clean2 <- gsub("[._-]", "", tolower(sample2))
    if (clean1 == clean2) return(0.9)
    
    # Substring matching
    if (grepl(tolower(sample1), tolower(sample2)) || grepl(tolower(sample2), tolower(sample1))) {
      return(0.85)
    }
    
    return(0.0)
  }
  
  # Validation function for DESeq2 compatibility
  validate_annotation_for_deseq2 <- function(annotation_df) {
    # Check basic structure
    if (is.null(annotation_df) || nrow(annotation_df) == 0) {
      return(list(valid = FALSE, message = "Empty annotation data"))
    }
    
    if (!"Condition" %in% colnames(annotation_df)) {
      return(list(valid = FALSE, message = "Missing 'Condition' column"))
    }
    
    # Check for sufficient non-NA conditions
    valid_conditions <- !is.na(annotation_df$Condition) & annotation_df$Condition != ""
    valid_count <- sum(valid_conditions)
    
    if (valid_count < 4) {  # Minimum 4 samples for meaningful analysis
      return(list(valid = FALSE, message = paste0("Insufficient samples with conditions: ", valid_count, " (minimum 4 required)")))
    }
    
    # Check group requirements
    condition_table <- table(annotation_df$Condition[valid_conditions])
    n_groups <- length(condition_table)
    min_group_size <- min(condition_table)
    
    if (n_groups < 2) {
      return(list(valid = FALSE, message = "Need at least 2 groups for differential analysis"))
    }
    
    if (min_group_size < 2) {
      return(list(
        valid = FALSE, 
        message = paste0("Group '", names(condition_table)[which.min(condition_table)], 
                        "' has only ", min_group_size, " sample(s). Each group needs ≥2 samples")
      ))
    }
    
    return(list(
      valid = TRUE,
      message = paste0("Validation passed: ", n_groups, " groups with ", paste(condition_table, collapse = ", "), " samples each")
    ))
  }
  
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
    # Strategy 1: Extract prefixes before underscore/hyphen/dot separators
    prefixes <- sub("^([^_\\-\\.]+)[_\\-\\.].*", "\\1", sample_names)
    
    # If no separators found, use the original names (remove trailing numbers)
    no_sep_idx <- prefixes == sample_names
    if (any(no_sep_idx)) {
      prefixes[no_sep_idx] <- sub("\\d+$", "", sample_names[no_sep_idx])
    }
    
    # Strategy 2: If all prefixes are the same or too simple, try smarter extraction
    if (length(unique(prefixes)) == 1) {
      # For cases like MC9_1, M1245_1, M242_1, MLM_1
      # Extract everything before the last number before separator
      smart_prefixes <- sub("^([A-Za-z]+\\d*)[_\\-\\.].*", "\\1", sample_names)
      
      # If that gives us better grouping, use it
      if (length(unique(smart_prefixes)) > length(unique(prefixes))) {
        prefixes <- smart_prefixes
      }
    }
    
    # Strategy 3: Final fallback - extract longest common prefix pattern
    if (length(unique(prefixes)) == 1) {
      # Try to find natural breaking points
      for (i in 2:max(nchar(sample_names))) {
        candidate_prefixes <- substr(sample_names, 1, i)
        # Stop when we start getting good separation
        if (length(unique(candidate_prefixes)) > 2 && length(unique(candidate_prefixes)) < length(sample_names)) {
          prefixes <- candidate_prefixes
          break
        }
      }
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
    
    # Create single radio button group for all suggestions
    suggestion_choices <- setNames(names(suggestions), 
                                 sapply(names(suggestions), function(name) {
                                   suggestion <- suggestions[[name]]
                                   paste0(suggestion$name, " (", round(suggestion$confidence * 100, 1), "% confidence)")
                                 }))
    
    tagList(
      # Single radio button group for suggestion selection
      wellPanel(
        style = "background-color: #f8f9fa; border: 2px solid #dee2e6;",
        h5("📋 Select a Pattern:"),
        radioButtons(
          ns("selected_suggestion"),
          NULL,
          choices = suggestion_choices,
          selected = names(suggestions)[1]  # Select first (highest confidence) by default
        )
      ),
      
      br(),
      
      # Dynamic preview of selected suggestion
      uiOutput(ns("selected_suggestion_preview"))
    )
  })
  
  # Render preview of selected suggestion
  output$selected_suggestion_preview <- renderUI({
    req(input$selected_suggestion, local_values$pattern_suggestions)
    
    suggestion <- local_values$pattern_suggestions[[input$selected_suggestion]]
    
    confidence_color <- if (suggestion$confidence > 0.7) "success" else 
                       if (suggestion$confidence > 0.4) "warning" else "danger"
    
    wellPanel(
      style = paste0("border-left: 6px solid ", 
                    switch(confidence_color, 
                           "success" = "#28a745",
                           "warning" = "#ffc107", 
                           "danger" = "#dc3545")),
      
      h5(paste0("📊 Preview: ", suggestion$name)),
      p(suggestion$description),
      
      div(
        class = paste0("alert alert-", confidence_color),
        style = "margin: 10px 0;",
        paste0("✨ Confidence Score: ", round(suggestion$confidence * 100, 1), "%")
      ),
      
      h6("Sample Group Assignments:"),
      
      # Show all sample assignments in a compact table format
      tags$div(
        style = "max-height: 200px; overflow-y: auto; background: #f8f9fa; padding: 10px; border-radius: 4px;",
        tags$table(
          class = "table table-sm table-striped",
          style = "margin: 0; font-size: 12px;",
          tags$thead(
            tags$tr(
              tags$th("Sample", style = "width: 60%;"),
              tags$th("Group", style = "width: 40%;")
            )
          ),
          tags$tbody(
            lapply(1:length(local_values$sample_names), function(i) {
              tags$tr(
                tags$td(local_values$sample_names[i]),
                tags$td(
                  tags$span(
                    suggestion$groups[i],
                    class = "badge badge-primary"
                  )
                )
              )
            })
          )
        )
      ),
      
      # Group summary
      br(),
      h6("Group Summary:"),
      renderTable({
        group_counts <- table(suggestion$groups)
        data.frame(
          Group = names(group_counts),
          Samples = as.integer(group_counts),
          stringsAsFactors = FALSE
        )
      }, bordered = TRUE, striped = TRUE, hover = TRUE, width = "100%")
    )
  })
  
  # Accept suggestion
  observeEvent(input$accept_suggestion, {
    req(input$selected_suggestion, local_values$pattern_suggestions)
    
    tryCatch({
      selected_name <- input$selected_suggestion
      suggestion <- local_values$pattern_suggestions[[selected_name]]
      
      if (is.null(suggestion)) {
        showNotification("❌ Selected suggestion not found", type = "error")
        return()
      }
      
      # Validate suggestion has required data
      if (is.null(suggestion$groups) || length(suggestion$groups) != length(local_values$sample_names)) {
        showNotification("❌ Invalid suggestion data", type = "error")
        return()
      }
      
      # Create comprehensive annotation data frame
      annotation_df <- data.frame(
        Sample = local_values$sample_names,
        Condition = suggestion$groups,
        Group_Size = as.numeric(table(suggestion$groups)[suggestion$groups]),
        Confidence = suggestion$confidence,
        Detection_Method = selected_name,
        stringsAsFactors = FALSE
      )
      
      # Validate minimum requirements for DESeq2
      group_counts <- table(annotation_df$Condition)
      min_group_size <- min(group_counts)
      n_groups <- length(group_counts)
      
      if (n_groups < 2) {
        showNotification("❌ Need at least 2 groups for analysis", type = "error")
        return()
      }
      
      if (min_group_size < 2) {
        showNotification(
          paste0("⚠️ Warning: Group '", names(group_counts)[which.min(group_counts)], 
                "' has only ", min_group_size, " sample(s). Minimum 2 samples per group recommended."),
          type = "warning",
          duration = 8
        )
      }
      
      # Apply the annotation
      local_values$current_annotation <- annotation_df
      values$annotation_data <- annotation_df
      
      # Success notification with detailed info
      showNotification(
        paste0("✅ Applied ", suggestion$name, " successfully! ",
               n_groups, " groups created with ", 
               round(suggestion$confidence * 100, 1), "% confidence"),
        type = "message",
        duration = 5
      )
      
      # Log success for debugging
      cat("Successfully applied suggestion:", selected_name, "\n")
      cat("Groups created:", paste(names(group_counts), collapse = ", "), "\n")
      cat("Sample counts:", paste(group_counts, collapse = ", "), "\n")
      
    }, error = function(e) {
      showNotification(
        paste0("❌ Error applying suggestion: ", e$message),
        type = "error"
      )
      cat("Error in accept_suggestion:", e$message, "\n")
    })
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
  
  # Status outputs for UI reactivity
  output$panno_status <- reactive({
    !is.null(values$panno_annotation)
  })
  outputOptions(output, "panno_status", suspendWhenHidden = FALSE)
  
  # pAnno information display
  output$panno_info <- renderUI({
    if (!is.null(values$panno_annotation)) {
      panno_data <- values$panno_annotation
      n_samples <- nrow(panno_data)
      n_groups <- length(unique(panno_data$Condition))
      
      div(
        class = "alert alert-success",
        h6("📋 pAnno File Detected"),
        p(paste0("Annotation file loaded with ", n_samples, " samples and ", n_groups, " groups."),
          style = "margin-bottom: 5px;"),
        p("Sample annotation will be applied automatically.", style = "margin: 0; font-size: 12px;")
      )
    }
  })
  
  # Return annotation data
  return(reactive({ local_values$current_annotation }))
}