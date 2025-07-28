# Enhanced Sample Annotation Module for Prairie Genomics Suite v5
# Multi-group support with smart pattern detection and batch effect handling
# Based on Emory DESeq2 methodology

# UI Function
enhancedSampleAnnotationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        h4("🧬 Multi-Group Sample Annotation"),
        
        wellPanel(
          h5("📊 Detection Summary"),
          
          # Sample detection results
          uiOutput(ns("sample_detection_summary")),
          
          br(),
          
          h5("🎯 Group Assignment Method"),
          radioButtons(
            ns("annotation_method"),
            NULL,
            choices = list(
              "🤖 Auto-detect patterns" = "auto",
              "📝 Manual assignment" = "manual",
              "📄 Upload clinical file" = "upload"
            ),
            selected = "auto"
          ),
          
          # Upload clinical data option
          conditionalPanel(
            condition = paste0("input['", ns("annotation_method"), "'] == 'upload'"),
            fileInput(
              ns("clinical_file"),
              "Upload Clinical/Metadata File",
              accept = c(".csv", ".tsv", ".txt", ".xlsx"),
              placeholder = "Optional: CSV/Excel with sample metadata"
            ),
            
            conditionalPanel(
              condition = paste0("output['", ns("has_clinical_file"), "']"),
              checkboxInput(
                ns("clinical_has_header"),
                "File has header row",
                value = TRUE
              ),
              
              uiOutput(ns("clinical_column_mapping"))
            )
          )
        ),
        
        # Batch effect detection
        wellPanel(
          h5("🔬 Batch Effect Detection"),
          
          uiOutput(ns("batch_detection_ui")),
          
          conditionalPanel(
            condition = paste0("output['", ns("batch_detected"), "']"),
            div(
              class = "alert alert-warning",
              style = "margin-top: 10px;",
              h6("⚠️ Potential Batch Effects Detected"),
              p("Consider batch correction in the analysis step.")
            )
          )
        )
      ),
      
      column(
        6,
        h4("📋 Sample Group Assignment"),
        
        # Auto-detection results
        conditionalPanel(
          condition = paste0("input['", ns("annotation_method"), "'] == 'auto'"),
          
          wellPanel(
            h5("🤖 Detected Patterns"),
            
            uiOutput(ns("auto_detection_results")),
            
            br(),
            
            h5("⚙️ Pattern Refinement"),
            
            fluidRow(
              column(
                6,
                numericInput(
                  ns("min_group_size"),
                  "Minimum group size:",
                  value = 2,
                  min = 1,
                  max = 10
                )
              ),
              column(
                6,
                selectInput(
                  ns("pattern_separator"),
                  "Pattern separator:",
                  choices = list(
                    "Underscore (_)" = "_",
                    "Dash (-)" = "-",
                    "Dot (.)" = "\\.",
                    "Space ( )" = " ",
                    "Custom" = "custom"
                  ),
                  selected = "_"
                )
              )
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("pattern_separator"), "'] == 'custom'"),
              textInput(
                ns("custom_separator"),
                "Custom separator (regex):",
                value = "_",
                placeholder = "Enter regex pattern"
              )
            ),
            
            actionButton(
              ns("rerun_detection"),
              "🔄 Re-run Detection",
              class = "btn-info"
            )
          )
        ),
        
        # Manual assignment interface
        conditionalPanel(
          condition = paste0("input['", ns("annotation_method"), "'] == 'manual'"),
          
          wellPanel(
            h5("📝 Manual Group Assignment"),
            
            uiOutput(ns("manual_assignment_ui"))
          )
        ),
        
        # Group management
        wellPanel(
          h5("🏷️ Group Management"),
          
          fluidRow(
            column(
              6,
              textInput(
                ns("new_group_name"),
                "Create new group:",
                placeholder = "Enter group name"
              )
            ),
            column(
              6,
              br(),
              actionButton(
                ns("add_group"),
                "➕ Add Group",
                class = "btn-success btn-sm"
              )
            )
          ),
          
          br(),
          
          # Group renaming
          uiOutput(ns("group_rename_ui"))
        )
      )
    ),
    
    hr(),
    
    # Final annotation table
    fluidRow(
      column(
        12,
        wellPanel(
          h4("📊 Final Sample Annotation"),
          
          # Validation status
          uiOutput(ns("validation_status")),
          
          br(),
          
          # Sample annotation table
          DT::dataTableOutput(ns("annotation_table")),
          
          br(),
          
          # Action buttons
          fluidRow(
            column(
              4,
              actionButton(
                ns("save_annotation"),
                "💾 Save Annotation",
                class = "btn-primary btn-lg",
                style = "width: 100%;"
              )
            ),
            column(
              4,
              downloadButton(
                ns("download_annotation"),
                "📥 Download Annotation",
                class = "btn-info btn-lg",
                style = "width: 100%;"
              )
            ),
            column(
              4,
              actionButton(
                ns("reset_annotation"),
                "🔄 Reset All",
                class = "btn-warning btn-lg",
                style = "width: 100%;"
              )
            )
          )
        )
      )
    ),
    
    # Status messages
    div(id = ns("status_messages"))
  )
}

# Server Function
enhancedSampleAnnotation <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    sample_names = NULL,
    detected_groups = NULL,
    current_annotation = NULL,
    batch_info = NULL,
    clinical_data = NULL,
    validation_results = NULL,
    status_messages = list()
  )
  
  # Initialize when expression data is available
  observe({
    req(values$expression_data)
    
    # Extract sample names
    sample_names <- colnames(values$expression_data)
    if (length(sample_names) == 0) return()
    
    local_values$sample_names <- sample_names
    cat("Initialized sample annotation with", length(sample_names), "samples\n")
    
    # Add small delay to ensure UI is ready
    Sys.sleep(0.1)
    
    # Run initial auto-detection only if we have valid sample names
    if (length(sample_names) >= 2) {
      auto_detect_groups()
    } else {
      add_status_message("⚠️ Need at least 2 samples for group detection", "warning")
    }
    
    # Detect potential batch effects
    detect_batch_effects()
  })
  
  # Enhanced auto-detection function with better error handling
  auto_detect_groups <- function() {
    req(local_values$sample_names)
    
    tryCatch({
      # Get separator pattern with proper defaults
      separator <- "_"  # Default separator
      if (!is.null(input$pattern_separator) && input$pattern_separator != "") {
        if (input$pattern_separator == "custom") {
          if (!is.null(input$custom_separator) && input$custom_separator != "") {
            separator <- input$custom_separator
          }
        } else {
          separator <- input$pattern_separator
        }
      }
      
      # Get minimum group size with proper default
      min_group_size <- 2  # Default minimum group size
      if (!is.null(input$min_group_size) && is.numeric(input$min_group_size)) {
        min_group_size <- max(1, as.integer(input$min_group_size))
      }
      
      cat("Auto-detecting groups with separator:", separator, "and min size:", min_group_size, "\n")
      
      # Apply detection with validation
      if (length(local_values$sample_names) < 2) {
        stop("Need at least 2 samples for group detection")
      }
      
      detection_results <- detect_sample_patterns_robust(
        local_values$sample_names,
        separator = separator,
        min_group_size = min_group_size
      )
      
      if (is.null(detection_results) || is.null(detection_results$groups)) {
        stop("Pattern detection returned invalid results")
      }
      
      local_values$detected_groups <- detection_results
      
      # Create initial annotation
      create_initial_annotation(detection_results)
      
      add_status_message("✅ Sample pattern detection completed", "success")
      cat("Detection completed successfully. Found", length(unique(detection_results$groups)), "groups\n")
      
    }, error = function(e) {
      add_status_message(paste("❌ Pattern detection failed:", e$message), "danger")
      cat("Auto-detection error:", e$message, "\n")
      
      # Create fallback annotation
      tryCatch({
        create_fallback_annotation()
        add_status_message("Created fallback annotation. Please use manual assignment.", "warning")
      }, error = function(e2) {
        cat("Fallback annotation failed:", e2$message, "\n")
      })
    })
  }
  
  # Robust and simplified pattern detection function
  detect_sample_patterns_robust <- function(sample_names, separator = "_", min_group_size = 2) {
    tryCatch({
      cat("Starting pattern detection for", length(sample_names), "samples\n")
      
      # Define common biological keywords
      bio_keywords <- c("Control", "Ctrl", "C", "Vehicle", "Veh", "Treatment", "Treat", 
                       "T", "Drug", "Compound", "B1", "TransB", "DN", "SM", "aN", 
                       "Activated", "Naive", "Low", "Med", "High")
      
      # Try different strategies in order of preference
      strategies <- list()
      
      # Strategy 1: Split by separator and try different positions
      if (separator != "" && separator != " ") {
        for (sample in sample_names[1:min(5, length(sample_names))]) {  # Test with first few samples
          parts <- unlist(strsplit(sample, separator, fixed = TRUE))
          if (length(parts) > 1) {
            # Try first part
            strategies[["first_part"]] <- sapply(sample_names, function(x) {
              parts <- unlist(strsplit(x, separator, fixed = TRUE))
              if (length(parts) > 0) parts[1] else x
            })
            
            # Try last part  
            strategies[["last_part"]] <- sapply(sample_names, function(x) {
              parts <- unlist(strsplit(x, separator, fixed = TRUE))
              if (length(parts) > 0) tail(parts, 1) else x
            })
            
            # Try second part if available
            if (length(parts) >= 2) {
              strategies[["second_part"]] <- sapply(sample_names, function(x) {
                parts <- unlist(strsplit(x, separator, fixed = TRUE))
                if (length(parts) >= 2) parts[2] else parts[1]
              })
            }
            break
          }
        }
      }
      
      # Strategy 2: Look for biological keywords
      bio_groups <- c()
      for (sample in sample_names) {
        found_keyword <- FALSE
        for (keyword in bio_keywords) {
          if (grepl(keyword, sample, ignore.case = TRUE)) {
            bio_groups <- c(bio_groups, keyword)
            found_keyword <- TRUE
            break
          }
        }
        if (!found_keyword) {
          bio_groups <- c(bio_groups, "Unknown")
        }
      }
      if (length(unique(bio_groups)) > 1) {
        strategies[["bio_keywords"]] <- bio_groups
      }
      
      # Strategy 3: Simple numeric/alphabetic patterns
      numeric_groups <- gsub("[^0-9]", "", sample_names)
      if (length(unique(numeric_groups[numeric_groups != ""])) > 1) {
        strategies[["numeric"]] <- ifelse(numeric_groups == "", "Group1", paste0("Group", numeric_groups))
      }
      
      # Strategy 4: Fallback - simple alternating groups
      if (length(strategies) == 0) {
        n_samples <- length(sample_names)
        strategies[["fallback"]] <- rep(c("Group1", "Group2"), length.out = n_samples)
      }
      
      # Evaluate strategies
      best_strategy <- NULL
      best_score <- -1
      best_result <- NULL
      
      for (strategy_name in names(strategies)) {
        groups <- strategies[[strategy_name]]
        group_counts <- table(groups)
        
        # Calculate score
        n_groups <- length(group_counts)
        min_size <- min(group_counts)
        
        score <- 0
        if (n_groups >= 2 && n_groups <= 10) score <- score + 2
        if (min_size >= min_group_size) score <- score + 3
        if (n_groups >= 2 && n_groups <= 5) score <- score + 1  # Prefer moderate number of groups
        
        cat("Strategy", strategy_name, "- Groups:", n_groups, "Min size:", min_size, "Score:", score, "\n")
        
        if (score > best_score) {
          best_score <- score
          best_strategy <- strategy_name
          best_result <- list(
            groups = groups,
            group_counts = group_counts,
            n_groups = n_groups,
            min_size = min_size
          )
        }
      }
      
      # Create final result
      if (is.null(best_result)) {
        stop("No valid grouping strategy found")
      }
      
      confidence <- min(1, best_score / 6)  # Normalize to 0-1
      
      result <- list(
        method = "robust_detection",
        strategy_used = best_strategy,
        groups = best_result$groups,
        group_counts = best_result$group_counts,
        confidence = confidence,
        separator = separator,
        n_groups = best_result$n_groups,
        min_size = best_result$min_size
      )
      
      cat("Best strategy:", best_strategy, "with confidence:", round(confidence, 2), "\n")
      return(result)
      
    }, error = function(e) {
      cat("Error in robust pattern detection:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Fallback annotation function
  create_fallback_annotation <- function() {
    req(local_values$sample_names)
    
    tryCatch({
      n_samples <- length(local_values$sample_names)
      
      # Create simple alternating groups
      if (n_samples >= 4) {
        # Create 2 groups for balanced design
        groups <- rep(c("Group1", "Group2"), length.out = n_samples)
      } else {
        # All samples in one group if too few
        groups <- rep("Group1", n_samples)
      }
      
      annotation_df <- data.frame(
        Sample = local_values$sample_names,
        Condition = groups,
        Group_Size = as.numeric(table(groups)[groups]),
        Confidence = 0.3,  # Low confidence for fallback
        Detection_Method = "fallback",
        stringsAsFactors = FALSE
      )
      
      local_values$current_annotation <- annotation_df
      
      # Store fallback detection info
      local_values$detected_groups <- list(
        method = "fallback",
        strategy_used = "alternating",
        groups = groups,
        group_counts = table(groups),
        confidence = 0.3,
        separator = "_"
      )
      
      cat("Fallback annotation created with", length(unique(groups)), "groups\n")
      
    }, error = function(e) {
      cat("Error creating fallback annotation:", e$message, "\n")
    })
  }
  
  # Create initial annotation from detection results
  create_initial_annotation <- function(detection_results) {
    req(detection_results, local_values$sample_names)
    
    annotation_df <- data.frame(
      Sample = local_values$sample_names,
      Condition = detection_results$groups,
      Group_Size = as.numeric(detection_results$group_counts[detection_results$groups]),
      Confidence = detection_results$confidence,
      Detection_Method = detection_results$strategy_used,
      stringsAsFactors = FALSE
    )
    
    local_values$current_annotation <- annotation_df
  }
  
  # Batch effect detection
  detect_batch_effects <- function() {
    req(local_values$sample_names)
    
    # Look for batch indicators in sample names
    batch_patterns <- c(
      "Batch", "batch", "Run", "run", "Seq", "seq", 
      "Lane", "lane", "Chip", "chip", "Plate", "plate"
    )
    
    batch_info <- list()
    
    for (pattern in batch_patterns) {
      # Check if pattern exists in sample names
      matches <- grep(pattern, local_values$sample_names, value = TRUE)
      if (length(matches) > 0) {
        # Extract batch identifiers
        batch_ids <- gsub(paste0(".*", pattern, "([0-9]+).*"), "\\1", matches, ignore.case = TRUE)
        if (length(unique(batch_ids)) > 1) {
          batch_info[[pattern]] <- list(
            samples = matches,
            batch_ids = batch_ids,
            n_batches = length(unique(batch_ids))
          )
        }
      }
    }
    
    local_values$batch_info <- batch_info
  }
  
  # Clinical file upload handling
  observeEvent(input$clinical_file, {
    req(input$clinical_file)
    
    tryCatch({
      # Read clinical file
      file_path <- input$clinical_file$datapath
      file_ext <- tools::file_ext(input$clinical_file$name)
      
      if (file_ext == "csv") {
        clinical_data <- readr::read_csv(file_path, col_names = input$clinical_has_header)
      } else if (file_ext %in% c("tsv", "txt")) {
        clinical_data <- readr::read_tsv(file_path, col_names = input$clinical_has_header)
      } else if (file_ext == "xlsx") {
        if (requireNamespace("readxl", quietly = TRUE)) {
          clinical_data <- readxl::read_excel(file_path, col_names = input$clinical_has_header)
        } else {
          stop("readxl package not available for Excel files")
        }
      } else {
        stop("Unsupported file format")
      }
      
      local_values$clinical_data <- clinical_data
      add_status_message("✅ Clinical file loaded successfully", "success")
      
    }, error = function(e) {
      add_status_message(paste("❌ Error loading clinical file:", e$message), "danger")
    })
  })
  
  # Re-run detection when parameters change
  observeEvent(input$rerun_detection, {
    auto_detect_groups()
  })
  
  # Handle manual assignment
  observeEvent(input$apply_manual, {
    req(local_values$sample_names)
    
    tryCatch({
      cat("Applying manual group assignments...\n")
      
      # Collect manual assignments
      manual_groups <- c()
      for (i in seq_along(local_values$sample_names)) {
        input_id <- paste0("manual_group_", i)
        group_value <- input[[input_id]]
        
        if (is.null(group_value) || group_value == "" || group_value == "Select group...") {
          manual_groups <- c(manual_groups, "Unassigned")
        } else {
          manual_groups <- c(manual_groups, group_value)
        }
      }
      
      # Validate manual assignments
      group_counts <- table(manual_groups)
      unassigned_count <- sum(manual_groups == "Unassigned")
      
      if (unassigned_count > 0) {
        add_status_message(paste("⚠️", unassigned_count, "samples remain unassigned"), "warning")
      }
      
      # Remove unassigned samples for analysis
      valid_samples <- manual_groups != "Unassigned"
      if (sum(valid_samples) < 2) {
        stop("Need at least 2 assigned samples")
      }
      
      valid_group_counts <- table(manual_groups[valid_samples])
      if (min(valid_group_counts) < 2) {
        add_status_message("⚠️ Some groups have less than 2 samples", "warning")
      }
      
      # Create annotation data frame
      annotation_df <- data.frame(
        Sample = local_values$sample_names,
        Condition = manual_groups,
        Group_Size = as.numeric(group_counts[manual_groups]),
        Confidence = 1.0,  # High confidence for manual assignment
        Detection_Method = "manual",
        stringsAsFactors = FALSE
      )
      
      local_values$current_annotation <- annotation_df
      
      # Update detected groups info
      local_values$detected_groups <- list(
        method = "manual",
        strategy_used = "user_assignment",
        groups = manual_groups,
        group_counts = group_counts,
        confidence = 1.0,
        separator = "manual"
      )
      
      add_status_message("✅ Manual assignments applied successfully", "success")
      cat("Manual assignment completed. Groups:", paste(names(valid_group_counts), collapse = ", "), "\n")
      
    }, error = function(e) {
      add_status_message(paste("❌ Manual assignment failed:", e$message), "danger")
      cat("Manual assignment error:", e$message, "\n")
    })
  })
  
  # Save annotation
  observeEvent(input$save_annotation, {
    req(local_values$current_annotation)
    
    # Validate annotation
    validation <- validate_annotation(local_values$current_annotation)
    
    if (validation$valid) {
      values$annotation_data <- local_values$current_annotation
      values$batch_data <- local_values$batch_info
      
      add_status_message("✅ Sample annotation saved successfully", "success")
    } else {
      add_status_message(paste("❌ Validation failed:", validation$message), "danger")
    }
  })
  
  # Validation function
  validate_annotation <- function(annotation_df) {
    if (is.null(annotation_df) || nrow(annotation_df) == 0) {
      return(list(valid = FALSE, message = "No annotation data"))
    }
    
    # Check for required columns
    required_cols <- c("Sample", "Condition")
    missing_cols <- setdiff(required_cols, colnames(annotation_df))
    if (length(missing_cols) > 0) {
      return(list(valid = FALSE, message = paste("Missing columns:", paste(missing_cols, collapse = ", "))))
    }
    
    # Check for minimum group sizes
    group_counts <- table(annotation_df$Condition)
    min_size <- min(group_counts)
    if (min_size < 2) {
      return(list(valid = FALSE, message = paste("Groups must have at least 2 samples. Found group with", min_size, "samples")))
    }
    
    # Check for at least 2 groups
    if (length(group_counts) < 2) {
      return(list(valid = FALSE, message = "At least 2 groups required for analysis"))
    }
    
    local_values$validation_results <- group_counts
    return(list(valid = TRUE, message = "Validation passed"))
  }
  
  # Helper function to add status messages
  add_status_message <- function(message, type = "info") {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    
    alert_class <- switch(
      type,
      "success" = "alert-success",
      "info" = "alert-info", 
      "warning" = "alert-warning",
      "danger" = "alert-danger",
      "alert-info"
    )
    
    new_message <- list(
      message = message,
      type = alert_class,
      timestamp = timestamp
    )
    
    local_values$status_messages <- append(local_values$status_messages, list(new_message))
    
    # Keep only last 5 messages
    if (length(local_values$status_messages) > 5) {
      local_values$status_messages <- tail(local_values$status_messages, 5)
    }
  }
  
  # UI Outputs
  
  # Sample detection summary
  output$sample_detection_summary <- renderUI({
    req(local_values$sample_names, local_values$detected_groups)
    
    n_samples <- length(local_values$sample_names)
    n_groups <- length(local_values$detected_groups$group_counts)
    confidence <- round(local_values$detected_groups$confidence * 100, 1)
    
    div(
      class = if (confidence > 70) "alert alert-success" else "alert alert-warning",
      h6("Sample Detection Results:"),
      tags$ul(
        tags$li(paste("Total samples:", n_samples)),
        tags$li(paste("Detected groups:", n_groups)),
        tags$li(paste("Confidence:", paste0(confidence, "%"))),
        tags$li(paste("Strategy:", local_values$detected_groups$strategy_used))
      )
    )
  })
  
  # Auto-detection results
  output$auto_detection_results <- renderUI({
    req(local_values$detected_groups)
    
    group_counts <- local_values$detected_groups$group_counts
    
    group_items <- lapply(names(group_counts), function(group) {
      tags$li(paste0(group, ": ", group_counts[[group]], " samples"))
    })
    
    tagList(
      tags$ul(group_items),
      if (local_values$detected_groups$confidence < 0.7) {
        div(
          class = "alert alert-warning",
          style = "margin-top: 10px;",
          "⚠️ Low confidence detection. Consider manual assignment or adjusting parameters."
        )
      }
    )
  })
  
  # Manual assignment UI
  output$manual_assignment_ui <- renderUI({
    req(local_values$sample_names)
    
    # Create group selection for each sample
    sample_inputs <- lapply(seq_along(local_values$sample_names), function(i) {
      sample_name <- local_values$sample_names[i]
      
      # Get current assignment if available
      current_group <- if (!is.null(local_values$current_annotation)) {
        local_values$current_annotation$Condition[i]
      } else {
        ""
      }
      
      # Available groups
      available_groups <- if (!is.null(local_values$detected_groups)) {
        names(local_values$detected_groups$group_counts)
      } else {
        c("Group1", "Group2")
      }
      
      fluidRow(
        column(
          6,
          tags$label(sample_name, style = "font-weight: normal;")
        ),
        column(
          6,
          selectInput(
            ns(paste0("manual_group_", i)),
            NULL,
            choices = c("Select group..." = "", available_groups),
            selected = if (current_group %in% available_groups) current_group else "",
            width = "100%"
          )
        )
      )
    })
    
    do.call(tagList, c(sample_inputs, list(
      br(),
      actionButton(
        ns("apply_manual"),
        "✅ Apply Manual Assignment",
        class = "btn-success"
      )
    )))
  })
  
  # Batch detection UI
  output$batch_detection_ui <- renderUI({
    req(local_values$batch_info)
    
    if (length(local_values$batch_info) == 0) {
      p("No batch effects detected in sample names", style = "color: green;")
    } else {
      batch_summaries <- lapply(names(local_values$batch_info), function(pattern) {
        info <- local_values$batch_info[[pattern]]
        p(paste0("Found ", info$n_batches, " ", pattern, " batches"))
      })
      
      do.call(tagList, batch_summaries)
    }
  })
  
  # Annotation table
  output$annotation_table <- DT::renderDataTable({
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
  
  # Validation status
  output$validation_status <- renderUI({
    req(local_values$current_annotation)
    
    validation <- validate_annotation(local_values$current_annotation)
    
    if (validation$valid) {
      div(
        class = "alert alert-success",
        h6("✅ Validation Passed"),
        p(validation$message),
        if (!is.null(local_values$validation_results)) {
          tagList(
            p("Group sizes:"),
            tags$ul(
              lapply(names(local_values$validation_results), function(group) {
                tags$li(paste0(group, ": ", local_values$validation_results[[group]], " samples"))
              })
            )
          )
        }
      )
    } else {
      div(
        class = "alert alert-danger",
        h6("❌ Validation Failed"),
        p(validation$message)
      )
    }
  })
  
  # Conditional outputs for UI reactivity
  output$has_clinical_file <- reactive({
    !is.null(local_values$clinical_data)
  })
  outputOptions(output, "has_clinical_file", suspendWhenHidden = FALSE)
  
  output$batch_detected <- reactive({
    !is.null(local_values$batch_info) && length(local_values$batch_info) > 0
  })
  outputOptions(output, "batch_detected", suspendWhenHidden = FALSE)
  
  # Download handler
  output$download_annotation <- downloadHandler(
    filename = function() {
      paste0("sample_annotation_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(local_values$current_annotation)) {
        write.csv(local_values$current_annotation, file, row.names = FALSE)
      }
    }
  )
  
  # Return the annotation data
  return(reactive({ local_values$current_annotation }))
}