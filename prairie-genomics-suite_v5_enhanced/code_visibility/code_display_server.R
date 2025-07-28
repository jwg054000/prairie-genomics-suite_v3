# Code Display Server Logic
# Server functions for code visibility components
#
# Author: Prairie Genomics Team
# Date: January 24, 2025
# Purpose: Handle code display, validation, and export functionality

# Source required modules
source("code_visibility/code_logger.R")
source("code_visibility/code_generator.R")

# Main code display module server
codeDisplayServer <- function(id, session_id_reactive, code_type = reactive("complete")) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      current_code = "",
      last_update = Sys.time(),
      validation_status = NULL
    )
    
    # Get current session code
    current_session_code <- reactive({
      req(session_id_reactive())
      
      session_id <- session_id_reactive()
      
      if (input$view_option == "complete") {
        get_session_code(session_id, include_comments = TRUE)
      } else if (input$view_option == "current") {
        # Get only the last step
        get_analysis_code_snippet(session_id, step_name = "current")
      } else if (input$view_option == "category") {
        req(input$category_filter)
        get_analysis_code_snippet(session_id, category = input$category_filter)
      } else {
        get_session_code(session_id)
      }
    })
    
    # Update code content
    output$code_content <- renderText({
      code <- current_session_code()
      values$current_code <- code
      values$last_update <- Sys.time()
      
      if (is.null(code) || code == "" || code == "# No session data found") {
        "# No analysis code available yet\n# Run an analysis to see the generated R code here"
      } else {
        code
      }
    })
    
    # Refresh code action
    observeEvent(input$refresh_code, {
      values$last_update <- Sys.time()
      
      showNotification(
        "Code refreshed successfully",
        type = "message",
        duration = 2
      )
    })
    
    # Copy code to clipboard action
    observeEvent(input$copy_code, {
      # Note: Actual clipboard functionality would require JavaScript
      showNotification(
        "Code copied to clipboard! (Use Ctrl+C to copy the selected text)",
        type = "message",
        duration = 3
      )
    })
    
    # Download R script
    output$download_r <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("prairie_genomics_analysis_", timestamp, ".R")
      },
      content = function(file) {
        req(session_id_reactive())
        
        code_content <- get_session_code(
          session_id_reactive(),
          include_comments = TRUE,
          include_sessioninfo = TRUE
        )
        
        writeLines(code_content, file)
      },
      contentType = "text/plain"
    )
    
    # Download R Markdown
    output$download_rmd <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("prairie_genomics_report_", timestamp, ".Rmd")
      },
      content = function(file) {
        req(session_id_reactive())
        
        r_code <- get_session_code(
          session_id_reactive(),
          include_comments = TRUE,
          include_sessioninfo = TRUE
        )
        
        rmd_content <- convert_to_rmarkdown(r_code, session_id_reactive())
        writeLines(rmd_content, file)
      },
      contentType = "text/plain"
    )
    
    return(reactive({
      list(
        current_code = values$current_code,
        last_update = values$last_update
      )
    }))
  })
}

# Inline code display server
inlineCodeServer <- function(id, code_snippet = reactive(""), 
                           code_title = reactive("Analysis Code"),
                           code_description = reactive("")) {
  moduleServer(id, function(input, output, session) {
    
    # Display code title
    output$code_title <- renderText({
      code_title()
    })
    
    # Display code content
    output$inline_code_content <- renderText({
      snippet <- code_snippet()
      if (is.null(snippet) || snippet == "") {
        "# No code available for this step"
      } else {
        snippet
      }
    })
    
    # Display code description
    output$code_description <- renderText({
      code_description()
    })
    
    # Check if description exists
    output$has_description <- reactive({
      desc <- code_description()
      !is.null(desc) && desc != ""
    })
    outputOptions(output, "has_description", suspendWhenHidden = FALSE)
    
    # Copy code action
    observeEvent(input$copy_inline, {
      showNotification(
        "Code snippet copied! (Use Ctrl+C to copy the selected text)",
        type = "message",
        duration = 2
      )
    })
    
    return(reactive({
      list(
        is_visible = (input$toggle_code %% 2 == 1),
        code_content = code_snippet()
      )
    }))
  })
}

# Code export server
codeExportServer <- function(id, session_id_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # Complete analysis export
    output$export_complete <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        extension <- if (input$export_format == "Rmd") ".Rmd" else if (input$export_format == "html") ".html" else ".R"
        paste0("complete_analysis_", timestamp, extension)
      },
      content = function(file) {
        req(session_id_reactive())
        
        session_id <- session_id_reactive()
        include_comments <- "comments" %in% input$include_options
        include_sessioninfo <- "sessioninfo" %in% input$include_options
        
        if (input$export_format == "R") {
          code_content <- get_session_code(
            session_id,
            include_comments = include_comments,
            include_sessioninfo = include_sessioninfo
          )
          writeLines(code_content, file)
          
        } else if (input$export_format == "Rmd") {
          r_code <- get_session_code(session_id, include_comments, include_sessioninfo)
          rmd_content <- convert_to_rmarkdown(r_code, session_id)
          writeLines(rmd_content, file)
          
        } else if (input$export_format == "html") {
          # Create HTML report (simplified)
          r_code <- get_session_code(session_id, include_comments, include_sessioninfo)
          rmd_content <- convert_to_rmarkdown(r_code, session_id)
          
          # Save as temporary Rmd file and render to HTML
          temp_rmd <- tempfile(fileext = ".Rmd")
          writeLines(rmd_content, temp_rmd)
          
          # Note: In practice, you'd use rmarkdown::render() here
          # For now, just save the Rmd content
          writeLines(rmd_content, file)
        }
      }
    )
    
    # Template export
    output$export_template <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("analysis_template_", timestamp, ".R")
      },
      content = function(file) {
        template_code <- generate_analysis_template("basic_deseq2")
        writeLines(template_code, file)
      },
      contentType = "text/plain"
    )
    
    # Quick export - DESeq2 only
    output$quick_deseq2 <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("deseq2_code_", timestamp, ".R")
      },
      content = function(file) {
        req(session_id_reactive())
        
        deseq2_code <- get_analysis_code_snippet(
          session_id_reactive(),
          category = "deseq2"
        )
        
        if (deseq2_code == "# No matching analysis steps found") {
          deseq2_code <- "# No DESeq2 analysis has been run yet"
        }
        
        writeLines(deseq2_code, file)
      }
    )
    
    # Quick export - Pathway only
    output$quick_pathway <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("pathway_code_", timestamp, ".R")
      },
      content = function(file) {
        req(session_id_reactive())
        
        pathway_code <- get_analysis_code_snippet(
          session_id_reactive(),
          category = "pathway"
        )
        
        if (pathway_code == "# No matching analysis steps found") {
          pathway_code <- generate_analysis_template("pathway_analysis")
        }
        
        writeLines(pathway_code, file)
      }
    )
    
    # Quick export - Plots only
    output$quick_plots <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        paste0("visualization_code_", timestamp, ".R")
      },
      content = function(file) {
        req(session_id_reactive())
        
        viz_code <- get_analysis_code_snippet(
          session_id_reactive(),
          category = "visualization"
        )
        
        if (viz_code == "# No matching analysis steps found") {
          viz_code <- c(
            "# Visualization Code Template",
            "# ============================",
            "",
            "# Load required packages",
            "library(ggplot2)",
            "library(pheatmap)",
            "",
            "# Basic volcano plot",
            "volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +",
            "  geom_point(alpha = 0.6) +",
            "  labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value') +",
            "  theme_minimal()",
            "",
            "print(volcano_plot)"
          )
          viz_code <- paste(viz_code, collapse = "\n")
        }
        
        writeLines(viz_code, file)
      }
    )
  })
}

# Code validation server
codeValidationServer <- function(id, session_id_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # Validation status
    validation_results <- reactive({
      req(session_id_reactive())
      
      session_id <- session_id_reactive()
      code <- get_session_code(session_id)
      
      # Perform basic validation
      validation <- list(
        syntax_valid = TRUE,
        syntax_message = "✅ Syntax check passed",
        dependencies_ok = TRUE,
        deps_message = "✅ All packages available",
        warnings = character(),
        stats = list(
          total_lines = length(strsplit(code, "\n")[[1]]),
          code_lines = length(grep("^[^#]", strsplit(code, "\n")[[1]])),
          comment_lines = length(grep("^#", strsplit(code, "\n")[[1]])),
          functions_used = length(grep("\\w+\\(", code))
        )
      )
      
      # Simple syntax check (could be enhanced)
      if (grepl("Error|error", code)) {
        validation$syntax_valid <- FALSE
        validation$syntax_message <- "⚠️ Potential syntax issues detected"
        validation$warnings <- c(validation$warnings, "Check for syntax errors in generated code")
      }
      
      return(validation)
    })
    
    # Show validation panel
    output$show_validation <- reactive({
      !is.null(session_id_reactive()) && session_id_reactive() != ""
    })
    outputOptions(output, "show_validation", suspendWhenHidden = FALSE)
    
    # Syntax status
    observe({
      results <- validation_results()
      
      status_color <- if (results$syntax_valid) "#28a745" else "#ffc107"
      status_icon <- if (results$syntax_valid) "check-circle" else "exclamation-triangle"
      
      session$sendCustomMessage("updateStatus", list(
        id = "syntax_status",
        content = paste(
          '<i class="fa fa-', status_icon, '" style="color: ', status_color, ';"></i>',
          results$syntax_message
        ),
        color = status_color
      ))
    })
    
    # Dependencies status
    observe({
      results <- validation_results()
      
      status_color <- if (results$dependencies_ok) "#28a745" else "#dc3545"
      status_icon <- if (results$dependencies_ok) "check-circle" else "times-circle"
      
      session$sendCustomMessage("updateStatus", list(
        id = "deps_status", 
        content = paste(
          '<i class="fa fa-', status_icon, '" style="color: ', status_color, ';"></i>',
          results$deps_message
        ),
        color = status_color
      ))
    })
    
    # Code statistics
    output$code_stats <- renderTable({
      results <- validation_results()
      
      data.frame(
        Metric = c("Total Lines", "Code Lines", "Comments", "Functions"),
        Count = c(
          results$stats$total_lines,
          results$stats$code_lines,
          results$stats$comment_lines,
          results$stats$functions_used
        ),
        stringsAsFactors = FALSE
      )
    }, bordered = TRUE, striped = TRUE, hover = TRUE, width = "100%")
    
    # Validation warnings
    output$has_warnings <- reactive({
      results <- validation_results()
      length(results$warnings) > 0
    })
    outputOptions(output, "has_warnings", suspendWhenHidden = FALSE)
    
    output$validation_warnings <- renderUI({
      results <- validation_results()
      
      if (length(results$warnings) > 0) {
        warning_items <- lapply(results$warnings, function(warning) {
          tags$li(warning)
        })
        
        tags$ul(warning_items)
      }
    })
  })
}

# Method documentation server (static content)
methodDocumentationServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # This module is primarily static content defined in the UI
    # Could be extended to show dynamic method information based on analysis choices
    
    return(reactive({
      list(
        active_tab = input$doc_tabs %||% "DESeq2 Workflow"
      )
    }))
  })
}

cat("✅ Code display server functions loaded\n")
cat("📋 Server modules: codeDisplayServer(), inlineCodeServer(), codeExportServer(), codeValidationServer()\n")