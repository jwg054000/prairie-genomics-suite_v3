# DESeq2 Analysis Module for Prairie Genomics Suite
# Handles differential expression analysis using DESeq2

# UI Function
deseq2AnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        4,
        h4("Analysis Settings"),
        
        wellPanel(
          h5("🎯 Comparison Settings"),
          
          uiOutput(ns("contrast_selection")),
          
          br(),
          
          h5("📊 Filtering Parameters"),
          
          numericInput(
            ns("padj_cutoff"),
            "Adjusted p-value cutoff:",
            value = 0.05,
            min = 0.001,
            max = 0.1,
            step = 0.001
          ),
          
          numericInput(
            ns("fc_cutoff"), 
            "Log2 fold change cutoff:",
            value = 1.0,
            min = 0.1,
            max = 5.0,
            step = 0.1
          ),
          
          br(),
          
          h5("⚙️ Advanced Options"),
          
          checkboxInput(
            ns("independent_filtering"),
            "Independent filtering",
            value = TRUE
          ),
          
          selectInput(
            ns("fit_type"),
            "Fit type:",
            choices = list(
              "parametric" = "parametric",
              "local" = "local", 
              "mean" = "mean"
            ),
            selected = "parametric"
          ),
          
          selectInput(
            ns("test"),
            "Statistical test:",
            choices = list(
              "Wald" = "Wald",
              "LRT" = "LRT"
            ),
            selected = "Wald"
          )
        )
      ),
      
      column(
        8,
        h4("Analysis Status"),
        
        # Analysis button
        conditionalPanel(
          condition = paste0("!output['", ns("analysis_running"), "']"),
          actionButton(
            ns("run_analysis"),
            "🚀 Run DESeq2 Analysis",
            class = "btn-primary btn-lg",
            style = "width: 100%; height: 60px; font-size: 18px;"
          )
        ),
        
        # Analysis progress
        conditionalPanel(
          condition = paste0("output['", ns("analysis_running"), "']"),
          
          h5("Running DESeq2 Analysis..."),
          # Simple progress display without shinyWidgets
          div(
            style = "background-color: #3c8dbc; color: white; padding: 8px; border-radius: 4px; text-align: center; margin: 10px 0;",
            id = ns("analysis_progress"),
            "Progress: 0%"
          ),
          
          div(id = ns("progress_text"), "Initializing...")
        ),
        
        br(), br(),
        
        # Results summary
        conditionalPanel(
          condition = paste0("output['", ns("show_results"), "']"),
          
          wellPanel(
            h4("📊 Analysis Results Summary"),
            
            fluidRow(
              column(
                3,
                valueBoxOutput(ns("total_genes_box"), width = 12)
              ),
              column(
                3,
                valueBoxOutput(ns("significant_genes_box"), width = 12)
              ),
              column(
                3,
                valueBoxOutput(ns("upregulated_box"), width = 12)
              ),
              column(
                3,
                valueBoxOutput(ns("downregulated_box"), width = 12)
              )
            ),
            
            br(),
            
            h5("Top Significant Genes:"),
            DT::dataTableOutput(ns("top_genes_table"))
          )
        )
      )
    ),
    
    br(),
    
    # Status messages
    div(id = ns("status_messages"))
  )
}

# Server Function
deseq2Analysis <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    dds_object = NULL,
    results_object = NULL,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    status_messages = list(),
    available_contrasts = NULL
  )
  
  # Update available contrasts when annotation data changes
  observe({
    if (!is.null(values$annotation_data)) {
      conditions <- unique(values$annotation_data$Condition)
      
      if (length(conditions) >= 2) {
        # Create all possible pairwise contrasts
        contrasts <- combn(conditions, 2, simplify = FALSE)
        contrast_names <- sapply(contrasts, function(x) paste(x[1], "vs", x[2]))
        names(contrasts) <- contrast_names
        
        local_values$available_contrasts <- contrasts
      }
    }
  })
  
  # Render contrast selection
  output$contrast_selection <- renderUI({
    req(local_values$available_contrasts)
    
    contrast_choices <- names(local_values$available_contrasts)
    
    selectInput(
      ns("selected_contrast"),
      "Select comparison:",
      choices = contrast_choices,
      selected = contrast_choices[1]
    )
  })
  
  # Run DESeq2 analysis
  observeEvent(input$run_analysis, {
    req(values$expression_data, values$annotation_data)
    req(input$selected_contrast)
    
    # Check if DESeq2 is available
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      add_status_message("❌ DESeq2 package not available. This feature requires Bioconductor installation.", "danger")
      return()
    }
    
    # Show progress
    local_values$analysis_running <- TRUE
    
    # Run analysis synchronously (simplified for debugging)
    result <- run_deseq2_analysis(
      expression_data = values$expression_data,
      annotation_data = values$annotation_data,
      contrast = local_values$available_contrasts[[input$selected_contrast]],
      padj_cutoff = input$padj_cutoff,
      fc_cutoff = input$fc_cutoff,
      independent_filtering = input$independent_filtering,
      fit_type = input$fit_type,
      test = input$test
    )
    
    if (result$success) {
      local_values$dds_object <- result$dds
      local_values$results_object <- result$results
      local_values$analysis_complete <- TRUE
      values$deseq2_results <- result$results_df
      values$filtered_data <- result$filtered_data
      
      add_status_message("✅ DESeq2 analysis completed successfully!", "success")
    } else {
      add_status_message(paste0("❌ Analysis failed: ", result$error), "danger")
    }
    
    # Hide progress, show button
    local_values$analysis_running <- FALSE
  })
  
  # DESeq2 analysis function
  run_deseq2_analysis <- function(expression_data, annotation_data, contrast, 
                                 padj_cutoff, fc_cutoff, independent_filtering, 
                                 fit_type, test) {
    tryCatch({
      # Update progress
      update_progress(10, "Preparing data...")
      
      # Prepare count matrix
      count_matrix <- as.matrix(expression_data)
      count_matrix <- round(count_matrix)
      count_matrix[count_matrix < 0] <- 0
      
      # Prepare sample data
      sample_data <- annotation_data
      rownames(sample_data) <- sample_data$Sample
      sample_data <- sample_data[colnames(count_matrix), , drop = FALSE]
      
      # Update progress
      update_progress(25, "Creating DESeq2 object...")
      
      # Create DESeqDataSet
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_data,
        design = ~ Condition
      )
      
      # Filter low count genes
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]
      
      # Update progress
      update_progress(50, "Running DESeq2...")
      
      # Run DESeq2
      dds <- DESeq(
        dds,
        fitType = fit_type,
        test = test
      )
      
      # Update progress
      update_progress(75, "Extracting results...")
      
      # Extract results
      res <- results(
        dds,
        contrast = c("Condition", contrast[1], contrast[2]),
        independentFiltering = independent_filtering
      )
      
      # Convert to data frame
      results_df <- as.data.frame(res)
      results_df$gene <- rownames(results_df)
      
      # Add significance flags
      results_df$significant <- !is.na(results_df$padj) & 
                              results_df$padj < padj_cutoff & 
                              abs(results_df$log2FoldChange) > fc_cutoff
      
      results_df$regulation <- ifelse(
        results_df$significant,
        ifelse(results_df$log2FoldChange > 0, "Up", "Down"),
        "NS"
      )
      
      # Sort by adjusted p-value
      results_df <- results_df[order(results_df$padj, na.last = TRUE), ]
      
      # Update progress
      update_progress(100, "Analysis complete!")
      
      return(list(
        success = TRUE,
        dds = dds,
        results = res,
        results_df = results_df,
        filtered_data = counts(dds, normalized = TRUE)
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        error = e$message
      ))
    })
  }
  
  # Helper function to update progress (simplified for shinyapps.io)
  update_progress <- function(value, text) {
    # Simple console output instead of complex progress updates
    cat(paste0("Progress: ", value, "% - ", text, "\n"))
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
  
  # Render status messages
  output[[paste0(ns(""), "status_messages")]] <- renderUI({
    if (length(local_values$status_messages) > 0) {
      message_divs <- lapply(local_values$status_messages, function(msg) {
        div(
          class = paste("alert", msg$type),
          style = "margin-bottom: 5px; padding: 8px 12px;",
          tags$small(paste0("[", msg$timestamp, "] ")),
          msg$message
        )
      })
      do.call(tagList, message_divs)
    }
  })
  
  # Control results visibility
  output$show_results <- reactive({
    local_values$analysis_complete && !is.null(values$deseq2_results)
  })
  outputOptions(output, "show_results", suspendWhenHidden = FALSE)
  
  # Control analysis running state
  output$analysis_running <- reactive({
    local_values$analysis_running
  })
  outputOptions(output, "analysis_running", suspendWhenHidden = FALSE)
  
  # Results summary boxes
  output$total_genes_box <- renderValueBox({
    total_genes <- if (!is.null(values$deseq2_results)) {
      nrow(values$deseq2_results)
    } else 0
    
    valueBox(
      value = formatC(total_genes, format = "d", big.mark = ","),
      subtitle = "Genes Tested",
      icon = icon("dna"),
      color = "blue",
      width = 12
    )
  })
  
  output$significant_genes_box <- renderValueBox({
    significant_genes <- if (!is.null(values$deseq2_results)) {
      sum(values$deseq2_results$significant, na.rm = TRUE)
    } else 0
    
    valueBox(
      value = formatC(significant_genes, format = "d", big.mark = ","),
      subtitle = "Significant",
      icon = icon("star"),
      color = "yellow",
      width = 12
    )
  })
  
  output$upregulated_box <- renderValueBox({
    upregulated <- if (!is.null(values$deseq2_results)) {
      sum(values$deseq2_results$regulation == "Up", na.rm = TRUE)
    } else 0
    
    valueBox(
      value = formatC(upregulated, format = "d", big.mark = ","),
      subtitle = "Upregulated",
      icon = icon("arrow-up"),
      color = "green",
      width = 12
    )
  })
  
  output$downregulated_box <- renderValueBox({
    downregulated <- if (!is.null(values$deseq2_results)) {
      sum(values$deseq2_results$regulation == "Down", na.rm = TRUE)
    } else 0
    
    valueBox(
      value = formatC(downregulated, format = "d", big.mark = ","),
      subtitle = "Downregulated", 
      icon = icon("arrow-down"),
      color = "red",
      width = 12
    )
  })
  
  # Top genes table
  output$top_genes_table <- DT::renderDataTable({
    req(values$deseq2_results)
    
    # Get top 20 significant genes
    top_genes <- values$deseq2_results[values$deseq2_results$significant == TRUE, ]
    top_genes <- head(top_genes, 20)
    
    # Select and format columns
    display_data <- data.frame(
      Gene = top_genes$gene,
      "Log2 FC" = round(top_genes$log2FoldChange, 3),
      "P-value" = formatC(top_genes$pvalue, format = "e", digits = 2),
      "Adj. P-value" = formatC(top_genes$padj, format = "e", digits = 2),
      "Regulation" = top_genes$regulation,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 10,
        dom = 'tip',
        columnDefs = list(
          list(className = 'dt-center', targets = c(1, 2, 3, 4)),
          list(
            targets = 4,
            render = JS(
              "function(data, type, row, meta) {",
              "  if(data == 'Up') {",
              "    return '<span class=\"badge badge-success\">Up</span>';",
              "  } else if(data == 'Down') {",
              "    return '<span class=\"badge badge-danger\">Down</span>';", 
              "  } else {",
              "    return '<span class=\"badge badge-secondary\">NS</span>';",
              "  }",
              "}"
            )
          )
        )
      ),
      escape = FALSE,
      rownames = FALSE
    )
  })
  
  # Return results
  return(reactive({ local_values$results_object }))
}