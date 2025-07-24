# Data Upload Module for Prairie Genomics Suite
# Handles file upload and data preprocessing

# UI Function
dataUploadUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        h4("Upload Expression Data"),
        
        fileInput(
          ns("expression_file"),
          "Choose Expression Data File",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          width = "100%"
        ),
        
        radioButtons(
          ns("file_type"),
          "File Type:",
          choices = if (exists("excel_support") && excel_support) {
            list(
              "Auto-detect" = "auto",
              "CSV (comma-separated)" = "csv", 
              "TSV (tab-separated)" = "tsv",
              "Excel" = "excel"
            )
          } else {
            list(
              "Auto-detect" = "auto",
              "CSV (comma-separated)" = "csv", 
              "TSV (tab-separated)" = "tsv"
            )
          },
          selected = "auto",
          inline = TRUE
        ),
        
        checkboxInput(
          ns("has_header"),
          "File has header row",
          value = TRUE
        ),
        
        checkboxInput(
          ns("has_rownames"),
          "First column contains gene names",
          value = TRUE
        )
      ),
      
      column(
        6,
        h4("Data Processing Options"),
        
        numericInput(
          ns("min_count"),
          "Minimum read count per gene:",
          value = 10,
          min = 1,
          max = 1000,
          step = 1
        ),
        
        numericInput(
          ns("min_samples"),
          "Minimum samples with counts:",
          value = 3,
          min = 1,
          max = 50,
          step = 1
        ),
        
        checkboxInput(
          ns("log_transform"),
          "Apply log2 transformation for preview",
          value = FALSE
        ),
        
        br(),
        
        actionButton(
          ns("process_data"),
          "Process Data",
          class = "btn-primary",
          style = "width: 100%;"
        )
      )
    ),
    
    br(),
    
    # Status messages
    conditionalPanel(
      condition = paste0("output['", ns("show_status"), "']"),
      div(id = ns("status_messages"))
    ),
    
    br(),
    
    # Data summary
    conditionalPanel(
      condition = paste0("output['", ns("show_summary"), "']"),
      fluidRow(
        column(
          3,
          valueBoxOutput(ns("n_genes"), width = 12)
        ),
        column(
          3, 
          valueBoxOutput(ns("n_samples"), width = 12)
        ),
        column(
          3,
          valueBoxOutput(ns("total_counts"), width = 12)
        ),
        column(
          3,
          valueBoxOutput(ns("median_counts"), width = 12)
        )
      )
    )
  )
}

# Server Function
dataUpload <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    raw_data = NULL,
    processed_data = NULL,
    file_info = NULL,
    status_messages = list()
  )
  
  # File upload observer
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    tryCatch({
      # Add status message
      add_status_message("📁 Reading file...", "info")
      
      file_path <- input$expression_file$datapath
      file_ext <- tools::file_ext(input$expression_file$name)
      
      # Determine file type
      file_type <- input$file_type
      if (file_type == "auto") {
        file_type <- switch(
          tolower(file_ext),
          "csv" = "csv",
          "tsv" = "tsv", 
          "txt" = "tsv",
          "xlsx" = "excel",
          "xls" = "excel",
          "csv"  # default
        )
      }
      
      # Read data based on file type
      if (file_type == "excel") {
        if (exists("excel_support") && excel_support && requireNamespace("readxl", quietly = TRUE)) {
          raw_data <- readxl::read_excel(
            file_path,
            col_names = input$has_header
          )
        } else {
          add_status_message("❌ Excel support not available. Please use CSV or TSV format.", "danger")
          return()
        }
      } else if (file_type == "csv") {
        raw_data <- readr::read_csv(
          file_path,
          col_names = input$has_header,
          show_col_types = FALSE
        )
      } else {  # tsv
        raw_data <- readr::read_tsv(
          file_path,
          col_names = input$has_header,
          show_col_types = FALSE
        )
      }
      
      # Handle gene names
      if (input$has_rownames) {
        gene_names <- raw_data[[1]]
        raw_data <- raw_data[, -1]
        raw_data <- as.data.frame(raw_data)
        rownames(raw_data) <- gene_names
      } else {
        raw_data <- as.data.frame(raw_data)
        rownames(raw_data) <- paste0("Gene_", 1:nrow(raw_data))
      }
      
      # Convert to numeric
      raw_data <- raw_data %>%
        mutate_all(~ as.numeric(as.character(.)))
      
      # Remove rows with all NA
      raw_data <- raw_data[rowSums(is.na(raw_data)) < ncol(raw_data), ]
      
      # Replace remaining NA with 0
      raw_data[is.na(raw_data)] <- 0
      
      local_values$raw_data <- raw_data
      local_values$file_info <- list(
        name = input$expression_file$name,
        size = input$expression_file$size,
        type = file_type
      )
      
      add_status_message(
        paste0("✅ Successfully loaded ", nrow(raw_data), " genes and ", 
               ncol(raw_data), " samples"), 
        "success"
      )
      
    }, error = function(e) {
      add_status_message(
        paste0("❌ Error reading file: ", e$message), 
        "danger"
      )
    })
  })
  
  # Process data when button is clicked
  observeEvent(input$process_data, {
    req(local_values$raw_data)
    
    tryCatch({
      add_status_message("🔄 Processing data...", "info")
      
      processed_data <- local_values$raw_data
      
      # Filter by minimum count
      if (input$min_count > 0) {
        # Count samples with sufficient reads per gene
        sufficient_samples <- rowSums(processed_data >= input$min_count)
        keep_genes <- sufficient_samples >= input$min_samples
        
        filtered_genes <- sum(!keep_genes)
        processed_data <- processed_data[keep_genes, ]
        
        if (filtered_genes > 0) {
          add_status_message(
            paste0("🧹 Filtered out ", filtered_genes, " low-count genes"), 
            "info"
          )
        }
      }
      
      # Ensure all values are positive for DESeq2
      processed_data[processed_data < 0] <- 0
      processed_data <- round(processed_data)
      
      local_values$processed_data <- processed_data
      values$expression_data <- processed_data
      
      add_status_message(
        paste0("✅ Data processing complete! Final dataset: ", 
               nrow(processed_data), " genes × ", ncol(processed_data), " samples"),
        "success"
      )
      
    }, error = function(e) {
      add_status_message(
        paste0("❌ Error processing data: ", e$message),
        "danger"
      )
    })
  })
  
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
  
  # Control visibility of status and summary sections
  output$show_status <- reactive({
    length(local_values$status_messages) > 0
  })
  outputOptions(output, "show_status", suspendWhenHidden = FALSE)
  
  output$show_summary <- reactive({
    !is.null(local_values$processed_data)
  })
  outputOptions(output, "show_summary", suspendWhenHidden = FALSE)
  
  # Summary value boxes
  output$n_genes <- renderValueBox({
    n_genes <- if (!is.null(local_values$processed_data)) nrow(local_values$processed_data) else 0
    
    valueBox(
      value = formatC(n_genes, format = "d", big.mark = ","),
      subtitle = "Genes",
      icon = icon("dna"),
      color = "green",
      width = 12
    )
  })
  
  output$n_samples <- renderValueBox({
    n_samples <- if (!is.null(local_values$processed_data)) ncol(local_values$processed_data) else 0
    
    valueBox(
      value = n_samples,
      subtitle = "Samples",
      icon = icon("vials"),
      color = "blue", 
      width = 12
    )
  })
  
  output$total_counts <- renderValueBox({
    total_counts <- if (!is.null(local_values$processed_data)) {
      formatC(sum(local_values$processed_data), format = "d", big.mark = ",")
    } else "0"
    
    valueBox(
      value = total_counts,
      subtitle = "Total Counts",
      icon = icon("calculator"),
      color = "purple",
      width = 12
    )
  })
  
  output$median_counts <- renderValueBox({
    median_counts <- if (!is.null(local_values$processed_data)) {
      formatC(median(colSums(local_values$processed_data)), format = "d", big.mark = ",")
    } else "0"
    
    valueBox(
      value = median_counts,
      subtitle = "Median Library Size",
      icon = icon("chart-line"),
      color = "orange",
      width = 12
    )
  })
  
  # Return processed data
  return(reactive({ local_values$processed_data }))
}