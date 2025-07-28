# Bulletproof Data Tester - Handles any file gracefully
# Maximum error protection and graceful degradation

library(shiny)
library(DT)

# Bulletproof file loading
load_expression_safe <- function(file_path) {
  tryCatch({
    # Try to read the file
    data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Basic validation
    if (nrow(data) == 0) {
      return(list(error = "File is empty"))
    }
    
    if (ncol(data) < 2) {
      return(list(error = "File must have at least 2 columns"))
    }
    
    # Check if first column looks like gene names
    first_col <- data[, 1]
    if (is.character(first_col) || is.factor(first_col)) {
      # First column is gene names
      gene_names <- as.character(first_col)
      
      # Get numeric columns
      numeric_cols <- data[, -1, drop = FALSE]
      
      # Try to convert to numeric
      for (i in 1:ncol(numeric_cols)) {
        if (!is.numeric(numeric_cols[, i])) {
          numeric_cols[, i] <- as.numeric(as.character(numeric_cols[, i]))
        }
      }
      
      # Check for conversion issues
      if (any(is.na(numeric_cols) & !is.na(data[, -1]))) {
        return(list(error = "Some data could not be converted to numbers"))
      }
      
      # Create matrix
      data_matrix <- as.matrix(numeric_cols)
      rownames(data_matrix) <- make.unique(gene_names)
      
      return(list(
        data = data_matrix,
        genes = nrow(data_matrix),
        samples = ncol(data_matrix),
        sample_names = colnames(data_matrix),
        gene_names = rownames(data_matrix)[1:min(5, nrow(data_matrix))]
      ))
      
    } else {
      # All columns are numeric
      data_matrix <- as.matrix(data)
      rownames(data_matrix) <- paste0("Gene_", 1:nrow(data_matrix))
      
      return(list(
        data = data_matrix,
        genes = nrow(data_matrix),
        samples = ncol(data_matrix),
        sample_names = colnames(data_matrix),
        gene_names = paste0("Gene_", 1:min(5, nrow(data_matrix)))
      ))
    }
    
  }, error = function(e) {
    return(list(error = paste("Failed to load file:", e$message)))
  })
}

load_metadata_safe <- function(file_path) {
  tryCatch({
    # Try to read the file
    data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Basic validation
    if (nrow(data) == 0) {
      return(list(error = "Metadata file is empty"))
    }
    
    if (ncol(data) == 0) {
      return(list(error = "Metadata file has no columns"))
    }
    
    return(list(
      data = data,
      samples = nrow(data),
      columns = colnames(data),
      sample_ids = if ("Sample_ID" %in% colnames(data)) data$Sample_ID else NULL
    ))
    
  }, error = function(e) {
    return(list(error = paste("Failed to load metadata:", e$message)))
  })
}

# Ultra-safe UI
ui <- fluidPage(
  titlePanel("🔧 Bulletproof Data Tester"),
  
  div(class = "alert alert-info",
      "This version handles any file gracefully and won't crash. Upload your files to test!"),
  
  br(),
  
  fluidRow(
    column(6,
      wellPanel(
        h3("📁 File Upload"),
        
        fileInput("expr_file", 
                 "Expression Matrix:",
                 accept = c(".csv", ".tsv", ".txt")),
        
        fileInput("meta_file", 
                 "Sample Metadata:",
                 accept = c(".csv", ".tsv", ".txt")),
        
        hr(),
        
        h4("📊 Status"),
        verbatimTextOutput("status", placeholder = TRUE)
      )
    ),
    
    column(6,
      wellPanel(
        h3("ℹ️ File Information"), 
        verbatimTextOutput("info", placeholder = TRUE)
      )
    )
  ),
  
  fluidRow(
    column(12,
      wellPanel(
        h3("👀 Data Preview"),
        
        conditionalPanel(
          condition = "output.has_expr_data",
          h4("Expression Data (First 5 genes, 5 samples)"),
          DT::dataTableOutput("expr_preview")
        ),
        
        conditionalPanel(
          condition = "output.has_meta_data", 
          h4("Metadata"),
          DT::dataTableOutput("meta_preview")
        ),
        
        conditionalPanel(
          condition = "!output.has_expr_data && !output.has_meta_data",
          div(class = "alert alert-info", "Upload files to see preview")
        )
      )
    )
  )
)

# Bulletproof server
server <- function(input, output, session) {
  
  # Safe reactive storage
  values <- reactiveValues(
    expr_result = NULL,
    meta_result = NULL,
    messages = c("Ready to upload files...")
  )
  
  # Expression file handler
  observeEvent(input$expr_file, {
    req(input$expr_file)
    
    values$messages <- c(values$messages, 
                        paste("📊 Loading expression file:", input$expr_file$name))
    
    # Load with bulletproof function
    result <- load_expression_safe(input$expr_file$datapath)
    values$expr_result <- result
    
    if (!is.null(result$error)) {
      values$messages <- c(values$messages, paste("❌", result$error))
    } else {
      values$messages <- c(values$messages, 
                          paste("✅ Expression data loaded:", result$genes, "genes ×", result$samples, "samples"))
    }
  })
  
  # Metadata file handler  
  observeEvent(input$meta_file, {
    req(input$meta_file)
    
    values$messages <- c(values$messages, 
                        paste("📋 Loading metadata file:", input$meta_file$name))
    
    # Load with bulletproof function
    result <- load_metadata_safe(input$meta_file$datapath)
    values$meta_result <- result
    
    if (!is.null(result$error)) {
      values$messages <- c(values$messages, paste("❌", result$error))
    } else {
      values$messages <- c(values$messages, 
                          paste("✅ Metadata loaded:", result$samples, "samples with", length(result$columns), "columns"))
    }
  })
  
  # Status display
  output$status <- renderText({
    paste(tail(values$messages, 10), collapse = "\n")
  })
  
  # Information display
  output$info <- renderText({
    info_lines <- c()
    
    # Expression info
    if (!is.null(values$expr_result) && is.null(values$expr_result$error)) {
      expr <- values$expr_result
      info_lines <- c(info_lines, 
                     "=== EXPRESSION MATRIX ===",
                     paste("Dimensions:", expr$genes, "genes ×", expr$samples, "samples"),
                     paste("Sample names:", paste(head(expr$sample_names, 3), collapse = ", "), "..."),
                     paste("First genes:", paste(head(expr$gene_names, 3), collapse = ", "), "..."),
                     "")
    }
    
    # Metadata info
    if (!is.null(values$meta_result) && is.null(values$meta_result$error)) {
      meta <- values$meta_result
      info_lines <- c(info_lines,
                     "=== METADATA ===", 
                     paste("Samples:", meta$samples),
                     paste("Columns:", paste(meta$columns, collapse = ", ")),
                     "")
      
      if (!is.null(meta$sample_ids)) {
        info_lines <- c(info_lines,
                       paste("Sample IDs:", paste(head(meta$sample_ids, 3), collapse = ", "), "..."))
      }
    }
    
    if (length(info_lines) == 0) {
      "No files loaded yet."
    } else {
      paste(info_lines, collapse = "\n")
    }
  })
  
  # Check if data exists for conditional panels
  output$has_expr_data <- reactive({
    !is.null(values$expr_result) && is.null(values$expr_result$error)
  })
  
  output$has_meta_data <- reactive({
    !is.null(values$meta_result) && is.null(values$meta_result$error)
  })
  
  outputOptions(output, "has_expr_data", suspendWhenHidden = FALSE)
  outputOptions(output, "has_meta_data", suspendWhenHidden = FALSE)
  
  # Expression preview table
  output$expr_preview <- DT::renderDataTable({
    req(values$expr_result)
    req(is.null(values$expr_result$error))
    
    # Show only first 5 genes and 5 samples for safety
    data_matrix <- values$expr_result$data
    preview_matrix <- data_matrix[1:min(5, nrow(data_matrix)), 
                                 1:min(5, ncol(data_matrix))]
    
    # Convert to data frame
    preview_df <- data.frame(
      Gene = rownames(preview_matrix),
      preview_matrix,
      check.names = FALSE
    )
    
    DT::datatable(preview_df,
                  options = list(dom = 't', scrollX = TRUE),
                  rownames = FALSE)
  })
  
  # Metadata preview table
  output$meta_preview <- DT::renderDataTable({
    req(values$meta_result)
    req(is.null(values$meta_result$error))
    
    DT::datatable(values$meta_result$data,
                  options = list(dom = 'tp', scrollX = TRUE, pageLength = 15),
                  rownames = FALSE)
  })
}

# Run the bulletproof app
shinyApp(ui = ui, server = server)