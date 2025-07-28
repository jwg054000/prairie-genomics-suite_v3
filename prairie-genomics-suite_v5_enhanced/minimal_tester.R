# Minimal Data Tester - Ultra-simple version for stability
# Just basic file upload and display

library(shiny)
library(DT)
library(readr)

# Minimal file loading functions
load_csv_simple <- function(file_path) {
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
    
    # If first column is character, assume it's gene names
    if (is.character(data[[1]])) {
      gene_names <- data[[1]]
      numeric_data <- data[, -1]
      data_matrix <- as.matrix(numeric_data)
      rownames(data_matrix) <- gene_names
      return(data_matrix)
    } else {
      return(as.matrix(data))
    }
  }, error = function(e) {
    stop(paste("Error loading file:", e$message))
  })
}

# Ultra-minimal UI
ui <- fluidPage(
  titlePanel("🧬 Minimal Data Tester"),
  
  br(),
  
  fluidRow(
    column(6,
      wellPanel(
        h3("📊 Upload Files"),
        
        fileInput("expr_file", "Expression Matrix (CSV):",
                 accept = ".csv"),
        
        fileInput("meta_file", "Sample Metadata (CSV):",
                 accept = ".csv"),
        
        hr(),
        
        verbatimTextOutput("status")
      )
    ),
    
    column(6,
      wellPanel(
        h3("📋 File Info"),
        verbatimTextOutput("file_info")
      )
    )
  ),
  
  fluidRow(
    column(12,
      wellPanel(
        h3("🔍 Data Preview"),
        tabsetPanel(
          tabPanel("Expression Data", 
                  br(),
                  DT::dataTableOutput("expr_table")),
          tabPanel("Metadata",
                  br(), 
                  DT::dataTableOutput("meta_table"))
        )
      )
    )
  )
)

# Ultra-minimal server
server <- function(input, output, session) {
  
  # Storage
  values <- reactiveValues(
    expr_data = NULL,
    meta_data = NULL,
    messages = c()
  )
  
  # Expression file
  observeEvent(input$expr_file, {
    req(input$expr_file)
    
    values$messages <- c(values$messages, paste("Loading expression file:", input$expr_file$name))
    
    tryCatch({
      expr_data <- load_csv_simple(input$expr_file$datapath)
      values$expr_data <- expr_data
      values$messages <- c(values$messages, 
                          paste("✅ Loaded:", nrow(expr_data), "genes ×", ncol(expr_data), "samples"))
      
    }, error = function(e) {
      values$messages <- c(values$messages, paste("❌ Error:", e$message))
    })
  })
  
  # Metadata file
  observeEvent(input$meta_file, {
    req(input$meta_file)
    
    values$messages <- c(values$messages, paste("Loading metadata file:", input$meta_file$name))
    
    tryCatch({
      meta_data <- read_csv(input$meta_file$datapath, show_col_types = FALSE)
      values$meta_data <- meta_data
      values$messages <- c(values$messages, 
                          paste("✅ Loaded:", nrow(meta_data), "samples with", ncol(meta_data), "columns"))
      
    }, error = function(e) {
      values$messages <- c(values$messages, paste("❌ Error:", e$message))
    })
  })
  
  # Status output
  output$status <- renderText({
    paste(values$messages, collapse = "\n")
  })
  
  # File info
  output$file_info <- renderText({
    info <- c()
    
    if (!is.null(values$expr_data)) {
      info <- c(info, paste("Expression Matrix:"))
      info <- c(info, paste("- Dimensions:", nrow(values$expr_data), "×", ncol(values$expr_data)))
      info <- c(info, paste("- First gene:", rownames(values$expr_data)[1]))
      info <- c(info, paste("- Sample names:", paste(head(colnames(values$expr_data)), collapse = ", ")))
      info <- c(info, "")
    }
    
    if (!is.null(values$meta_data)) {
      info <- c(info, paste("Sample Metadata:"))
      info <- c(info, paste("- Samples:", nrow(values$meta_data)))
      info <- c(info, paste("- Columns:", paste(colnames(values$meta_data), collapse = ", ")))
      
      if ("Sample_ID" %in% colnames(values$meta_data)) {
        info <- c(info, paste("- Sample IDs:", paste(head(values$meta_data$Sample_ID), collapse = ", ")))
      }
    }
    
    paste(info, collapse = "\n")
  })
  
  # Expression table
  output$expr_table <- DT::renderDataTable({
    req(values$expr_data)
    
    # Show only first 10 genes and 6 samples to avoid memory issues
    preview_data <- values$expr_data[1:min(10, nrow(values$expr_data)), 
                                    1:min(6, ncol(values$expr_data))]
    
    # Convert to data frame with gene names
    df <- data.frame(
      Gene = rownames(preview_data),
      preview_data,
      check.names = FALSE
    )
    
    DT::datatable(df, 
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  # Metadata table
  output$meta_table <- DT::renderDataTable({
    req(values$meta_data)
    
    DT::datatable(values$meta_data,
                  options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE)
  })
}

# Run app
shinyApp(ui = ui, server = server)