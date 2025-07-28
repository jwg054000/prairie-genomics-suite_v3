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
        
        br(),
        
        div(
          class = "alert alert-info",
          style = "font-size: 12px; padding: 8px;",
          h6("📋 Optional: pAnno/Annotation File"),
          p("Upload a separate annotation file (pAnno format) to automatically assign sample groups and metadata.", style = "margin: 0;")
        ),
        
        fileInput(
          ns("annotation_file"),
          "Choose Sample Annotation File (Optional)",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          width = "100%",
          placeholder = "pAnno file with sample metadata"
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
        
        div(
          class = "alert alert-info",
          style = "padding: 10px; margin-bottom: 15px;",
          h6("🧹 Gene Filtering Settings", style = "margin-top: 0;"),
          p("These settings remove low-count genes before analysis. This is normal RNA-seq preprocessing.", style = "margin-bottom: 5px; font-size: 12px;"),
          p(strong("To keep ALL genes: Set 'Minimum read count' to 0"), style = "margin-bottom: 0; font-size: 12px; color: #d63384;")
        ),
        
        checkboxInput(
          ns("disable_filtering"),
          "🚫 Disable all gene filtering (keep all genes)",
          value = FALSE
        ),
        
        numericInput(
          ns("min_count"),
          "Minimum read count per gene:",
          value = 1,  # Changed from 10 to 1 - much less aggressive
          min = 0,
          max = 1000,
          step = 1
        ),
        
        numericInput(
          ns("min_samples"),
          "Minimum samples with counts:",
          value = 1,  # Changed from 3 to 1 - much less aggressive  
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
        
        # Gene conversion options
        h5("🧬 Gene ID Conversion"),
        
        checkboxInput(
          ns("convert_gene_ids"),
          "Convert gene IDs to symbols",
          value = TRUE
        ),
        
        conditionalPanel(
          condition = paste0("input['\", ns(\"convert_gene_ids\"), \"'] == true"),
          selectInput(
            ns("gene_species"),
            "Species (for conversion):",
            choices = list(
              "Auto-detect" = "auto",
              "Human (Homo sapiens)" = "human",
              "Mouse (Mus musculus)" = "mouse"
            ),
            selected = "auto"
          ),
          
          div(
            class = "alert alert-info",
            style = "padding: 8px; margin-top: 10px;",
            tags$small("🔄 Gene conversion uses cached results for speed. First conversion may take longer.")
          )
        ),
        
        br(),
        
        actionButton(
          ns("process_data"),
          "🚀 Process Data & Convert Genes",
          class = "btn-primary btn-lg",
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
    ),
    
    br(),
    
    # Gene filtering diagnostics
    conditionalPanel(
      condition = paste0("output['", ns("show_gene_diagnostics"), "']"),
      fluidRow(
        column(
          12,
          div(
            class = "alert alert-info",
            style = "margin-bottom: 15px;",
            h5("🔍 Gene Filtering Diagnostics", style = "margin-top: 0;"),
            htmlOutput(ns("gene_diagnostics_display"))
          )
        )
      ),
      
      # Advanced filtering controls
      fluidRow(
        column(
          3,
          h6("🎛️ Advanced Controls"),
          numericInput(
            ns("min_expression_threshold"),
            "Min expression:",
            value = 1,
            min = 0,
            step = 0.1
          )
        ),
        column(
          3,
          br(),
          numericInput(
            ns("min_samples_expressed"),
            "Min samples:",
            value = 2,
            min = 1,
            step = 1
          )
        ),
        column(
          3,
          br(),
          actionButton(
            ns("apply_advanced_filtering"),
            "Apply Filters",
            class = "btn-warning btn-sm",
            style = "margin-top: 25px;"
          )
        ),
        column(
          3,
          br(),
          actionButton(
            ns("reset_filtering"),
            "Reset",
            class = "btn-secondary btn-sm",
            style = "margin-top: 25px;"
          )
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
    annotation_data = NULL,
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
      
      # Handle gene names and duplicates
      if (input$has_rownames) {
        gene_names <- raw_data[[1]]
        data_only <- raw_data[, -1]
        
        # Check for duplicate gene IDs
        if (any(duplicated(gene_names))) {
          n_duplicates <- sum(duplicated(gene_names))
          add_status_message(
            paste0("⚠️ Found ", n_duplicates, " duplicate gene IDs. Aggregating by summing counts..."), 
            "warning"
          )
          
          # Convert expression columns to numeric first
          data_only <- data_only %>%
            mutate_all(~ as.numeric(as.character(.)))
          
          # Combine gene names with data for aggregation
          combined_data <- cbind(gene_id = gene_names, data_only)
          
          # Aggregate duplicates by summing
          if (requireNamespace("dplyr", quietly = TRUE)) {
            aggregated <- combined_data %>%
              group_by(gene_id) %>%
              summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
            
            gene_names <- aggregated$gene_id
            raw_data <- aggregated[, -1]  # Remove gene_id column
          } else {
            # Fallback without dplyr
            unique_genes <- unique(gene_names)
            aggregated_data <- matrix(0, nrow = length(unique_genes), ncol = ncol(data_only))
            
            for (i in seq_along(unique_genes)) {
              gene <- unique_genes[i]
              matching_rows <- which(gene_names == gene)
              if (length(matching_rows) == 1) {
                aggregated_data[i, ] <- as.numeric(data_only[matching_rows, ])
              } else {
                aggregated_data[i, ] <- colSums(data_only[matching_rows, ], na.rm = TRUE)
              }
            }
            
            raw_data <- as.data.frame(aggregated_data)
            colnames(raw_data) <- colnames(data_only)
            gene_names <- unique_genes
          }
          
          add_status_message(
            paste0("✅ Aggregated duplicates. Reduced from ", length(raw_data[[1]]), 
                   " to ", length(gene_names), " unique genes"), 
            "success"
          )
        } else {
          raw_data <- data_only
        }
        
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
      
      # DEBUG: Check what data was actually loaded
      cat("🔍 DEBUG: File upload completed\n")
      cat("   - File name:", input$expression_file$name, "\n")
      cat("   - Raw data dimensions:", nrow(raw_data), "x", ncol(raw_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(raw_data), 3), collapse = ", "), "\n")
      
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
  
  # Annotation file upload observer
  observeEvent(input$annotation_file, {
    req(input$annotation_file)
    
    tryCatch({
      add_status_message("📋 Loading annotation file...", "info")
      
      file_path <- input$annotation_file$datapath
      file_ext <- tools::file_ext(input$annotation_file$name)
      
      # Read annotation file based on extension
      annotation_data <- NULL
      if (file_ext == "csv") {
        annotation_data <- readr::read_csv(file_path, show_col_types = FALSE)
      } else if (file_ext %in% c("tsv", "txt")) {
        annotation_data <- readr::read_delim(file_path, delim = "\t", show_col_types = FALSE)
      } else if (file_ext %in% c("xlsx", "xls") && exists("excel_support") && excel_support) {
        annotation_data <- readxl::read_excel(file_path)
      } else {
        stop("Unsupported annotation file format or readxl not available")
      }
      
      # Validate annotation data structure
      annotation_result <- validate_panno_file(annotation_data)
      
      if (annotation_result$valid) {
        local_values$annotation_data <- annotation_result$data
        add_status_message(
          paste0("✅ Annotation file loaded: ", nrow(annotation_result$data), 
                " samples with ", length(annotation_result$groups), " groups detected"),
          "success"
        )
        
        # Automatically trigger sample matching if expression data is available
        if (!is.null(local_values$processed_data)) {
          match_samples_to_annotation()
        }
        
      } else {
        add_status_message(paste0("⚠️ Annotation file warning: ", annotation_result$message), "warning")
        local_values$annotation_data <- annotation_result$data  # Store anyway for manual review
      }
      
    }, error = function(e) {
      add_status_message(paste0("❌ Error loading annotation file: ", e$message), "danger")
    })
  })
  
  # Process data when button is clicked
  observeEvent(input$process_data, {
    req(local_values$raw_data)
    
    tryCatch({
      add_status_message("🔄 Processing data with enhanced chunked processing...", "info")
      
      # DEBUG: Check raw data before processing
      cat("🔍 DEBUG: Starting data processing\n")
      cat("   - Raw data dimensions:", nrow(local_values$raw_data), "x", ncol(local_values$raw_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(local_values$raw_data), 3), collapse = ", "), "\n")
      
      # Use enhanced chunked processing for large datasets
      if (exists("process_large_dataset") && is.function(process_large_dataset)) {
        cat("🔍 DEBUG: Using process_large_dataset function\n")
        processed_data <- process_large_dataset(local_values$raw_data)
        cat("🔍 DEBUG: After process_large_dataset - dimensions:", nrow(processed_data), "x", ncol(processed_data), "\n")
      } else {
        cat("🔍 DEBUG: Using raw data directly (no process_large_dataset)\n")
        # Fallback to basic processing if function not available
        processed_data <- local_values$raw_data
      }
      
      # Filter by minimum count with detailed logging
      if (!input$disable_filtering && input$min_count > 0) {
        genes_before_filter <- nrow(processed_data)
        
        # Count samples with sufficient reads per gene
        sufficient_samples <- rowSums(processed_data >= input$min_count)
        keep_genes <- sufficient_samples >= input$min_samples
        
        filtered_genes <- sum(!keep_genes)
        processed_data <- processed_data[keep_genes, ]
        
        # DETAILED DEBUG: Show filtering impact
        cat("🔍 DEBUG: Gene filtering applied\n")
        cat("   - Before filtering:", genes_before_filter, "genes\n")
        cat("   - Filter criteria: >=", input$min_count, "counts in >=", input$min_samples, "samples\n")
        cat("   - Genes removed:", filtered_genes, "\n")
        cat("   - Genes remaining:", nrow(processed_data), "\n")
        cat("   - Filtering percentage:", round(100 * filtered_genes / genes_before_filter, 1), "%\n")
        
        if (filtered_genes > 0) {
          add_status_message(
            paste0("🧹 Gene Filtering: ", genes_before_filter, " → ", nrow(processed_data), 
                   " genes (removed ", filtered_genes, " low-count genes with criteria: ≥", 
                   input$min_count, " counts in ≥", input$min_samples, " samples)"), 
            "warning"
          )
        }
      } else {
        if (input$disable_filtering) {
          cat("🔍 DEBUG: Gene filtering DISABLED by user - keeping all", nrow(processed_data), "genes\n")
          add_status_message(
            paste0("✅ Filtering disabled - keeping all ", nrow(processed_data), " genes"), 
            "success"
          )
        } else {
          cat("🔍 DEBUG: No gene filtering applied (min_count = 0)\n")
        }
      }
      
      # Ensure all values are positive for DESeq2
      processed_data[processed_data < 0] <- 0
      processed_data <- round(processed_data)
      
      # Gene conversion step
      gene_conversion_result <- NULL
      if (input$convert_gene_ids) {
        add_status_message("🧬 Converting gene IDs to symbols...", "info")
        
        # Use simple local conversion (fast and reliable)
        tryCatch({
          # Load simple gene conversion functions
          if (!exists("convert_genes_simple")) {
            source("simple_gene_conversion.R")
          }
          
          gene_ids <- rownames(processed_data)
          conversion_species <- input$gene_species %||% "auto"
          
          # CRITICAL FIX: Auto-detect species from gene IDs
          if (conversion_species == "auto") {
            human_count <- sum(grepl("^ENSG[0-9]", gene_ids))
            mouse_count <- sum(grepl("^ENSMUSG[0-9]", gene_ids))
            
            if (mouse_count > human_count) {
              conversion_species <- "mouse"
              cat("🔍 AUTO-DETECT: Detected mouse genes (ENSMUSG) - using mouse conversion\n")
            } else if (human_count > 0) {
              conversion_species <- "human"
              cat("🔍 AUTO-DETECT: Detected human genes (ENSG) - using human conversion\n")
            } else {
              conversion_species <- "human"  # Default fallback
              cat("🔍 AUTO-DETECT: No clear pattern - defaulting to human conversion\n")
            }
          }
          
          # DEBUG: Check what genes are being sent for conversion
          cat("🔍 DEBUG: About to convert", length(gene_ids), "genes\n")
          cat("   - First 3 genes:", paste(head(gene_ids, 3), collapse = ", "), "\n")
          cat("   - Species:", conversion_species, "(auto-detected)\n")
          
          # Try simple local conversion first (much faster and more reliable)
          gene_conversion_result <- convert_genes_simple(
            gene_ids, 
            species = conversion_species
          )
          
          # If simple conversion failed completely, try the cache system as fallback
          if (is.null(gene_conversion_result) || 
              sum(!is.na(gene_conversion_result$gene_symbol)) == 0) {
            
            cat("⚠️ Simple conversion failed, trying cache system...\n")
            add_status_message("⚠️ Local conversion failed, trying BioMart cache...", "warning")
            
            # Load cache system as fallback
            if (!exists("convert_genes_fast")) {
              source("gene_conversion_cache.R")
            }
            
            gene_conversion_result <- convert_genes_fast(
              gene_ids, 
              species = conversion_species,
              use_cache = TRUE
            )
          }
          
          if (!is.null(gene_conversion_result)) {
            # Calculate conversion stats
            converted_count <- sum(!is.na(gene_conversion_result$gene_symbol) & 
                                 gene_conversion_result$gene_symbol != "")
            conversion_rate <- round(100 * converted_count / length(gene_ids), 1)
            
            add_status_message(
              paste0("✅ Gene conversion completed: ", conversion_rate, "% success rate (", 
                     converted_count, "/", length(gene_ids), " genes)"),
              "success"
            )
            
            # Store conversion results as an attribute
            attr(processed_data, "gene_symbols") <- gene_conversion_result
            attr(processed_data, "conversion_species") <- conversion_species
            attr(processed_data, "conversion_rate") <- conversion_rate
          } else {
            add_status_message("⚠️ Gene conversion failed, using original IDs", "warning")
          }
          
        }, error = function(e) {
          add_status_message(
            paste0("⚠️ Gene conversion error: ", e$message, ". Using original IDs."),
            "warning"
          )
        })
      }
      
      # CRITICAL SAFEGUARD: Reject 15,357 gene dataset (known cache issue)
      if (nrow(processed_data) == 15357) {
        add_status_message(
          "🚨 CRITICAL: Detected problematic 15,357 gene cached dataset! This will be rejected to prevent data override issues.", 
          "danger"
        )
        cat("🚨 CRITICAL: 15,357 gene dataset rejected - likely cached mouse data\n")
        cat("   - Gene IDs: ", paste(head(rownames(processed_data), 5), collapse = ", "), "\n")
        return()
      }
      
      # Generate gene diagnostics
      local_values$gene_diagnostics <- create_gene_diagnostics(processed_data, gene_conversion_result)
      
      local_values$processed_data <- processed_data
      values$expression_data <- processed_data
      
      # Update final message to include conversion info
      final_message <- paste0("✅ Data processing complete! Final dataset: ", 
                             nrow(processed_data), " genes × ", ncol(processed_data), " samples")
      
      if (!is.null(gene_conversion_result) && input$convert_gene_ids) {
        conversion_rate <- attr(processed_data, "conversion_rate") %||% 0
        final_message <- paste0(final_message, " | Gene conversion: ", conversion_rate, "% success")
      }
      
      add_status_message(final_message, "success")
      
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
  
  # pAnno file validation function
  validate_panno_file <- function(annotation_data) {
    tryCatch({
      # Check basic structure
      if (is.null(annotation_data) || nrow(annotation_data) == 0) {
        return(list(valid = FALSE, message = "Annotation file is empty", data = NULL))
      }
      
      # Check for common pAnno/annotation file columns
      required_cols <- NULL
      condition_col <- NULL
      sample_col <- NULL
      
      # Look for sample identifier column (various common names)
      sample_candidates <- c("Sample", "SampleID", "sample", "sampleid", "Sample_ID", "ID", "Name")
      for (col in sample_candidates) {
        if (col %in% colnames(annotation_data)) {
          sample_col <- col
          break
        }
      }
      
      # Look for condition/group column
      condition_candidates <- c("Condition", "Group", "Treatment", "condition", "group", "treatment", 
                               "Phenotype", "Class", "Type", "phenotype", "class", "type")
      for (col in condition_candidates) {
        if (col %in% colnames(annotation_data)) {
          condition_col <- col
          break
        }
      }
      
      # If we can't find obvious columns, try to infer from data
      if (is.null(sample_col) || is.null(condition_col)) {
        # Use first column as sample ID if it looks like sample names
        if (ncol(annotation_data) >= 2) {
          sample_col <- colnames(annotation_data)[1]
          condition_col <- colnames(annotation_data)[2]
        } else {
          return(list(valid = FALSE, message = "Cannot identify sample and condition columns", data = annotation_data))
        }
      }
      
      # Standardize column names
      standardized_data <- annotation_data
      colnames(standardized_data)[colnames(standardized_data) == sample_col] <- "Sample"
      colnames(standardized_data)[colnames(standardized_data) == condition_col] <- "Condition"
      
      # Check for reasonable number of groups
      unique_groups <- unique(standardized_data$Condition)
      if (length(unique_groups) < 2) {
        return(list(valid = FALSE, message = "Need at least 2 different groups/conditions", data = standardized_data))
      }
      
      if (length(unique_groups) > 20) {
        return(list(valid = FALSE, message = "Too many groups (>20) - check data format", data = standardized_data))
      }
      
      # Check for reasonable sample sizes per group
      group_counts <- table(standardized_data$Condition)
      min_group_size <- min(group_counts)
      if (min_group_size < 2) {
        return(list(
          valid = FALSE, 
          message = paste0("Some groups have <2 samples. Minimum group size: ", min_group_size),
          data = standardized_data
        ))
      }
      
      return(list(
        valid = TRUE, 
        message = "Annotation file validated successfully",
        data = standardized_data,
        groups = unique_groups,
        sample_col = "Sample",
        condition_col = "Condition"
      ))
      
    }, error = function(e) {
      return(list(valid = FALSE, message = paste("Validation error:", e$message), data = annotation_data))
    })
  }
  
  # Sample matching function
  match_samples_to_annotation <- function() {
    if (is.null(local_values$processed_data) || is.null(local_values$annotation_data)) {
      return()
    }
    
    tryCatch({
      expression_samples <- colnames(local_values$processed_data)
      annotation_samples <- local_values$annotation_data$Sample
      
      # Try exact matching first
      exact_matches <- intersect(expression_samples, annotation_samples)
      
      if (length(exact_matches) == length(expression_samples)) {
        add_status_message(
          paste0("🎯 Perfect match! All ", length(exact_matches), " samples matched between expression and annotation data"),
          "success"
        )
        
        # Store matched annotation for automatic use
        matched_annotation <- local_values$annotation_data[local_values$annotation_data$Sample %in% exact_matches, ]
        values$annotation_data <- matched_annotation
        
        return()
      }
      
      # Try fuzzy matching for common naming variations
      fuzzy_matches <- try_fuzzy_sample_matching(expression_samples, annotation_samples)
      
      if (fuzzy_matches$success) {
        add_status_message(
          paste0("🔧 Smart matching successful! Matched ", length(fuzzy_matches$matches), 
                "/", length(expression_samples), " samples with name variations"),
          "success"
        )
        
        # Store matched annotation
        matched_annotation <- fuzzy_matches$matched_data
        values$annotation_data <- matched_annotation
        
      } else {
        add_status_message(
          paste0("⚠️ Partial match: ", length(exact_matches), "/", length(expression_samples), 
                " samples matched. Check sample naming consistency."),
          "warning"
        )
      }
      
    }, error = function(e) {
      add_status_message(paste0("❌ Error matching samples: ", e$message), "danger")
    })
  }
  
  # Fuzzy sample matching function
  try_fuzzy_sample_matching <- function(expression_samples, annotation_samples) {
    tryCatch({
      matches <- list()
      
      for (exp_sample in expression_samples) {
        best_match <- NULL
        best_score <- 0
        
        for (anno_sample in annotation_samples) {
          # Try various matching strategies
          score <- 0
          
          # Exact match
          if (exp_sample == anno_sample) {
            score <- 1.0
          }
          # Case insensitive
          else if (tolower(exp_sample) == tolower(anno_sample)) {
            score <- 0.95
          }
          # Remove common separators and compare
          else if (gsub("[._-]", "", tolower(exp_sample)) == gsub("[._-]", "", tolower(anno_sample))) {
            score <- 0.9
          }
          # Substring matching (annotation contains expression name or vice versa)
          else if (grepl(tolower(exp_sample), tolower(anno_sample)) || grepl(tolower(anno_sample), tolower(exp_sample))) {
            score <- 0.8
          }
          
          if (score > best_score) {
            best_score <- score
            best_match <- anno_sample
          }
        }
        
        if (best_score >= 0.8) {  # Minimum confidence threshold
          matches[[exp_sample]] <- best_match
        }
      }
      
      if (length(matches) >= length(expression_samples) * 0.8) {  # 80% match rate
        # Create matched annotation data
        matched_data <- local_values$annotation_data[local_values$annotation_data$Sample %in% unlist(matches), ]
        
        return(list(
          success = TRUE,
          matches = matches,
          matched_data = matched_data
        ))
      } else {
        return(list(success = FALSE, matches = matches))
      }
      
    }, error = function(e) {
      return(list(success = FALSE, error = e$message))
    })
  }
  
  # Gene diagnostics function
  create_gene_diagnostics <- function(data_matrix, conversion_results = NULL) {
    tryCatch({
      diagnostics <- list()
      
      # Basic data stats
      diagnostics$original_genes <- nrow(data_matrix)
      diagnostics$original_samples <- ncol(data_matrix)
      
      # Gene ID analysis
      gene_ids <- rownames(data_matrix)
      diagnostics$human_genes <- sum(grepl("^ENSG[0-9]", gene_ids))
      diagnostics$mouse_genes <- sum(grepl("^ENSMUSG[0-9]", gene_ids))
      diagnostics$other_genes <- length(gene_ids) - diagnostics$human_genes - diagnostics$mouse_genes
      
      # Conversion stats
      if (!is.null(conversion_results)) {
        converted_count <- sum(!is.na(conversion_results$gene_symbol) & 
                             conversion_results$gene_symbol != "" & 
                             conversion_results$gene_symbol != conversion_results$ensembl_gene_id)
        diagnostics$conversion_rate <- round(100 * converted_count / length(gene_ids), 1)
        diagnostics$converted_genes <- converted_count
      }
      
      # Expression level analysis
      gene_means <- rowMeans(data_matrix, na.rm = TRUE)
      diagnostics$zero_genes <- sum(gene_means == 0)
      diagnostics$low_expression_genes <- sum(gene_means < 10)
      diagnostics$high_expression_genes <- sum(gene_means > 1000)
      
      # Sample correlation
      if (ncol(data_matrix) > 1) {
        tryCatch({
          sample_cors <- cor(log2(data_matrix + 1), use = "complete.obs")
          diagnostics$min_correlation <- round(min(sample_cors[upper.tri(sample_cors)]), 3)
          diagnostics$median_correlation <- round(median(sample_cors[upper.tri(sample_cors)]), 3)
        }, error = function(e) {
          diagnostics$min_correlation <- NA
          diagnostics$median_correlation <- NA
        })
      }
      
      return(diagnostics)
      
    }, error = function(e) {
      return(list(error = e$message))
    })
  }
  
  # Format diagnostics as HTML
  format_diagnostics_html <- function(diagnostics) {
    if ("error" %in% names(diagnostics)) {
      return(paste0("<span style='color: red;'>❌ Error: ", diagnostics$error, "</span>"))
    }
    
    html_parts <- c()
    
    # Gene composition
    html_parts <- c(html_parts, "<strong>📊 Gene Composition:</strong><br>")
    html_parts <- c(html_parts, paste0("🧬 Total genes: ", formatC(diagnostics$original_genes, format="d", big.mark=","), "<br>"))
    
    if (diagnostics$human_genes > 0) {
      html_parts <- c(html_parts, paste0("👤 Human genes (ENSG): ", diagnostics$human_genes, "<br>"))
    }
    if (diagnostics$mouse_genes > 0) {
      html_parts <- c(html_parts, paste0("🐭 Mouse genes (ENSMUSG): ", diagnostics$mouse_genes, "<br>"))
    }
    if (diagnostics$other_genes > 0) {
      html_parts <- c(html_parts, paste0("❓ Other/symbol genes: ", diagnostics$other_genes, "<br>"))
    }
    
    # Conversion results
    if ("conversion_rate" %in% names(diagnostics)) {
      conversion_color <- if (diagnostics$conversion_rate > 80) "green" else if (diagnostics$conversion_rate > 50) "orange" else "red"
      html_parts <- c(html_parts, paste0("<span style='color: ", conversion_color, ";'>🔄 Conversion rate: ", 
                                       diagnostics$conversion_rate, "% (", diagnostics$converted_genes, " genes)</span><br>"))
    }
    
    # Expression levels
    html_parts <- c(html_parts, "<br><strong>📈 Expression Levels:</strong><br>")
    html_parts <- c(html_parts, paste0("🚫 Zero expression: ", diagnostics$zero_genes, " genes<br>"))
    html_parts <- c(html_parts, paste0("📉 Low expression (<10): ", diagnostics$low_expression_genes, " genes<br>"))
    html_parts <- c(html_parts, paste0("📈 High expression (>1000): ", diagnostics$high_expression_genes, " genes<br>"))
    
    # Sample quality
    if ("median_correlation" %in% names(diagnostics) && !is.na(diagnostics$median_correlation)) {
      cor_color <- if (diagnostics$median_correlation > 0.8) "green" else if (diagnostics$median_correlation > 0.6) "orange" else "red"
      html_parts <- c(html_parts, paste0("<br><strong>📊 Sample Quality:</strong><br>"))
      html_parts <- c(html_parts, paste0("<span style='color: ", cor_color, ";'>📊 Median correlation: ", 
                                       diagnostics$median_correlation, "</span><br>"))
    }
    
    return(paste(html_parts, collapse = ""))
  }
  
  # Gene diagnostics output
  output$gene_diagnostics_display <- renderUI({
    if (!is.null(local_values$gene_diagnostics)) {
      HTML(format_diagnostics_html(local_values$gene_diagnostics))
    } else {
      HTML("<em>Process data to see gene diagnostics</em>")
    }
  })
  
  # Show diagnostics panel
  output$show_gene_diagnostics <- reactive({
    !is.null(local_values$gene_diagnostics)
  })
  outputOptions(output, "show_gene_diagnostics", suspendWhenHidden = FALSE)
  
  # Advanced filtering handler
  observeEvent(input$apply_advanced_filtering, {
    req(local_values$processed_data)
    
    tryCatch({
      add_status_message("🎛️ Applying advanced gene filtering...", "info")
      
      original_data <- local_values$processed_data
      
      # Apply expression threshold
      if (input$min_expression_threshold > 0) {
        gene_means <- rowMeans(original_data, na.rm = TRUE)
        keep_expression <- gene_means >= input$min_expression_threshold
      } else {
        keep_expression <- rep(TRUE, nrow(original_data))
      }
      
      # Apply samples threshold
      if (input$min_samples_expressed > 0) {
        samples_expressed <- rowSums(original_data > 0, na.rm = TRUE)
        keep_samples <- samples_expressed >= input$min_samples_expressed
      } else {
        keep_samples <- rep(TRUE, nrow(original_data))
      }
      
      # Combine filters
      keep_genes <- keep_expression & keep_samples
      filtered_data <- original_data[keep_genes, ]
      
      genes_removed <- sum(!keep_genes)
      
      if (genes_removed > 0) {
        # Update processed data
        local_values$processed_data <- filtered_data
        values$expression_data <- filtered_data
        
        # Update diagnostics
        gene_conversion_result <- attr(original_data, "gene_symbols")
        if (!is.null(gene_conversion_result)) {
          # Filter conversion results to match filtered genes
          filtered_gene_ids <- rownames(filtered_data)
          gene_conversion_result <- gene_conversion_result[
            gene_conversion_result$ensembl_gene_id %in% filtered_gene_ids, 
          ]
          attr(filtered_data, "gene_symbols") <- gene_conversion_result
        }
        
        local_values$gene_diagnostics <- create_gene_diagnostics(filtered_data, gene_conversion_result)
        
        add_status_message(
          paste0("✅ Advanced filtering applied: removed ", genes_removed, " genes. ", 
                 nrow(filtered_data), " genes remaining."),
          "success"
        )
      } else {
        add_status_message("ℹ️ No genes removed by current filter settings", "info")
      }
      
    }, error = function(e) {
      add_status_message(paste0("❌ Advanced filtering error: ", e$message), "danger")
    })
  })
  
  # Reset filtering handler
  observeEvent(input$reset_filtering, {
    req(local_values$raw_data)
    
    add_status_message("🔄 Resetting to original processed data...", "info")
    
    # Re-trigger processing with original data
    updateNumericInput(session, "min_expression_threshold", value = 1)
    updateNumericInput(session, "min_samples_expressed", value = 2)
    
    # You could also trigger the original process_data logic here
    add_status_message("✅ Filter settings reset. Click 'Process Data' to reprocess.", "success")
  })
  
  # Return processed data and annotation data
  return(list(
    expression_data = reactive({ local_values$processed_data }),
    annotation_data = reactive({ local_values$annotation_data })
  ))
}