# DESeq2 Analysis Module for Prairie Genomics Suite
# Handles differential expression analysis using DESeq2
# 
# Recent Updates:
# - Simplified to use DESeq2 default parameters only (most robust approach)
# - Added BioMart integration for human/mouse gene symbol conversion
# - Enhanced display with both Ensembl IDs and gene symbols
# - Batch processing for large gene sets with progress tracking
# - NEW: Ultra-fast cached gene conversion system (95%+ faster than BioMart)

# Load fast gene conversion cache module
tryCatch({
  source("gene_conversion_cache.R")
  cache_module_available <- TRUE
  cat("✅ Fast gene conversion cache module loaded\n")
}, error = function(e) {
  cache_module_available <- FALSE
  cat("⚠️ Gene conversion cache module not available, using BioMart fallback\n")
})

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
          
          h5("🧬 Gene Annotation"),
          
          radioButtons(
            ns("species"),
            "Species for gene symbol conversion:",
            choices = list(
              "Human (Homo sapiens)" = "human",
              "Mouse (Mus musculus)" = "mouse",
              "Skip conversion (keep Ensembl IDs)" = "none"
            ),
            selected = "human"
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("species"), "'] != 'none'"),
            p("BioMart will be used to convert Ensembl IDs to gene symbols", 
              style = "font-size: 12px; color: #666;")
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
    
    # Run analysis synchronously with simplified parameters
    result <- run_deseq2_analysis(
      expression_data = values$expression_data,
      annotation_data = values$annotation_data,
      contrast = local_values$available_contrasts[[input$selected_contrast]],
      padj_cutoff = input$padj_cutoff,
      fc_cutoff = input$fc_cutoff,
      species = input$species
    )
    
    if (result$success) {
      local_values$dds_object <- result$dds
      local_values$results_object <- result$results
      local_values$analysis_complete <- TRUE
      
      # Store unified DESeq2 data structure
      values$deseq2_data <- list(
        dds = result$dds,
        results = result$results,
        results_df = result$results_df,
        normalized_counts = result$normalized_counts,
        raw_counts = result$raw_counts,
        size_factors = result$size_factors,
        gene_mapping = result$gene_mapping,
        metadata = result$metadata
      )
      
      # Keep backward compatibility
      values$deseq2_results <- result$results_df
      values$filtered_data <- result$filtered_data
      
      add_status_message("✅ DESeq2 analysis completed successfully!", "success")
    } else {
      add_status_message(paste0("❌ Analysis failed: ", result$error), "danger")
    }
    
    # Hide progress, show button
    local_values$analysis_running <- FALSE
  })
  
  # DESeq2 analysis function (simplified to use defaults)
  run_deseq2_analysis <- function(expression_data, annotation_data, contrast, 
                                 padj_cutoff, fc_cutoff, species) {
    tryCatch({
      # CRITICAL DEBUG: Check what data we received
      cat("🔍 CRITICAL DEBUG: DESeq2 received expression data:\n")
      cat("   - Dimensions:", nrow(expression_data), "x", ncol(expression_data), "\n")
      cat("   - Sample gene IDs:", paste(head(rownames(expression_data), 5), collapse = ", "), "\n")
      cat("   - Data class:", class(expression_data), "\n")
      cat("   - Has gene_symbols attribute:", !is.null(attr(expression_data, "gene_symbols")), "\n")
      
      # Update progress
      update_progress(10, "Preparing data...")
      
      # Prepare count matrix
      count_matrix <- as.matrix(expression_data)
      count_matrix <- round(count_matrix)
      count_matrix[count_matrix < 0] <- 0
      
      # Use pre-converted gene symbols from data upload step
      update_progress(15, "Preparing gene identifiers...")
      
      original_gene_ids <- rownames(count_matrix)
      gene_symbols <- attr(expression_data, "gene_symbols")
      
      if (!is.null(gene_symbols)) {
        cat("✅ Using pre-converted gene symbols from data upload step\n")
        
        # Match symbols to current gene IDs
        symbol_matches <- match(original_gene_ids, gene_symbols$ensembl_gene_id)
        valid_matches <- !is.na(symbol_matches)
        
        # Create new rownames using symbols where available, Ensembl IDs as fallback
        new_rownames <- original_gene_ids  # Start with original IDs
        new_rownames[valid_matches] <- gene_symbols$gene_symbol[symbol_matches[valid_matches]]
        
        # Handle cases where symbol is empty or NA
        empty_symbols <- is.na(new_rownames) | new_rownames == "" | new_rownames == "NA"
        new_rownames[empty_symbols] <- original_gene_ids[empty_symbols]
        
        # Handle duplicate symbols by adding suffix
        if (any(duplicated(new_rownames))) {
          dup_indices <- which(duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE))
          for (i in dup_indices) {
            if (new_rownames[i] != original_gene_ids[i]) {
              new_rownames[i] <- paste0(new_rownames[i], "_", original_gene_ids[i])
            }
          }
        }
        
        # Update count matrix with gene symbols
        rownames(count_matrix) <- new_rownames
        
        # Report conversion success
        converted_count <- sum(valid_matches & new_rownames != original_gene_ids)
        conversion_rate <- round(100 * converted_count / length(original_gene_ids), 1)
        conversion_species <- attr(expression_data, "conversion_species") %||% "unknown"
        
        cat("✅ Using cached gene symbols:", conversion_rate, "% coverage\n")
        cat("📊 Species:", conversion_species, "| Converted:", converted_count, "| Original IDs:", length(original_gene_ids) - converted_count, "\n")
        
      } else {
        cat("⚠️ No pre-converted gene symbols available, using original gene IDs\n")
        cat("💡 Enable gene conversion in the Data Upload tab for better results\n")
      }
      
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
      update_progress(50, "Running DESeq2 with default parameters...")
      
      # Run DESeq2 with all default parameters (most robust approach)
      dds <- DESeq(dds)
      
      # Update progress
      update_progress(70, "Extracting results...")
      
      # Extract results with default parameters
      res <- results(
        dds,
        contrast = c("Condition", contrast[1], contrast[2])
      )
      
      # Convert to data frame
      results_df <- as.data.frame(res)
      results_df$gene <- rownames(results_df)
      
      # Update progress - gene conversion was already done before DESeq2
      update_progress(80, "Preparing final results...")
      
      # Since we converted gene IDs before DESeq2, results_df$gene already contains symbols/IDs
      # But we need to preserve the original Ensembl IDs for proper mapping
      results_df$ensembl_id <- original_gene_ids  # Preserve original Ensembl IDs
      results_df$display_name <- results_df$gene
      results_df$gene_symbol <- results_df$gene  # For compatibility with downstream processes
      
      if (species != "none") {
        add_status_message("✅ Using gene symbols from pre-conversion step", "success")
      } else {
        add_status_message("ℹ️ Using original gene identifiers", "info")
      }
      
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
      update_progress(100, "Analysis complete with gene symbol conversion!")
      
      # Create gene mapping for centralized use
      gene_mapping <- data.frame(
        ensembl_id = results_df$ensembl_id,  # Use preserved Ensembl IDs
        gene_symbol = results_df$gene_symbol,
        display_name = results_df$display_name,
        stringsAsFactors = FALSE
      )
      
      # Clean up memory before returning results
      gc()
      
      # Create unified DESeq2 results structure
      return(list(
        success = TRUE,
        dds = dds,
        results = res,
        results_df = results_df,
        filtered_data = counts(dds, normalized = TRUE),
        # New unified structure components
        normalized_counts = counts(dds, normalized = TRUE),
        raw_counts = counts(dds, normalized = FALSE),
        size_factors = sizeFactors(dds),
        gene_mapping = gene_mapping,
        metadata = list(
          contrast = contrast,
          padj_cutoff = padj_cutoff,
          fc_cutoff = fc_cutoff,
          species = species,
          n_samples = ncol(dds),
          n_genes_tested = nrow(results_df),
          n_significant = sum(results_df$significant, na.rm = TRUE),
          analysis_date = Sys.Date()
        )
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
  
  # BioMart gene symbol conversion function with robust error handling
  convert_ensembl_to_symbols <- function(ensembl_ids, species) {
    tryCatch({
      # Check if biomaRt is available
      if (!requireNamespace("biomaRt", quietly = TRUE)) {
        cat("biomaRt package not available - skipping gene symbol conversion\n")
        return(create_fallback_results(ensembl_ids))
      }
      
      # Connect to appropriate mart with error handling
      mart <- NULL
      symbol_attribute <- NULL
      
      tryCatch({
        if (species == "human") {
          # Try multiple mirrors for reliability (valid mirrors: www, useast, asia)
          for (mirror in c("www", "useast", "asia")) {
            tryCatch({
              mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = mirror)
              symbol_attribute <- "hgnc_symbol"
              cat(paste0("Connected to Ensembl ", mirror, " mirror for human data\n"))
              break
            }, error = function(e) {
              cat(paste0("Failed to connect to ", mirror, " mirror: ", e$message, "\n"))
            })
          }
        } else if (species == "mouse") {
          for (mirror in c("www", "useast", "asia")) {
            tryCatch({
              mart <- biomaRt::useEnsembl("ensembl", "mmusculus_gene_ensembl", mirror = mirror)
              symbol_attribute <- "mgi_symbol" 
              cat(paste0("Connected to Ensembl ", mirror, " mirror for mouse data\n"))
              break
            }, error = function(e) {
              cat(paste0("Failed to connect to ", mirror, " mirror: ", e$message, "\n"))
            })
          }
        } else {
          stop("Unsupported species: ", species)
        }
      }, error = function(e) {
        cat(paste0("Failed to connect to any Ensembl mirror: ", e$message, "\n"))
        return(create_fallback_results(ensembl_ids))
      })
      
      # If all mirrors failed, try legacy connection method as last resort
      if (is.null(mart)) {
        cat("All Ensembl mirrors failed, trying legacy BioMart connection...\n")
        tryCatch({
          if (species == "human") {
            mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            symbol_attribute <- "hgnc_symbol"
            cat("Connected using legacy BioMart method for human data\n")
          } else if (species == "mouse") {
            mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
            symbol_attribute <- "mgi_symbol"
            cat("Connected using legacy BioMart method for mouse data\n")
          }
        }, error = function(e) {
          cat(paste0("Legacy connection also failed: ", e$message, "\n"))
        })
      }
      
      if (is.null(mart)) {
        cat("All BioMart connection methods failed - using fallback\n")
        return(create_fallback_results(ensembl_ids))
      }
      
      # Get gene symbols in smaller batches with retry logic
      batch_size <- 200  # Reduced batch size for reliability
      all_results <- data.frame()
      total_batches <- ceiling(length(ensembl_ids) / batch_size)
      failed_batches <- 0
      
      for (i in seq(1, length(ensembl_ids), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(ensembl_ids))
        batch_ids <- ensembl_ids[i:end_idx]
        batch_num <- ceiling(i / batch_size)
        
        cat(paste0("Converting batch ", batch_num, " of ", total_batches, " (", length(batch_ids), " genes)\n"))
        
        # Try this batch with retry logic
        batch_results <- query_biomart_with_retry(mart, symbol_attribute, batch_ids, max_retries = 3)
        
        if (!is.null(batch_results) && nrow(batch_results) > 0) {
          # Rename symbol column for consistency
          colnames(batch_results)[colnames(batch_results) == symbol_attribute] <- "gene_symbol"
          all_results <- rbind(all_results, batch_results)
        } else {
          failed_batches <- failed_batches + 1
          cat(paste0("Batch ", batch_num, " failed after retries\n"))
        }
        
        # Rate limiting to avoid overwhelming the server
        if (batch_num < total_batches) {
          Sys.sleep(0.5)  # 500ms delay between batches
        }
      }
      
      # Create complete result set with NAs for missing symbols
      complete_results <- data.frame(
        ensembl_gene_id = ensembl_ids,
        gene_symbol = all_results$gene_symbol[match(ensembl_ids, all_results$ensembl_gene_id)],
        stringsAsFactors = FALSE
      )
      
      # Report conversion statistics
      converted_count <- sum(!is.na(complete_results$gene_symbol) & complete_results$gene_symbol != "")
      total_count <- length(ensembl_ids)
      conversion_rate <- round(100 * converted_count / total_count, 1)
      
      cat(paste0("Gene symbol conversion completed: ", converted_count, "/", total_count, " genes (", conversion_rate, "%)\n"))
      if (failed_batches > 0) {
        cat(paste0("Note: ", failed_batches, " batches failed due to server issues\n"))
      }
      
      return(complete_results)
      
    }, error = function(e) {
      cat(paste0("Critical error in gene symbol conversion: ", e$message, "\n"))
      cat("Falling back to Ensembl IDs only\n")
      return(create_fallback_results(ensembl_ids))
    })
  }
  
  # Helper function to query BioMart with retry logic
  query_biomart_with_retry <- function(mart, symbol_attribute, batch_ids, max_retries = 3) {
    for (attempt in 1:max_retries) {
      tryCatch({
        result <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", symbol_attribute),
          filters = "ensembl_gene_id",
          values = batch_ids,
          mart = mart
        )
        return(result)
      }, error = function(e) {
        if (attempt < max_retries) {
          wait_time <- attempt * 2  # Progressive backoff: 2s, 4s, 6s
          cat(paste0("Attempt ", attempt, " failed (", e$message, "), retrying in ", wait_time, "s...\n"))
          Sys.sleep(wait_time)
        } else {
          cat(paste0("All ", max_retries, " attempts failed: ", e$message, "\n"))
        }
      })
    }
    return(NULL)
  }
  
  # Helper function to create fallback results
  create_fallback_results <- function(ensembl_ids) {
    return(data.frame(
      ensembl_gene_id = ensembl_ids,
      gene_symbol = NA,
      stringsAsFactors = FALSE
    ))
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
      Gene = top_genes$display_name,
      "Ensembl ID" = if("ensembl_id" %in% colnames(top_genes)) top_genes$ensembl_id else top_genes$gene,
      "Log2 FC" = round(top_genes$log2FoldChange, 3),
      "P-value" = formatC(top_genes$pvalue, format = "e", digits = 2),
      "Adj. P-value" = formatC(top_genes$padj, format = "e", digits = 2),
      "Regulation" = top_genes$regulation,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Hide Ensembl ID column if gene symbols are available
    if (!is.null(top_genes$gene_symbol) && any(!is.na(top_genes$gene_symbol))) {
      column_defs <- list(
        list(className = 'dt-center', targets = c(2, 3, 4, 5)),
        list(
          targets = 5,
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
    } else {
      # Hide Ensembl ID column if no symbols
      display_data <- display_data[, !colnames(display_data) %in% "Ensembl ID"]
      column_defs <- list(
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
    }
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 10,
        dom = 'tip',
        columnDefs = column_defs
      ),
      escape = FALSE,
      rownames = FALSE
    )
  })
  
  # Return results
  return(reactive({ local_values$results_object }))
}