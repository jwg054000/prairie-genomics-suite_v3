# Phase 1 - Async Handlers for Prairie Genomics Suite
# Promises-based async processing to eliminate UI blocking
# 
# This module implements async handlers for computationally intensive tasks:
# - DESeq2 differential expression analysis
# - Gene symbol conversion and annotation
# - Pathway enrichment analysis
# - Large dataset processing and visualization

library(promises)
library(future)
library(callr)

# Configure async processing
plan(multisession, workers = 2)  # Use 2 background processes

# Async DESeq2 Analysis Handler
async_deseq2_analysis <- function(expression_data, annotation_data, contrast, 
                                 padj_cutoff = 0.05, fc_cutoff = 1.0, species = "human") {
  
  # Create promise for async DESeq2 processing
  future_promise({
    # Load required libraries in background process
    library(DESeq2)
    library(biomaRt)
    
    # Progress callback for UI updates
    progress <- function(message, value = NULL) {
      list(message = message, value = value, timestamp = Sys.time())
    }
    
    # Async DESeq2 analysis implementation
    tryCatch({
      # Step 1: Data preparation
      progress("Preparing count matrix...", 10)
      count_matrix <- as.matrix(expression_data)
      count_matrix <- round(count_matrix)
      count_matrix[count_matrix < 0] <- 0
      
      # Step 2: Sample data preparation
      progress("Preparing sample metadata...", 20)
      sample_data <- annotation_data
      rownames(sample_data) <- sample_data$Sample
      sample_data <- sample_data[colnames(count_matrix), , drop = FALSE]
      
      # Step 3: Create DESeqDataSet
      progress("Creating DESeq2 object...", 30)
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_data,
        design = ~ Condition
      )
      
      # Step 4: Filter low count genes
      progress("Filtering low count genes...", 40)
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]
      
      # Step 5: Run DESeq2 analysis
      progress("Running DESeq2 differential expression...", 60)
      dds <- DESeq(dds)
      
      # Step 6: Extract results
      progress("Extracting results...", 80)
      res <- results(dds, contrast = c("Condition", contrast[1], contrast[2]))
      
      # Step 7: Format results
      progress("Formatting results and gene conversion...", 90)
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
      
      progress("Analysis complete!", 100)
      
      # Return comprehensive results structure
      return(list(
        success = TRUE,
        dds = dds,
        results = res,
        results_df = results_df,
        normalized_counts = counts(dds, normalized = TRUE),
        raw_counts = counts(dds, normalized = FALSE),
        size_factors = sizeFactors(dds),
        metadata = list(
          contrast = contrast,
          padj_cutoff = padj_cutoff,
          fc_cutoff = fc_cutoff,
          species = species,
          n_samples = ncol(dds),
          n_genes_tested = nrow(results_df),
          n_significant = sum(results_df$significant, na.rm = TRUE),
          analysis_date = Sys.Date()
        ),
        progress_log = "DESeq2 analysis completed successfully"
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        error = paste("DESeq2 analysis failed:", e$message),
        progress_log = paste("Error occurred:", e$message)
      ))
    })
  })
}

# Async Gene Conversion Handler
async_gene_conversion <- function(gene_ids, species = "human", id_type = "ensembl") {
  
  future_promise({
    library(biomaRt)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    
    tryCatch({
      # Progress tracking
      progress <- function(msg, val) list(message = msg, value = val, timestamp = Sys.time())
      
      progress("Initializing gene conversion...", 10)
      
      # Use offline annotation databases first for speed
      if (species == "human" && id_type == "ensembl") {
        progress("Using org.Hs.eg.db for human gene conversion...", 30)
        
        # Convert Ensembl to Symbol using offline database
        symbols <- tryCatch({
          mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", 
                keytype = "ENSEMBL", multiVals = "first")
        }, error = function(e) NULL)
        
        if (!is.null(symbols)) {
          progress("Offline conversion successful", 80)
          
          results <- data.frame(
            original_id = gene_ids,
            gene_symbol = symbols,
            conversion_method = "org.Hs.eg.db",
            stringsAsFactors = FALSE
          )
          
          # Calculate success rate
          success_rate <- round(100 * sum(!is.na(results$gene_symbol)) / nrow(results), 1)
          
          progress(paste("Gene conversion complete:", success_rate, "% success"), 100)
          
          return(list(
            success = TRUE,
            results = results,
            success_rate = success_rate,
            method = "offline_annotation_db"
          ))
        }
      }
      
      # Fallback to BioMart if offline method fails
      progress("Falling back to BioMart online conversion...", 40)
      
      # BioMart connection with retry logic
      mart <- NULL
      for (mirror in c("www", "useast", "asia")) {
        tryCatch({
          if (species == "human") {
            mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = mirror)
            symbol_attr <- "hgnc_symbol"
          } else if (species == "mouse") {
            mart <- useEnsembl("ensembl", "mmusculus_gene_ensembl", mirror = mirror)
            symbol_attr <- "mgi_symbol"
          }
          break
        }, error = function(e) NULL)
      }
      
      if (is.null(mart)) {
        return(list(
          success = FALSE,
          error = "Failed to connect to BioMart servers",
          method = "biomart_fallback"
        ))
      }
      
      progress("Connected to BioMart, converting genes...", 60)
      
      # Batch conversion with progress tracking
      batch_size <- 200
      all_results <- data.frame()
      total_batches <- ceiling(length(gene_ids) / batch_size)
      
      for (i in seq(1, length(gene_ids), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(gene_ids))
        batch_ids <- gene_ids[i:end_idx]
        batch_num <- ceiling(i / batch_size)
        
        batch_progress <- 60 + (30 * batch_num / total_batches)
        progress(paste("Converting batch", batch_num, "of", total_batches), batch_progress)
        
        batch_result <- tryCatch({
          getBM(attributes = c("ensembl_gene_id", symbol_attr),
                filters = "ensembl_gene_id",
                values = batch_ids,
                mart = mart)
        }, error = function(e) NULL)
        
        if (!is.null(batch_result) && nrow(batch_result) > 0) {
          colnames(batch_result)[2] <- "gene_symbol"
          all_results <- rbind(all_results, batch_result)
        }
        
        Sys.sleep(0.3)  # Rate limiting
      }
      
      # Create complete results with NAs for missing
      complete_results <- data.frame(
        original_id = gene_ids,
        gene_symbol = all_results$gene_symbol[match(gene_ids, all_results$ensembl_gene_id)],
        conversion_method = "biomart",
        stringsAsFactors = FALSE
      )
      
      success_rate <- round(100 * sum(!is.na(complete_results$gene_symbol)) / nrow(complete_results), 1)
      progress(paste("BioMart conversion complete:", success_rate, "% success"), 100)
      
      return(list(
        success = TRUE,
        results = complete_results,
        success_rate = success_rate,
        method = "biomart_online"
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        error = paste("Gene conversion failed:", e$message),
        method = "failed"
      ))
    })
  })
}

# Async Pathway Analysis Handler
async_pathway_analysis <- function(deseq2_results, analysis_type = "GO", species = "human") {
  
  future_promise({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(DOSE)
    library(enrichplot)
    
    tryCatch({
      progress <- function(msg, val) list(message = msg, value = val, timestamp = Sys.time())
      
      progress("Preparing gene lists for pathway analysis...", 10)
      
      # Prepare significant genes
      sig_genes <- deseq2_results[deseq2_results$significant == TRUE, ]
      
      if (nrow(sig_genes) < 10) {
        return(list(
          success = FALSE,
          error = "Insufficient significant genes for pathway analysis (minimum 10 required)",
          n_genes = nrow(sig_genes)
        ))
      }
      
      progress(paste("Found", nrow(sig_genes), "significant genes"), 20)
      
      # Convert gene IDs if needed
      gene_list <- sig_genes$gene
      
      # Determine organism database
      if (species == "human") {
        orgdb <- org.Hs.eg.db
        organism <- "hsa"
      } else if (species == "mouse") {
        orgdb <- org.Mm.eg.db
        organism <- "mmu"
      } else {
        return(list(success = FALSE, error = "Unsupported species for pathway analysis"))
      }
      
      progress(paste("Running", analysis_type, "pathway analysis..."), 50)
      
      # Run pathway analysis based on type
      if (analysis_type == "GO") {
        enrich_result <- enrichGO(
          gene = gene_list,
          OrgDb = orgdb,
          ont = "BP",  # Biological Process
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
      } else if (analysis_type == "KEGG") {
        enrich_result <- enrichKEGG(
          gene = gene_list,
          organism = organism,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH"
        )
      }
      
      progress("Processing pathway results...", 80)
      
      if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
        return(list(
          success = FALSE,
          error = "No significant pathways found",
          n_input_genes = length(gene_list)
        ))
      }
      
      # Format results
      pathway_results <- as.data.frame(enrich_result@result)
      n_pathways <- nrow(pathway_results)
      
      progress(paste("Found", n_pathways, "enriched pathways"), 100)
      
      return(list(
        success = TRUE,
        results = pathway_results,
        enrichment_object = enrich_result,
        analysis_type = analysis_type,
        species = species,
        n_input_genes = length(gene_list),
        n_pathways = n_pathways,
        analysis_date = Sys.Date()
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        error = paste("Pathway analysis failed:", e$message)
      ))
    })
  })
}

# Async Large Dataset Processing Handler
async_large_dataset_processing <- function(file_path, chunk_size = 5000) {
  
  future_promise({
    library(readr)
    library(dplyr)
    
    tryCatch({
      progress <- function(msg, val) list(message = msg, value = val, timestamp = Sys.time())
      
      progress("Starting large dataset processing...", 5)
      
      # Read file in chunks to handle large datasets
      file_info <- file.info(file_path)
      file_size_mb <- round(file_info$size / 1024^2, 2)
      
      progress(paste("Processing", file_size_mb, "MB file in chunks..."), 10)
      
      # Read first chunk to determine structure
      first_chunk <- read_csv(file_path, n_max = 100, show_col_types = FALSE)
      n_cols <- ncol(first_chunk)
      col_names <- colnames(first_chunk)
      
      progress(paste("Detected", n_cols, "columns"), 20)
      
      # Process file in chunks
      all_data <- NULL
      chunk_num <- 0
      total_rows <- 0
      
      # Estimate total rows for progress tracking
      estimated_rows <- nrow(read_csv(file_path, col_select = 1, show_col_types = FALSE))
      
      # Read and process in chunks
      con <- file(file_path, "r")
      header <- readLines(con, n = 1)
      close(con)
      
      chunk_start <- 1
      while (TRUE) {
        chunk_data <- tryCatch({
          read_csv(file_path, skip = chunk_start, n_max = chunk_size, 
                  col_names = col_names, show_col_types = FALSE)
        }, error = function(e) NULL)
        
        if (is.null(chunk_data) || nrow(chunk_data) == 0) break
        
        chunk_num <- chunk_num + 1
        total_rows <- total_rows + nrow(chunk_data)
        
        # Progress based on estimated rows
        progress_val <- min(95, 20 + (70 * total_rows / estimated_rows))
        progress(paste("Processed chunk", chunk_num, "-", total_rows, "rows"), progress_val)
        
        # Data quality checks and processing
        chunk_data <- chunk_data %>%
          # Remove completely empty rows
          filter(rowSums(is.na(.)) < ncol(.)) %>%
          # Basic data cleaning
          mutate_if(is.character, ~na_if(., "")) %>%
          mutate_if(is.character, ~na_if(., "NA"))
        
        # Combine chunks
        if (is.null(all_data)) {
          all_data <- chunk_data
        } else {
          all_data <- bind_rows(all_data, chunk_data)
        }
        
        chunk_start <- chunk_start + chunk_size + 1
        
        # Memory management
        if (chunk_num %% 10 == 0) {
          gc()  # Garbage collection every 10 chunks
        }
      }
      
      progress("Finalizing processed data...", 100)
      
      # Final data summary
      summary_stats <- list(
        total_rows = nrow(all_data),
        total_cols = ncol(all_data),
        chunks_processed = chunk_num,
        file_size_mb = file_size_mb,
        memory_usage_mb = round(object.size(all_data) / 1024^2, 2),
        processing_time = Sys.time()
      )
      
      return(list(
        success = TRUE,
        data = all_data,
        summary = summary_stats,
        column_names = colnames(all_data)
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        error = paste("Large dataset processing failed:", e$message)
      ))
    })
  })
}

# Async Progress Tracking Utility
create_async_progress_handler <- function(session, output_id) {
  function(promise_obj, success_callback = NULL, error_callback = NULL) {
    
    # Set up progress tracking
    progress_value <- reactiveVal(0)
    progress_message <- reactiveVal("Initializing...")
    
    # Render progress UI
    output[[output_id]] <- renderUI({
      div(
        class = "progress-container",
        style = "margin: 20px 0;",
        
        div(
          class = "progress-bar-container",
          style = "background-color: #f0f0f0; border-radius: 10px; height: 25px; overflow: hidden;",
          
          div(
            class = "progress-bar",
            style = paste0(
              "background: linear-gradient(90deg, #28a745, #20c997); ",
              "height: 100%; transition: width 0.3s ease; color: white; ",
              "display: flex; align-items: center; justify-content: center; ",
              "font-weight: bold; width: ", progress_value(), "%;"
            ),
            paste0(progress_value(), "%")
          )
        ),
        
        div(
          class = "progress-message",
          style = "margin-top: 10px; color: #666; font-size: 14px;",
          progress_message()
        )
      )
    })
    
    # Handle promise resolution
    promise_obj %>%
      then(
        onFulfilled = function(result) {
          progress_value(100)
          progress_message("Complete!")
          
          if (!is.null(success_callback)) {
            success_callback(result)
          }
          
          return(result)
        },
        onRejected = function(error) {
          progress_value(0)
          progress_message(paste("Error:", error$message))
          
          if (!is.null(error_callback)) {
            error_callback(error)
          }
          
          return(error)
        }
      )
  }
}

# Export async handlers
list(
  async_deseq2_analysis = async_deseq2_analysis,
  async_gene_conversion = async_gene_conversion,
  async_pathway_analysis = async_pathway_analysis,
  async_large_dataset_processing = async_large_dataset_processing,
  create_async_progress_handler = create_async_progress_handler
)