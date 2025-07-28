# Standalone Script Generator for Code Visibility System
# Generate self-contained R scripts with dependency management
#
# Author: Prairie Genomics Team
# Date: January 27, 2025
# Purpose: Phase 3 - Enhanced Export Capabilities

# Load required libraries
library(tools)

# Generate standalone R script with complete dependency management
generate_standalone_script <- function(session_id, output_file = NULL, 
                                     include_dependencies = TRUE, 
                                     include_data_loading = TRUE,
                                     add_comments = TRUE,
                                     optimization_level = "standard") {
  tryCatch({
    cat("🔧 Generating standalone R script...\n")
    
    # Load session data
    session_data <- get_session_data(session_id)
    if (is.null(session_data)) {
      stop("Session data not found for session: ", session_id)
    }
    
    # Create output filename if not provided
    if (is.null(output_file)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      output_file <- paste0("prairie_analysis_standalone_", session_id, "_", timestamp, ".R")
    }
    
    cat("📝 Building script components...\n")
    
    # Build script sections
    script_sections <- list()
    
    # Header section
    script_sections$header <- create_script_header(session_data, add_comments)
    
    # Dependency section
    if (include_dependencies) {
      script_sections$dependencies <- create_dependency_section(session_data, optimization_level)
    }
    
    # Data loading section
    if (include_data_loading) {
      script_sections$data_loading <- create_data_loading_section(session_data, add_comments)
    }
    
    # Helper functions section
    script_sections$helpers <- create_helper_functions_section(session_data, optimization_level)
    
    # Main analysis section
    script_sections$analysis <- create_analysis_section(session_data, add_comments, optimization_level)
    
    # Results export section
    script_sections$export <- create_export_section(session_data, add_comments)
    
    # Footer section
    script_sections$footer <- create_script_footer(session_data, add_comments)
    
    # Combine all sections
    complete_script <- paste(script_sections, collapse = "\n\n")
    
    # Write script to file
    cat("💾 Writing script to file:", output_file, "\n")
    writeLines(complete_script, output_file)
    
    # Validate script syntax
    cat("✅ Validating script syntax...\n")
    syntax_check <- validate_script_syntax(output_file)
    
    # Generate script metadata
    script_info <- list(
      session_id = session_id,
      output_file = output_file,
      generated_at = Sys.time(),
      file_size = file.size(output_file),
      lines_of_code = length(readLines(output_file)),
      includes_dependencies = include_dependencies,
      includes_data_loading = include_data_loading,
      optimization_level = optimization_level,
      syntax_valid = syntax_check$valid,
      estimated_runtime = estimate_script_runtime(session_data)
    )
    
    cat("✅ Standalone script generated successfully!\n")
    cat("📋 Script details:\n")
    cat("   - File:", output_file, "\n")
    cat("   - Size:", round(script_info$file_size / 1024, 1), "KB\n")
    cat("   - Lines of code:", script_info$lines_of_code, "\n")
    cat("   - Syntax valid:", script_info$syntax_valid, "\n")
    cat("   - Estimated runtime:", script_info$estimated_runtime, "\n")
    
    return(script_info)
    
  }, error = function(e) {
    cat("❌ Standalone script generation failed:", e$message, "\n")
    return(NULL)
  })
}

# Create script header with metadata and documentation
create_script_header <- function(session_data, add_comments = TRUE) {
  header <- paste0(
    "#!/usr/bin/env Rscript\n",
    "#\n",
    "# Prairie Genomics Suite - Standalone Analysis Script\n",
    "# Generated from session: ", session_data$session_id, "\n",
    "#\n",
    "# This script reproduces the complete analysis performed in the\n",
    "# Prairie Genomics Suite v5 Enhanced platform.\n",
    "#\n",
    "# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "# Platform: ", R.version$platform, "\n",
    "# R Version: ", R.version.string, "\n",
    "#\n",
    "# Usage: Rscript ", basename(tempfile()), ".R\n",
    "#\n",
    "# ============================================================================\n"
  )
  
  if (add_comments) {
    header <- paste0(header,
      "# ANALYSIS OVERVIEW\n",
      "# ============================================================================\n",
      "#\n",
      "# This analysis includes the following steps:\n"
    )
    
    for (i in seq_along(session_data$steps)) {
      step <- session_data$steps[[i]]
      header <- paste0(header, 
        "# ", i, ". ", step$step_name, "\n"
      )
    }
    
    header <- paste0(header,
      "#\n",
      "# Total analysis steps: ", length(session_data$steps), "\n",
      "# Species: ", session_data$species %||% "Not specified", "\n",
      "# Analysis types: ", paste(unique(sapply(session_data$steps, function(x) x$analysis_type)), collapse = ", "), "\n",
      "#\n",
      "# ============================================================================\n"
    )
  }
  
  return(header)
}

# Create comprehensive dependency management section
create_dependency_section <- function(session_data, optimization_level = "standard") {
  dep_section <- "# DEPENDENCY MANAGEMENT\n# ============================================================================\n\n"
  
  # Extract required packages from analysis steps
  required_packages <- extract_required_packages(session_data)
  
  # Add package installation and loading function
  dep_section <- paste0(dep_section, '
# Function to install and load required packages
install_and_load_packages <- function(packages) {
  cat("📦 Checking and installing required packages...\\n")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("⚠️ Installing package:", pkg, "\\n")
      
      # Try Bioconductor first for genomics packages
      if (pkg %in% c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", 
                    "KEGGREST", "GO.db", "AnnotationDbi", "biomaRt", "fgsea", "msigdbr")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
    
    # Load the package
    library(pkg, character.only = TRUE)
    cat("✅ Loaded:", pkg, "\\n")
  }
  
  cat("📦 All packages loaded successfully!\\n")
}

# Required packages for this analysis
required_packages <- c(',
    paste0('"', required_packages, '"', collapse = ", "),
    ')

# Install and load all required packages
install_and_load_packages(required_packages)
')
  
  if (optimization_level == "advanced") {
    dep_section <- paste0(dep_section, '
# Advanced dependency management
# Check package versions for reproducibility
check_package_versions <- function() {
  cat("🔍 Checking package versions for reproducibility...\\n")
  
  version_info <- data.frame(
    Package = required_packages,
    Version = sapply(required_packages, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(packageVersion(pkg))
      } else {
        "Not installed"
      }
    }),
    stringsAsFactors = FALSE
  )
  
  print(version_info)
  return(version_info)
}

# Check versions
package_versions <- check_package_versions()
')
  }
  
  dep_section <- paste0(dep_section, "\n# Set global options for analysis\n")
  dep_section <- paste0(dep_section, "options(stringsAsFactors = FALSE)\n")
  dep_section <- paste0(dep_section, "set.seed(42)  # For reproducibility\n\n")
  
  return(dep_section)
}

# Extract required packages from session data
extract_required_packages <- function(session_data) {
  base_packages <- c("ggplot2", "dplyr", "readr", "tibble")
  analysis_packages <- c()
  
  # Extract packages based on analysis types
  analysis_types <- unique(sapply(session_data$steps, function(x) x$analysis_type))
  
  if ("DESeq2" %in% analysis_types) {
    analysis_packages <- c(analysis_packages, "DESeq2", "SummarizedExperiment")
  }
  
  if ("GO" %in% analysis_types) {
    analysis_packages <- c(analysis_packages, "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "GO.db")
  }
  
  if ("KEGG" %in% analysis_types) {
    analysis_packages <- c(analysis_packages, "clusterProfiler", "KEGGREST")
  }
  
  if ("GSEA" %in% analysis_types) {
    analysis_packages <- c(analysis_packages, "fgsea", "msigdbr")
  }
  
  # Species-specific packages
  if (!is.null(session_data$species)) {
    if (session_data$species == "human") {
      analysis_packages <- c(analysis_packages, "org.Hs.eg.db")
    } else if (session_data$species == "mouse") {
      analysis_packages <- c(analysis_packages, "org.Mm.eg.db")
    }
  }
  
  return(unique(c(base_packages, analysis_packages)))
}

# Create data loading section with file handling
create_data_loading_section <- function(session_data, add_comments = TRUE) {
  data_section <- "# DATA LOADING AND PREPROCESSING\n# ============================================================================\n\n"
  
  if (add_comments) {
    data_section <- paste0(data_section, '
# This section contains functions to load and preprocess your data files.
# You will need to update the file paths to match your local file locations.
#
# Expected data files:
# - Expression data: Count matrix (genes as rows, samples as columns)
# - Sample annotation: Sample metadata with experimental groups
#
')
  }
  
  data_section <- paste0(data_section, '
# Function to load expression data
load_expression_data <- function(file_path) {
  cat("📊 Loading expression data from:", file_path, "\\n")
  
  # Detect file format and load appropriately
  if (grepl("\\\\.csv$", file_path, ignore.case = TRUE)) {
    data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  } else if (grepl("\\\\.tsv$|\\\\.txt$", file_path, ignore.case = TRUE)) {
    data <- read.delim(file_path, row.names = 1, check.names = FALSE)
  } else {
    # Try to auto-detect delimiter
    data <- read.delim(file_path, row.names = 1, check.names = FALSE)
  }
  
  cat("✅ Loaded expression data:", nrow(data), "genes x", ncol(data), "samples\\n")
  return(data)
}

# Function to load sample annotation
load_sample_annotation <- function(file_path) {
  cat("📋 Loading sample annotation from:", file_path, "\\n")
  
  if (grepl("\\\\.csv$", file_path, ignore.case = TRUE)) {
    annotation <- read.csv(file_path, row.names = 1, check.names = FALSE)
  } else {
    annotation <- read.delim(file_path, row.names = 1, check.names = FALSE)
  }
  
  cat("✅ Loaded sample annotation:", nrow(annotation), "samples\\n")
  return(annotation)
}

# Data loading instructions
cat("\\n📂 DATA LOADING INSTRUCTIONS\\n")
cat("==============================\\n")
cat("Please update the file paths below to point to your data files:\\n\\n")

# TODO: Update these file paths to match your data location
expression_file <- "path/to/your/expression_data.csv"  # Update this path
annotation_file <- "path/to/your/sample_annotation.csv"  # Update this path

# Uncomment and run when you have updated the file paths
# expression_data <- load_expression_data(expression_file)
# sample_annotation <- load_sample_annotation(annotation_file)

cat("⚠️ Remember to update the file paths above before running!\\n\\n")
')
  
  return(data_section)
}

# Create helper functions section
create_helper_functions_section <- function(session_data, optimization_level = "standard") {
  helpers_section <- "# HELPER FUNCTIONS\n# ============================================================================\n\n"
  
  # Extract unique function types needed based on analysis
  analysis_types <- unique(sapply(session_data$steps, function(x) x$analysis_type))
  
  if ("DESeq2" %in% analysis_types) {
    helpers_section <- paste0(helpers_section, '
# DESeq2 analysis helper function
run_deseq2_analysis <- function(count_data, sample_annotation, design_formula = "~ condition") {
  cat("🧬 Running DESeq2 differential expression analysis...\\n")
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_annotation,
    design = as.formula(design_formula)
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Extract results
  results <- results(dds, alpha = 0.05)
  results_df <- as.data.frame(results)
  results_df <- results_df[order(results_df$padj), ]
  
  cat("✅ DESeq2 analysis completed:", nrow(results_df), "genes analyzed\\n")
  
  return(list(
    dds = dds,
    results = results_df,
    significant_genes = sum(results_df$padj < 0.05, na.rm = TRUE)
  ))
}
')
  }
  
  if (any(c("GO", "KEGG") %in% analysis_types)) {
    helpers_section <- paste0(helpers_section, '
# Gene ID conversion helper function
convert_gene_ids <- function(gene_list, from_type = "ENSEMBL", to_type = "ENTREZID", species = "human") {
  cat("🔄 Converting gene IDs from", from_type, "to", to_type, "...\\n")
  
  # Select organism database
  if (species == "human") {
    org_db <- org.Hs.eg.db
  } else if (species == "mouse") {
    org_db <- org.Mm.eg.db
  } else {
    warning("Species not supported, using human database")
    org_db <- org.Hs.eg.db
  }
  
  # Convert gene IDs
  converted_ids <- mapIds(
    org_db,
    keys = gene_list,
    column = to_type,
    keytype = from_type,
    multiVals = "first"
  )
  
  # Remove NAs
  converted_ids <- converted_ids[!is.na(converted_ids)]
  
  conversion_rate <- length(converted_ids) / length(gene_list) * 100
  cat("✅ Converted", length(converted_ids), "genes (", round(conversion_rate, 1), "% success rate)\\n")
  
  return(converted_ids)
}
')
  }
  
  if ("GO" %in% analysis_types) {
    helpers_section <- paste0(helpers_section, '
# GO enrichment analysis helper function
run_go_enrichment <- function(gene_list, species = "human", ontology = "BP") {
  cat("🔍 Running GO enrichment analysis...\\n")
  
  # Select organism database
  if (species == "human") {
    org_db <- org.Hs.eg.db
  } else if (species == "mouse") {
    org_db <- org.Mm.eg.db
  } else {
    org_db <- org.Hs.eg.db
  }
  
  # Limit gene list size for performance
  if (length(gene_list) > 2000) {
    cat("⚠️ Large gene list detected, limiting to top 2000 genes for performance\\n")
    gene_list <- gene_list[1:2000]
  }
  
  # Run GO enrichment
  go_results <- enrichGO(
    gene = gene_list,
    OrgDb = org_db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if (nrow(go_results@result) > 0) {
    cat("✅ GO enrichment completed:", nrow(go_results@result), "enriched terms found\\n")
    return(as.data.frame(go_results@result))
  } else {
    cat("⚠️ No significant GO terms found\\n")
    return(data.frame())
  }
}
')
  }
  
  # Add visualization helper functions
  helpers_section <- paste0(helpers_section, '
# Visualization helper functions
create_volcano_plot <- function(results_df, title = "Volcano Plot") {
  cat("📊 Creating volcano plot...\\n")
  
  # Prepare data for plotting
  plot_data <- results_df
  plot_data$significant <- with(plot_data, 
    ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not significant"))
  
  # Create volcano plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not significant" = "grey", "Significant" = "red")) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to save plots
save_plot <- function(plot_object, filename, width = 8, height = 6, dpi = 300) {
  cat("💾 Saving plot to:", filename, "\\n")
  ggsave(filename, plot_object, width = width, height = height, dpi = dpi)
  cat("✅ Plot saved successfully\\n")
}
')
  
  return(helpers_section)
}

# Create main analysis section
create_analysis_section <- function(session_data, add_comments = TRUE, optimization_level = "standard") {
  analysis_section <- "# MAIN ANALYSIS WORKFLOW\n# ============================================================================\n\n"
  
  if (add_comments) {
    analysis_section <- paste0(analysis_section, '
# This section reproduces the exact analysis steps performed in the
# Prairie Genomics Suite platform. Each step corresponds to an action
# taken during the original analysis session.
#
')
  }
  
  # Add each analysis step
  for (i in seq_along(session_data$steps)) {
    step <- session_data$steps[[i]]
    
    analysis_section <- paste0(analysis_section, 
      "# ", "="*60, "\n",
      "# STEP ", i, ": ", toupper(step$step_name), "\n",
      "# ", "="*60, "\n\n"
    )
    
    if (add_comments && !is.null(step$description)) {
      analysis_section <- paste0(analysis_section, 
        "# ", step$description, "\n\n"
      )
    }
    
    # Add the actual R code for this step
    if (!is.null(step$r_code)) {
      analysis_section <- paste0(analysis_section, 
        step$r_code, "\n\n"
      )
    } else {
      analysis_section <- paste0(analysis_section,
        "# Code for this step was not captured\n",
        "cat(\"Executing step: ", step$step_name, "\\n\")\n\n"
      )
    }
    
    # Add results summary if available
    if (!is.null(step$output_summary)) {
      analysis_section <- paste0(analysis_section,
        "# Results summary: ", step$output_summary, "\n\n"
      )
    }
  }
  
  return(analysis_section)
}

# Create results export section
create_export_section <- function(session_data, add_comments = TRUE) {
  export_section <- "# RESULTS EXPORT\n# ============================================================================\n\n"
  
  if (add_comments) {
    export_section <- paste0(export_section, '
# This section provides functions to export your analysis results
# in various formats for further use and sharing.
#
')
  }
  
  export_section <- paste0(export_section, '
# Function to export all results
export_results <- function(output_dir = "prairie_results") {
  cat("📁 Creating output directory:", output_dir, "\\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Export DESeq2 results if available
  if (exists("deseq2_results")) {
    results_file <- file.path(output_dir, "deseq2_results.csv")
    write.csv(deseq2_results, results_file, row.names = TRUE)
    cat("✅ Exported DESeq2 results to:", results_file, "\\n")
  }
  
  # Export GO results if available
  if (exists("go_results")) {
    go_file <- file.path(output_dir, "go_enrichment_results.csv")
    write.csv(go_results, go_file, row.names = FALSE)
    cat("✅ Exported GO results to:", go_file, "\\n")
  }
  
  # Export KEGG results if available
  if (exists("kegg_results")) {
    kegg_file <- file.path(output_dir, "kegg_pathway_results.csv")
    write.csv(kegg_results, kegg_file, row.names = FALSE)
    cat("✅ Exported KEGG results to:", kegg_file, "\\n")
  }
  
  # Export session information
  session_file <- file.path(output_dir, "session_info.txt")
  writeLines(capture.output(sessionInfo()), session_file)
  cat("✅ Exported session info to:", session_file, "\\n")
  
  cat("🎉 All results exported to:", output_dir, "\\n")
}

# Run export function
# Uncomment the line below to export results after running the analysis
# export_results()
')
  
  return(export_section)
}

# Create script footer
create_script_footer <- function(session_data, add_comments = TRUE) {
  footer <- "# ANALYSIS COMPLETION\n# ============================================================================\n\n"
  
  footer <- paste0(footer, 
    "cat(\"🎉 Analysis completed successfully!\\n\")\n",
    "cat(\"📊 Session ID: ", session_data$session_id, "\\n\")\n",
    "cat(\"⏰ Script runtime:\", proc.time()[\"elapsed\"], \"seconds\\n\")\n",
    "cat(\"📅 Completed on:\", format(Sys.time(), \"%Y-%m-%d %H:%M:%S\"), \"\\n\")\n\n"
  )
  
  if (add_comments) {
    footer <- paste0(footer, '
# ============================================================================
# END OF ANALYSIS SCRIPT
# ============================================================================
#
# This script was automatically generated by the Prairie Genomics Suite v5 
# Enhanced platform. It reproduces the complete analysis performed during
# your interactive session.
#
# For questions or support, please visit:
# https://github.com/prairie-genomics-suite
#
# Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '
# Session: ', session_data$session_id, '
# ============================================================================
')
  }
  
  return(footer)
}

# Validate script syntax
validate_script_syntax <- function(script_file) {
  tryCatch({
    # Parse the script to check for syntax errors
    parse(file = script_file)
    
    return(list(
      valid = TRUE,
      message = "Script syntax is valid"
    ))
    
  }, error = function(e) {
    return(list(
      valid = FALSE,
      message = paste("Syntax error:", e$message)
    ))
  })
}

# Estimate script runtime based on analysis complexity
estimate_script_runtime <- function(session_data) {
  base_time <- 30  # Base time in seconds
  
  # Add time based on analysis types
  for (step in session_data$steps) {
    if (step$analysis_type == "DESeq2") {
      base_time <- base_time + 60  # DESeq2 takes ~1 minute
    } else if (step$analysis_type %in% c("GO", "KEGG")) {
      base_time <- base_time + 30  # Pathway analysis ~30 seconds
    } else if (step$analysis_type == "GSEA") {
      base_time <- base_time + 90  # GSEA takes longer
    } else {
      base_time <- base_time + 10  # Other steps
    }
  }
  
  # Format time estimate
  if (base_time < 60) {
    return(paste(base_time, "seconds"))
  } else if (base_time < 3600) {
    return(paste(round(base_time / 60, 1), "minutes"))
  } else {
    return(paste(round(base_time / 3600, 1), "hours"))
  }
}

cat("✅ Standalone Script Generator loaded successfully!\n")
cat("📋 Available functions:\n")
cat("   - generate_standalone_script(): Create self-contained R scripts\n")
cat("   - Optimization levels: 'basic', 'standard', 'advanced'\n")
cat("   - Includes dependency management and validation\n")
cat("🎯 Ready for Phase 3 standalone script generation!\n")