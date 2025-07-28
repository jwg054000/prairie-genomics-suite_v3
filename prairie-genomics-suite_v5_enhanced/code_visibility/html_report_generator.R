# HTML Report Generator for Code Visibility System
# Generate comprehensive HTML reports with embedded plots and interactive elements
#
# Author: Prairie Genomics Team
# Date: January 27, 2025
# Purpose: Phase 3 - Enhanced Export Capabilities

# Load required libraries
if (!requireNamespace("htmltools", quietly = TRUE)) {
  cat("⚠️ Installing htmltools package...\n")
  install.packages("htmltools")
}

if (!requireNamespace("plotly", quietly = TRUE)) {
  cat("⚠️ Installing plotly package...\n") 
  install.packages("plotly")
}

library(htmltools)
library(plotly)
library(ggplot2)

# Generate comprehensive HTML report with code, results, and plots
generate_html_report <- function(session_id, output_file = NULL, include_plots = TRUE, 
                                interactive_plots = TRUE, theme = "default") {
  tryCatch({
    cat("📄 Generating comprehensive HTML report...\n")
    
    # Load session data
    session_data <- get_session_data(session_id)
    if (is.null(session_data)) {
      stop("Session data not found for session: ", session_id)
    }
    
    # Create output filename if not provided
    if (is.null(output_file)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      output_file <- paste0("prairie_genomics_report_", session_id, "_", timestamp, ".html")
    }
    
    cat("📊 Building report sections...\n")
    
    # Build HTML report sections
    html_content <- list(
      create_html_header(session_data, theme),
      create_html_summary(session_data),
      create_html_methods(session_data),
      create_html_code_section(session_data),
      create_html_results_section(session_data, include_plots, interactive_plots),
      create_html_appendix(session_data),
      create_html_footer()
    )
    
    # Combine all sections
    full_html <- do.call(tagList, html_content)
    
    # Save HTML file
    cat("💾 Saving HTML report to:", output_file, "\n")
    save_html(full_html, file = output_file, selfcontained = TRUE)
    
    # Generate report metadata
    report_info <- list(
      session_id = session_id,
      output_file = output_file,
      generated_at = Sys.time(),
      file_size = file.size(output_file),
      sections = c("Summary", "Methods", "Code", "Results", "Appendix"),
      includes_plots = include_plots,
      interactive = interactive_plots,
      theme = theme
    )
    
    cat("✅ HTML report generated successfully!\n")
    cat("📋 Report details:\n")
    cat("   - File:", output_file, "\n")
    cat("   - Size:", round(report_info$file_size / 1024, 1), "KB\n")
    cat("   - Sections:", length(report_info$sections), "\n")
    cat("   - Interactive plots:", interactive_plots, "\n")
    
    return(report_info)
    
  }, error = function(e) {
    cat("❌ HTML report generation failed:", e$message, "\n")
    return(NULL)
  })
}

# Create HTML header with metadata and CSS
create_html_header <- function(session_data, theme = "default") {
  # Define CSS styles based on theme
  css_styles <- switch(theme,
    "default" = "
      body { font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; line-height: 1.6; color: #333; }
      .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }
      .section { margin: 30px 0; padding: 20px; border-left: 4px solid #667eea; background: #f8f9fa; border-radius: 5px; }
      .code-block { background: #f4f4f4; border: 1px solid #ddd; border-radius: 5px; padding: 15px; font-family: 'Courier New', monospace; overflow-x: auto; }
      .method-box { background: #e3f2fd; border-left: 4px solid #2196f3; padding: 15px; margin: 10px 0; border-radius: 5px; }
      .result-box { background: #f3e5f5; border-left: 4px solid #9c27b0; padding: 15px; margin: 10px 0; border-radius: 5px; }
      .plot-container { text-align: center; margin: 20px 0; padding: 15px; background: white; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
      h1, h2, h3 { color: #333; }
      .timestamp { color: #666; font-size: 0.9em; }
      .code-title { background: #667eea; color: white; padding: 8px 15px; margin: 0 0 0 0; border-radius: 5px 5px 0 0; font-weight: bold; }
    ",
    "scientific" = "
      body { font-family: 'Times New Roman', serif; margin: 40px; line-height: 1.8; color: #000; background: white; }
      .header { background: #2c3e50; color: white; padding: 30px; margin-bottom: 30px; text-align: center; }
      .section { margin: 30px 0; padding: 20px; border: 1px solid #bdc3c7; background: #fff; }
      .code-block { background: #ecf0f1; border: 1px solid #bdc3c7; padding: 15px; font-family: 'Courier New', monospace; font-size: 0.9em; }
      .method-box { background: #fff; border: 2px solid #3498db; padding: 15px; margin: 10px 0; }
      .result-box { background: #fff; border: 2px solid #e74c3c; padding: 15px; margin: 10px 0; }
      h1, h2, h3 { color: #2c3e50; text-align: center; }
    ",
    "modern" = "
      body { font-family: 'Helvetica Neue', Arial, sans-serif; margin: 20px; line-height: 1.6; color: #2c3e50; background: #f8f9fa; }
      .header { background: linear-gradient(45deg, #ff6b6b, #4ecdc4); color: white; padding: 40px; border-radius: 15px; margin-bottom: 30px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); }
      .section { margin: 30px 0; padding: 25px; background: white; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }
      .code-block { background: #1e1e1e; color: #f8f8f2; border-radius: 8px; padding: 20px; font-family: 'Fira Code', monospace; overflow-x: auto; }
      .method-box { background: linear-gradient(90deg, #667eea, #764ba2); color: white; padding: 20px; margin: 15px 0; border-radius: 10px; }
      .result-box { background: linear-gradient(90deg, #f093fb, #f5576c); color: white; padding: 20px; margin: 15px 0; border-radius: 10px; }
    "
  )
  
  # Create header content
  header_content <- div(
    class = "header",
    h1("Prairie Genomics Suite - Analysis Report"),
    p(paste("Session ID:", session_data$session_id)),
    p(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), class = "timestamp"),
    if (!is.null(session_data$user_info)) {
      p(paste("Analysis by:", session_data$user_info$name %||% "User"))
    }
  )
  
  # Combine head and header
  tagList(
    tags$head(
      tags$title("Prairie Genomics Analysis Report"),
      tags$style(HTML(css_styles)),
      tags$meta(charset = "UTF-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0"),
      tags$meta(name = "generator", content = "Prairie Genomics Suite v5 Enhanced")
    ),
    header_content
  )
}

# Create analysis summary section
create_html_summary <- function(session_data) {
  div(
    class = "section",
    h2("📊 Analysis Summary"),
    
    div(
      class = "method-box",
      h4("📋 Analysis Overview"),
      p("This report contains the complete computational analysis performed using the Prairie Genomics Suite v5 Enhanced platform."),
      tags$ul(
        tags$li(paste("Total analysis steps:", length(session_data$steps))),
        tags$li(paste("Analysis types:", paste(unique(sapply(session_data$steps, function(x) x$analysis_type)), collapse = ", "))),
        tags$li(paste("Species analyzed:", session_data$species %||% "Not specified")),
        tags$li(paste("Session duration:", format_duration(session_data$end_time, session_data$start_time)))
      )
    ),
    
    if (length(session_data$steps) > 0) {
      div(
        h4("🔬 Analysis Workflow"),
        tags$ol(
          lapply(session_data$steps, function(step) {
            tags$li(
              strong(step$step_name),
              " - ",
              step$description %||% "Analysis step completed",
              span(paste(" (", format(step$timestamp, "%H:%M:%S"), ")"), class = "timestamp")
            )
          })
        )
      )
    }
  )
}

# Create methods and documentation section
create_html_methods <- function(session_data) {
  div(
    class = "section", 
    h2("🔬 Methods and Documentation"),
    
    div(
      class = "method-box",
      h4("📚 Statistical Methods"),
      p("This analysis employed established bioinformatics methods with appropriate statistical corrections:"),
      tags$ul(
        if (any(sapply(session_data$steps, function(x) x$analysis_type == "DESeq2"))) {
          tags$li("Differential Expression: DESeq2 with Benjamini-Hochberg FDR correction (Love et al., 2014)")
        },
        if (any(sapply(session_data$steps, function(x) x$analysis_type == "GO"))) {
          tags$li("Gene Ontology Enrichment: Over-representation analysis with hypergeometric test")
        },
        if (any(sapply(session_data$steps, function(x) x$analysis_type == "KEGG"))) {
          tags$li("KEGG Pathway Analysis: Enrichment testing with multiple testing correction")
        },
        if (any(sapply(session_data$steps, function(x) x$analysis_type == "GSEA"))) {
          tags$li("Gene Set Enrichment Analysis: Pre-ranked GSEA with permutation testing")
        }
      )
    ),
    
    div(
      class = "method-box",
      h4("💻 Computational Environment"),
      tags$ul(
        tags$li(paste("R Version:", R.version.string)),
        tags$li(paste("Platform:", R.version$platform)),
        tags$li("Key Packages: DESeq2, clusterProfiler, ggplot2, plotly"),
        tags$li("Analysis Platform: Prairie Genomics Suite v5 Enhanced")
      )
    ),
    
    div(
      class = "method-box",
      h4("📖 References"),
      tags$ul(
        tags$li("Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550."),
        tags$li("Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284-287."),
        tags$li("Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS 102, 15545-15550.")
      )
    )
  )
}

# Create comprehensive code section
create_html_code_section <- function(session_data) {
  div(
    class = "section",
    h2("💻 Complete Analysis Code"),
    
    p("The following R code can be used to reproduce this entire analysis:"),
    
    div(
      class = "code-title",
      "🔧 Complete Reproducible R Script"
    ),
    div(
      class = "code-block",
      pre(
        generate_complete_script(session_data$session_id)
      )
    ),
    
    h3("📋 Step-by-Step Code Breakdown"),
    
    lapply(session_data$steps, function(step) {
      div(
        h4(paste("Step", step$step_number, ":", step$step_name)),
        if (!is.null(step$description)) {
          p(step$description)
        },
        div(
          class = "code-title",
          paste("💻", step$step_name, "Code")
        ),
        div(
          class = "code-block",
          pre(step$r_code %||% "# Code not available for this step")
        ),
        if (!is.null(step$output_summary)) {
          div(
            class = "result-box",
            h5("📊 Results Summary"),
            p(step$output_summary)
          )
        }
      )
    })
  )
}

# Create results section with embedded plots
create_html_results_section <- function(session_data, include_plots = TRUE, interactive_plots = TRUE) {
  div(
    class = "section",
    h2("📈 Analysis Results"),
    
    if (include_plots) {
      div(
        h3("📊 Generated Visualizations"),
        p("The following plots were generated during the analysis:"),
        
        # Embed plots if they exist
        lapply(session_data$steps, function(step) {
          if (!is.null(step$plots) && length(step$plots) > 0) {
            div(
              h4(paste("Plots from", step$step_name)),
              lapply(step$plots, function(plot_info) {
                create_embedded_plot(plot_info, interactive_plots)
              })
            )
          }
        })
      )
    },
    
    # Summary statistics
    div(
      h3("📋 Analysis Statistics"),
      create_results_summary_table(session_data)
    ),
    
    # Export information
    div(
      class = "method-box",
      h4("💾 Data Export"),
      p("All results can be exported in multiple formats:"),
      tags$ul(
        tags$li("Raw results: CSV/TSV format"),
        tags$li("Plots: PNG, PDF, SVG formats"),
        tags$li("Code: R script (.R) and R Markdown (.Rmd)"),
        tags$li("Complete analysis: HTML report (this document)")
      )
    )
  )
}

# Create embedded plot with plotly support
create_embedded_plot <- function(plot_info, interactive = TRUE) {
  div(
    class = "plot-container",
    h5(plot_info$title %||% "Analysis Plot"),
    
    if (interactive && !is.null(plot_info$ggplot_object)) {
      tryCatch({
        # Convert ggplot to plotly for interactivity
        plotly_plot <- ggplotly(plot_info$ggplot_object)
        plotly_plot
      }, error = function(e) {
        # Fallback to static image
        if (!is.null(plot_info$file_path) && file.exists(plot_info$file_path)) {
          tags$img(src = plot_info$file_path, style = "max-width: 100%; height: auto;")
        } else {
          div(
            class = "result-box",
            p("⚠️ Plot not available - ", plot_info$title %||% "Unknown plot")
          )
        }
      })
    } else if (!is.null(plot_info$file_path) && file.exists(plot_info$file_path)) {
      tags$img(src = plot_info$file_path, style = "max-width: 100%; height: auto;")
    } else {
      div(
        class = "result-box",
        p("📊 Plot generated: ", plot_info$title %||% "Analysis visualization")
      )
    }
  )
}

# Create results summary table
create_results_summary_table <- function(session_data) {
  # Extract key results from analysis steps
  results_summary <- data.frame(
    Step = sapply(session_data$steps, function(x) x$step_name),
    Type = sapply(session_data$steps, function(x) x$analysis_type %||% "Unknown"),
    Status = sapply(session_data$steps, function(x) if(x$success) "✅ Success" else "❌ Failed"),
    Duration = sapply(session_data$steps, function(x) {
      if (!is.null(x$execution_time)) {
        paste(round(x$execution_time, 2), "seconds")
      } else {
        "Unknown"
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Convert to HTML table
  tags$table(
    style = "width: 100%; border-collapse: collapse; margin: 20px 0;",
    tags$thead(
      tags$tr(
        lapply(names(results_summary), function(col) {
          tags$th(col, style = "border: 1px solid #ddd; padding: 12px; background: #f2f2f2; text-align: left;")
        })
      )
    ),
    tags$tbody(
      lapply(1:nrow(results_summary), function(i) {
        tags$tr(
          lapply(results_summary[i, ], function(cell) {
            tags$td(cell, style = "border: 1px solid #ddd; padding: 12px;")
          })
        )
      })
    )
  )
}

# Create appendix section
create_html_appendix <- function(session_data) {
  div(
    class = "section",
    h2("📎 Appendix"),
    
    div(
      class = "method-box",
      h4("🔧 Session Information"),
      div(
        class = "code-block",
        pre(
          paste(
            "Session ID:", session_data$session_id,
            "\nStart Time:", format(session_data$start_time),
            "\nEnd Time:", format(session_data$end_time %||% Sys.time()),
            "\nTotal Steps:", length(session_data$steps),
            "\nPlatform:", R.version$platform,
            "\nR Version:", R.version.string
          )
        )
      )
    ),
    
    div(
      class = "method-box", 
      h4("📦 Package Versions"),
      div(
        class = "code-block",
        pre(
          capture.output(sessionInfo())
        )
      )
    ),
    
    div(
      class = "method-box",
      h4("ℹ️ About Prairie Genomics Suite"),
      p("This analysis was performed using the Prairie Genomics Suite v5 Enhanced, a comprehensive platform for genomics data analysis featuring:"),
      tags$ul(
        tags$li("Complete code visibility and reproducibility"),
        tags$li("Integrated differential expression analysis"),
        tags$li("Comprehensive pathway analysis"),
        tags$li("Interactive visualizations"),
        tags$li("Multiple export formats"),
        tags$li("Educational documentation")
      ),
      p("For more information, visit: https://github.com/prairie-genomics-suite")
    )
  )
}

# Create footer
create_html_footer <- function() {
  div(
    style = "text-align: center; margin-top: 50px; padding: 20px; border-top: 1px solid #ddd; color: #666;",
    p("Generated by Prairie Genomics Suite v5 Enhanced"),
    p(paste("Report created on:", format(Sys.time(), "%Y-%m-%d at %H:%M:%S"))),
    p("© 2025 Prairie Genomics Team - Open Source Bioinformatics Platform")
  )
}

# Helper function to format duration
format_duration <- function(end_time, start_time) {
  if (is.null(end_time) || is.null(start_time)) return("Unknown")
  
  duration <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  if (duration < 1) {
    paste(round(duration * 60, 1), "seconds")
  } else if (duration < 60) {
    paste(round(duration, 1), "minutes") 
  } else {
    hours <- floor(duration / 60)
    mins <- round(duration %% 60)
    paste(hours, "hours", mins, "minutes")
  }
}

# Enhanced function to generate reports with different templates
generate_themed_html_report <- function(session_id, template = "default", 
                                       output_dir = "reports", custom_css = NULL) {
  tryCatch({
    cat("🎨 Generating themed HTML report...\n")
    cat("📋 Template:", template, "\n")
    
    # Create output directory
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Generate filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(output_dir, paste0("prairie_report_", template, "_", timestamp, ".html"))
    
    # Generate report with theme
    report_info <- generate_html_report(
      session_id = session_id,
      output_file = output_file,
      include_plots = TRUE,
      interactive_plots = (template %in% c("default", "modern")),
      theme = template
    )
    
    return(report_info)
    
  }, error = function(e) {
    cat("❌ Themed report generation failed:", e$message, "\n")
    return(NULL)
  })
}

cat("✅ HTML Report Generator loaded successfully!\n")
cat("📋 Available functions:\n")
cat("   - generate_html_report(): Create comprehensive HTML reports\n")
cat("   - generate_themed_html_report(): Create themed reports\n") 
cat("   - Templates: 'default', 'scientific', 'modern'\n")
cat("🎯 Ready for Phase 3 HTML report generation!\n")