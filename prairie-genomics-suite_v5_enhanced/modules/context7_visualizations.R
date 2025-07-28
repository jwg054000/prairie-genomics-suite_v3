# Context7-Inspired Advanced Visualizations for Prairie Genomics Suite v5
# R/Plotly implementation of Context7 best practices for genomics visualization
# Features accessibility-first design, publication quality, and advanced interactivity

# Load required libraries
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  plotly, ggplot2, dplyr, RColorBrewer, 
  pheatmap, car, rgl, ggrepel, DT
)

# UI Function
context7VisualizationsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        3,
        h4("🎨 Context7 Visualizations"),
        
        wellPanel(
          h5("📊 Plot Selection"),
          
          selectInput(
            ns("plot_type"),
            "Choose visualization:",
            choices = list(
              "🌋 Enhanced Volcano Plot" = "volcano",
              "🔥 Interactive Heatmap" = "heatmap", 
              "📈 Advanced PCA Plot" = "pca",
              "📊 MA Plot" = "ma_plot",
              "📉 Expression Profile" = "expression_profile"
            ),
            selected = "volcano"
          ),
          
          # Data source selection
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] != 'pca'"),
            
            uiOutput(ns("comparison_selection"))
          )
        ),
        
        # Plot-specific controls
        wellPanel(
          h5("⚙️ Plot Controls"),
          
          # Volcano plot controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
            
            h6("🎯 Significance Thresholds"),
            
            sliderInput(
              ns("volcano_fc_cutoff"),
              "Log₂ Fold Change Cutoff:",
              min = 0, max = 3, value = 1.0, step = 0.1
            ),
            
            sliderInput(
              ns("volcano_p_cutoff"),
              "Adj. P-value Cutoff:",
              min = 0.001, max = 0.1, value = 0.05, step = 0.001
            ),
            
            h6("🎨 Appearance"),
            
            selectInput(
              ns("volcano_color_scheme"),
              "Color scheme:",
              choices = list(
                "Accessibility Friendly" = "accessible",
                "Publication Quality" = "publication",
                "High Contrast" = "high_contrast",
                "Emory Style" = "emory"
              ),
              selected = "accessible"
            ),
            
            checkboxInput(
              ns("volcano_show_labels"),
              "Show gene labels",
              value = TRUE
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("volcano_show_labels"), "']"),
              numericInput(
                ns("volcano_n_labels"),
                "Number of labels:",
                value = 10,
                min = 0,
                max = 50
              )
            ),
            
            h6("🔍 Gene Highlighting"),
            
            textInput(
              ns("volcano_highlight_genes"),
              "Highlight genes (comma-separated):",
              placeholder = "e.g., TP53, BRCA1, MYC"
            ),
            
            checkboxInput(
              ns("volcano_enable_selection"),
              "Enable interactive selection",
              value = TRUE
            )
          ),
          
          # PCA plot controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'pca'"),
            
            h6("📐 Dimensions"),
            
            radioButtons(
              ns("pca_dimensions"),
              "Plot dimensions:",
              choices = list(
                "2D Plot" = "2d",
                "3D Plot" = "3d"
              ),
              selected = "2d"
            ),
            
            h6("🎨 Styling"),
            
            checkboxInput(
              ns("pca_show_ellipses"),
              "Show confidence ellipses",
              value = TRUE
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("pca_show_ellipses"), "']"),
              numericInput(
                ns("pca_confidence_level"),
                "Confidence level:",
                value = 0.95,
                min = 0.5,
                max = 0.99,
                step = 0.05
              )
            ),
            
            checkboxInput(
              ns("pca_show_labels"),
              "Show sample labels",
              value = TRUE
            ),
            
            selectInput(
              ns("pca_color_by"),
              "Color samples by:",
              choices = NULL  # Will be populated dynamically
            )
          ),
          
          # Heatmap controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
            
            h6("🧬 Gene Selection"),
            
            numericInput(
              ns("heatmap_n_genes"),
              "Number of top genes:",
              value = 50,
              min = 10,
              max = 200
            ),
            
            selectInput(
              ns("heatmap_gene_selection"),
              "Gene selection method:",
              choices = list(
                "Most significant" = "significant",
                "Highest variance" = "variance",
                "Custom list" = "custom"
              ),
              selected = "significant"
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("heatmap_gene_selection"), "'] == 'custom'"),
              textAreaInput(
                ns("heatmap_custom_genes"),
                "Gene list (one per line):",
                placeholder = "TP53\nBRCA1\nMYC\n...",
                height = "100px"
              )
            ),
            
            h6("🎨 Appearance"),
            
            selectInput(
              ns("heatmap_color_scheme"),
              "Color palette:",
              choices = list(
                "Red-Blue" = "RdBu",
                "Red-White-Blue" = "RdBu_r", 
                "Viridis" = "viridis",
                "Plasma" = "plasma"
              ),
              selected = "RdBu_r"
            ),
            
            checkboxInput(
              ns("heatmap_scale_rows"),
              "Scale by rows (Z-score)",
              value = TRUE
            ),
            
            checkboxInput(
              ns("heatmap_cluster_rows"),
              "Cluster genes",
              value = TRUE
            ),
            
            checkboxInput(
              ns("heatmap_cluster_cols"),
              "Cluster samples",
              value = TRUE
            )
          )
        ),
        
        # Export controls
        wellPanel(
          h5("📥 Export Options"),
          
          selectInput(
            ns("export_format"),
            "Export format:",
            choices = list(
              "PNG (High Resolution)" = "png",
              "PDF (Vector)" = "pdf",
              "SVG (Vector)" = "svg",
              "HTML (Interactive)" = "html"
            ),
            selected = "png"
          ),
          
          fluidRow(
            column(
              6,
              numericInput(
                ns("export_width"),
                "Width (px):",
                value = 800,
                min = 400,
                max = 3000
              )
            ),
            column(
              6,
              numericInput(
                ns("export_height"),
                "Height (px):",
                value = 600,
                min = 300,
                max = 2000
              )
            )
          ),
          
          downloadButton(
            ns("download_plot"),
            "📥 Download Plot",
            class = "btn-success",
            style = "width: 100%; margin-top: 10px;"
          )
        )
      ),
      
      column(
        9,
        # Main plot area
        wellPanel(
          style = "min-height: 600px;",
          
          # Plot title and controls
          fluidRow(
            column(
              8,
              h4(textOutput(ns("plot_title")))
            ),
            column(
              4,
              div(
                style = "text-align: right; padding-top: 10px;",
                actionButton(
                  ns("refresh_plot"),
                  "🔄 Refresh Plot",
                  class = "btn-info btn-sm"
                ),
                actionButton(
                  ns("fullscreen_plot"),
                  "🔍 Fullscreen",
                  class = "btn-secondary btn-sm"
                )
              )
            )
          ),
          
          hr(),
          
          # Plot output
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
            plotlyOutput(ns("volcano_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'pca'"),
            plotlyOutput(ns("pca_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
            plotOutput(ns("heatmap_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'ma_plot'"),
            plotlyOutput(ns("ma_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'expression_profile'"),
            plotlyOutput(ns("expression_profile_plot"), height = "600px")
          )
        ),
        
        # Plot information and statistics
        wellPanel(
          h5("📊 Plot Information"),
          
          fluidRow(
            column(
              6,
              h6("📈 Statistics"),
              uiOutput(ns("plot_statistics"))
            ),
            column(
              6,
              h6("🔍 Selected Points"),
              uiOutput(ns("selected_points_info"))
            )
          )
        )
      )
    )
  )
}

# Server Function
context7Visualizations <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values for this module
  local_values <- reactiveValues(
    current_plot = NULL,
    selected_points = NULL,
    plot_data = NULL,
    available_comparisons = NULL
  )
  
  # Color schemes (Context7-inspired, accessibility-first)
  color_schemes <- list(
    accessible = list(
      significant_up = "#E69F00",      # Orange - colorblind safe
      significant_down = "#56B4E9",    # Sky blue - colorblind safe  
      non_significant = "#999999",     # Neutral gray
      highlight = "#CC79A7",           # Purple-pink
      background = "#FFFFFF",
      grid = "#F0F0F0"
    ),
    publication = list(
      significant_up = "#D62728",      # Nature red
      significant_down = "#1F77B4",    # Nature blue
      non_significant = "#7F7F7F",     # Gray
      highlight = "#FF7F0E",           # Orange
      background = "#FFFFFF",
      grid = "#E5E5E5"
    ),
    high_contrast = list(
      significant_up = "#FF0000",      # Pure red
      significant_down = "#0000FF",    # Pure blue
      non_significant = "#808080",     # Medium gray
      highlight = "#FFFF00",           # Yellow
      background = "#FFFFFF", 
      grid = "#CCCCCC"
    ),
    emory = list(
      significant_up = "#fc8d62",      # Emory orange
      significant_down = "#66c2a5",    # Emory teal
      non_significant = "#b3b3b3",     # Light gray
      highlight = "#8da0cb",           # Emory purple
      background = "#FFFFFF",
      grid = "#f0f0f0"
    )
  )
  
  # Update available comparisons when results change
  observe({
    if (!is.null(values$deseq2_results)) {
      # Single comparison case
      local_values$available_comparisons <- list("Current Analysis" = values$deseq2_results)
    } else if (!is.null(values$comparison_results)) {
      # Multiple comparisons case
      local_values$available_comparisons <- values$comparison_results
    }
  })
  
  # Update PCA color options
  observe({
    req(values$annotation_data)
    
    color_options <- colnames(values$annotation_data)
    # Remove non-informative columns
    color_options <- color_options[!color_options %in% c("Sample", "Group_Size", "Detection_Method")]
    
    updateSelectInput(session, "pca_color_by",
                     choices = color_options,
                     selected = "Condition")
  })
  
  # Update comparison selection
  output$comparison_selection <- renderUI({
    req(local_values$available_comparisons)
    
    comparison_names <- names(local_values$available_comparisons)
    
    selectInput(
      ns("selected_comparison"),
      "Select comparison:",
      choices = comparison_names,
      selected = comparison_names[1]
    )
  })
  
  # Enhanced volcano plot (Context7-inspired)
  output$volcano_plot <- renderPlotly({
    req(input$selected_comparison, local_values$available_comparisons)
    
    # Get selected comparison data
    comparison_data <- local_values$available_comparisons[[input$selected_comparison]]
    
    if (is.null(comparison_data)) return(NULL)
    
    # Extract results dataframe
    if ("results_df" %in% names(comparison_data)) {
      results_df <- comparison_data$results_df
    } else {
      results_df <- comparison_data
    }
    
    # Create enhanced volcano plot
    create_context7_volcano_plot(
      results_df = results_df,
      fc_cutoff = input$volcano_fc_cutoff,
      p_cutoff = input$volcano_p_cutoff,
      color_scheme = input$volcano_color_scheme,
      show_labels = input$volcano_show_labels,
      n_labels = input$volcano_n_labels,
      highlight_genes = parse_gene_list(input$volcano_highlight_genes),
      enable_selection = input$volcano_enable_selection
    )
  })
  
  # Enhanced PCA plot
  output$pca_plot <- renderPlotly({
    req(values$filtered_data, values$annotation_data, input$pca_color_by)
    
    if (input$pca_dimensions == "2d") {
      create_context7_pca_2d(
        expression_data = values$filtered_data,
        metadata = values$annotation_data,
        color_by = input$pca_color_by,
        show_ellipses = input$pca_show_ellipses,
        confidence_level = input$pca_confidence_level,
        show_labels = input$pca_show_labels
      )
    } else {
      create_context7_pca_3d(
        expression_data = values$filtered_data,
        metadata = values$annotation_data,
        color_by = input$pca_color_by,
        show_ellipses = input$pca_show_ellipses,
        confidence_level = input$pca_confidence_level
      )
    }
  })
  
  # Enhanced heatmap
  output$heatmap_plot <- renderPlot({
    req(input$selected_comparison, local_values$available_comparisons, values$filtered_data)
    
    # Get comparison data
    comparison_data <- local_values$available_comparisons[[input$selected_comparison]]
    
    create_context7_heatmap(
      expression_data = values$filtered_data,
      results_data = comparison_data,
      annotation_data = values$annotation_data,
      n_genes = input$heatmap_n_genes,
      gene_selection = input$heatmap_gene_selection,
      custom_genes = if (input$heatmap_gene_selection == "custom") {
        strsplit(input$heatmap_custom_genes, "\n")[[1]]
      } else NULL,
      color_scheme = input$heatmap_color_scheme,
      scale_rows = input$heatmap_scale_rows,
      cluster_rows = input$heatmap_cluster_rows,
      cluster_cols = input$heatmap_cluster_cols
    )
  })
  
  # Plot title
  output$plot_title <- renderText({
    plot_titles <- list(
      volcano = "🌋 Enhanced Volcano Plot",
      heatmap = "🔥 Interactive Expression Heatmap",
      pca = paste("📈 Advanced PCA Plot -", input$pca_dimensions),
      ma_plot = "📊 MA Plot",
      expression_profile = "📉 Expression Profile Plot"
    )
    
    title <- plot_titles[[input$plot_type]] %||% "Visualization"
    
    if (input$plot_type != "pca" && !is.null(input$selected_comparison)) {
      title <- paste(title, "-", input$selected_comparison)
    }
    
    return(title)
  })
  
  # Plot statistics
  output$plot_statistics <- renderUI({
    req(input$plot_type)
    
    if (input$plot_type == "volcano" && !is.null(local_values$available_comparisons)) {
      comparison_data <- local_values$available_comparisons[[input$selected_comparison]]
      
      if (!is.null(comparison_data)) {
        results_df <- if ("results_df" %in% names(comparison_data)) {
          comparison_data$results_df
        } else {
          comparison_data
        }
        
        n_total <- nrow(results_df)
        n_sig <- sum(results_df$significant, na.rm = TRUE)
        n_up <- sum(results_df$regulation == "Up", na.rm = TRUE)
        n_down <- sum(results_df$regulation == "Down", na.rm = TRUE)
        
        tagList(
          p(paste("Total genes:", formatC(n_total, format = "d", big.mark = ","))),
          p(paste("Significant:", formatC(n_sig, format = "d", big.mark = ","), 
                  paste0("(", round(n_sig/n_total*100, 1), "%)"))),
          p(paste("Up-regulated:", formatC(n_up, format = "d", big.mark = ","))),
          p(paste("Down-regulated:", formatC(n_down, format = "d", big.mark = ",")))
        )
      }
    } else if (input$plot_type == "pca" && !is.null(values$filtered_data)) {
      n_genes <- nrow(values$filtered_data)
      n_samples <- ncol(values$filtered_data)
      
      tagList(
        p(paste("Genes analyzed:", formatC(n_genes, format = "d", big.mark = ","))),
        p(paste("Samples:", n_samples))
      )
    }
  })
  
  # Refresh plot
  observeEvent(input$refresh_plot, {
    # Force reactive recalculation
    local_values$plot_data <- Sys.time()
  })
}

# Context7-enhanced volcano plot function
create_context7_volcano_plot <- function(results_df, fc_cutoff = 1.0, p_cutoff = 0.05,
                                        color_scheme = "accessible", show_labels = TRUE,
                                        n_labels = 10, highlight_genes = NULL,
                                        enable_selection = TRUE) {
  
  # Get color scheme
  colors <- color_schemes[[color_scheme]] %||% color_schemes$accessible
  
  # Prepare data
  df <- results_df %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      significance_category = case_when(
        padj <= p_cutoff & log2FoldChange >= fc_cutoff ~ "significant_up",
        padj <= p_cutoff & log2FoldChange <= -fc_cutoff ~ "significant_down",
        TRUE ~ "non_significant"
      ),
      fold_change_linear = 2^abs(log2FoldChange)
    )
  
  # Create base plot
  p <- plot_ly(
    data = df,
    x = ~log2FoldChange,
    y = ~neg_log10_padj,
    color = ~significance_category,
    colors = c(
      "non_significant" = colors$non_significant,
      "significant_up" = colors$significant_up,
      "significant_down" = colors$significant_down
    ),
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = ~ifelse(significance_category == "non_significant", 4, 6),
      opacity = ~ifelse(significance_category == "non_significant", 0.6, 0.8),
      line = list(width = 0.5, color = 'rgba(255,255,255,0.5)')
    ),
    text = ~gene,
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "Log₂ Fold Change: %{x:.3f}<br>",
      "-Log₁₀(Adj P-value): %{y:.2f}<br>",
      "Linear FC: %{customdata:.2f}x<br>",
      "<extra></extra>"
    ),
    customdata = ~fold_change_linear,
    name = ""
  ) %>%
    layout(
      title = list(
        text = "",
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Log₂ Fold Change",
        showgrid = TRUE,
        gridcolor = colors$grid,
        showline = TRUE,
        linecolor = 'black',
        mirror = TRUE,
        zeroline = TRUE,
        zerolinecolor = colors$grid,
        zerolinewidth = 2
      ),
      yaxis = list(
        title = "-Log₁₀(Adjusted P-value)",
        showgrid = TRUE,
        gridcolor = colors$grid,
        showline = TRUE,
        linecolor = 'black',
        mirror = TRUE
      ),
      plot_bgcolor = colors$background,
      paper_bgcolor = colors$background,
      showlegend = FALSE,
      hovermode = 'closest'
    )
  
  # Add significance threshold lines
  p <- p %>%
    add_hline(
      y = -log10(p_cutoff),
      line = list(dash = "dash", color = colors$significant_up, width = 2),
      opacity = 0.7
    ) %>%
    add_vline(
      x = fc_cutoff,
      line = list(dash = "dash", color = colors$significant_up, width = 2),
      opacity = 0.7
    ) %>%
    add_vline(
      x = -fc_cutoff,
      line = list(dash = "dash", color = colors$significant_down, width = 2),
      opacity = 0.7
    )
  
  # Add gene labels for top significant genes
  if (show_labels && n_labels > 0) {
    sig_genes <- df %>%
      dplyr::filter(significance_category != "non_significant") %>%
      dplyr::arrange(padj) %>%
      head(n_labels)
    
    if (nrow(sig_genes) > 0) {
      p <- p %>%
        add_trace(
          data = sig_genes,
          x = ~log2FoldChange,
          y = ~neg_log10_padj,
          text = ~gene,
          mode = 'text',
          textposition = 'top center',
          textfont = list(size = 10, color = colors$significant_up),
          showlegend = FALSE,
          hoverinfo = 'skip'
        )
    }
  }
  
  # Highlight specific genes
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    highlight_df <- df %>%
      dplyr::filter(toupper(gene) %in% toupper(highlight_genes))
    
    if (nrow(highlight_df) > 0) {
      p <- p %>%
        add_trace(
          data = highlight_df,
          x = ~log2FoldChange,
          y = ~neg_log10_padj,
          text = ~gene,
          mode = 'markers+text',
          marker = list(
            size = 10,
            color = colors$highlight,
            symbol = 'star',
            line = list(width = 2, color = colors$background)
          ),
          textposition = 'top center',
          textfont = list(size = 10, color = colors$highlight),
          name = "Highlighted Genes",
          showlegend = TRUE,
          hovertemplate = paste(
            "<b>🌟 %{text}</b><br>",
            "Log₂ Fold Change: %{x:.3f}<br>",
            "-Log₁₀(Adj P-value): %{y:.2f}<br>",
            "<extra></extra>"
          )
        )
    }
  }
  
  # Configure plot for publication quality
  p <- p %>%
    config(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c('lasso2d', 'select2d'),
      toImageButtonOptions = list(
        format = 'png',
        filename = 'volcano_plot',
        height = 600,
        width = 800,
        scale = 2
      )
    )
  
  return(p)
}

# Context7-enhanced 2D PCA plot
create_context7_pca_2d <- function(expression_data, metadata, color_by = "Condition",
                                  show_ellipses = TRUE, confidence_level = 0.95,
                                  show_labels = TRUE) {
  
  # Prepare data for PCA
  log_data <- log2(expression_data + 1)
  
  # Select most variable genes
  gene_vars <- apply(log_data, 1, var, na.rm = TRUE)
  top_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, nrow(log_data))]
  
  # Perform PCA
  pca_result <- prcomp(t(log_data[top_genes, ]), center = TRUE, scale. = TRUE)
  
  # Create plot data
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = colnames(expression_data)
  )
  
  # Add metadata
  pca_data <- merge(pca_data, metadata, by = "Sample", all.x = TRUE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Create plot
  p <- plot_ly(
    data = pca_data,
    x = ~PC1,
    y = ~PC2,
    color = as.formula(paste0("~", color_by)),
    type = 'scatter',
    mode = 'markers',
    marker = list(size = 10, opacity = 0.8),
    text = ~Sample,
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "PC1: %{x:.2f}<br>",
      "PC2: %{y:.2f}<br>",
      paste0(color_by, ": %{color}<br>"),
      "<extra></extra>"
    )
  ) %>%
    layout(
      title = "",
      xaxis = list(
        title = paste0("PC1 (", round(var_explained[1], 1), "% variance)")
      ),
      yaxis = list(
        title = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
      ),
      plot_bgcolor = 'white',
      showlegend = TRUE
    )
  
  # Add confidence ellipses
  if (show_ellipses && requireNamespace("car", quietly = TRUE)) {
    groups <- unique(pca_data[[color_by]])
    
    for (group in groups) {
      group_data <- pca_data[pca_data[[color_by]] == group, ]
      
      if (nrow(group_data) >= 4) {
        tryCatch({
          # Calculate confidence ellipse
          ellipse_coords <- car::ellipse(
            center = c(mean(group_data$PC1), mean(group_data$PC2)),
            shape = cov(cbind(group_data$PC1, group_data$PC2)),
            radius = sqrt(qchisq(confidence_level, 2)),
            segments = 100
          )
          
          p <- p %>%
            add_trace(
              x = ellipse_coords[, 1],
              y = ellipse_coords[, 2],
              mode = 'lines',
              line = list(dash = 'dash', width = 2),
              showlegend = FALSE,
              hoverinfo = 'skip',
              name = paste("Ellipse", group)
            )
        }, error = function(e) {
          # Skip ellipse if calculation fails
        })
      }
    }
  }
  
  return(p)
}

# Context7-enhanced 3D PCA plot (based on Emory methodology)
create_context7_pca_3d <- function(expression_data, metadata, color_by = "Condition",
                                  show_ellipses = TRUE, confidence_level = 0.95) {
  
  # Prepare data for PCA (same as 2D)
  log_data <- log2(expression_data + 1)
  gene_vars <- apply(log_data, 1, var, na.rm = TRUE)
  top_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, nrow(log_data))]
  
  # Perform PCA
  pca_result <- prcomp(t(log_data[top_genes, ]), center = TRUE, scale. = TRUE)
  
  # Create plot data with 3 components
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2], 
    PC3 = pca_result$x[, 3],
    Sample = colnames(expression_data)
  )
  
  # Add metadata
  pca_data <- merge(pca_data, metadata, by = "Sample", all.x = TRUE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:3] * 100
  
  # Create 3D plot
  p <- plot_ly(
    data = pca_data,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = as.formula(paste0("~", color_by)),
    type = 'scatter3d',
    mode = 'markers',
    marker = list(size = 8, opacity = 0.8),
    text = ~Sample,
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "PC1: %{x:.2f}<br>",
      "PC2: %{y:.2f}<br>",
      "PC3: %{z:.2f}<br>",
      paste0(color_by, ": %{color}<br>"),
      "<extra></extra>"
    )
  ) %>%
    layout(
      title = "",
      scene = list(
        xaxis = list(title = paste0("PC1 (", round(var_explained[1], 1), "% variance)")),
        yaxis = list(title = paste0("PC2 (", round(var_explained[2], 1), "% variance)")),
        zaxis = list(title = paste0("PC3 (", round(var_explained[3], 1), "% variance)"))
      )
    )
  
  # Add 3D confidence ellipses (simplified version)
  if (show_ellipses && requireNamespace("car", quietly = TRUE)) {
    # Implementation would add 3D ellipses using car::ellipse3d
    # This is complex and would require rgl integration
  }
  
  return(p)
}

# Context7-enhanced heatmap
create_context7_heatmap <- function(expression_data, results_data, annotation_data,
                                   n_genes = 50, gene_selection = "significant",
                                   custom_genes = NULL, color_scheme = "RdBu_r",
                                   scale_rows = TRUE, cluster_rows = TRUE,
                                   cluster_cols = TRUE) {
  
  # Select genes based on method
  if (gene_selection == "custom" && !is.null(custom_genes)) {
    selected_genes <- intersect(custom_genes, rownames(expression_data))
  } else if (gene_selection == "significant") {
    # Get most significant genes from results
    if ("results_df" %in% names(results_data)) {
      sig_genes <- results_data$results_df %>%
        dplyr::filter(significant == TRUE) %>%
        dplyr::arrange(padj) %>%
        head(n_genes)
      selected_genes <- intersect(sig_genes$gene, rownames(expression_data))
    } else {
      selected_genes <- rownames(expression_data)[1:min(n_genes, nrow(expression_data))]
    }
  } else {
    # Highest variance genes
    gene_vars <- apply(expression_data, 1, var, na.rm = TRUE)
    top_var_genes <- order(gene_vars, decreasing = TRUE)[1:n_genes]
    selected_genes <- rownames(expression_data)[top_var_genes]
  }
  
  # Subset expression data
  heatmap_data <- expression_data[selected_genes, , drop = FALSE]
  
  # Create annotation for samples
  sample_annotation <- annotation_data %>%
    dplyr::select(Sample, Condition) %>%
    tibble::column_to_rownames("Sample")
  
  # Ensure sample order matches
  sample_annotation <- sample_annotation[colnames(heatmap_data), , drop = FALSE]
  
  # Create color annotations
  condition_colors <- RColorBrewer::brewer.pal(
    max(3, length(unique(sample_annotation$Condition))), 
    "Set2"
  )
  names(condition_colors) <- unique(sample_annotation$Condition)
  
  annotation_colors <- list(Condition = condition_colors)
  
  # Create heatmap
  pheatmap::pheatmap(
    heatmap_data,
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    color = if (color_scheme == "RdBu_r") {
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    } else {
      RColorBrewer::brewer.pal(11, color_scheme)
    },
    scale = if (scale_rows) "row" else "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 8,
    cellwidth = 12,
    cellheight = 8
  )
}

# Helper function to parse gene list
parse_gene_list <- function(gene_string) {
  if (is.null(gene_string) || gene_string == "") {
    return(NULL)
  }
  
  # Split by comma and clean
  genes <- strsplit(gene_string, ",")[[1]]
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  
  return(genes)
}