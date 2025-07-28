# Visualization Module for Prairie Genomics Suite
# Handles interactive plots and visualizations
# 
# Fixed: Added ggrepel package handling with fallbacks to prevent 
# "could not find function 'geom_text_repel'" errors

# Helper function to check if ggrepel is available
check_ggrepel <- function() {
  if (exists("ggrepel_available") && ggrepel_available) {
    return(requireNamespace("ggrepel", quietly = TRUE))
  }
  return(FALSE)
}

# UI Function
visualizationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        3,
        wellPanel(
          h4("🎨 Plot Controls"),
          
          # Plot type selection
          radioButtons(
            ns("plot_type"),
            "Select Visualization:",
            choices = list(
              "🌋 Volcano Plot" = "volcano",
              "🔥 Heatmap" = "heatmap", 
              "📊 PCA Plot" = "pca",
              "📈 MA Plot" = "ma"
            ),
            selected = "volcano"
          ),
          
          hr(),
          
          # Volcano plot controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
            
            h5("Volcano Plot Settings"),
            
            numericInput(
              ns("volcano_padj"),
              "P-value cutoff:",
              value = 0.05,
              min = 0.001,
              max = 0.1,
              step = 0.001
            ),
            
            numericInput(
              ns("volcano_fc"),
              "Fold change cutoff:",
              value = 1.0,
              min = 0.1,
              max = 5.0,
              step = 0.1
            ),
            
            selectInput(
              ns("volcano_color"),
              "Color scheme:",
              choices = list(
                "Default" = "default",
                "Colorblind friendly" = "colorblind",
                "Publication" = "publication"
              ),
              selected = "default"
            ),
            
            textInput(
              ns("highlight_genes"),
              "Highlight genes (comma-separated):",
              placeholder = "e.g., TP53, BRCA1, MYC"
            ),
            
            checkboxInput(
              ns("show_labels"),
              "Show gene labels",
              value = TRUE
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("show_labels"), "']"),
              numericInput(
                ns("n_labels"),
                "Number of labels:",
                value = 10,
                min = 0,
                max = 50,
                step = 1
              )
            )
          ),
          
          # Heatmap controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
            
            h5("Heatmap Settings"),
            
            numericInput(
              ns("heatmap_top_n"),
              "Top N genes:",
              value = 50,
              min = 10,
              max = 500,
              step = 10
            ),
            
            selectInput(
              ns("heatmap_clustering"),
              "Clustering:",
              choices = list(
                "Both" = "both",
                "Rows only" = "row",
                "Columns only" = "column",
                "None" = "none"
              ),
              selected = "both"
            ),
            
            selectInput(
              ns("heatmap_scale"),
              "Scale data:",
              choices = list(
                "By row (Z-score)" = "row",
                "By column" = "column", 
                "None" = "none"
              ),
              selected = "row"
            ),
            
            selectInput(
              ns("heatmap_color"),
              "Color palette:",
              choices = list(
                "Blue-Red" = "RdBu",
                "Blue-White-Red" = "RdBu_r",
                "Viridis" = "viridis",
                "Plasma" = "plasma"
              ),
              selected = "RdBu_r"
            )
          ),
          
          # PCA controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'pca'"),
            
            h5("PCA Settings"),
            
            selectInput(
              ns("pca_dims"),
              "Dimensions:",
              choices = list(
                "PC1 vs PC2" = "1_2",
                "PC1 vs PC3" = "1_3",
                "PC2 vs PC3" = "2_3"
              ),
              selected = "1_2"
            ),
            
            checkboxInput(
              ns("pca_labels"),
              "Show sample labels",
              value = TRUE
            ),
            
            checkboxInput(
              ns("pca_ellipses"),
              "Show confidence ellipses",
              value = TRUE
            )
          ),
          
          hr(),
          
          # Export controls
          h5("📥 Export"),
          
          selectInput(
            ns("export_format"),
            "Format:",
            choices = list(
              "PNG" = "png",
              "PDF" = "pdf",
              "SVG" = "svg"
            ),
            selected = "png"
          ),
          
          numericInput(
            ns("export_width"),
            "Width (inches):",
            value = 8,
            min = 4,
            max = 20,
            step = 1
          ),
          
          numericInput(
            ns("export_height"),
            "Height (inches):",
            value = 6,
            min = 4,
            max = 20,
            step = 1
          ),
          
          downloadButton(
            ns("download_plot"),
            "Download Plot",
            class = "btn-primary",
            style = "width: 100%;"
          )
        )
      ),
      
      column(
        9,
        # Plot display area
        div(
          style = "border: 1px solid #ddd; border-radius: 5px; padding: 15px; background: white;",
          
          # Plot title
          h4(textOutput(ns("plot_title")), style = "text-align: center; margin-bottom: 20px;"),
          
          # Plot output
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
            plotlyOutput(ns("volcano_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
            plotOutput(ns("heatmap_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'pca'"),
            plotlyOutput(ns("pca_plot"), height = "600px")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'ma'"),
            plotlyOutput(ns("ma_plot"), height = "600px")
          )
        ),
        
        br(),
        
        # Plot statistics
        wellPanel(
          h5("📊 Plot Statistics"),
          fluidRow(
            column(3, valueBoxOutput(ns("stat1"), width = 12)),
            column(3, valueBoxOutput(ns("stat2"), width = 12)),
            column(3, valueBoxOutput(ns("stat3"), width = 12)),
            column(3, valueBoxOutput(ns("stat4"), width = 12))
          )
        )
      )
    )
  )
}

# Server Function
visualization <- function(input, output, session, values) {
  ns <- session$ns
  
  # Reactive values
  local_values <- reactiveValues(
    current_plot = NULL,
    plot_data = NULL
  )
  
  # Plot title
  output$plot_title <- renderText({
    switch(
      input$plot_type,
      "volcano" = "Interactive Volcano Plot",
      "heatmap" = "Expression Heatmap",
      "pca" = "Principal Component Analysis",
      "ma" = "MA Plot"
    )
  })
  
  # Volcano plot
  output$volcano_plot <- renderPlotly({
    req(values$deseq2_results)
    
    results_df <- values$deseq2_results
    
    # Apply current thresholds
    results_df$significant <- !is.na(results_df$padj) & 
                            results_df$padj < input$volcano_padj & 
                            abs(results_df$log2FoldChange) > input$volcano_fc
    
    # Color points
    colors <- get_volcano_colors(input$volcano_color)
    
    results_df$color <- ifelse(
      results_df$significant,
      ifelse(results_df$log2FoldChange > 0, colors$up, colors$down),
      colors$ns
    )
    
    # Create base plot
    p <- ggplot(results_df, aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = I(color),
      text = paste0(
        "Gene: ", ifelse(!is.na(display_name), display_name, gene), "<br>",
        ifelse(!is.na(gene_symbol) & gene_symbol != display_name, paste0("Ensembl: ", gene, "<br>"), ""),
        "Log2 FC: ", round(log2FoldChange, 3), "<br>",
        "Adj. P-value: ", formatC(padj, format = "e", digits = 2), "<br>",
        "P-value: ", formatC(pvalue, format = "e", digits = 2)
      )
    )) +
      geom_point(alpha = 0.7, size = 1) +
      
      # Add threshold lines
      geom_hline(yintercept = -log10(input$volcano_padj), 
                linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = c(-input$volcano_fc, input$volcano_fc), 
                linetype = "dashed", color = "gray50") +
      
      labs(
        x = "Log2 Fold Change",
        y = "-Log10(Adjusted P-value)"
      ) +
      
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    
    # Add gene labels if requested
    if (input$show_labels && input$n_labels > 0) {
      top_genes <- results_df[results_df$significant == TRUE, ]
      top_genes <- head(top_genes[order(top_genes$padj), ], input$n_labels)
      
      if (nrow(top_genes) > 0) {
        # Check if ggrepel is available
        if (check_ggrepel()) {
          p <- p + ggrepel::geom_text_repel(
            data = top_genes,
            aes(label = ifelse(!is.na(display_name), display_name, gene)),
            size = 3,
            color = "black",
            min.segment.length = 0.1,
            max.overlaps = 20
          )
        } else {
          # Fallback to basic geom_text if ggrepel is not available
          p <- p + geom_text(
            data = top_genes,
            aes(label = ifelse(!is.na(display_name), display_name, gene)),
            size = 3,
            color = "black",
            vjust = -0.5,
            hjust = 0.5
          )
          # Show notification once per session
          if (!exists("ggrepel_warning_shown") || !ggrepel_warning_shown) {
            showNotification(
              "⚠️ ggrepel package not available. Using basic text labels (may overlap).",
              type = "warning",
              duration = 5
            )
            assign("ggrepel_warning_shown", TRUE, envir = .GlobalEnv)
          }
        }
      }
    }
    
    # Highlight specific genes
    if (!is.null(input$highlight_genes) && input$highlight_genes != "") {
      highlight_list <- trimws(strsplit(input$highlight_genes, ",")[[1]])
      # Search in both gene IDs and symbols
      highlight_data <- results_df[
        results_df$gene %in% highlight_list | 
        (!is.na(results_df$gene_symbol) & results_df$gene_symbol %in% highlight_list) |
        (!is.na(results_df$display_name) & results_df$display_name %in% highlight_list), ]
      
      if (nrow(highlight_data) > 0) {
        p <- p + 
          geom_point(
            data = highlight_data,
            aes(x = log2FoldChange, y = -log10(padj)),
            color = "purple",
            size = 3,
            shape = 21,
            fill = "yellow",
            stroke = 2
          ) +
          # Add labels for highlighted genes
          if (check_ggrepel()) {
            ggrepel::geom_text_repel(
              data = highlight_data,
              aes(label = ifelse(!is.na(display_name), display_name, gene)),
              color = "purple",
              fontface = "bold",
              size = 4
            )
          } else {
            geom_text(
              data = highlight_data,
              aes(label = ifelse(!is.na(display_name), display_name, gene)),
              color = "purple",
              fontface = "bold",
              size = 4,
              vjust = -0.5,
              hjust = 0.5
            )
          }
      }
    }
    
    local_values$current_plot <- p
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        showlegend = FALSE,
        hovermode = "closest"
      )
  })
  
  # Heatmap plot
  output$heatmap_plot <- renderPlot({
    req(values$deseq2_results)
    
    # Use optimized data access if available
    if (!is.null(values$deseq2_data)) {
      results_df <- values$deseq2_data$results_df
      norm_counts <- values$deseq2_data$normalized_counts
    } else {
      results_df <- values$deseq2_results
      norm_counts <- values$filtered_data
    }
    
    req(norm_counts)
    
    # Get top significant genes
    sig_genes <- results_df[results_df$significant == TRUE, ]
    sig_genes <- head(sig_genes[order(sig_genes$padj), ], input$heatmap_top_n)
    
    if (nrow(sig_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "No significant genes found", cex = 1.5)
      return()
    }
    
    # Get expression data for these genes
    expr_data <- norm_counts[sig_genes$gene, ]
    
    # Use display names for row names if available
    if (!is.null(sig_genes$display_name)) {
      display_names <- ifelse(!is.na(sig_genes$display_name), sig_genes$display_name, sig_genes$gene)
      rownames(expr_data) <- display_names
    }
    
    # Transform data if requested
    if (input$heatmap_scale == "row") {
      expr_data <- t(scale(t(expr_data)))
    } else if (input$heatmap_scale == "column") {
      expr_data <- scale(expr_data)
    }
    
    # Set up clustering
    cluster_rows <- input$heatmap_clustering %in% c("both", "row")
    cluster_cols <- input$heatmap_clustering %in% c("both", "column")
    
    # Create annotation for samples
    annotation_col <- data.frame(
      Condition = values$annotation_data$Condition[match(colnames(expr_data), values$annotation_data$Sample)]
    )
    rownames(annotation_col) <- colnames(expr_data)
    
    # Create heatmap
    pheatmap(
      expr_data,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_rownames = nrow(expr_data) <= 50,
      show_colnames = TRUE,
      annotation_col = annotation_col,
      color = get_heatmap_colors(input$heatmap_color),
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10
    )
  })
  
  # PCA plot
  output$pca_plot <- renderPlotly({
    req(values$annotation_data)
    
    # Use optimized data access if available
    if (!is.null(values$deseq2_data)) {
      norm_counts <- values$deseq2_data$normalized_counts
    } else {
      norm_counts <- values$filtered_data
    }
    
    req(norm_counts)
    
    # Perform PCA
    pca_data <- t(norm_counts)
    pca_result <- prcomp(pca_data, scale. = TRUE)
    
    # Get variance explained
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
    
    # Extract dimensions
    dims <- as.integer(strsplit(input$pca_dims, "_")[[1]])
    pc_x <- pca_result$x[, dims[1]]
    pc_y <- pca_result$x[, dims[2]]
    
    # Create plot data
    plot_data <- data.frame(
      Sample = rownames(pca_result$x),
      PC_X = pc_x,
      PC_Y = pc_y,
      Condition = values$annotation_data$Condition[match(rownames(pca_result$x), values$annotation_data$Sample)]
    )
    
    # Create plot
    p <- ggplot(plot_data, aes(
      x = PC_X,
      y = PC_Y,
      color = Condition,
      text = paste0("Sample: ", Sample, "<br>Condition: ", Condition)
    )) +
      geom_point(size = 3, alpha = 0.8) +
      
      labs(
        x = paste0("PC", dims[1], " (", var_explained[dims[1]], "% variance)"),
        y = paste0("PC", dims[2], " (", var_explained[dims[2]], "% variance)")
      ) +
      
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank()
      )
    
    # Add confidence ellipses
    if (input$pca_ellipses) {
      p <- p + stat_ellipse(aes(color = Condition), alpha = 0.3)
    }
    
    # Add sample labels
    if (input$pca_labels) {
      if (check_ggrepel()) {
        p <- p + ggrepel::geom_text_repel(
          aes(label = Sample),
          size = 3,
          show.legend = FALSE
        )
      } else {
        p <- p + geom_text(
          aes(label = Sample),
          size = 3,
          show.legend = FALSE,
          vjust = -0.5,
          hjust = 0.5
        )
      }
    }
    
    local_values$current_plot <- p
    
    ggplotly(p, tooltip = "text")
  })
  
  # MA plot
  output$ma_plot <- renderPlotly({
    req(values$deseq2_results)
    
    results_df <- values$deseq2_results
    results_df$significant <- !is.na(results_df$padj) & 
                            results_df$padj < 0.05 & 
                            abs(results_df$log2FoldChange) > 1
    
    # Calculate A values (average expression)
    results_df$A <- log10(results_df$baseMean + 1)
    results_df$M <- results_df$log2FoldChange
    
    # Color points
    results_df$color <- ifelse(
      results_df$significant,
      ifelse(results_df$M > 0, "red", "blue"),
      "gray"
    )
    
    p <- ggplot(results_df, aes(
      x = A,
      y = M,
      color = I(color),
      text = paste0(
        "Gene: ", ifelse(!is.na(display_name), display_name, gene), "<br>",
        ifelse(!is.na(gene_symbol) & gene_symbol != display_name, paste0("Ensembl: ", gene, "<br>"), ""),
        "Log2 FC: ", round(M, 3), "<br>",
        "Mean Expression: ", round(baseMean, 1), "<br>",
        "Adj. P-value: ", formatC(padj, format = "e", digits = 2)
      )
    )) +
      geom_point(alpha = 0.6, size = 1) +
      
      labs(
        x = "Log10(Mean Expression)",
        y = "Log2 Fold Change"
      ) +
      
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    
    local_values$current_plot <- p
    
    ggplotly(p, tooltip = "text")
  })
  
  # Helper function for volcano plot colors
  get_volcano_colors <- function(scheme) {
    switch(
      scheme,
      "default" = list(up = "red", down = "blue", ns = "gray"),
      "colorblind" = list(up = "#E69F00", down = "#56B4E9", ns = "#999999"),
      "publication" = list(up = "#D62728", down = "#1F77B4", ns = "#7F7F7F")
    )
  }
  
  # Helper function for heatmap colors
  get_heatmap_colors <- function(palette) {
    switch(
      palette,
      "RdBu" = colorRampPalette(brewer.pal(11, "RdBu"))(100),
      "RdBu_r" = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
      "viridis" = viridis::viridis(100),
      "plasma" = viridis::plasma(100)
    )
  }
  
  # Plot statistics
  output$stat1 <- renderValueBox({
    if (input$plot_type == "volcano" && !is.null(values$deseq2_results)) {
      sig_count <- sum(values$deseq2_results$padj < input$volcano_padj & 
                      abs(values$deseq2_results$log2FoldChange) > input$volcano_fc, 
                      na.rm = TRUE)
      valueBox(
        value = formatC(sig_count, format = "d", big.mark = ","),
        subtitle = "Significant Genes",
        icon = icon("star"),
        color = "yellow"
      )
    } else if (input$plot_type == "heatmap") {
      valueBox(
        value = input$heatmap_top_n,
        subtitle = "Genes Shown",
        icon = icon("fire"),
        color = "orange"
      )
    } else {
      valueBox(value = "0", subtitle = "Stat 1", icon = icon("info"))
    }
  })
  
  output$stat2 <- renderValueBox({
    if (input$plot_type == "volcano" && !is.null(values$deseq2_results)) {
      up_count <- sum(values$deseq2_results$padj < input$volcano_padj & 
                     values$deseq2_results$log2FoldChange > input$volcano_fc, 
                     na.rm = TRUE)
      valueBox(
        value = formatC(up_count, format = "d", big.mark = ","),
        subtitle = "Upregulated",
        icon = icon("arrow-up"),
        color = "green"
      )
    } else {
      valueBox(value = "0", subtitle = "Stat 2", icon = icon("info"))
    }
  })
  
  output$stat3 <- renderValueBox({
    if (input$plot_type == "volcano" && !is.null(values$deseq2_results)) {
      down_count <- sum(values$deseq2_results$padj < input$volcano_padj & 
                       values$deseq2_results$log2FoldChange < -input$volcano_fc, 
                       na.rm = TRUE)
      valueBox(
        value = formatC(down_count, format = "d", big.mark = ","),
        subtitle = "Downregulated",
        icon = icon("arrow-down"),
        color = "red"
      )
    } else {
      valueBox(value = "0", subtitle = "Stat 3", icon = icon("info"))
    }
  })
  
  output$stat4 <- renderValueBox({
    if (!is.null(values$deseq2_results)) {
      total_count <- nrow(values$deseq2_results)
      valueBox(
        value = formatC(total_count, format = "d", big.mark = ","),
        subtitle = "Total Genes",
        icon = icon("dna"),
        color = "blue"
      )
    } else {
      valueBox(value = "0", subtitle = "Stat 4", icon = icon("info"))
    }
  })
  
  # Download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(input$plot_type, "_plot_", Sys.Date(), ".", input$export_format)
    },
    content = function(file) {
      if (!is.null(local_values$current_plot)) {
        ggsave(
          file,
          plot = local_values$current_plot,
          width = input$export_width,
          height = input$export_height,
          dpi = 300,
          device = input$export_format
        )
      }
    }
  )
}