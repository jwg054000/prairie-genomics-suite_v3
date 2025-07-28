# Code Display UI Components
# User interface elements for displaying and exporting R code
#
# Author: Prairie Genomics Team
# Date: January 24, 2025  
# Purpose: Provide interactive code visibility throughout the application

library(shiny)
library(shinydashboard)

# Main code display module UI
codeDisplayUI <- function(id, title = "Analysis Code", height = "400px") {
  ns <- NS(id)
  
  tagList(
    wellPanel(
      style = "background-color: #f8f9fa; border: 1px solid #dee2e6;",
      
      fluidRow(
        column(
          8,
          h4(
            icon("code"), 
            title,
            style = "margin-top: 0; color: #495057;"
          )
        ),
        column(
          4,
          div(
            style = "text-align: right; margin-top: 5px;",
            
            # View options dropdown
            selectInput(
              ns("view_option"),
              label = NULL,
              choices = list(
                "Complete Script" = "complete",
                "Current Step Only" = "current",
                "By Category" = "category"
              ),
              selected = "complete",
              width = "150px"
            )
          )
        )
      ),
      
      # Category filter (conditional)
      conditionalPanel(
        condition = paste0("input['", ns("view_option"), "'] == 'category'"),
        
        fluidRow(
          column(
            6,
            selectInput(
              ns("category_filter"),
              "Select Category:",
              choices = list(
                "Data Upload" = "data_upload",
                "Sample Annotation" = "sample_annotation", 
                "DESeq2 Analysis" = "deseq2",
                "Pathway Analysis" = "pathway",
                "Visualization" = "visualization"
              ),
              selected = "deseq2",
              width = "100%"
            )
          )
        )
      ),
      
      # Code display area
      div(
        style = paste0("height: ", height, "; overflow-y: auto; background-color: #ffffff; border: 1px solid #ced4da; border-radius: 4px;"),
        
        # Code content
        verbatimTextOutput(
          ns("code_content"),
          placeholder = FALSE
        )
      ),
      
      br(),
      
      # Action buttons
      fluidRow(
        column(
          6,
          actionButton(
            ns("copy_code"),
            "📋 Copy Code",
            class = "btn-primary btn-sm",
            style = "margin-right: 10px;"
          ),
          
          actionButton(
            ns("refresh_code"),
            "🔄 Refresh",
            class = "btn-secondary btn-sm"
          )
        ),
        column(
          6,
          div(
            style = "text-align: right;",
            
            downloadButton(
              ns("download_r"),
              "📄 Download .R",
              class = "btn-success btn-sm",
              style = "margin-right: 5px;"
            ),
            
            downloadButton(
              ns("download_rmd"),
              "📊 Download .Rmd", 
              class = "btn-info btn-sm"
            )
          )
        )
      )
    )
  )
}

# Inline code snippet UI (for embedding in other tabs)
inlineCodeUI <- function(id, button_text = "View Code", button_class = "btn-outline-primary btn-sm") {
  ns <- NS(id)
  
  tagList(
    # Toggle button
    actionButton(
      ns("toggle_code"),
      paste(icon("code"), button_text),
      class = button_class,
      style = "margin-bottom: 10px;"
    ),
    
    # Collapsible code panel
    conditionalPanel(
      condition = paste0("input['", ns("toggle_code"), "'] % 2 == 1"),
      
      wellPanel(
        style = "background-color: #f8f9fa; border-left: 4px solid #007bff; margin-top: 10px;",
        
        fluidRow(
          column(
            8,
            h5(
              icon("code"),
              textOutput(ns("code_title"), inline = TRUE),
              style = "margin: 0; color: #495057;"
            )
          ),
          column(
            4,
            div(
              style = "text-align: right;",
              actionButton(
                ns("copy_inline"),
                "📋 Copy",
                class = "btn-outline-secondary btn-xs"
              )
            )
          )
        ),
        
        # Code display
        div(
          style = "margin-top: 10px; max-height: 300px; overflow-y: auto;",
          
          verbatimTextOutput(
            ns("inline_code_content")
          )
        ),
        
        # Description
        conditionalPanel(
          condition = paste0("output['", ns("has_description"), "']"),
          
          div(
            style = "margin-top: 10px; padding: 8px; background-color: #e9ecef; border-radius: 4px; font-size: 0.9em;",
            
            strong("Description: "),
            textOutput(ns("code_description"), inline = TRUE)
          )
        )
      )
    )
  )
}

# Code export panel UI
codeExportUI <- function(id) {
  ns <- NS(id)
  
  wellPanel(
    h4(icon("download"), "Export Analysis Code"),
    
    p("Download the complete R code for your analysis to reproduce results independently."),
    
    fluidRow(
      column(
        6,
        
        h5("Export Format:"),
        radioButtons(
          ns("export_format"),
          label = NULL,
          choices = list(
            "R Script (.R)" = "R",
            "R Markdown (.Rmd)" = "Rmd",
            "Complete Report (HTML)" = "html"
          ),
          selected = "R"
        )
      ),
      column(
        6,
        
        h5("Include Options:"),
        checkboxGroupInput(
          ns("include_options"),
          label = NULL,
          choices = list(
            "Comments and descriptions" = "comments",
            "Session information" = "sessioninfo", 
            "Package versions" = "versions",
            "Execution timing" = "timing"
          ),
          selected = c("comments", "sessioninfo", "versions")
        )
      )
    ),
    
    br(),
    
    # Export buttons
    fluidRow(
      column(
        6,
        downloadButton(
          ns("export_complete"),
          "📦 Download Complete Analysis",
          class = "btn-primary",
          style = "width: 100%;"
        )
      ),
      column(
        6,
        downloadButton(
          ns("export_template"),
          "📋 Download Template",
          class = "btn-secondary",
          style = "width: 100%;"
        )
      )
    ),
    
    br(),
    
    # Quick access buttons
    h5("Quick Downloads:"),
    fluidRow(
      column(4,
        downloadButton(
          ns("quick_deseq2"),
          "DESeq2 Only",
          class = "btn-outline-primary btn-sm",
          style = "width: 100%;"
        )
      ),
      column(4,
        downloadButton(
          ns("quick_pathway"),
          "Pathway Only", 
          class = "btn-outline-success btn-sm",
          style = "width: 100%;"
        )
      ),
      column(4,
        downloadButton(
          ns("quick_plots"),
          "Plots Only",
          class = "btn-outline-info btn-sm", 
          style = "width: 100%;"
        )
      )
    )
  )
}

# Code validation display UI
codeValidationUI <- function(id) {
  ns <- NS(id)
  
  div(
    id = ns("validation_panel"),
    style = "margin-top: 15px;",
    
    conditionalPanel(
      condition = paste0("output['", ns("show_validation"), "']"),
      
      wellPanel(
        style = "border-left: 4px solid #28a745; background-color: #f8fff8;",
        
        h5(
          icon("check-circle", style = "color: #28a745;"),
          "Code Validation",
          style = "margin-top: 0; color: #155724;"
        ),
        
        fluidRow(
          column(
            6,
            
            h6("Syntax Check:"),
            div(
              id = ns("syntax_status"),
              style = "padding: 5px; border-radius: 3px; margin-bottom: 10px;",
              
              icon("spinner", class = "fa-spin"),
              " Checking syntax..."
            ),
            
            h6("Dependencies:"),
            div(
              id = ns("deps_status"),
              style = "padding: 5px; border-radius: 3px;",
              
              icon("spinner", class = "fa-spin"),
              " Verifying packages..."
            )
          ),
          column(
            6,
            
            h6("Code Statistics:"),
            tableOutput(ns("code_stats")),
            
            conditionalPanel(
              condition = paste0("output['", ns("has_warnings"), "']"),
              
              h6("Warnings:"),
              div(
                style = "background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 8px; border-radius: 4px; font-size: 0.9em;",
                
                uiOutput(ns("validation_warnings"))
              )
            )
          )
        )
      )
    )
  )
}

# Method documentation UI
methodDocumentationUI <- function(id) {
  ns <- NS(id)
  
  wellPanel(
    style = "background-color: #f8f9fa;",
    
    h4(icon("book"), "Analysis Methods Documentation"),
    
    tabsetPanel(
      id = ns("doc_tabs"),
      
      tabPanel(
        "DESeq2 Workflow",
        br(),
        div(
          h5("Differential Expression Analysis"),
          p("The DESeq2 analysis follows the standard workflow described in Love et al. (2014):"),
          
          tags$ol(
            tags$li(strong("Data Preparation:"), " Count matrix normalization and size factor estimation"),
            tags$li(strong("Dispersion Estimation:"), " Gene-wise and fitted dispersions calculated"),
            tags$li(strong("Statistical Testing:"), " Negative binomial GLM with Wald test"),
            tags$li(strong("Multiple Testing Correction:"), " Benjamini-Hochberg FDR adjustment")
          ),
          
          h6("Key Parameters:"),
          tags$ul(
            tags$li("Design formula: Specified by user based on experimental design"),
            tags$li("Significance threshold: padj < 0.05 (adjustable)"),
            tags$li("Fold change threshold: |log2FC| > 1 (adjustable)")
          ),
          
          h6("References:"),
          tags$ul(
            tags$li("Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550")
          )
        )
      ),
      
      tabPanel(
        "Pathway Analysis",
        br(),
        div(
          h5("Gene Set Enrichment Methods"),
          
          h6("Gene Ontology (GO) Analysis:"),
          p("Over-representation analysis using the hypergeometric test to identify enriched biological processes, molecular functions, and cellular components."),
          
          h6("KEGG Pathway Analysis:"),
          p("Pathway enrichment analysis using curated KEGG pathway gene sets with Fisher's exact test."),
          
          h6("Gene Set Enrichment Analysis (GSEA):"),
          p("Rank-based analysis using the Kolmogorov-Smirnov statistic to identify coordinated changes in predefined gene sets."),
          
          tags$ul(
            tags$li("Ranking metric: Signed p-value (sign(log2FC) × -log10(pvalue))"),
            tags$li("Gene set size: 15-500 genes"),
            tags$li("Permutation: 1000 permutations (fgseaMultilevel)")
          ),
          
          h6("References:"),
          tags$ul(
            tags$li("Subramanian, A. et al. (2005) Gene set enrichment analysis: a knowledge-based approach. PNAS 102:15545-15550"),
            tags$li("Yu, G. et al. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16:284-287")
          )
        )
      ),
      
      tabPanel(
        "Statistical Methods",
        br(),
        div(
          h5("Statistical Approaches Used"),
          
          h6("Multiple Testing Correction:"),
          p("All p-values are adjusted for multiple testing using the Benjamini-Hochberg method to control the false discovery rate (FDR)."),
          
          h6("Effect Size Considerations:"),
          p("In addition to statistical significance (padj < 0.05), biological significance is assessed using fold change thresholds (typically |log2FC| > 1)."),
          
          h6("Sample Size and Power:"),
          p("DESeq2 provides robust results with small sample sizes (n ≥ 3 per group) through empirical Bayes shrinkage of dispersion estimates."),
          
          h6("Assumptions:"),
          tags$ul(
            tags$li("Count data follows negative binomial distribution"),
            tags$li("Independent observations within groups"),
            tags$li("Adequate sequencing depth across samples")
          )
        )
      )
    )
  )
}

cat("✅ Code display UI components loaded\n")
cat("📋 UI modules: codeDisplayUI(), inlineCodeUI(), codeExportUI(), methodDocumentationUI()\n")