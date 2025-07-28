# Phase 1 Demo App - Modern UI Components Showcase
# This app demonstrates the new Phase 1 components in action
# Run this to see the visual improvements before integrating

library(shiny)
library(shinydashboard)
library(DT)

# Load Phase 1 components
source("phase1/components/modern_ui_components.R")
ui_components <- source("phase1/components/modern_ui_components.R")$value

# Demo UI
ui <- fluidPage(
  title = "🧬 Prairie Genomics Suite - Phase 1 Demo",
  
  # Include modern CSS and JavaScript
  includeCSS("www/css/modern_components.css"),
  includeScript("www/js/modern_interactions.js"),
  
  # Custom styling for demo
  tags$style(HTML("
    body { background-color: var(--gray-50); font-family: var(--font-family-sans); }
    .demo-header { text-align: center; padding: 2rem 0; background: white; margin-bottom: 2rem; box-shadow: var(--shadow-sm); }
    .demo-section { margin-bottom: 3rem; }
    .demo-title { color: var(--gray-900); font-size: 1.5rem; font-weight: 600; margin-bottom: 1rem; }
  ")),
  
  # Demo Header
  div(
    class = "demo-header",
    h1("🧬 Prairie Genomics Suite", style = "color: var(--primary-600); margin-bottom: 0.5rem;"),
    h2("Phase 1 Modern UI Components Demo", style = "color: var(--gray-600); font-size: 1.25rem; font-weight: 400;"),
    p("Experience the enhanced user interface with modern design components", 
      style = "color: var(--gray-500); margin-top: 1rem;")
  ),
  
  div(
    class = "container-modern",
    
    # Stats Cards Section
    div(
      class = "demo-section",
      h3("📊 Statistics Cards", class = "demo-title"),
      fluidRow(
        ui_components$stats_card(
          title = "Total Genes",
          value = "15,357",
          subtitle = "Analyzed",
          icon = "🧬",
          color = "primary",
          change = 5.2,
          width = 3
        ),
        ui_components$stats_card(
          title = "Significant",
          value = "1,247",
          subtitle = "Differentially expressed",
          icon = "⭐",
          color = "success",
          change = 12.8,
          width = 3
        ),
        ui_components$stats_card(
          title = "Upregulated",
          value = "623",
          subtitle = "Increased expression",
          icon = "📈",
          color = "success",
          change = -2.1,
          width = 3
        ),
        ui_components$stats_card(
          title = "Downregulated",
          value = "624",
          subtitle = "Decreased expression",
          icon = "📉",
          color = "error",
          change = 8.7,
          width = 3
        )
      )
    ),
    
    # Modern Cards Section
    div(
      class = "demo-section",
      h3("🃏 Modern Card Components", class = "demo-title"),
      fluidRow(
        ui_components$modern_card(
          title = "DESeq2 Analysis",
          subtitle = "Differential expression analysis results",
          body = div(
            p("Analysis completed successfully with 15,357 genes tested.", 
              style = "margin-bottom: 1rem;"),
            ui_components$modern_progress(
              "demo_progress_1",
              label = "Analysis Progress",
              percentage = 100,
              message = "Complete! Found 1,247 significant genes."
            )
          ),
          footer = div(
            tags$button(
              class = "btn-modern btn-modern-primary btn-modern-sm",
              "View Results"
            ),
            tags$button(
              class = "btn-modern btn-modern-secondary btn-modern-sm",
              style = "margin-left: 0.5rem;",
              "Export Data"
            )
          ),
          width = 6
        ),
        
        ui_components$modern_card(
          title = "Pathway Analysis",
          subtitle = "Gene ontology and pathway enrichment",
          body = div(
            p("Pathway analysis identified enriched biological processes.", 
              style = "margin-bottom: 1rem;"),
            ui_components$modern_progress(
              "demo_progress_2",
              label = "GO Analysis",
              percentage = 75,
              message = "Processing 2,847 pathways..."
            )
          ),
          footer = div(
            tags$button(
              class = "btn-modern btn-modern-success btn-modern-sm",
              "View Pathways"
            )
          ),
          width = 6
        )
      )
    ),
    
    # Form Components Section
    div(
      class = "demo-section",
      h3("📝 Enhanced Form Components", class = "demo-title"),
      fluidRow(
        ui_components$modern_card(
          title = "Analysis Parameters",
          body = div(
            fluidRow(
              column(6,
                ui_components$modern_input(
                  "demo_padj",
                  "Adjusted P-value Cutoff",
                  value = "0.05",
                  placeholder = "Enter p-value (e.g., 0.05)",
                  type = "number",
                  help_text = "FDR-corrected p-value threshold for significance"
                )
              ),
              column(6,
                ui_components$modern_input(
                  "demo_fc",
                  "Log2 Fold Change Cutoff",
                  value = "1.0",
                  placeholder = "Enter fold change (e.g., 1.0)",
                  type = "number",
                  help_text = "Minimum absolute log2 fold change"
                )
              )
            ),
            fluidRow(
              column(6,
                ui_components$modern_select(
                  "demo_species",
                  "Species",
                  choices = list(
                    "Human (Homo sapiens)" = "human",
                    "Mouse (Mus musculus)" = "mouse",
                    "Rat (Rattus norvegicus)" = "rat"
                  ),
                  selected = "human",
                  help_text = "Select organism for gene annotation"
                )
              ),
              column(6,
                ui_components$modern_select(
                  "demo_contrast",
                  "Comparison",
                  choices = list(
                    "Treatment vs Control" = "treatment_control",
                    "Time Point 2 vs 1" = "time2_time1",
                    "Condition A vs B" = "condA_condB"
                  ),
                  help_text = "Select experimental comparison"
                )
              )
            )
          ),
          footer = div(
            tags$button(
              class = "btn-modern btn-modern-primary",
              "🚀 Run Analysis"
            ),
            tags$button(
              class = "btn-modern btn-modern-ghost",
              style = "margin-left: 0.5rem;",
              "Reset Parameters"
            )
          ),
          width = 8
        ),
        
        # Alerts Section
        div(
          class = "col-md-4",
          ui_components$modern_alert(
            title = "Analysis Complete",
            message = "DESeq2 analysis finished successfully. 1,247 genes found to be significantly differentially expressed.",
            type = "success",
            dismissible = TRUE
          ),
          
          ui_components$modern_alert(
            title = "Processing",
            message = "Pathway analysis is currently running. This may take a few minutes.",
            type = "info",
            dismissible = FALSE
          ),
          
          ui_components$modern_alert(
            title = "High Memory Usage",
            message = "Large dataset detected. Consider reducing dataset size if performance issues occur.",
            type = "warning",
            dismissible = TRUE
          )
        )
      )
    ),
    
    # Data Table Section
    div(
      class = "demo-section",
      h3("📋 Enhanced Data Tables", class = "demo-title"),
      ui_components$modern_card(
        title = "Top Significant Genes",
        subtitle = "Interactive results table with modern styling",
        body = div(
          # Create demo data
          DT::dataTableOutput("demo_table")
        ),
        width = 12
      )
    ),
    
    # Authentication Demo Section
    div(
      class = "demo-section",
      h3("🔐 Authentication Components (Optional)", class = "demo-title"),
      fluidRow(
        ui_components$modern_card(
          title = "User Authentication",
          body = ui_components$auth_login_form("demo_auth"),
          width = 6
        ),
        ui_components$modern_card(
          title = "Create Account",
          body = ui_components$auth_register_form("demo_auth"),
          width = 6
        )
      )
    )
  ),
  
  # Footer
  div(
    style = "text-align: center; padding: 2rem; color: var(--gray-500); border-top: 1px solid var(--gray-200); margin-top: 3rem; background: white;",
    p("🧬 Prairie Genomics Suite - Phase 1 Enhanced UI Components"),
    p("Built with modern design principles and enhanced user experience", style = "font-size: 0.875rem;")
  )
)

# Demo Server
server <- function(input, output, session) {
  
  # Demo data table
  output$demo_table <- DT::renderDataTable({
    # Create realistic demo data
    demo_data <- data.frame(
      Gene = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "PTEN", "RB1"),
      `Gene Symbol` = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "PTEN", "RB1"),
      `Ensembl ID` = c("ENSG00000012048", "ENSG00000139618", "ENSG00000141510", 
                      "ENSG00000146648", "ENSG00000136997", "ENSG00000133703",
                      "ENSG00000121879", "ENSG00000142208", "ENSG00000171862", "ENSG00000139687"),
      `Log2 FC` = c(2.3, -1.8, 3.1, -2.7, 1.9, 2.8, -1.5, 2.1, -3.2, 1.7),
      `P-value` = c(0.001, 0.003, 0.0001, 0.002, 0.008, 0.0005, 0.012, 0.004, 0.0002, 0.015),
      `Adj P-value` = c(0.02, 0.04, 0.005, 0.03, 0.08, 0.01, 0.09, 0.06, 0.008, 0.12),
      Regulation = c("Up", "Down", "Up", "Down", "Up", "Up", "Down", "Up", "Down", "Up"),
      check.names = FALSE
    )
    
    # Format the data with badges for regulation
    demo_data$Regulation <- sapply(demo_data$Regulation, function(reg) {
      if (reg == "Up") {
        '<span class="badge-modern badge-modern-success">Up</span>'
      } else if (reg == "Down") {
        '<span class="badge-modern badge-modern-error">Down</span>'
      } else {
        '<span class="badge-modern badge-modern-gray">NS</span>'
      }
    })
    
    # Use modern data table
    ui_components$modern_datatable(
      demo_data,
      options = list(
        pageLength = 5,
        dom = "tip",  # Simplified for demo
        columnDefs = list(
          list(className = "dt-center", targets = c(3, 4, 5, 6)),
          list(width = "100px", targets = 6)
        )
      )
    )
  })
  
  # Add some interactivity for demo purposes
  observeEvent(input$demo_auth_email, {
    if (!is.null(input$demo_auth_email) && input$demo_auth_email != "") {
      showNotification(
        "Demo: Email field updated",
        type = "message",
        duration = 2
      )
    }
  })
}

# Run the demo app
shinyApp(ui = ui, server = server)