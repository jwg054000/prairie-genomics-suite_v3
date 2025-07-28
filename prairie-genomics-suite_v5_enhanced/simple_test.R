# Simple test to verify R Shiny setup
# Tests core functionality without complex dependencies

cat("🧬 Simple R Shiny Test\n")
cat("====================\n")

# Test core packages only
core_packages <- c("shiny", "shinydashboard", "DT", "ggplot2", "dplyr")

for (pkg in core_packages) {
  cat(paste0("Loading ", pkg, "... "))
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("OK\n")
  } else {
    cat("FAILED\n")
    quit(status = 1)
  }
}

# Create simple UI
ui <- dashboardPage(
  dashboardHeader(title = "Prairie Genomics Test"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Test", tabName = "test")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "test",
        h2("Prairie Genomics Suite - R Shiny"),
        p("✅ Basic Shiny functionality working!"),
        fluidRow(
          box(
            title = "Test Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("test_table")
          )
        )
      )
    )
  )
)

# Simple server
server <- function(input, output, session) {
  output$test_table <- DT::renderDataTable({
    test_data <- data.frame(
      Gene = paste0("Gene_", 1:5),
      Expression = c(100, 150, 200, 75, 300),
      Condition = c("Control", "Treatment", "Control", "Treatment", "Control")
    )
    
    DT::datatable(test_data, options = list(pageLength = 10))
  })
}

cat("✅ UI and server created successfully\n")
cat("🚀 Ready to run with: shiny::runApp(list(ui = ui, server = server))\n")

# Test the app object creation
app <- list(ui = ui, server = server)
cat("✅ App object created\n")
cat("✅ All tests passed!\n")