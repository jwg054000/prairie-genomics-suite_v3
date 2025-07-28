# Phase 1 - Modern UI Components for Prairie Genomics Suite
# Reusable UI components with modern design and enhanced functionality
# 
# This module provides:
# - Modern card layouts with animations
# - Enhanced form components with validation
# - Progress indicators and loading states
# - Alert and notification systems
# - Interactive data tables and visualizations

library(shiny)
library(shinydashboard)
library(DT)

# ===========================================
# MODERN CARD COMPONENTS
# ===========================================

# Modern card with header, body, and footer
modern_card <- function(title = NULL, subtitle = NULL, body = NULL, footer = NULL,
                       status = "primary", width = 12, collapsible = FALSE,
                       id = NULL, class = "card-modern") {
  
  card_id <- id %||% paste0("card_", sample(1000:9999, 1))
  
  tagList(
    tags$div(
      class = paste("col-md", width),
      tags$div(
        id = card_id,
        class = class,
        `data-status` = status,
        
        # Card Header
        if (!is.null(title) || !is.null(subtitle)) {
          tags$div(
            class = "card-modern-header",
            if (!is.null(title)) {
              tags$h3(class = "card-modern-title", title)
            },
            if (!is.null(subtitle)) {
              tags$p(class = "card-modern-subtitle", subtitle)
            },
            if (collapsible) {
              tags$button(
                class = "btn-modern btn-modern-ghost btn-modern-sm card-collapse-btn",
                `data-target` = paste0("#", card_id, "_body"),
                "⌄"
              )
            }
          )
        },
        
        # Card Body
        tags$div(
          id = paste0(card_id, "_body"),
          class = "card-modern-body",
          body
        ),
        
        # Card Footer
        if (!is.null(footer)) {
          tags$div(
            class = "card-modern-footer",
            footer
          )
        }
      )
    )
  )
}

# Statistics card with number display
stats_card <- function(title, value, subtitle = NULL, icon = NULL, 
                      color = "primary", change = NULL, width = 3) {
  
  change_element <- if (!is.null(change)) {
    change_class <- if (change >= 0) "text-success" else "text-danger"
    change_icon <- if (change >= 0) "↗" else "↘"
    
    tags$div(
      class = paste("stats-change", change_class),
      style = "font-size: 0.875rem; margin-top: 0.5rem;",
      change_icon, " ", abs(change), "%"
    )
  }
  
  modern_card(
    title = NULL,
    body = tags$div(
      class = "stats-card-content",
      style = "text-align: center; padding: 1rem 0;",
      
      if (!is.null(icon)) {
        tags$div(
          class = paste0("stats-icon text-", color),
          style = "font-size: 2rem; margin-bottom: 0.5rem;",
          icon
        )
      },
      
      tags$div(
        class = "stats-value",
        style = paste0("font-size: 2.5rem; font-weight: 700; color: var(--", color, "-600); margin-bottom: 0.25rem;"),
        value
      ),
      
      tags$div(
        class = "stats-title",
        style = "font-size: 1rem; font-weight: 600; color: var(--gray-700); margin-bottom: 0.25rem;",
        title
      ),
      
      if (!is.null(subtitle)) {
        tags$div(
          class = "stats-subtitle",
          style = "font-size: 0.875rem; color: var(--gray-500);",
          subtitle
        )
      },
      
      change_element
    ),
    width = width,
    class = "card-modern stats-card"
  )
}

# ===========================================
# MODERN FORM COMPONENTS
# ===========================================

# Enhanced form input with modern styling
modern_input <- function(inputId, label, value = "", placeholder = NULL,
                        type = "text", help_text = NULL, required = FALSE,
                        validation = NULL, width = NULL) {
  
  input_class <- "form-input-modern"
  if (!is.null(validation)) {
    input_class <- paste(input_class, validation)
  }
  
  container_style <- if (!is.null(width)) {
    paste0("width: ", width, ";")
  } else ""
  
  tags$div(
    class = "form-group-modern",
    style = container_style,
    
    tags$label(
      class = "form-label-modern",
      `for` = inputId,
      label,
      if (required) tags$span(class = "text-danger", " *")
    ),
    
    tags$input(
      id = inputId,
      name = inputId,
      type = type,
      class = input_class,
      value = value,
      placeholder = placeholder,
      required = if (required) NA else NULL
    ),
    
    if (!is.null(help_text)) {
      tags$div(
        class = "form-help-text",
        style = "font-size: 0.75rem; color: var(--gray-500); margin-top: 0.25rem;",
        help_text
      )
    }
  )
}

# Enhanced select input with modern styling
modern_select <- function(inputId, label, choices, selected = NULL,
                         multiple = FALSE, help_text = NULL, required = FALSE,
                         width = NULL) {
  
  container_style <- if (!is.null(width)) {
    paste0("width: ", width, ";")
  } else ""
  
  tags$div(
    class = "form-group-modern",
    style = container_style,
    
    tags$label(
      class = "form-label-modern",
      `for` = inputId,
      label,
      if (required) tags$span(class = "text-danger", " *")
    ),
    
    tags$select(
      id = inputId,
      name = inputId,
      class = "form-select-modern",
      multiple = if (multiple) NA else NULL,
      required = if (required) NA else NULL,
      
      lapply(names(choices), function(name) {
        value <- choices[[name]]
        tags$option(
          value = value,
          selected = if (!is.null(selected) && value %in% selected) NA else NULL,
          name
        )
      })
    ),
    
    if (!is.null(help_text)) {
      tags$div(
        class = "form-help-text",
        style = "font-size: 0.75rem; color: var(--gray-500); margin-top: 0.25rem;",
        help_text
      )
    }
  )
}

# ===========================================
# PROGRESS AND LOADING COMPONENTS
# ===========================================

# Modern progress bar with label and animation
modern_progress <- function(elementId, label = "Progress", percentage = 0,
                           message = "", striped = TRUE, animated = TRUE) {
  
  bar_class <- "progress-modern-bar"
  if (striped) bar_class <- paste(bar_class, "progress-striped")
  if (animated) bar_class <- paste(bar_class, "progress-animated")
  
  tags$div(
    id = elementId,
    class = "progress-modern-container",
    
    tags$div(
      class = "progress-modern-label",
      tags$span(label),
      tags$span(paste0(round(percentage), "%"))
    ),
    
    tags$div(
      class = "progress-modern",
      tags$div(
        class = bar_class,
        style = paste0("width: ", percentage, "%;"),
        role = "progressbar",
        `aria-valuenow` = percentage,
        `aria-valuemin` = "0",
        `aria-valuemax` = "100"
      )
    ),
    
    if (message != "") {
      tags$div(
        class = "progress-modern-message",
        style = "font-size: 0.875rem; color: var(--gray-600); margin-top: 0.5rem;",
        message
      )
    }
  )
}

# Loading overlay for containers
loading_overlay <- function(message = "Loading...") {
  tags$div(
    class = "loading-modern",
    style = "position: absolute; top: 0; left: 0; right: 0; bottom: 0; 
             background: rgba(255, 255, 255, 0.9); z-index: 1000;
             display: flex; flex-direction: column; align-items: center; 
             justify-content: center;",
    
    tags$div(class = "loading-modern-spinner"),
    tags$div(class = "loading-modern-text", message)
  )
}

# ===========================================
# ALERT AND NOTIFICATION COMPONENTS
# ===========================================

# Modern alert with icon and actions
modern_alert <- function(message, type = "info", title = NULL, 
                        dismissible = TRUE, actions = NULL) {
  
  icons <- list(
    success = "✅",
    warning = "⚠️",
    error = "❌",
    info = "ℹ️"
  )
  
  tags$div(
    class = paste0("alert-modern alert-modern-", type),
    role = "alert",
    
    tags$div(
      class = "alert-modern-icon",
      icons[[type]]
    ),
    
    tags$div(
      class = "alert-modern-content",
      
      if (!is.null(title)) {
        tags$div(class = "alert-modern-title", title)
      },
      
      tags$div(class = "alert-modern-message", message)
    ),
    
    if (!is.null(actions)) {
      tags$div(
        class = "alert-modern-actions",
        style = "margin-left: auto; display: flex; gap: 0.5rem;",
        actions
      )
    },
    
    if (dismissible) {
      tags$button(
        type = "button",
        class = "alert-close",
        style = "background: none; border: none; font-size: 1.25rem; 
                 cursor: pointer; color: inherit; opacity: 0.7; margin-left: auto;",
        onclick = "this.parentElement.remove()",
        "×"
      )
    }
  )
}

# Notification system (client-side)
send_notification <- function(session, message, type = "info", duration = 5000) {
  session$sendCustomMessage(
    type = "showNotification",
    message = list(
      message = message,
      type = type,
      duration = duration
    )
  )
}

# ===========================================
# ENHANCED DATA TABLE COMPONENTS
# ===========================================

# Modern data table with enhanced features
modern_datatable <- function(data, options = list(), selection = "none",
                           extensions = c("Buttons", "ColReorder", "Responsive"),
                           dom = "Bfrtip", buttons = c("copy", "csv", "excel"),
                           pageLength = 25) {
  
  # Default options with modern styling
  default_options <- list(
    pageLength = pageLength,
    dom = dom,
    buttons = list(
      list(extend = "copy", className = "btn-modern btn-modern-sm"),
      list(extend = "csv", className = "btn-modern btn-modern-sm"),
      list(extend = "excel", className = "btn-modern btn-modern-sm")
    ),
    columnDefs = list(
      list(className = "dt-center", targets = "_all")
    ),
    responsive = TRUE,
    colReorder = TRUE,
    scrollX = TRUE,
    autoWidth = TRUE,
    language = list(
      search = "Search:",
      lengthMenu = "Show _MENU_ entries",
      info = "Showing _START_ to _END_ of _TOTAL_ entries",
      paginate = list(
        first = "First",
        last = "Last",
        `next` = "Next",
        previous = "Previous"
      )
    )
  )
  
  # Merge with user options
  final_options <- modifyList(default_options, options)
  
  DT::datatable(
    data,
    options = final_options,
    selection = selection,
    extensions = extensions,
    escape = FALSE,
    rownames = FALSE,
    class = "table-modern display nowrap"
  )
}

# ===========================================
# AUTHENTICATION COMPONENTS
# ===========================================

# Login form with modern styling
auth_login_form <- function(inputId_prefix = "auth") {
  tags$form(
    id = "signin-form",
    class = "auth-form-modern",
    
    tags$div(
      class = "auth-header",
      style = "text-align: center; margin-bottom: 2rem;",
      tags$h2("Sign In", style = "color: var(--gray-900); margin-bottom: 0.5rem;"),
      tags$p("Access your Prairie Genomics projects", style = "color: var(--gray-600);")
    ),
    
    modern_input(
      paste0(inputId_prefix, "_email"),
      "Email Address",
      type = "email",
      placeholder = "Enter your email",
      required = TRUE
    ),
    
    modern_input(
      paste0(inputId_prefix, "_password"),
      "Password",
      type = "password",
      placeholder = "Enter your password",
      required = TRUE
    ),
    
    tags$div(
      class = "form-actions",
      style = "margin-top: 1.5rem;",
      
      tags$button(
        type = "submit",
        class = "btn-modern btn-modern-primary btn-modern-lg",
        style = "width: 100%; margin-bottom: 1rem;",
        "Sign In"
      ),
      
      tags$div(
        style = "text-align: center;",
        tags$a(
          href = "#signup",
          style = "color: var(--primary-600); text-decoration: none;",
          "Don't have an account? Sign up"
        )
      )
    )
  )
}

# Registration form with modern styling
auth_register_form <- function(inputId_prefix = "auth") {
  tags$form(
    id = "signup-form",
    class = "auth-form-modern",
    
    tags$div(
      class = "auth-header",
      style = "text-align: center; margin-bottom: 2rem;",
      tags$h2("Create Account", style = "color: var(--gray-900); margin-bottom: 0.5rem;"),
      tags$p("Start your genomics analysis journey", style = "color: var(--gray-600);")
    ),
    
    modern_input(
      paste0(inputId_prefix, "_display_name"),
      "Full Name",
      placeholder = "Enter your full name",
      required = TRUE
    ),
    
    modern_input(
      paste0(inputId_prefix, "_email_signup"),
      "Email Address",
      type = "email",
      placeholder = "Enter your email",
      required = TRUE
    ),
    
    modern_input(
      paste0(inputId_prefix, "_password_signup"),
      "Password",
      type = "password",
      placeholder = "Choose a password (min. 6 characters)",
      required = TRUE,
      help_text = "Must be at least 6 characters long"
    ),
    
    modern_input(
      paste0(inputId_prefix, "_confirm_password"),
      "Confirm Password",
      type = "password",
      placeholder = "Confirm your password",
      required = TRUE
    ),
    
    tags$div(
      class = "form-actions",
      style = "margin-top: 1.5rem;",
      
      tags$button(
        type = "submit",
        class = "btn-modern btn-modern-success btn-modern-lg",
        style = "width: 100%; margin-bottom: 1rem;",
        "Create Account"
      ),
      
      tags$div(
        style = "text-align: center;",
        tags$a(
          href = "#signin",
          style = "color: var(--primary-600); text-decoration: none;",
          "Already have an account? Sign in"
        )
      )
    )
  )
}

# ===========================================
# STEP INDICATOR COMPONENT
# ===========================================

# Multi-step process indicator
step_indicator <- function(steps, current_step = 1, completed_steps = NULL) {
  
  tags$div(
    class = "step-indicator-container",
    style = "margin: 2rem 0;",
    
    lapply(seq_along(steps), function(i) {
      step_class <- "step-indicator"
      
      if (!is.null(completed_steps) && i %in% completed_steps) {
        step_class <- paste(step_class, "step-completed")
      } else if (i == current_step) {
        step_class <- paste(step_class, "step-active")
      } else {
        step_class <- paste(step_class, "step-pending")
      }
      
      tags$div(
        class = step_class,
        
        tags$div(
          class = "step-number",
          style = "width: 2rem; height: 2rem; border-radius: 50%; 
                   display: flex; align-items: center; justify-content: center;
                   font-weight: 600; margin-right: 1rem;",
          i
        ),
        
        tags$div(
          class = "step-content",
          tags$div(class = "step-title", style = "font-weight: 600;", steps[[i]]$title),
          if (!is.null(steps[[i]]$description)) {
            tags$div(
              class = "step-description",
              style = "font-size: 0.875rem; opacity: 0.8;",
              steps[[i]]$description
            )
          }
        )
      )
    })
  )
}

# ===========================================
# EXPORT COMPONENTS
# ===========================================

# Export all components
list(
  modern_card = modern_card,
  stats_card = stats_card,
  modern_input = modern_input,
  modern_select = modern_select,
  modern_progress = modern_progress,
  loading_overlay = loading_overlay,
  modern_alert = modern_alert,
  send_notification = send_notification,
  modern_datatable = modern_datatable,
  auth_login_form = auth_login_form,
  auth_register_form = auth_register_form,
  step_indicator = step_indicator
)