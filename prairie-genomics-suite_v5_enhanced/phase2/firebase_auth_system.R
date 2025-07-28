# Phase 2 - Firebase Authentication System
# Complete user authentication with registration, login, and session management
# 
# This module provides:
# - User registration and login
# - Secure session management  
# - Password reset functionality
# - User profile management
# - Access control and permissions

library(httr)
library(jsonlite)
library(shiny)
library(digest)

# ===========================================
# FIREBASE CONFIGURATION
# ===========================================

# Firebase configuration (replace with your project details)
FIREBASE_CONFIG <- list(
  apiKey = "your-api-key-here",
  authDomain = "prairie-genomics.firebaseapp.com",
  projectId = "prairie-genomics",
  storageBucket = "prairie-genomics.appspot.com",
  messagingSenderId = "123456789",
  appId = "1:123456789:web:abcdef123456",
  databaseURL = "https://prairie-genomics-default-rtdb.firebaseio.com/"
)

# Firebase Auth REST API endpoints
FIREBASE_AUTH_ENDPOINTS <- list(
  signup = paste0("https://identitytoolkit.googleapis.com/v1/accounts:signUp?key=", FIREBASE_CONFIG$apiKey),
  signin = paste0("https://identitytoolkit.googleapis.com/v1/accounts:signInWithPassword?key=", FIREBASE_CONFIG$apiKey),
  refresh = paste0("https://securetoken.googleapis.com/v1/token?key=", FIREBASE_CONFIG$apiKey),
  reset_password = paste0("https://identitytoolkit.googleapis.com/v1/accounts:sendOobCode?key=", FIREBASE_CONFIG$apiKey),
  verify_email = paste0("https://identitytoolkit.googleapis.com/v1/accounts:sendOobCode?key=", FIREBASE_CONFIG$apiKey),
  user_info = paste0("https://identitytoolkit.googleapis.com/v1/accounts:lookup?key=", FIREBASE_CONFIG$apiKey)
)

# ===========================================
# FIREBASE AUTHENTICATION CORE FUNCTIONS
# ===========================================

# Firebase user registration
firebase_register_user <- function(email, password, display_name = NULL) {
  
  # Validate input
  if (is.null(email) || email == "" || !grepl("@", email)) {
    return(list(success = FALSE, error = "Invalid email address"))
  }
  
  if (is.null(password) || nchar(password) < 6) {
    return(list(success = FALSE, error = "Password must be at least 6 characters"))
  }
  
  # Prepare request body
  request_body <- list(
    email = email,
    password = password,
    returnSecureToken = TRUE
  )
  
  # Add display name if provided
  if (!is.null(display_name) && display_name != "") {
    request_body$displayName <- display_name
  }
  
  tryCatch({
    # Make registration request
    response <- POST(
      url = FIREBASE_AUTH_ENDPOINTS$signup,
      body = toJSON(request_body, auto_unbox = TRUE),
      add_headers("Content-Type" = "application/json"),
      timeout(30)
    )
    
    # Parse response  
    response_content <- content(response, "text", encoding = "UTF-8")
    response_data <- fromJSON(response_content)
    
    if (response$status_code == 200) {
      # Successful registration
      user_data <- list(
        success = TRUE,
        user_id = response_data$localId,
        email = response_data$email,
        display_name = response_data$displayName %||% email,
        id_token = response_data$idToken,
        refresh_token = response_data$refreshToken,
        expires_in = as.numeric(response_data$expiresIn),
        registered_at = Sys.time(),
        email_verified = FALSE
      )
      
      return(user_data)
      
    } else {
      # Registration failed
      error_message <- if (!is.null(response_data$error$message)) {
        response_data$error$message
      } else {
        "Registration failed"
      }
      
      return(list(success = FALSE, error = error_message))
    }
    
  }, error = function(e) {
    return(list(success = FALSE, error = paste("Network error:", e$message)))
  })
}

# Firebase user login
firebase_login_user <- function(email, password) {
  
  # Validate input
  if (is.null(email) || email == "") {
    return(list(success = FALSE, error = "Email is required"))
  }
  
  if (is.null(password) || password == "") {
    return(list(success = FALSE, error = "Password is required"))
  }
  
  # Prepare request body
  request_body <- list(
    email = email,
    password = password,
    returnSecureToken = TRUE
  )
  
  tryCatch({
    # Make login request
    response <- POST(
      url = FIREBASE_AUTH_ENDPOINTS$signin,
      body = toJSON(request_body, auto_unbox = TRUE),
      add_headers("Content-Type" = "application/json"),
      timeout(30)
    )
    
    # Parse response
    response_content <- content(response, "text", encoding = "UTF-8")
    response_data <- fromJSON(response_content)
    
    if (response$status_code == 200) {
      # Successful login
      user_data <- list(
        success = TRUE,
        user_id = response_data$localId,
        email = response_data$email,
        display_name = response_data$displayName %||% email,
        id_token = response_data$idToken,
        refresh_token = response_data$refreshToken,
        expires_in = as.numeric(response_data$expiresIn),
        logged_in_at = Sys.time(),
        email_verified = response_data$emailVerified %||% FALSE
      )
      
      return(user_data)
      
    } else {
      # Login failed
      error_message <- if (!is.null(response_data$error$message)) {
        switch(response_data$error$message,
               "INVALID_EMAIL" = "Invalid email address",
               "EMAIL_NOT_FOUND" = "No account found with this email",
               "INVALID_PASSWORD" = "Incorrect password",
               "USER_DISABLED" = "This account has been disabled",
               "TOO_MANY_ATTEMPTS_TRY_LATER" = "Too many failed attempts. Try again later.",
               response_data$error$message
        )
      } else {
        "Login failed"
      }
      
      return(list(success = FALSE, error = error_message))
    }
    
  }, error = function(e) {
    return(list(success = FALSE, error = paste("Network error:", e$message)))
  })
}

# Refresh Firebase ID token
firebase_refresh_token <- function(refresh_token) {
  
  if (is.null(refresh_token) || refresh_token == "") {
    return(list(success = FALSE, error = "Refresh token is required"))
  }
  
  request_body <- list(
    grant_type = "refresh_token",
    refresh_token = refresh_token
  )
  
  tryCatch({
    response <- POST(
      url = FIREBASE_AUTH_ENDPOINTS$refresh,
      body = toJSON(request_body, auto_unbox = TRUE),
      add_headers("Content-Type" = "application/json"),
      timeout(30)
    )
    
    response_content <- content(response, "text", encoding = "UTF-8")
    response_data <- fromJSON(response_content)
    
    if (response$status_code == 200) {
      return(list(
        success = TRUE,
        id_token = response_data$id_token,
        refresh_token = response_data$refresh_token,
        expires_in = as.numeric(response_data$expires_in),
        user_id = response_data$user_id
      ))
    } else {
      return(list(success = FALSE, error = "Token refresh failed"))
    }
    
  }, error = function(e) {
    return(list(success = FALSE, error = paste("Token refresh error:", e$message)))
  })
}

# Send password reset email
firebase_reset_password <- function(email) {
  
  if (is.null(email) || email == "") {
    return(list(success = FALSE, error = "Email is required"))
  }
  
  request_body <- list(
    requestType = "PASSWORD_RESET",
    email = email
  )
  
  tryCatch({
    response <- POST(
      url = FIREBASE_AUTH_ENDPOINTS$reset_password,
      body = toJSON(request_body, auto_unbox = TRUE),
      add_headers("Content-Type" = "application/json"),
      timeout(30)
    )
    
    if (response$status_code == 200) {
      return(list(success = TRUE, message = "Password reset email sent"))
    } else {
      return(list(success = FALSE, error = "Failed to send reset email"))
    }
    
  }, error = function(e) {
    return(list(success = FALSE, error = paste("Reset password error:", e$message)))
  })
}

# ===========================================  
# USER SESSION MANAGEMENT
# ===========================================

# Create user session manager
create_user_session_manager <- function() {
  
  # Session state
  session_data <- reactiveValues(
    is_authenticated = FALSE,
    user = NULL,
    token_expires_at = NULL,
    last_activity = NULL
  )
  
  # Session methods
  session_manager <- list(
    
    # Login user and create session
    login = function(email, password) {
      result <- firebase_login_user(email, password)
      
      if (result$success) {
        # Store user session
        session_data$is_authenticated <- TRUE
        session_data$user <- result
        session_data$token_expires_at <- Sys.time() + result$expires_in
        session_data$last_activity <- Sys.time()
        
        # Store in browser session storage (via JavaScript)
        session$sendCustomMessage(
          type = "storeUserSession",
          message = list(
            user_id = result$user_id,
            email = result$email,
            display_name = result$display_name,
            id_token = result$id_token,
            expires_at = as.numeric(session_data$token_expires_at)
          )
        )
        
        showNotification(
          paste("Welcome back,", result$display_name %||% result$email, "!"),
          type = "success",
          duration = 3
        )
      }
      
      return(result)
    },
    
    # Register new user
    register = function(email, password, display_name = NULL) {
      result <- firebase_register_user(email, password, display_name)
      
      if (result$success) {
        # Automatically log in after registration
        session_data$is_authenticated <- TRUE
        session_data$user <- result
        session_data$token_expires_at <- Sys.time() + result$expires_in
        session_data$last_activity <- Sys.time()
        
        showNotification(
          paste("Welcome to Prairie Genomics Suite,", result$display_name %||% result$email, "!"),
          type = "success",
          duration = 5
        )
      }
      
      return(result)
    },
    
    # Logout user
    logout = function() {
      session_data$is_authenticated <- FALSE
      session_data$user <- NULL
      session_data$token_expires_at <- NULL
      session_data$last_activity <- NULL
      
      # Clear browser session storage
      session$sendCustomMessage(
        type = "clearUserSession",
        message = list()
      )
      
      showNotification("You have been logged out", type = "message", duration = 3)
    },
    
    # Check if user is authenticated
    is_authenticated = function() {
      if (!session_data$is_authenticated) return(FALSE)
      
      # Check token expiration
      if (!is.null(session_data$token_expires_at) && 
          Sys.time() > session_data$token_expires_at) {
        # Token expired, try to refresh
        if (!is.null(session_data$user$refresh_token)) {
          refresh_result <- firebase_refresh_token(session_data$user$refresh_token)
          
          if (refresh_result$success) {
            # Update session with new token
            session_data$user$id_token <- refresh_result$id_token
            session_data$token_expires_at <- Sys.time() + refresh_result$expires_in
            return(TRUE)
          }
        }
        
        # Refresh failed, logout
        session_manager$logout()
        return(FALSE)
      }
      
      return(TRUE)
    },
    
    # Get current user
    get_user = function() {
      if (session_manager$is_authenticated()) {
        return(session_data$user)
      }
      return(NULL)
    },
    
    # Update last activity
    update_activity = function() {
      if (session_data$is_authenticated) {
        session_data$last_activity <- Sys.time()
      }
    },
    
    # Get session data (reactive)
    get_session_data = function() {
      return(session_data)
    }
  )
  
  return(session_manager)
}

# ===========================================
# AUTHENTICATION UI COMPONENTS
# ===========================================

# Create authentication UI
create_auth_ui <- function(session_manager) {
  
  # Get session data reactively
  session_data <- session_manager$get_session_data()
  
  # Authentication status UI
  auth_status_ui <- reactive({
    if (session_data$is_authenticated) {
      user <- session_data$user
      div(
        class = "user-info-panel",
        style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
        
        div(
          style = "display: flex; align-items: center; justify-content: space-between;",
          
          div(
            h4(style = "margin: 0; color: #28a745;", "✅ Signed In"),
            p(style = "margin: 5px 0 0 0; color: #666;",
              paste("Welcome,", user$display_name %||% user$email))
          ),
          
          actionButton(
            "logout_btn",
            "Sign Out",
            class = "btn btn-outline-secondary btn-sm",
            style = "margin-left: 15px;"
          )
        )
      )
    } else {
      div(
        class = "auth-panel",
        
        # Tab navigation
        div(
          class = "auth-tabs",
          style = "margin-bottom: 20px;",
          
          div(
            class = "btn-group btn-group-justified",
            role = "group",
            
            actionButton(
              "show_login_tab",
              "Sign In",
              class = "btn btn-primary active",
              style = "border-radius: 4px 0 0 4px;"
            ),
            
            actionButton(
              "show_register_tab", 
              "Create Account",
              class = "btn btn-outline-primary",
              style = "border-radius: 0 4px 4px 0;"
            )
          )
        ),
        
        # Login form
        conditionalPanel(
          condition = "input.show_login_tab % 2 == 1 || input.show_login_tab == 0",
          
          div(
            id = "login_form",
            
            h3("Sign In", style = "color: #333; margin-bottom: 20px;"),
            
            div(
              class = "form-group",
              textInput(
                "login_email",
                "Email Address",
                placeholder = "Enter your email",
                width = "100%"
              )
            ),
            
            div(
              class = "form-group",
              passwordInput(
                "login_password",
                "Password", 
                placeholder = "Enter your password",
                width = "100%"
              )
            ),
            
            div(
              class = "form-actions",
              style = "margin-top: 20px;",
              
              actionButton(
                "login_submit",
                "Sign In",
                class = "btn btn-primary btn-lg",
                style = "width: 100%; margin-bottom: 10px;"
              ),
              
              div(
                style = "text-align: center;",
                actionLink(
                  "forgot_password_link",
                  "Forgot your password?",
                  style = "color: #007bff; font-size: 14px;"
                )
              )
            )
          )
        ),
        
        # Registration form
        conditionalPanel(
          condition = "input.show_register_tab % 2 == 1",
          
          div(
            id = "register_form",
            
            h3("Create Account", style = "color: #333; margin-bottom: 20px;"),
            
            div(
              class = "form-group",
              textInput(
                "register_name",
                "Full Name",
                placeholder = "Enter your full name",
                width = "100%"
              )
            ),
            
            div(
              class = "form-group", 
              textInput(
                "register_email",
                "Email Address",
                placeholder = "Enter your email",
                width = "100%"
              )
            ),
            
            div(
              class = "form-group",
              passwordInput(
                "register_password",
                "Password",
                placeholder = "Choose a password (min. 6 characters)",
                width = "100%"
              )
            ),
            
            div(
              class = "form-group",
              passwordInput(
                "register_confirm_password",
                "Confirm Password",
                placeholder = "Confirm your password", 
                width = "100%"
              )
            ),
            
            div(
              class = "form-actions",
              style = "margin-top: 20px;",
              
              actionButton(
                "register_submit",
                "Create Account",
                class = "btn btn-success btn-lg",
                style = "width: 100%;"
              )
            )
          )
        )
      )
    }
  })
  
  return(auth_status_ui)
}

# ===========================================
# EXPORT FUNCTIONS
# ===========================================

list(
  firebase_register_user = firebase_register_user,
  firebase_login_user = firebase_login_user,
  firebase_refresh_token = firebase_refresh_token,
  firebase_reset_password = firebase_reset_password,
  create_user_session_manager = create_user_session_manager,
  create_auth_ui = create_auth_ui,
  FIREBASE_CONFIG = FIREBASE_CONFIG
)