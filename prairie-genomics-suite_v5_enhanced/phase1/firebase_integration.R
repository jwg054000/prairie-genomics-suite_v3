# Phase 1 - Firebase Integration for Prairie Genomics Suite
# Authentication, project management, and real-time features
# 
# This module provides:
# - Firebase Authentication (Google, email/password)
# - Firestore database integration for project storage
# - Real-time collaboration features
# - User session management

library(httr)
library(jsonlite)
library(digest)

# Firebase Configuration
FIREBASE_CONFIG <- list(
  apiKey = Sys.getenv("FIREBASE_API_KEY", ""),
  authDomain = Sys.getenv("FIREBASE_AUTH_DOMAIN", "prairie-genomics.firebaseapp.com"),
  projectId = Sys.getenv("FIREBASE_PROJECT_ID", "prairie-genomics"),
  storageBucket = Sys.getenv("FIREBASE_STORAGE_BUCKET", "prairie-genomics.appspot.com"),
  messagingSenderId = Sys.getenv("FIREBASE_MESSAGING_SENDER_ID", ""),
  appId = Sys.getenv("FIREBASE_APP_ID", "")
)

# Firebase Authentication Helper
firebase_auth <- function() {
  
  # Initialize authentication state
  auth_state <- reactiveValues(
    user = NULL,
    token = NULL,
    authenticated = FALSE,
    loading = FALSE,
    error = NULL
  )
  
  # Sign in with email and password
  sign_in_email_password <- function(email, password) {
    auth_state$loading <- TRUE
    auth_state$error <- NULL
    
    tryCatch({
      # Firebase REST API for authentication
      auth_url <- paste0(
        "https://identitytoolkit.googleapis.com/v1/accounts:signInWithPassword?key=",
        FIREBASE_CONFIG$apiKey
      )
      
      response <- POST(
        auth_url,
        body = list(
          email = email,
          password = password,
          returnSecureToken = TRUE
        ),
        encode = "json",
        add_headers("Content-Type" = "application/json")
      )
      
      if (status_code(response) == 200) {
        auth_data <- content(response, "parsed")
        
        # Store user data
        auth_state$user <- list(
          uid = auth_data$localId,
          email = auth_data$email,
          displayName = auth_data$displayName %||% "User",
          emailVerified = as.logical(auth_data$emailVerified %||% FALSE)
        )
        
        auth_state$token <- auth_data$idToken
        auth_state$authenticated <- TRUE
        auth_state$loading <- FALSE
        
        # Store token in session storage (client-side)
        session$sendCustomMessage(
          type = "store_auth_token",
          message = list(token = auth_data$idToken, user = auth_state$user)
        )
        
        return(list(success = TRUE, user = auth_state$user))
        
      } else {
        error_data <- content(response, "parsed")
        auth_state$error <- error_data$error$message %||% "Authentication failed"
        auth_state$loading <- FALSE
        
        return(list(success = FALSE, error = auth_state$error))
      }
      
    }, error = function(e) {
      auth_state$error <- paste("Authentication error:", e$message)
      auth_state$loading <- FALSE
      return(list(success = FALSE, error = auth_state$error))
    })
  }
  
  # Sign up with email and password
  sign_up_email_password <- function(email, password, display_name = NULL) {
    auth_state$loading <- TRUE
    auth_state$error <- NULL
    
    tryCatch({
      # Firebase REST API for registration
      signup_url <- paste0(
        "https://identitytoolkit.googleapis.com/v1/accounts:signUp?key=",
        FIREBASE_CONFIG$apiKey
      )
      
      response <- POST(
        signup_url,
        body = list(
          email = email,
          password = password,
          returnSecureToken = TRUE
        ),
        encode = "json",
        add_headers("Content-Type" = "application/json")
      )
      
      if (status_code(response) == 200) {
        auth_data <- content(response, "parsed")
        
        # Update profile with display name if provided
        if (!is.null(display_name)) {
          update_profile(auth_data$idToken, display_name)
        }
        
        # Store user data
        auth_state$user <- list(
          uid = auth_data$localId,
          email = auth_data$email,
          displayName = display_name %||% "User",
          emailVerified = FALSE
        )
        
        auth_state$token <- auth_data$idToken
        auth_state$authenticated <- TRUE
        auth_state$loading <- FALSE
        
        return(list(success = TRUE, user = auth_state$user))
        
      } else {
        error_data <- content(response, "parsed")
        auth_state$error <- error_data$error$message %||% "Registration failed"
        auth_state$loading <- FALSE
        
        return(list(success = FALSE, error = auth_state$error))
      }
      
    }, error = function(e) {
      auth_state$error <- paste("Registration error:", e$message)
      auth_state$loading <- FALSE
      return(list(success = FALSE, error = auth_state$error))
    })
  }
  
  # Update user profile
  update_profile <- function(token, display_name) {
    tryCatch({
      profile_url <- paste0(
        "https://identitytoolkit.googleapis.com/v1/accounts:update?key=",
        FIREBASE_CONFIG$apiKey
      )
      
      POST(
        profile_url,
        body = list(
          idToken = token,
          displayName = display_name,
          returnSecureToken = FALSE
        ),
        encode = "json",
        add_headers("Content-Type" = "application/json")
      )
    }, error = function(e) {
      cat("Profile update failed:", e$message, "\n")
    })
  }
  
  # Sign out
  sign_out <- function() {
    auth_state$user <- NULL
    auth_state$token <- NULL
    auth_state$authenticated <- FALSE
    auth_state$error <- NULL
    
    # Clear client-side storage
    session$sendCustomMessage(type = "clear_auth_token", message = list())
    
    return(list(success = TRUE))
  }
  
  # Verify token and refresh if needed
  verify_token <- function(token) {
    tryCatch({
      verify_url <- paste0(
        "https://identitytoolkit.googleapis.com/v1/accounts:lookup?key=",
        FIREBASE_CONFIG$apiKey
      )
      
      response <- POST(
        verify_url,
        body = list(idToken = token),
        encode = "json",
        add_headers("Content-Type" = "application/json")
      )
      
      if (status_code(response) == 200) {
        user_data <- content(response, "parsed")$users[[1]]
        
        auth_state$user <- list(
          uid = user_data$localId,
          email = user_data$email,
          displayName = user_data$displayName %||% "User",
          emailVerified = as.logical(user_data$emailVerified %||% FALSE)
        )
        
        auth_state$token <- token
        auth_state$authenticated <- TRUE
        
        return(list(success = TRUE, user = auth_state$user))
      } else {
        return(list(success = FALSE, error = "Invalid token"))
      }
      
    }, error = function(e) {
      return(list(success = FALSE, error = e$message))
    })
  }
  
  # Return authentication interface
  list(
    state = auth_state,
    sign_in = sign_in_email_password,
    sign_up = sign_up_email_password,
    sign_out = sign_out,
    verify_token = verify_token
  )
}

# Firestore Database Integration
firestore_db <- function() {
  
  # Get Firestore document
  get_document <- function(collection, document_id, token) {
    tryCatch({
      doc_url <- paste0(
        "https://firestore.googleapis.com/v1/projects/",
        FIREBASE_CONFIG$projectId,
        "/databases/(default)/documents/",
        collection, "/", document_id
      )
      
      response <- GET(
        doc_url,
        add_headers(
          "Authorization" = paste("Bearer", token),
          "Content-Type" = "application/json"
        )
      )
      
      if (status_code(response) == 200) {
        doc_data <- content(response, "parsed")
        return(list(success = TRUE, data = parse_firestore_fields(doc_data$fields)))
      } else {
        return(list(success = FALSE, error = "Document not found"))
      }
      
    }, error = function(e) {
      return(list(success = FALSE, error = e$message))
    })
  }
  
  # Set Firestore document
  set_document <- function(collection, document_id, data, token) {
    tryCatch({
      doc_url <- paste0(
        "https://firestore.googleapis.com/v1/projects/",
        FIREBASE_CONFIG$projectId,
        "/databases/(default)/documents/",
        collection, "/", document_id
      )
      
      # Convert R data to Firestore format
      firestore_data <- list(fields = format_firestore_fields(data))
      
      response <- PATCH(
        doc_url,
        body = firestore_data,
        encode = "json",
        add_headers(
          "Authorization" = paste("Bearer", token),
          "Content-Type" = "application/json"
        )
      )
      
      if (status_code(response) %in% c(200, 201)) {
        return(list(success = TRUE, id = document_id))
      } else {
        error_data <- content(response, "parsed")
        return(list(success = FALSE, error = error_data$error$message))
      }
      
    }, error = function(e) {
      return(list(success = FALSE, error = e$message))
    })
  }
  
  # Query collection
  query_collection <- function(collection, filters = NULL, token, limit = 50) {
    tryCatch({
      query_url <- paste0(
        "https://firestore.googleapis.com/v1/projects/",
        FIREBASE_CONFIG$projectId,
        "/databases/(default)/documents:runQuery"
      )
      
      # Build structured query
      query_body <- list(
        structuredQuery = list(
          from = list(list(collectionId = collection)),
          limit = list(value = limit)
        )
      )
      
      # Add filters if provided
      if (!is.null(filters)) {
        query_body$structuredQuery$where <- build_firestore_filters(filters)
      }
      
      response <- POST(
        query_url,
        body = query_body,
        encode = "json",
        add_headers(
          "Authorization" = paste("Bearer", token),
          "Content-Type" = "application/json"
        )
      )
      
      if (status_code(response) == 200) {
        query_results <- content(response, "parsed")
        
        # Parse results
        documents <- lapply(query_results, function(result) {
          if (!is.null(result$document)) {
            return(list(
              id = basename(result$document$name),
              data = parse_firestore_fields(result$document$fields)
            ))
          }
          return(NULL)
        })
        
        # Remove NULL results
        documents <- documents[!sapply(documents, is.null)]
        
        return(list(success = TRUE, documents = documents))
        
      } else {
        return(list(success = FALSE, error = "Query failed"))
      }
      
    }, error = function(e) {
      return(list(success = FALSE, error = e$message))
    })
  }
  
  # Helper function to format R data for Firestore
  format_firestore_fields <- function(data) {
    fields <- list()
    
    for (name in names(data)) {
      value <- data[[name]]
      
      if (is.character(value)) {
        fields[[name]] <- list(stringValue = value)
      } else if (is.numeric(value)) {
        if (is.integer(value)) {
          fields[[name]] <- list(integerValue = as.character(value))
        } else {
          fields[[name]] <- list(doubleValue = value)
        }
      } else if (is.logical(value)) {
        fields[[name]] <- list(booleanValue = value)
      } else if (is.list(value)) {
        fields[[name]] <- list(mapValue = list(fields = format_firestore_fields(value)))
      } else {
        # Convert to string as fallback
        fields[[name]] <- list(stringValue = as.character(value))
      }
    }
    
    return(fields)
  }
  
  # Helper function to parse Firestore fields to R data
  parse_firestore_fields <- function(fields) {
    if (is.null(fields)) return(list())
    
    data <- list()
    
    for (name in names(fields)) {
      field <- fields[[name]]
      
      if (!is.null(field$stringValue)) {
        data[[name]] <- field$stringValue
      } else if (!is.null(field$integerValue)) {
        data[[name]] <- as.integer(field$integerValue)
      } else if (!is.null(field$doubleValue)) {
        data[[name]] <- field$doubleValue
      } else if (!is.null(field$booleanValue)) {
        data[[name]] <- field$booleanValue
      } else if (!is.null(field$mapValue)) {
        data[[name]] <- parse_firestore_fields(field$mapValue$fields)
      } else {
        data[[name]] <- NULL
      }
    }
    
    return(data)
  }
  
  # Return database interface
  list(
    get = get_document,
    set = set_document,
    query = query_collection
  )
}

# Project Management Integration
project_manager <- function(auth, db) {
  
  # Create new project
  create_project <- function(project_name, description = "", user_uid, token) {
    project_id <- paste0("project_", digest(paste(project_name, user_uid, Sys.time()), "md5"))
    
    project_data <- list(
      name = project_name,
      description = description,
      owner_uid = user_uid,
      created_at = as.character(Sys.time()),
      updated_at = as.character(Sys.time()),
      collaborators = list(user_uid),
      datasets = list(),
      analyses = list(),
      status = "active"
    )
    
    result <- db$set("projects", project_id, project_data, token)
    
    if (result$success) {
      return(list(success = TRUE, project_id = project_id, data = project_data))
    } else {
      return(result)
    }
  }
  
  # Get user projects
  get_user_projects <- function(user_uid, token) {
    # Query projects where user is owner or collaborator
    filters <- list(
      list(field = "owner_uid", op = "EQUAL", value = user_uid)
    )
    
    result <- db$query("projects", filters, token, limit = 100)
    
    if (result$success) {
      return(list(success = TRUE, projects = result$documents))
    } else {
      return(result)
    }
  }
  
  # Update project
  update_project <- function(project_id, updates, token) {
    # Get existing project first
    existing <- db$get("projects", project_id, token)
    
    if (!existing$success) {
      return(existing)
    }
    
    # Merge updates
    updated_data <- c(existing$data, updates)
    updated_data$updated_at <- as.character(Sys.time())
    
    return(db$set("projects", project_id, updated_data, token))
  }
  
  # Save analysis results to project
  save_analysis_results <- function(project_id, analysis_type, results, token) {
    analysis_id <- paste0("analysis_", digest(paste(analysis_type, Sys.time()), "md5"))
    
    analysis_data <- list(
      id = analysis_id,
      type = analysis_type,
      results = results,
      created_at = as.character(Sys.time()),
      project_id = project_id
    )
    
    # Save analysis as separate document
    analysis_result <- db$set("analyses", analysis_id, analysis_data, token)
    
    if (analysis_result$success) {
      # Update project to reference this analysis
      project_result <- db$get("projects", project_id, token)
      
      if (project_result$success) {
        existing_analyses <- project_result$data$analyses %||% list()
        updated_analyses <- c(existing_analyses, analysis_id)
        
        update_result <- update_project(project_id, list(analyses = updated_analyses), token)
        
        if (update_result$success) {
          return(list(success = TRUE, analysis_id = analysis_id))
        } else {
          return(update_result)
        }
      }
    }
    
    return(analysis_result)
  }
  
  # Return project management interface
  list(
    create = create_project,
    get_user_projects = get_user_projects,
    update = update_project,
    save_analysis = save_analysis_results
  )
}

# Export Firebase integration
list(
  config = FIREBASE_CONFIG,
  auth = firebase_auth,
  db = firestore_db,
  projects = project_manager
)