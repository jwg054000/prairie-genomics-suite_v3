# 🔧 Prairie Genomics Suite - Error Logging System Guide

## 📋 Overview

The Prairie Genomics Suite now includes a comprehensive error logging system designed to make debugging and bug fixing much easier. This system automatically captures errors, user actions, system information, and performance metrics.

## 🗂️ File Structure

```
prairie-genomics-suite_v5_enhanced/
├── logging_system.R          # Core logging infrastructure
├── debug_logs_viewer.R       # Tools for viewing and analyzing logs
├── logs/                     # Automatically created log directory
│   ├── prairie_genomics_2025-01-28.log    # Daily log files
│   ├── prairie_genomics_2025-01-27.log
│   └── ...
└── ERROR_LOGGING_GUIDE.md    # This guide
```

## 🚀 Quick Start

### Viewing Recent Logs
```r
# Load the debug viewer
source("debug_logs_viewer.R")

# View last 50 log entries
view_recent_logs(50)

# Show error summary for today
show_error_summary(1)

# Search for specific issues
search_logs("KEGG")
search_logs("pathway analysis failed")
```

### Quick Debugging Commands
```r
# Show recent errors only
debug_recent_errors()

# Show pathway analysis issues
debug_pathway_issues()

# Show user actions that led to errors
debug_user_actions()
```

## 📊 Log Levels

The system uses 4 log levels:

1. **ERROR** - Critical failures that prevent functionality
2. **WARNING** - Issues that don't stop execution but may cause problems
3. **INFO** - General information about app operations
4. **DEBUG** - Detailed information for development

## 🔍 What Gets Logged

### User Actions
Every user interaction is logged with context:
```
[2025-01-28 14:23:15] session_20250128_142301_1234 [INFO] USER_ACTION: User performed: PATHWAY_ANALYSIS_CLICKED
    Details: {
      "analysis_type": "GO",
      "species": "mouse",
      "deseq2_available": true,
      "deseq2_genes": 15357
    }
```

### Errors with Stack Traces
Complete error information for debugging:
```
[2025-01-28 14:25:30] session_20250128_142301_1234 [ERROR] PATHWAY_ANALYSIS: Pathway analysis failed
    Details: {
      "error_message": "object 'gene_list' not found",
      "analysis_type": "KEGG",
      "species": "mouse",
      "input_parameters": {
        "padj_cutoff": 0.05,
        "fc_cutoff": 1.5
      }
    }
    Stack Trace:
      run_pathway_analysis()
      enrichKEGG()
```

### Performance Metrics
Analysis timing and resource usage:
```
[2025-01-28 14:26:45] session_20250128_142301_1234 [INFO] PERFORMANCE: Analysis performance: GO
    Details: {
      "duration_seconds": 23.45,
      "input_gene_count": 400,
      "output_pathway_count": 67,
      "genes_per_second": 17.06,
      "memory_used_mb": 145.2
    }
```

### System Information
Captured at startup and with errors:
```
[2025-01-28 14:20:00] session_20250128_142301_1234 [INFO] SYSTEM_INFO: System information captured
    Details: {
      "r_version": "R version 4.3.0",
      "platform": "darwin",
      "os": "Darwin",
      "user": "researcher",
      "loaded_packages": ["shiny", "ggplot2", "DESeq2"]
    }
```

## 🛠️ Using the Logging System

### For Developers - Adding Logging to New Functions

```r
# Add error logging to any function
my_analysis_function <- function(data) {
  with_error_logging("MY_ANALYSIS", "Custom analysis function", {
    # Your analysis code here
    result <- complex_calculation(data)
    return(result)
  })
}

# Log user actions
log_user_action("BUTTON_CLICKED", list(
  button_id = "run_analysis",
  user_inputs = list(species = "mouse", threshold = 0.05)
))

# Log specific errors
log_error("DATA_VALIDATION", "Invalid gene count", list(
  expected_min = 100,
  actual_count = 50,
  file_name = "uploaded_data.csv"
))

# Log performance
start_time <- Sys.time()
# ... do analysis ...
end_time <- Sys.time()
log_analysis_performance("GO_ANALYSIS", 
                        as.numeric(difftime(end_time, start_time, units = "secs")),
                        gene_list, 
                        nrow(results))
```

### For Bug Fixing - Common Debugging Workflows

#### 1. User Reports "Pathway Analysis Failed"
```r
# Step 1: Load debugging tools
source("debug_logs_viewer.R")

# Step 2: Search for pathway issues
debug_pathway_issues()

# Step 3: Look for the user's session
search_logs("PATHWAY_ANALYSIS_CLICKED")

# Step 4: Check what preceded the error
view_recent_logs(200)
```

#### 2. App Crashes or Unexpected Behavior
```r
# Check recent errors
debug_recent_errors()

# Get error summary
show_error_summary(2)  # Last 2 days

# Search for specific error patterns
search_logs("Error in")
search_logs("could not find function")
```

#### 3. Performance Issues
```r
# Look for slow analyses
search_logs("PERFORMANCE")

# Check memory usage patterns
search_logs("memory_mb")

# Find timeout issues
search_logs("timeout")
```

## 📁 Log File Management

### Automatic Features
- **Daily Rotation**: New log file created each day
- **Size Limits**: Files automatically rotated when they exceed 10MB
- **File Limits**: Only keeps last 5 log files automatically

### Manual Maintenance
```r
# Clean old logs (keeps last 7 days)
clean_old_logs(7)

# View current log configuration
get_log_config()
```

## 🎯 Common Error Patterns and Solutions

### 1. Missing Gene Conversion
**Log Pattern**: `WARNING: Genes still appear to be Entrez IDs`
**Solution**: Check org.Mm.eg.db package installation

### 2. Pathway Analysis Inflation
**Log Pattern**: `Found 4038 enriched pathways`
**Solution**: Check stringent filtering parameters (qvalueCutoff)

### 3. Memory Issues
**Log Pattern**: `memory_used_mb: 950+`
**Solution**: Reduce gene list size or enable chunked processing

### 4. Species Detection Failures
**Log Pattern**: `Auto-detect species failed`
**Solution**: Verify gene ID format (ENSG vs ENSMUSG)

## 🔧 Configuration

### Changing Log Level
```r
# Show more detailed logs (DEBUG level)
LOG_CONFIG$current_level <- 4

# Show only errors and warnings
LOG_CONFIG$current_level <- 2
```

### Disabling Logging
```r
# Temporarily disable
LOG_CONFIG$enabled <- FALSE

# Re-enable
LOG_CONFIG$enabled <- TRUE
initialize_logging()
```

## 📊 Log Analysis Tips

### Finding Patterns
```r
# Most common errors
show_error_summary(7)

# User workflow analysis
search_logs("USER_ACTION")

# Performance bottlenecks
search_logs("duration_seconds.*[2-9][0-9]")  # > 20 seconds
```

### Session Tracking
Each user session gets a unique ID, making it easy to track a user's entire workflow:
```r
# Find all actions in a specific session
search_logs("session_20250128_142301_1234")
```

## 🚨 Emergency Debugging

If the app is completely broken:

1. **Check if logging system loaded**:
   ```bash
   grep "Logging system initialized" logs/prairie_genomics_*.log
   ```

2. **Find the last error before crash**:
   ```r
   view_recent_logs(50)
   ```

3. **Check system resources**:
   ```r
   search_logs("memory_mb")
   search_logs("SYSTEM_INFO")
   ```

## 📝 Best Practices

### For Developers
1. **Always use `with_error_logging()`** for critical functions
2. **Log user actions** at the start of event handlers
3. **Include relevant context** in error details
4. **Use appropriate log levels** (ERROR for failures, INFO for normal operations)

### For Bug Fixing
1. **Start with error summary** to understand scope
2. **Use session IDs** to track user workflows
3. **Look for patterns** in error categories
4. **Check performance logs** for resource issues

## 🎉 Success!

With this logging system, debugging should be much easier! The logs provide:
- **Complete context** for every error
- **User workflow tracking** to reproduce issues
- **Performance metrics** to identify bottlenecks
- **System information** for environment issues

Happy debugging! 🐛→✅