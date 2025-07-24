# 🔧 Prairie Genomics Suite R Shiny - Troubleshooting Guide

## 🚀 Quick Test Workflow

### 1. Test Data Upload
1. **Launch app**: `R -f run_app.R`
2. **Go to**: 📁 Data Upload tab
3. **Upload**: Select `test_data.csv` (provided)
4. **Settings**: Keep defaults (CSV, header row, gene names in first column)
5. **Click**: "Process Data" button
6. **Expected**: Green success message and data preview table

### 2. Test Sample Annotation  
1. **Go to**: 🧬 Sample Annotation tab
2. **Click**: "🔍 Detect Sample Patterns" 
3. **Expected**: Pattern suggestions appear with confidence scores
4. **Click**: "✅ Accept Suggestion" for highest confidence pattern
5. **Expected**: Green success message and annotation table appears

### 3. Test DESeq2 Analysis
1. **Go to**: 🚀 DESeq2 Analysis tab
2. **Settings**: Keep default parameters
3. **Click**: "🚀 Run DESeq2 Analysis"
4. **Expected**: Analysis completes and results summary appears

## ❌ Common Issues & Solutions

### Issue: "Invalid color" errors
**Fix**: ✅ **FIXED** - Updated all valueBox colors to valid Shiny dashboard colors

### Issue: "No package called shinyjs"
**Fix**: ✅ **FIXED** - Removed shinyjs dependencies, using standard Shiny conditionalPanel

### Issue: showNotification type errors
**Fix**: ✅ **FIXED** - Updated all notification types to valid values (default, message, warning, error)

### Issue: Process Data button not working
**Symptoms**: Button click has no effect or throws error
**Fix**: ✅ **FIXED** - Fixed readr namespace issues and data processing logic

### Issue: readxl package missing
**Solution**: 
```r
install.packages("readxl")
```
Or use CSV/TSV files instead of Excel

### Issue: DESeq2 package missing
**Solution**:
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## 🧪 Testing Commands

### Quick Package Check
```r
# Test all required packages
source("simple_test.R")
```

### Install Missing Packages
```r
# Run installation script
source("install.R")
```

### Check App Structure
```r
# Verify all files exist
file.exists(c("app.R", "modules/data_upload.R", "modules/sample_annotation.R", 
              "modules/deseq2_analysis.R", "modules/visualization.R"))
```

## 📊 Test Data Information

**File**: `test_data.csv`
- **Format**: CSV with header
- **Dimensions**: 10 genes × 6 samples
- **Groups**: 3 Control samples + 3 Treatment samples
- **Gene IDs**: ENSG format (standard Ensembl)
- **Expression**: Simulated count data

**Sample Names Pattern**:
- Control_1, Control_2, Control_3
- Treatment_1, Treatment_2, Treatment_3

**Expected Detection**: Prefix-based grouping with high confidence (>80%)

## 🔍 Debug Mode

### Enable Detailed Logging
```r
# Add to start of app.R
options(shiny.trace = TRUE)
options(shiny.error = browser)
```

### Check Console Output
- Look for error messages in R console
- Check browser developer tools (F12) for JavaScript errors
- Monitor network tab for failed requests

## 📞 Getting Help

### Error Reporting
When reporting issues, include:
1. **R version**: `R.version.string`
2. **Package versions**: Run `simple_test.R`
3. **Error message**: Copy full error from console
4. **Steps to reproduce**: What you clicked/uploaded
5. **Browser**: Chrome/Firefox/Safari and version

### Common Fixes
1. **Restart R session**: Sometimes package conflicts resolve with fresh start
2. **Clear browser cache**: Hard refresh (Ctrl+F5) can fix display issues  
3. **Update packages**: `update.packages(ask = FALSE)`
4. **Check working directory**: Must be in `prairie-genomics-suite_shiny/`

## ✅ Success Indicators

**Data Upload Success**:
- ✅ Green "Data processing complete" message
- ✅ Data preview table appears
- ✅ Summary boxes show correct gene/sample counts

**Sample Annotation Success**:
- ✅ Pattern detection finds groupings
- ✅ Confidence scores displayed
- ✅ Annotation table shows all samples assigned

**DESeq2 Analysis Success**:
- ✅ Analysis completes without errors
- ✅ Results summary shows significant genes
- ✅ Value boxes display statistics

**Visualization Success**:
- ✅ Interactive plots load without errors
- ✅ Plot controls respond to changes
- ✅ Export functionality works

---

**🧬 Prairie Genomics Suite R Shiny v1.0 - Making genomics analysis accessible!**