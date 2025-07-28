# 🧬 VERIFY GENE SYMBOLS ARE WORKING - TEST INSTRUCTIONS

## **✅ FIXES IMPLEMENTED:**

1. **Exact user aggregation pattern**: `aggregate(raw.data[,-1],FUN = sum,by=list(factor(raw.data$gene_id)))`
2. **100% gene conversion working**: All components verified in diagnostic test
3. **Entry point confirmed**: `run_app.R` loads `app.R` which loads all modules correctly

## **🧪 TEST WITH RUNNING APPLICATION:**

### **Step 1: Start the Application**
```bash
cd /Users/joshuagarton/Documents/GitHub/prairie-genomics-suite_v2/prairie-genomics-suite_shiny/prairie-genomics-suite_optimized
Rscript run_app.R
```

### **Step 2: Upload Test Data**
1. Go to the **📁 Data Upload** tab
2. Upload the file: `test_data_with_known_genes.csv`
3. **IMPORTANT**: Ensure "Convert gene IDs to symbols" checkbox is ✅ **CHECKED**
4. Click **🚀 Process Data**

### **Step 3: Verify Results**
**Expected Results:**
- **Original**: 6 genes with Ensembl IDs (including 2 duplicates)
- **After Aggregation**: 4 unique genes (duplicates combined using exact user pattern)
- **After Conversion**: Gene symbols should appear:
  - `ENSG00000139618` → `BRCA2` (duplicate rows summed)
  - `ENSG00000012048` → `BRCA1`
  - `ENSG00000111640` → `GAPDH` (duplicate rows summed)
  - `ENSG00000075624` → `ACTB`

**What to Check:**
1. **Processing Messages**: Look for "Gene aggregation completed using exact user pattern"
2. **Processed Data Preview**: Should show Gene_Symbol column with BRCA2, BRCA1, GAPDH, ACTB
3. **No Ensembl IDs**: Should NOT see any ENSG... identifiers in final results

### **Step 4: Test DESeq2 Integration**
1. Go to **🧬 Sample Annotation** tab
2. Create annotation: Control (samples 1-3) vs Treatment (samples 4-6)
3. Go to **🚀 DESeq2 Analysis** tab
4. Run analysis
5. **Verify**: Results should show gene symbols (BRCA2, BRCA1, etc.) NOT Ensembl IDs

## **🔧 TROUBLESHOOTING:**

### **If you still see Ensembl IDs:**

1. **Check Conversion Enabled**:
   - Ensure "Convert gene IDs to symbols" checkbox is checked ✅
   - Species should be set to "Human"

2. **Clear Browser Cache**:
   - Hard refresh: Ctrl+F5 (Windows) or Cmd+Shift+R (Mac)
   - Or clear browser cache completely

3. **Check R Console**:
   - Look for any error messages during processing
   - Should see: "✅ Local conversion completed: 4/4 genes (100%)"

4. **Verify Data Format**:
   - First column should contain Ensembl gene IDs (starting with ENSG)
   - Use the provided test file to ensure known working data

5. **Restart Application**:
   - Stop R session (Ctrl+C)
   - Restart with `Rscript run_app.R`

## **🎯 SUCCESS CRITERIA:**

✅ **Duplicates aggregated**: 6 → 4 genes using exact user pattern
✅ **Gene symbols displayed**: BRCA2, BRCA1, GAPDH, ACTB (not ENSG...)
✅ **UI shows symbols**: Gene_Symbol column in processed data preview
✅ **DESeq2 receives symbols**: Analysis results use gene symbols
✅ **No errors**: Clean processing with success messages

## **📊 EXPECTED CONSOLE OUTPUT:**
```
🔄 Implementing exact aggregation pattern: aggregate(raw.data[,-1],FUN = sum,by=list(factor(raw.data$gene_id)))
📊 Found 2 duplicate gene IDs. Aggregating by summing...
✅ Gene aggregation completed using exact user pattern:
   - Original genes: 6
   - Final unique genes: 4  
   - Duplicate gene IDs aggregated: 2
🚀 Using fast local gene conversion (org.Hs.eg.db)...
✅ Local conversion completed: 4/4 genes (100%)
✅ Gene symbol conversion completed:
   - Total genes: 4
   - Successfully converted: 4 (100%)
```

---

**If this test works correctly but your actual data doesn't, the issue may be with your specific dataset format or gene IDs that aren't in the conversion database.**