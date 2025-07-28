# 🐭🧬 SPECIES MISMATCH FIX - Human vs Mouse Gene Conversion

## **🔍 PROBLEM IDENTIFIED:**
You mentioned "are you accidentally trying to convert to human vs mouse genes?" - this suggests your data contains **mouse genes** but the app is set to convert them as **human genes**.

## **🎯 QUICK SOLUTION:**

### **Step 1: Identify Your Gene Type**
Look at the first column of your data file. Your gene IDs should look like one of these:

**🧬 Human Ensembl IDs:**
```
ENSG00000139618
ENSG00000012048  
ENSG00000111640
```

**🐭 Mouse Ensembl IDs:**
```
ENSMUSG00000051951
ENSMUSG00000020122
ENSMUSG00000057666
```

### **Step 2: Set Correct Species in UI**
1. Start the app: `Rscript run_app.R`
2. Go to **📁 Data Upload** tab
3. Upload your data file
4. **CRITICAL**: In the "Species" dropdown, select:
   - **"Human"** if your IDs start with `ENSG`
   - **"Mouse"** if your IDs start with `ENSMUSG`

### **Step 3: Verify Conversion Settings**
- ✅ Ensure "Convert gene IDs to symbols" checkbox is **CHECKED**
- ✅ Ensure species matches your gene ID format
- ✅ Click "🚀 Process Data"

## **🧪 TEST FILES PROVIDED:**

### **For Human Genes:**
Use `test_data_with_known_genes.csv`:
- Contains: `ENSG00000139618, ENSG00000012048, ENSG00000111640, ENSG00000075624`
- Should convert to: `BRCA2, BRCA1, GAPDH, ACTB`
- Set species to: **Human**

### **For Mouse Genes:**
Use `test_mouse_genes.csv`:
- Contains: `ENSMUSG00000051951, ENSMUSG00000020122, ENSMUSG00000057666, ENSMUSG00000030636`
- Should convert to mouse gene symbols
- Set species to: **Mouse**

## **🔧 TROUBLESHOOTING:**

### **If you see this error:**
```
❌ Local conversion failed: None of the keys entered are valid keys for 'ENSEMBL'
```
**Solution:** You have mouse genes but species is set to human (or vice versa)

### **If conversion rate is 0%:**
1. **Check species setting** - most common cause
2. **Check gene ID format** - ensure they're valid Ensembl IDs
3. **Check internet connection** - needed for biomaRt fallback
4. **Verify gene IDs are current** - old/deprecated IDs may not convert

### **If you still see Ensembl IDs instead of gene symbols:**
1. **Clear browser cache** - hard refresh (Ctrl+F5)
2. **Check processing messages** - look for conversion success messages
3. **Verify checkbox is checked** - gene conversion must be enabled
4. **Restart app** - stop R session and restart with `Rscript run_app.R`

## **📊 EXPECTED CONSOLE OUTPUT:**
```bash
# For correct species:
✅ org.Hs.eg.db loaded for fast human gene conversion
🚀 Using fast local gene conversion (org.Hs.eg.db) for human genes...
✅ Local human conversion completed: 4/4 genes (100%)

# For wrong species:
❌ Local conversion failed: None of the keys entered are valid keys for 'ENSEMBL'
```

## **🎯 SUMMARY:**
The gene conversion IS working - the issue is likely a **species mismatch**:
- Mouse genes (ENSMUSG...) being processed as Human
- Human genes (ENSG...) being processed as Mouse

**Set the correct species in the UI dropdown and the conversion should work perfectly.**