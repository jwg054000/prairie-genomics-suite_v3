# 🧬 Gene ID to Symbol Conversion Guide

## Overview

The Prairie Genomics Suite now includes automatic gene ID to symbol conversion in the Data Upload tab. This feature converts Ensembl gene IDs (like `ENSG00000000003`) to human-readable gene symbols (like `TSPAN6`).

## ✨ Features

- **Automatic Detection**: Recognizes Ensembl gene IDs and determines if conversion is appropriate
- **Species Support**: Human and mouse gene conversion
- **Robust Error Handling**: Graceful fallback when biomaRt is unavailable
- **Chunked Processing**: Handles large datasets efficiently
- **Progress Tracking**: Real-time conversion progress updates
- **Smart Caching**: Reduces redundant network requests

## 🎯 How to Use

### 1. Upload Your Data
- Navigate to the **Data Upload** tab
- Select your RNA-seq expression file (CSV, TSV, or Excel)
- Configure file format options (separator, header)

### 2. Configure Gene Conversion
- **Enable Conversion**: Check "Convert gene IDs to symbols" (enabled by default)
- **Select Species**: Choose "Human" or "Mouse" 
- The system will automatically detect if your genes are Ensembl IDs

### 3. Process Data
- Click "🚀 Process Data"
- Watch the progress bar for conversion status
- Conversion happens automatically after data validation

## 📊 What Happens During Conversion

1. **Detection**: System checks if ≥50% of genes match Ensembl ID pattern
2. **Connection**: Connects to Ensembl biomaRt database
3. **Chunked Processing**: Converts genes in batches of 1,000 to avoid timeouts
4. **Mapping**: Maps gene IDs to official gene symbols
5. **Fallback**: Uses original IDs when symbols aren't available
6. **Statistics**: Reports conversion success rate

## 🎮 Example Workflow

```
Upload: genes.csv with ENSG00000000003, ENSG00000000005, ...
↓
Enable conversion: ✅ Human species
↓
Process: 🚀 Process Data
↓
Result: TSPAN6, TNMD, ... (with 85% conversion rate)
```

## ⚙️ Configuration Options

### Enable/Disable Conversion
- **Enabled**: Attempts to convert Ensembl IDs to symbols
- **Disabled**: Keeps original gene IDs unchanged

### Species Selection
- **Human**: Uses `hsapiens_gene_ensembl` dataset
- **Mouse**: Uses `mmusculus_gene_ensembl` dataset

## 🛡️ Error Handling

### Network Issues
- **Connection Timeout**: Falls back to original gene IDs
- **Service Unavailable**: Continues with IDs, shows warning
- **Partial Failure**: Reports conversion rate, mixed symbols/IDs

### Data Issues
- **Non-Ensembl IDs**: Automatically skips conversion
- **Mixed ID Types**: Attempts conversion, reports success rate
- **Duplicate Symbols**: Adds suffixes (`GENE_1`, `GENE_2`)

## 📈 Performance

### Speed
- **Small datasets** (<1,000 genes): ~10-30 seconds
- **Large datasets** (>10,000 genes): ~2-5 minutes
- **Timeout protection**: 30-second limit per chunk

### Memory
- **Efficient chunking**: Processes 1,000 genes at a time
- **Smart cleanup**: Removes temporary data after conversion
- **Memory monitoring**: Tracks usage during conversion

## 🔧 Troubleshooting

### Common Issues

#### "Most genes don't appear to be Ensembl IDs"
- **Cause**: Your gene IDs don't match Ensembl format
- **Solution**: This is normal for non-Ensembl data; conversion is skipped automatically

#### "biomaRt connection failed"
- **Cause**: Internet connection or Ensembl server issues
- **Solution**: The app continues with original IDs; try again later

#### "Conversion rate: 0%"
- **Cause**: Gene IDs not found in selected species database
- **Solution**: Check species selection or verify gene ID format

### Performance Tips

1. **Stable Internet**: Ensure reliable connection for biomaRt
2. **Species Match**: Select correct species for your data
3. **Batch Processing**: Large datasets are automatically chunked
4. **Memory Management**: Close other applications for large conversions

## 📋 Supported Gene ID Formats

### ✅ Supported
- `ENSG00000000003` (Human Ensembl)
- `ENSG00000000003.14` (With version numbers)
- `ENSMUSG00000000001` (Mouse Ensembl)

### ❌ Not Supported (but app continues normally)
- Gene symbols already (`TP53`, `BRCA1`)
- RefSeq IDs (`NM_000546`)
- Entrez IDs (`7157`)
- Custom identifiers

## 🎯 Benefits

### For Users
- **Readable Results**: Gene symbols instead of cryptic IDs
- **Better Interpretation**: Easier to understand biological meaning
- **Downstream Analysis**: Compatible with pathway analysis tools
- **Publication Ready**: Standard gene nomenclature

### For Analysis
- **Standardized Names**: Consistent gene naming across analyses
- **Pathway Mapping**: Better compatibility with gene set databases
- **Literature Search**: Easier to find relevant publications
- **Cross-Platform**: Consistent with other genomics tools

## 📞 Support

### If Conversion Fails
- Check internet connection
- Verify gene ID format
- Try different species selection
- Disable conversion to proceed with original IDs

### If Performance is Slow
- Reduce dataset size for testing
- Check network speed
- Consider running during off-peak hours
- Use local biomaRt installation (advanced users)

---

**🧬 Gene conversion makes your genomics data more interpretable and analysis-ready!**