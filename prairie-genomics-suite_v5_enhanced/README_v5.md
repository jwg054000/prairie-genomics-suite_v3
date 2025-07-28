# 🚀 Prairie Genomics Suite v5 Enhanced

**Advanced Multi-Group RNA-seq Analysis Platform with Context7 Visualizations**

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-1.7+-green.svg)](https://shiny.rstudio.com/)
[![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-orange.svg)](https://bioconductor.org/packages/DESeq2/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## ✨ What's New in v5 Enhanced

Prairie Genomics Suite v5 represents a major advancement in genomics analysis platforms, incorporating advanced methodologies from leading research institutions and Context7-inspired visualization best practices.

### 🔬 **Scientific Enhancements**
- **Multi-group experimental designs** supporting complex comparisons (GroupA, GroupB, GroupC, etc.)
- **Advanced DESeq2 methodology** with enhanced statistical rigor
- **Comprehensive batch effect correction** using ComBat-seq and limma methods
- **All pairwise comparison analysis** with automatic contrast generation
- **Enhanced gene filtering** with customizable thresholds
- **Log fold change shrinkage** for improved accuracy

### 🎨 **Context7-Inspired Visualizations**
- **Accessibility-first design** with colorblind-friendly palettes
- **Publication-quality exports** in multiple formats (PNG, PDF, SVG)
- **Interactive volcano plots** with gene search and highlighting
- **Advanced PCA plots** with 2D/3D options and confidence ellipses
- **Smart heatmaps** with intelligent gene selection and clustering
- **Performance optimization** for large genomics datasets

### 🚀 **User Experience Improvements**
- **Smart pattern detection** for automatic sample annotation
- **Batch effect diagnostics** with visual feedback
- **Real-time progress tracking** throughout analysis pipeline
- **Enhanced error handling** with helpful guidance
- **Comprehensive export system** for results and visualizations

---

## 📋 Quick Start Guide

### 1. **Installation & Setup**

```r
# Install required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  shiny, shinydashboard, DT, plotly, ggplot2, dplyr, readr,
  DESeq2, sva, MASS, pheatmap, RColorBrewer, car, ggrepel,
  readxl, logger
)

# Run the enhanced application
shiny::runApp("app_enhanced.R")
```

### 2. **Data Upload**
- Upload RNA-seq count matrix (genes × samples)
- Supported formats: CSV, TSV, Excel
- Use clear sample naming patterns (e.g., `Group_Sample_Rep`)

### 3. **Sample Annotation**
- **Auto-detection**: Let the system detect patterns automatically
- **Manual assignment**: Drag-and-drop interface for custom grouping
- **Batch detection**: Automatic identification of potential batch effects

### 4. **Advanced Analysis**
- **Multi-group DESeq2**: All pairwise comparisons with statistical rigor
- **Batch correction**: Apply ComBat-seq or limma methods if needed
- **Quality control**: Comprehensive filtering and validation

### 5. **Context7 Visualizations**
- **Enhanced volcano plots**: Interactive with gene highlighting
- **3D PCA plots**: With confidence ellipses and transparency
- **Smart heatmaps**: Automatic gene selection and clustering

### 6. **Export Results**
- **Statistical results**: All comparisons with significance testing
- **Publication plots**: High-resolution, vector graphics
- **Comprehensive reports**: HTML and PDF formats

---

## 🏗️ Architecture Overview

```
prairie-genomics-suite_v5_enhanced/
├── app_enhanced.R                          # Main enhanced application
├── modules/
│   ├── enhanced_sample_annotation.R        # Multi-group annotation with batch detection
│   ├── enhanced_deseq2_analysis.R         # Advanced DESeq2 with Emory methodology
│   └── context7_visualizations.R          # Context7-inspired plotting system
├── utils/
│   └── batch_correction.R                 # ComBat-seq and limma utilities
├── test_data.csv                          # Simple test dataset
└── README_v5.md                           # This documentation
```

### **Key Components**

#### 🧬 **Enhanced Sample Annotation Module**
- **Smart pattern detection** with multiple strategies
- **Multi-group support** for complex experimental designs
- **Batch effect detection** from sample naming patterns
- **Confidence scoring** for automatic annotations
- **Manual override** with intuitive interface

#### 🚀 **Advanced DESeq2 Analysis Module**
- **Emory methodology** implementation with proper statistical practices
- **Multiple comparison support** with all pairwise contrasts
- **Batch correction integration** using ComBat-seq and limma
- **Enhanced filtering** with customizable parameters
- **LFC shrinkage** for improved fold change estimates

#### 🎨 **Context7 Visualization Module**
- **Accessibility-first** color schemes (colorblind-safe)
- **Interactive features** with advanced hover templates
- **Performance optimization** for large datasets
- **Publication quality** with vector graphics export
- **Real-time customization** with live parameter updates

---

## 🔬 Scientific Methodology

### **DESeq2 Analysis Pipeline**

Our enhanced DESeq2 implementation follows established best practices and includes:

1. **Data Preprocessing**
   ```r
   # Gene filtering (standard methodology)
   keep1 <- rowSums(counts) >= 10        # Total count threshold
   keep2 <- rowSums(counts > 0) >= 3     # Minimum samples expressed
   ```

2. **Batch Effect Correction**
   ```r
   # ComBat-seq for RNA-seq data
   corrected_counts <- ComBat_seq(counts, batch = batch_vector, group = conditions)
   
   # Alternative: limma removeBatchEffect  
   corrected_log <- removeBatchEffect(log_counts, batch = batch_vector, design = design)
   ```

3. **Statistical Analysis**
   ```r
   # Enhanced DESeq2 workflow
   dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ Condition)
   dds <- DESeq(dds, fitType = "parametric", test = "Wald")
   
   # Log fold change shrinkage
   res_lfc <- lfcShrink(dds, contrast = contrast, type = "normal")
   ```

4. **Multiple Comparisons**
   - All pairwise comparisons automatically generated
   - Proper multiple testing correction (Benjamini-Hochberg)
   - Effect size estimation with confidence intervals

### **Visualization Best Practices**

Following **Context7 principles** for genomics visualization:

- **Accessibility**: Colorblind-friendly palettes with high contrast
- **Clarity**: Intelligent labeling and hover information
- **Interactivity**: Real-time parameter adjustment
- **Performance**: Optimized for large genomics datasets
- **Publication Quality**: Vector graphics with proper typography

---

## 📊 Example Workflows

### **Workflow 1: Simple Two-Group Comparison**
```r
# 1. Upload data: Control vs Treatment (3 replicates each)
# 2. Auto-detection identifies groups with 95% confidence
# 3. Run standard DESeq2 analysis
# 4. Generate volcano plot with top 10 gene labels
# 5. Export results and publication-quality figures
```

### **Workflow 2: Multi-Group Analysis Example**
```r
# 1. Upload data: GroupA, GroupB, GroupC, GroupD, GroupE (5 samples each)
# 2. Smart detection identifies 5 groups from naming patterns
# 3. Generate 10 pairwise comparisons automatically
# 4. Apply batch correction if detected
# 5. Create comparison-specific heatmaps and volcano plots
# 6. Export comprehensive results with overlap analysis
```

### **Workflow 3: Time Course with Batch Effects**
```r
# 1. Upload data: Multiple timepoints with batch information
# 2. Batch effect detection identifies 3 sequencing runs
# 3. Apply ComBat-seq correction
# 4. Analyze temporal patterns with multiple contrasts
# 5. Generate 3D PCA with confidence ellipses
# 6. Export interactive HTML report
```

---

## 🎯 Key Features Deep Dive

### **🧬 Multi-Group Sample Annotation**

The enhanced annotation system supports complex experimental designs:

- **Pattern Detection Strategies**:
  - Last part of sample name (Sample_Group_Rep → Group)
  - First part identification (Group_Sample_Rep → Group)  
  - Biological keyword matching (Control, Treatment, etc.)
  - Custom separator patterns with regex support

- **Confidence Scoring**:
  - Group balance assessment
  - Minimum size validation
  - Pattern consistency evaluation
  - Alternative strategy comparison

- **Batch Effect Detection**:
  - Sample name pattern analysis
  - PCA-based batch association testing
  - Hierarchical clustering evaluation
  - Visual diagnostics with before/after plots

### **🚀 Advanced DESeq2 Implementation**

Our implementation follows rigorous statistical practices:

- **Enhanced Filtering**:
  - Customizable count thresholds
  - Sample-based expression requirements
  - Variance-based gene selection
  - Interactive parameter adjustment

- **Statistical Rigor**:
  - Proper dispersion estimation
  - Multiple testing correction
  - Effect size shrinkage
  - Confidence interval calculation

- **Batch Correction**:
  - ComBat-seq for count data
  - limma removeBatchEffect for continuous data
  - Design matrix integration
  - Diagnostic visualizations

### **🎨 Context7-Enhanced Visualizations**

Publication-quality plots with advanced interactivity:

- **Volcano Plots**:
  - Gene search and highlighting
  - Dynamic threshold adjustment
  - Advanced hover information
  - Multiple color schemes
  - Publication export options

- **PCA Plots**:
  - 2D and 3D visualization options
  - Confidence ellipses with transparency
  - Sample metadata integration
  - Interactive exploration

- **Heatmaps**:
  - Smart gene selection (significance, variance, custom)
  - Hierarchical clustering options
  - Sample annotation tracks
  - Multiple color palettes

---

## 💻 System Requirements

### **Minimum Requirements**
- **R**: Version 4.0 or higher
- **RAM**: 4GB (8GB+ recommended for large datasets)
- **Storage**: 1GB free space
- **Browser**: Modern browser with JavaScript enabled

### **Recommended Setup**
- **R**: Version 4.3+
- **RAM**: 16GB+ for datasets >50,000 genes
- **CPU**: Multi-core processor for parallel processing
- **Storage**: SSD for improved performance

### **Package Dependencies**

#### Core Packages
```r
# UI and visualization
shiny (>= 1.7.0)
shinydashboard (>= 0.7.0)  
DT (>= 0.18)
plotly (>= 4.10.0)
ggplot2 (>= 3.3.0)

# Data manipulation
dplyr (>= 1.0.0)
readr (>= 2.0.0)
```

#### Statistical Analysis
```r
# DESeq2 and batch correction
DESeq2 (>= 1.30.0)      # Bioconductor
sva (>= 3.40.0)         # ComBat-seq
limma (>= 3.46.0)       # Alternative batch correction

# Advanced statistics
MASS (>= 7.3-50)
car (>= 3.0-10)         # Confidence ellipses
```

#### Visualization Enhancement
```r
# Advanced plotting
pheatmap (>= 1.0.12)    # Enhanced heatmaps
RColorBrewer (>= 1.1-2) # Color palettes
ggrepel (>= 0.9.0)      # Intelligent labeling
rgl (>= 0.107.0)        # 3D visualization
```

#### Optional Packages
```r
readxl (>= 1.3.0)       # Excel file support
logger (>= 0.2.0)       # Enhanced logging
cluster (>= 2.1.0)      # Clustering metrics
```

---

## 📈 Performance Optimization

### **Large Dataset Handling**

The v5 enhanced version includes several performance optimizations:

- **Smart Sampling**: Automatic data subsampling for visualization
- **Progressive Loading**: Incremental data processing
- **Memory Management**: Efficient memory usage patterns
- **Caching System**: Results caching for repeated operations

### **Performance Benchmarks**

| Dataset Size | Processing Time | Memory Usage | Recommendation |
|--------------|----------------|--------------|----------------|
| < 10,000 genes | < 1 minute | < 2GB | Standard settings |
| 10,000-30,000 genes | 2-5 minutes | 4-8GB | Enable optimization |
| > 30,000 genes | 5-15 minutes | 8-16GB | Use sampling options |

### **Optimization Settings**

```r
# For large datasets, adjust these parameters:
- Max points for full interactivity: 10,000
- Max points for basic visualization: 50,000  
- Sampling threshold: 100,000
- Memory monitoring: Automatic garbage collection
```

---

## 🔧 Troubleshooting

### **Common Issues and Solutions**

#### **Installation Problems**
```r
# Issue: Bioconductor packages not installing
# Solution:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "sva", "limma"))
```

#### **Memory Issues**
```r
# Issue: Out of memory errors
# Solutions:
1. Increase R memory limit: memory.limit(16000)  # Windows
2. Enable data sampling in visualization settings
3. Process data in chunks for very large datasets
4. Close other applications to free RAM
```

#### **Data Format Issues**
```r
# Issue: Expression data not recognized
# Solutions:
1. Ensure genes are in rows, samples in columns
2. Use numeric data only (no text in expression values)
3. Check for missing values and handle appropriately
4. Verify sample names match annotation data
```

#### **Analysis Failures**
```r
# Issue: DESeq2 analysis fails
# Common causes and solutions:
1. Insufficient replicates: Need ≥2 per group (≥3 recommended)
2. All-zero genes: Automatically filtered in preprocessing
3. Design matrix problems: Check sample annotation format
4. Memory constraints: Use filtering to reduce gene set
```

### **Getting Help**

1. **Built-in Help**: Check the Help & Documentation tab
2. **System Status**: Monitor package availability and memory usage
3. **Log Files**: Review analysis logs for detailed error messages
4. **Community Support**: Visit our GitHub issues page
5. **Contact**: Email support with detailed error descriptions

---

## 📊 Data Format Specifications

### **Expression Data Format**

Your RNA-seq count matrix should follow these specifications:

```
Gene_ID,Sample1_Group1,Sample2_Group1,Sample1_Group2,Sample2_Group2
Gene_1,150,132,87,95
Gene_2,45,52,234,198
Gene_3,1203,1087,1456,1234
...
```

**Requirements**:
- **First column**: Gene identifiers (unique)
- **Subsequent columns**: Raw count data (integers)
- **Headers**: Clear sample names with group patterns
- **Format**: CSV, TSV, or Excel (.xlsx)

### **Sample Naming Conventions**

For optimal automatic detection, use consistent naming patterns:

**Recommended Patterns**:
```
# Pattern 1: Group_Sample_Replicate
Control_Sample_1, Control_Sample_2, Treatment_Sample_1, Treatment_Sample_2

# Pattern 2: Group_Replicate  
Control_1, Control_2, Control_3, Treatment_1, Treatment_2, Treatment_3

# Pattern 3: Sample_Group (Multi-group style)
GroupA_1, GroupA_2, GroupB_1, GroupB_2, GroupC_1, GroupC_2, GroupD_1, GroupD_2, GroupE_1, GroupE_2
```

**Batch Information** (if applicable):
```
# Include batch identifiers in sample names
Control_Batch1_1, Control_Batch1_2, Control_Batch2_1, Control_Batch2_2
```

### **Optional Clinical Data**

You can optionally upload a separate clinical/metadata file:

```
Sample,Condition,Batch,TimePoint,Treatment
Sample1_Control,Control,Batch1,0h,Vehicle
Sample2_Control,Control,Batch1,0h,Vehicle  
Sample1_Treatment,Treatment,Batch2,24h,Drug
Sample2_Treatment,Treatment,Batch2,24h,Drug
```

---

## 🎓 Educational Resources

### **Getting Started with RNA-seq Analysis**

1. **Introduction to DESeq2**: [Bioconductor Vignette](https://bioconductor.org/packages/DESeq2/)
2. **Batch Effect Correction**: [ComBat-seq Documentation](https://bioconductor.org/packages/sva/)
3. **Experimental Design**: Best practices for RNA-seq studies
4. **Statistical Interpretation**: Understanding p-values and fold changes

### **Advanced Topics**

1. **Multi-group Comparisons**: Statistical considerations for complex designs
2. **Batch Effect Detection**: Methods and diagnostic approaches
3. **Visualization Best Practices**: Creating publication-quality figures
4. **Performance Optimization**: Handling large genomics datasets

### **Example Studies**

The methodology implemented in v5 has been successfully used in:
- **Cell differentiation studies**
- **Immune cell profiling** projects
- **Multi-condition treatment** comparisons
- **Time-course expression** analyses

---

## 🤝 Contributing

We welcome contributions to Prairie Genomics Suite v5 Enhanced!

### **How to Contribute**

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### **Development Guidelines**

- Follow **R style conventions** (lintr package)
- Include **comprehensive documentation** for new features
- Add **unit tests** for statistical functions
- Ensure **backwards compatibility** when possible
- Test with **multiple dataset sizes** and types

### **Priority Areas for Contribution**

- Additional **batch correction methods**
- New **visualization types** (pathway plots, network diagrams)
- **Performance optimizations** for very large datasets
- **Additional export formats** (PowerPoint, LaTeX)
- **Integration** with other genomics platforms

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### **Citation**

If you use Prairie Genomics Suite v5 Enhanced in your research, please cite:

```
Prairie Genomics Suite v5 Enhanced: Advanced Multi-Group RNA-seq Analysis Platform
with Context7 Visualizations and Batch Effect Correction
GitHub Repository: https://github.com/your-username/prairie-genomics-suite-v5
```

---

## 🙏 Acknowledgments

### **Scientific Methodology**
- **Research community** for DESeq2 best practices and multi-group analysis workflows
- **Bioconductor Core Team** for DESeq2, sva, and limma packages
- **Context7 Team** for visualization accessibility and design principles

### **Open Source Community**
- **R Core Team** and **Shiny developers**
- **Plotly team** for interactive visualization capabilities
- **All contributors** to the genomics analysis ecosystem

### **Testing and Feedback**
- **Beta testers** from genomics research communities
- **Accessibility consultants** for colorblind-friendly design
- **Performance testing** with large-scale genomics datasets

---

## 📞 Support & Contact

### **Technical Support**
- **GitHub Issues**: [Report bugs and request features](https://github.com/your-username/prairie-genomics-suite-v5/issues)
- **Documentation**: Comprehensive help available in-app
- **Community Forum**: Connect with other users

### **Academic Collaboration**
- **Research partnerships**: Contact for collaborative projects
- **Method development**: Contribute new statistical approaches
- **Publication support**: Assistance with methods sections

### **Commercial Licensing**
- **Enterprise solutions**: Scalable deployment options
- **Custom development**: Tailored analysis pipelines
- **Training workshops**: On-site and virtual training available

---

**🧬 Transform your RNA-seq data into publishable insights with Prairie Genomics Suite v5 Enhanced!**

*Last updated: January 2025 - v5.0.0 Enhanced Release*