# 🧬 Prairie Genomics Suite - R Shiny Implementation

An interactive R Shiny application for RNA-seq differential expression analysis with enhanced visualizations and user-friendly workflows.

## ✨ Features

### 📁 **Data Upload & Processing**
- Support for CSV, TSV, and Excel files
- Automatic file format detection
- Intelligent data validation and preprocessing
- Low-count gene filtering with customizable parameters

### 🧬 **Smart Sample Annotation**  
- **Automatic pattern detection** from sample names
- **Confidence-based suggestions** for experimental groupings
- **Manual annotation interface** with drag-and-drop functionality
- **Real-time validation** of experimental design

### 🚀 **DESeq2 Analysis**
- **One-click analysis** with optimal default parameters
- **Advanced parameter customization** for experienced users  
- **Real-time progress tracking** with detailed status updates
- **Comprehensive result summaries** with statistical metrics

### 📊 **Interactive Visualizations**
- **🌋 Volcano plots** with dynamic thresholds and gene highlighting
- **🔥 Expression heatmaps** with hierarchical clustering
- **📈 PCA plots** with confidence ellipses and metadata overlay
- **📊 MA plots** for quality assessment
- **Publication-ready exports** (PNG, PDF, SVG)

### 💫 **Enhanced User Experience**
- **Progress tracking** across the entire workflow
- **Smart navigation** with workflow-aware tab enabling
- **Contextual help** and guided parameter selection
- **Real-time statistics** and data summaries

## 🚀 Quick Start

### Prerequisites
- R (≥ 4.0.0)
- RStudio (recommended)

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/your-username/prairie-genomics-suite_v2.git
cd prairie-genomics-suite_v2/prairie-genomics-suite_shiny
```

2. **Install R dependencies**
```r
# Run the installation script
source("install.R")
```

3. **Launch the application**
```r
# In R or RStudio
shiny::runApp("app.R")
```

Or from command line:
```bash
R -e "shiny::runApp('app.R')"
```

## 📖 User Guide

### 1. Data Upload
1. Navigate to the **📁 Data Upload** tab
2. Choose your expression data file (genes × samples)
3. **Optional**: Upload a pAnno/annotation file for automatic sample grouping
   - Supports CSV, TSV, Excel formats
   - Should contain sample names and condition/group assignments
   - Will automatically match and apply to your expression data
4. Configure file reading options:
   - File type (auto-detected)
   - Header row presence
   - Gene names in first column
5. Set filtering parameters:
   - Minimum read count per gene
   - Minimum samples with counts
6. Click **Process Data** to validate and prepare your data

### 2. Sample Annotation
1. Go to **🧬 Sample Annotation** tab
2. **If pAnno file uploaded**: Annotation is applied automatically!
   - Sample groups are assigned based on your annotation file
   - Review the automatic assignments in the summary table
   - Proceed directly to DESeq2 analysis
3. **Without pAnno file - Automatic Detection**: Click **🔍 Detect Sample Patterns**
   - Review suggested groupings with confidence scores
   - Select the best suggestion and click **✅ Accept Suggestion**
4. **Manual Setup**: If automatic detection fails
   - Click **✏️ Manual Setup**
   - Assign each sample to experimental groups
   - Add custom group names as needed
   - Save your annotation

### 3. DESeq2 Analysis
1. Open **🚀 DESeq2 Analysis** tab
2. **Quick Analysis** (recommended):
   - Select comparison of interest
   - Use default parameters for optimal results
   - Click **🚀 Run DESeq2 Analysis**
3. **Advanced Options**:
   - Adjust p-value and fold change cutoffs
   - Modify statistical test parameters
   - Configure filtering options

### 4. Explore Results
1. Visit **📊 Visualizations** tab
2. **Volcano Plot**:
   - Adjust significance thresholds in real-time
   - Highlight genes of interest
   - Export publication-ready figures
3. **Heatmap**:
   - Customize gene selection and clustering
   - Apply different scaling methods
   - Choose color schemes
4. **PCA Plot**:
   - Explore sample relationships
   - Assess batch effects
   - Validate experimental design

### 5. Export Results
1. Go to **📋 Results Export** tab
2. Download options:
   - Complete DESeq2 results table
   - Significant genes only
   - Filtered expression data
   - High-resolution plots

## 🔧 Technical Details

### Architecture
```
prairie-genomics-suite_shiny/
├── app.R                      # Main Shiny application
├── modules/                   # Modular UI/Server components
│   ├── data_upload.R         # File upload and preprocessing
│   ├── sample_annotation.R   # Sample grouping and design
│   ├── deseq2_analysis.R     # Statistical analysis pipeline
│   └── visualization.R       # Interactive plots and charts
├── install.R                 # Dependency installation script
└── README.md                 # This documentation
```

### Key Technologies
- **Shiny**: Interactive web application framework
- **DESeq2**: Robust differential expression analysis
- **Plotly**: Interactive visualizations with zoom/pan
- **DT**: Enhanced data tables with search/filter
- **ggplot2**: Publication-quality static plots

### Performance Optimizations
- **Modular architecture** for maintainability and performance
- **Reactive programming** for efficient data flow
- **Progressive loading** of large datasets
- **Smart caching** of computation-intensive results
- **Background processing** for DESeq2 analysis

## 📊 Example Workflow

### Sample Data Format

#### Expression Data Format
Your expression data should be formatted as:
```
Gene        Sample1    Sample2    Sample3    Sample4
ENSG0001    150        200        175        225
ENSG0002    50         75         60         80
ENSG0003    300        280        320        290
...
```

#### pAnno/Annotation File Format (Optional)
For automatic sample annotation, provide a separate file with sample metadata:
```
Sample      Condition   Treatment   Batch
Sample1     Control     Vehicle     1
Sample2     Control     Vehicle     1  
Sample3     Treatment   Drug_A      2
Sample4     Treatment   Drug_A      2
```

**pAnno File Requirements:**
- **First column**: Sample names (must match expression data column names)
- **Second column**: Condition/Group assignments (e.g., Control, Treatment)
- **Additional columns**: Optional metadata (Treatment, Batch, TimePoint, etc.)
- **Supported formats**: CSV, TSV, Excel (.xlsx)
- **File naming**: Any filename ending in common annotation terms (pAnno, annotation, metadata, clinical, etc.)

### Expected Results
- **Volcano Plot**: Interactive exploration of differential expression
- **Heatmap**: Clustered visualization of top significant genes  
- **PCA Plot**: Sample relationships and quality assessment
- **Results Table**: Comprehensive statistics for all genes

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📝 Citation

If you use Prairie Genomics Suite in your research, please cite:

```
Prairie Genomics Suite: An Interactive R Shiny Platform for RNA-seq Analysis
[Your Name et al.]
[Journal/Preprint]
[Year]
```

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/your-username/prairie-genomics-suite_v2/issues)
- **Documentation**: See inline help and tooltips
- **Contact**: [your-email@domain.com]

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **DESeq2 team** for the robust statistical framework
- **Shiny team** for the excellent web application platform
- **Bioconductor community** for comprehensive genomics tools
- **R community** for the amazing ecosystem of packages

---

**🧬 Prairie Genomics Suite - Making genomics analysis accessible to everyone!**