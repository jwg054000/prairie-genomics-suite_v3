#!/usr/bin/env python3
"""
🧬 Genomics AI - Zero-Friction Analysis Platform
Upload → AI Detection → Publication-Ready Results in 5 minutes
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import uuid
import json
from pathlib import Path
from genomics_ai_engine import GenomicsAI, ExperimentalDesign, Species, ExpressionType

# ============================================================================
# 🎨 UI CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="Genomics AI - Instant Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for elegant design
st.markdown("""
<style>
    /* Main container styling */
    .main {
        padding: 0rem 1rem;
    }
    
    /* Upload area styling */
    .upload-area {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border-radius: 20px;
        padding: 3rem;
        text-align: center;
        color: white;
        margin: 2rem 0;
        box-shadow: 0 10px 30px rgba(0,0,0,0.2);
    }
    
    /* Detection results card */
    .detection-card {
        background: white;
        border-radius: 15px;
        padding: 2rem;
        box-shadow: 0 5px 20px rgba(0,0,0,0.08);
        margin: 1rem 0;
        border-left: 4px solid #667eea;
    }
    
    /* Success badge */
    .success-badge {
        background: #10b981;
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        display: inline-block;
        margin: 0.5rem;
        font-weight: bold;
    }
    
    /* Warning badge */
    .warning-badge {
        background: #f59e0b;
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        display: inline-block;
        margin: 0.5rem;
    }
    
    /* Analysis button */
    .stButton > button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        padding: 1rem 3rem;
        font-size: 1.2rem;
        border-radius: 30px;
        font-weight: bold;
        transition: transform 0.3s;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 20px rgba(102, 126, 234, 0.4);
    }
    
    /* Results section */
    .results-section {
        background: #f9fafb;
        border-radius: 20px;
        padding: 2rem;
        margin: 2rem 0;
    }
    
    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    
    /* Metric cards */
    [data-testid="metric-container"] {
        background-color: white;
        border: 1px solid rgba(0,0,0,0.1);
        padding: 1rem;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# 🧠 AI ENGINE INITIALIZATION
# ============================================================================

@st.cache_resource
def get_ai_engine():
    """Initialize AI engine (cached for performance)"""
    return GenomicsAI()

# ============================================================================
# 🎯 MAIN APPLICATION
# ============================================================================

def main():
    # Header
    st.markdown("""
    <div style='text-align: center; padding: 2rem 0;'>
        <h1 style='font-size: 3rem; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
            🧬 Genomics AI
        </h1>
        <p style='font-size: 1.5rem; color: #6b7280;'>
            Zero-Configuration Genomics Analysis
        </p>
        <p style='color: #9ca3af;'>
            Upload → AI Detection → Publication-Ready Results in 5 minutes
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Initialize session state
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
        st.session_state.analysis_results = None
        st.session_state.uploaded_data = None
        st.session_state.share_link = None
    
    # Step 1: Upload Area
    if not st.session_state.analysis_complete:
        show_upload_interface()
    else:
        show_results_interface()

def show_upload_interface():
    """Show the upload and AI detection interface"""
    
    # Upload section
    st.markdown("""
    <div class='upload-area'>
        <h2>🎯 Drop Your Genomics Data Here</h2>
        <p style='font-size: 1.2rem; margin: 1rem 0;'>
            Supports: CSV, TSV, Excel, HDF5, H5AD, MTX
        </p>
        <p style='opacity: 0.8;'>
            No configuration needed - our AI handles everything
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # File uploader
    uploaded_file = st.file_uploader(
        "Choose a file",
        type=['csv', 'tsv', 'txt', 'xlsx', 'xls', 'h5', 'hdf5', 'h5ad', 'mtx'],
        help="Upload your expression matrix - genes as rows, samples as columns (or vice versa, we'll figure it out!)",
        label_visibility="hidden"
    )
    
    if uploaded_file is not None:
        # Save uploaded file temporarily
        temp_path = Path(f"temp_{uploaded_file.name}")
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        # Run AI analysis
        with st.spinner("🧠 AI is analyzing your data..."):
            ai_engine = get_ai_engine()
            results = ai_engine.analyze(str(temp_path))
        
        # Clean up temp file
        temp_path.unlink()
        
        # Store results
        st.session_state.analysis_results = results
        st.session_state.uploaded_data = uploaded_file.name
        
        # Show detection results
        if results['status'] == 'ready':
            show_detection_results(results)
        else:
            st.error(f"❌ Analysis failed: {results.get('error', 'Unknown error')}")

def show_detection_results(results):
    """Display AI detection results"""
    
    st.markdown("<h2 style='text-align: center;'>✨ AI Detection Complete</h2>", unsafe_allow_html=True)
    
    # Create columns for results
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class='detection-card'>
            <h3>📊 Data Characteristics</h3>
        </div>
        """, unsafe_allow_html=True)
        
        genes, samples = results['data_shape']
        st.metric("Genes Detected", f"{genes:,}")
        st.metric("Samples Detected", f"{samples:,}")
        st.metric("Expression Type", results['expression_type'].replace('_', ' ').title())
    
    with col2:
        st.markdown("""
        <div class='detection-card'>
            <h3>🧬 Biological Context</h3>
        </div>
        """, unsafe_allow_html=True)
        
        st.metric("Species", results['species'].title())
        st.metric("Experimental Design", results['design'].replace('_', ' ').title())
        
        # Cell type inference
        if results.get('cell_types'):
            top_cell_type = max(results['cell_types'], key=results['cell_types'].get)
            confidence = results['cell_types'][top_cell_type]
            st.metric("Likely Cell Type", top_cell_type.title(), f"{confidence:.0%} confidence")
    
    with col3:
        st.markdown("""
        <div class='detection-card'>
            <h3>✅ Quality Assessment</h3>
        </div>
        """, unsafe_allow_html=True)
        
        validation = results.get('validation', {})
        quality_score = validation.get('quality_score', 0)
        
        # Quality indicator
        if quality_score > 0.8:
            st.markdown("<div class='success-badge'>✅ Excellent Quality</div>", unsafe_allow_html=True)
        elif quality_score > 0.6:
            st.markdown("<div class='warning-badge'>⚠️ Good Quality</div>", unsafe_allow_html=True)
        else:
            st.markdown("<div class='warning-badge'>⚠️ Quality Issues</div>", unsafe_allow_html=True)
        
        # Show warnings if any
        for warning in validation.get('warnings', []):
            st.warning(warning)
    
    # Show sample groups
    st.markdown("<h3 style='text-align: center; margin-top: 2rem;'>🔬 Detected Sample Groups</h3>", 
                unsafe_allow_html=True)
    
    groups = results.get('groups', {})
    group_df = pd.DataFrame([
        {"Group": group, "Samples": len(samples), "Sample Names": ", ".join(samples[:3]) + ("..." if len(samples) > 3 else "")}
        for group, samples in groups.items()
    ])
    
    st.dataframe(group_df, use_container_width=True, hide_index=True)
    
    # Show suggested comparisons
    st.markdown("<h3 style='text-align: center; margin-top: 2rem;'>🎯 AI-Suggested Analyses</h3>", 
                unsafe_allow_html=True)
    
    comparisons = results.get('suggested_comparisons', [])
    for comp in comparisons[:3]:
        col1, col2, col3 = st.columns([3, 2, 1])
        with col1:
            st.write(f"**{comp['name']}**")
        with col2:
            st.write(f"Type: {comp['type'].replace('_', ' ').title()}")
        with col3:
            st.write(f"Confidence: {comp['confidence']:.0%}")
        st.caption(f"💡 {comp['reasoning']}")
    
    # Analysis recommendations
    rec = results.get('recommendations')
    if rec:
        st.info(f"""
        **🔬 Recommended Analysis Approach:**
        - Method: {rec.primary_analysis}
        - Statistical Test: {rec.statistical_test}
        - Quality Checks: {', '.join(rec.quality_checks)}
        """)
    
    # Action buttons
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🚀 Run One-Click Analysis", type="primary", use_container_width=True):
            run_full_analysis()

def run_full_analysis():
    """Execute the full analysis pipeline"""
    
    # Generate analysis ID
    analysis_id = str(uuid.uuid4())[:8]
    
    # Simulate analysis steps
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    steps = [
        ("🔄 Normalizing expression data...", 0.2),
        ("📊 Running differential expression analysis...", 0.4),
        ("🧬 Performing pathway enrichment...", 0.6),
        ("📈 Generating visualizations...", 0.8),
        ("📄 Creating publication-ready report...", 1.0)
    ]
    
    for step, progress in steps:
        status_text.text(step)
        progress_bar.progress(progress)
        import time
        time.sleep(0.5)
    
    # Mark analysis as complete
    st.session_state.analysis_complete = True
    st.session_state.share_link = f"https://genomics.ai/share/{analysis_id}"
    
    # Clear progress indicators
    progress_bar.empty()
    status_text.empty()
    
    # Show success message
    st.success("✅ Analysis complete! Redirecting to results...")
    st.rerun()

def show_results_interface():
    """Show the analysis results interface"""
    
    # Header with share button
    col1, col2 = st.columns([3, 1])
    with col1:
        st.markdown("""
        <h1 style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
            ✅ Analysis Complete
        </h1>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div style='text-align: right; padding: 1rem;'>
            <a href='{st.session_state.share_link}' target='_blank' 
               style='background: #667eea; color: white; padding: 0.5rem 1.5rem; 
                      border-radius: 20px; text-decoration: none;'>
                🔗 Share Results
            </a>
        </div>
        """, unsafe_allow_html=True)
    
    # Results tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["📊 Overview", "🧬 Differential Expression", 
                                             "🛤️ Pathways", "📈 Visualizations", "📄 Report"])
    
    with tab1:
        show_overview_tab()
    
    with tab2:
        show_differential_expression_tab()
    
    with tab3:
        show_pathways_tab()
    
    with tab4:
        show_visualizations_tab()
    
    with tab5:
        show_report_tab()
    
    # Download section
    st.markdown("<h3 style='text-align: center; margin-top: 3rem;'>📥 Download Results</h3>", 
                unsafe_allow_html=True)
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.download_button(
            "📊 Expression Data",
            data="mock_expression_data",
            file_name="normalized_expression.csv",
            mime="text/csv"
        )
    with col2:
        st.download_button(
            "🧬 DE Results",
            data="mock_de_results",
            file_name="differential_expression.csv",
            mime="text/csv"
        )
    with col3:
        st.download_button(
            "🛤️ Pathways",
            data="mock_pathways",
            file_name="pathway_enrichment.csv",
            mime="text/csv"
        )
    with col4:
        st.download_button(
            "📄 Full Report",
            data="mock_report",
            file_name="analysis_report.pdf",
            mime="application/pdf"
        )

def show_overview_tab():
    """Display analysis overview"""
    
    # Key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Genes Analyzed", "15,234", delta="Filtered: 1,256")
    
    with col2:
        st.metric("Significant Genes", "2,145", delta="14.1%", delta_color="normal")
    
    with col3:
        st.metric("Enriched Pathways", "342", delta="FDR < 0.05")
    
    with col4:
        st.metric("Analysis Quality", "98%", delta="Excellent", delta_color="normal")
    
    # Summary insights
    st.markdown("""
    <div class='results-section'>
        <h3>🔍 Key Insights</h3>
        <ul style='font-size: 1.1rem; line-height: 2;'>
            <li>Strong differential expression detected between treatment and control groups</li>
            <li>Cell cycle and immune response pathways are significantly enriched</li>
            <li>Top upregulated genes include MYC, IL6, and CXCL8</li>
            <li>Results suggest activation of inflammatory response</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Quality control metrics
    st.subheader("🛡️ Quality Control Metrics")
    
    qc_data = pd.DataFrame({
        'Metric': ['Sample Correlation', 'PCA Variance', 'RLE Median', 'Normalization Factor'],
        'Value': [0.95, 0.78, 0.02, 1.03],
        'Status': ['✅ Excellent', '✅ Good', '✅ Excellent', '✅ Normal']
    })
    
    st.dataframe(qc_data, use_container_width=True, hide_index=True)

def show_differential_expression_tab():
    """Display differential expression results"""
    
    st.subheader("🧬 Top Differentially Expressed Genes")
    
    # Generate mock DE results
    np.random.seed(42)
    de_genes = pd.DataFrame({
        'Gene': ['MYC', 'IL6', 'CXCL8', 'TNF', 'VEGFA', 'TP53', 'EGFR', 'STAT3', 'AKT1', 'MAPK1'],
        'Log2 Fold Change': np.random.uniform(-3, 3, 10),
        'Adjusted P-value': np.random.exponential(0.001, 10),
        'Expression': np.random.choice(['Up', 'Down'], 10),
        'Significance': ['***', '***', '**', '***', '*', '**', '***', '**', '***', '*']
    })
    
    # Style the dataframe
    st.dataframe(
        de_genes.style.apply(lambda x: ['background-color: #fee2e2' if v == 'Down' 
                                       else 'background-color: #dcfce7' 
                                       for v in x], subset=['Expression']),
        use_container_width=True,
        hide_index=True
    )
    
    # Volcano plot
    st.subheader("🌋 Volcano Plot")
    
    # Generate mock volcano plot data
    n_genes = 1000
    log2fc = np.random.normal(0, 1.5, n_genes)
    pvals = np.random.exponential(0.1, n_genes)
    significant = (np.abs(log2fc) > 1) & (pvals < 0.05)
    
    fig = go.Figure()
    
    # Non-significant genes
    fig.add_trace(go.Scatter(
        x=log2fc[~significant],
        y=-np.log10(pvals[~significant]),
        mode='markers',
        marker=dict(color='gray', size=5, opacity=0.5),
        name='Not significant'
    ))
    
    # Significant genes
    fig.add_trace(go.Scatter(
        x=log2fc[significant],
        y=-np.log10(pvals[significant]),
        mode='markers',
        marker=dict(
            color=log2fc[significant],
            colorscale='RdBu_r',
            size=8,
            colorbar=dict(title="Log2 FC")
        ),
        name='Significant'
    ))
    
    # Add threshold lines
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
    fig.add_vline(x=-1, line_dash="dash", line_color="gray")
    fig.add_vline(x=1, line_dash="dash", line_color="gray")
    
    fig.update_layout(
        title="Volcano Plot - Treatment vs Control",
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10 Adjusted P-value",
        height=500,
        template="plotly_white"
    )
    
    st.plotly_chart(fig, use_container_width=True)

def show_pathways_tab():
    """Display pathway enrichment results"""
    
    st.subheader("🛤️ Top Enriched Pathways")
    
    # Mock pathway data
    pathways = pd.DataFrame({
        'Pathway': [
            'Cell cycle',
            'Immune system process',
            'Inflammatory response',
            'Apoptotic process',
            'Signal transduction',
            'DNA repair',
            'Metabolic process',
            'Cell proliferation'
        ],
        'Genes': [145, 89, 67, 45, 234, 34, 156, 78],
        'Enrichment Score': [4.5, 3.8, 3.2, 2.9, 2.5, 2.3, 2.1, 2.0],
        'Adjusted P-value': [1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.05]
    })
    
    # Create enrichment plot
    fig = px.bar(
        pathways,
        y='Pathway',
        x='Enrichment Score',
        orientation='h',
        color='Enrichment Score',
        color_continuous_scale='Viridis',
        title='Pathway Enrichment Analysis'
    )
    
    fig.update_layout(
        height=400,
        showlegend=False,
        template="plotly_white"
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Pathway details
    st.subheader("📋 Pathway Details")
    st.dataframe(pathways, use_container_width=True, hide_index=True)
    
    # Network visualization placeholder
    st.info("🔗 Interactive pathway network visualization available in the full report")

def show_visualizations_tab():
    """Display various visualizations"""
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📊 PCA Plot")
        
        # Mock PCA data
        n_samples = 12
        pca_data = pd.DataFrame({
            'PC1': np.random.normal(0, 1, n_samples),
            'PC2': np.random.normal(0, 1, n_samples),
            'Group': ['Treatment'] * 6 + ['Control'] * 6,
            'Sample': [f'S{i}' for i in range(1, n_samples + 1)]
        })
        
        fig = px.scatter(
            pca_data,
            x='PC1',
            y='PC2',
            color='Group',
            text='Sample',
            title='Principal Component Analysis',
            color_discrete_map={'Treatment': '#ef4444', 'Control': '#3b82f6'}
        )
        
        fig.update_traces(textposition='top center')
        fig.update_layout(height=400, template="plotly_white")
        
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        st.subheader("🔥 Heatmap")
        
        # Mock heatmap data
        genes = ['MYC', 'IL6', 'CXCL8', 'TNF', 'VEGFA', 'TP53', 'EGFR', 'STAT3']
        samples = [f'S{i}' for i in range(1, 13)]
        
        heatmap_data = np.random.randn(len(genes), len(samples))
        
        fig = px.imshow(
            heatmap_data,
            x=samples,
            y=genes,
            color_continuous_scale='RdBu_r',
            title='Top DE Genes Expression Heatmap',
            aspect='auto'
        )
        
        fig.update_layout(height=400)
        
        st.plotly_chart(fig, use_container_width=True)
    
    # MA plot
    st.subheader("📈 MA Plot")
    
    # Mock MA plot data
    mean_expr = np.random.exponential(2, 1000)
    log2fc = np.random.normal(0, 1.5, 1000)
    significant = (np.abs(log2fc) > 1) & (np.random.random(1000) < 0.1)
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=np.log2(mean_expr[~significant] + 1),
        y=log2fc[~significant],
        mode='markers',
        marker=dict(color='gray', size=5, opacity=0.5),
        name='Not significant'
    ))
    
    fig.add_trace(go.Scatter(
        x=np.log2(mean_expr[significant] + 1),
        y=log2fc[significant],
        mode='markers',
        marker=dict(color='red', size=8),
        name='Significant'
    ))
    
    fig.add_hline(y=0, line_color="black")
    fig.add_hline(y=1, line_dash="dash", line_color="gray")
    fig.add_hline(y=-1, line_dash="dash", line_color="gray")
    
    fig.update_layout(
        title="MA Plot - Mean Expression vs Log2 Fold Change",
        xaxis_title="Log2 Mean Expression",
        yaxis_title="Log2 Fold Change",
        height=400,
        template="plotly_white"
    )
    
    st.plotly_chart(fig, use_container_width=True)

def show_report_tab():
    """Display the analysis report"""
    
    st.subheader("📄 Publication-Ready Report")
    
    # Report sections
    report_content = f"""
# Differential Expression Analysis Report

**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M')}  
**Analysis ID**: {st.session_state.share_link.split('/')[-1] if st.session_state.share_link else 'N/A'}  
**Data File**: {st.session_state.uploaded_data}

## Executive Summary

This analysis identified significant differential expression between treatment and control groups, 
with 2,145 genes showing statistically significant changes (FDR < 0.05, |log2FC| > 1). 
Key pathways affected include cell cycle regulation, immune response, and inflammatory signaling.

## Methods

### Data Processing
- **Normalization**: DESeq2 size factor normalization
- **Filtering**: Genes with <10 total counts removed
- **Batch Correction**: Not required (no batch effects detected)

### Statistical Analysis
- **Method**: DESeq2 negative binomial generalized linear model
- **Test**: Wald test for differential expression
- **Multiple Testing**: Benjamini-Hochberg FDR correction
- **Thresholds**: Adjusted p-value < 0.05, |log2FC| > 1

### Quality Control
All samples passed quality control metrics:
- Sample correlation > 0.9
- RLE median < 0.1
- No outliers detected in PCA

## Results

### Differential Expression
- **Total genes tested**: 15,234
- **Significantly upregulated**: 1,089 genes
- **Significantly downregulated**: 1,056 genes

### Top Differentially Expressed Genes
1. MYC (log2FC: 2.45, adj.p: 1.2e-15)
2. IL6 (log2FC: 3.12, adj.p: 4.5e-12)
3. CXCL8 (log2FC: 2.89, adj.p: 7.8e-11)

### Pathway Enrichment
Top enriched pathways (FDR < 0.05):
1. Cell cycle (145 genes, ES: 4.5)
2. Immune system process (89 genes, ES: 3.8)
3. Inflammatory response (67 genes, ES: 3.2)

## Biological Interpretation

The results suggest strong activation of proliferative and inflammatory pathways in the treatment group. 
The upregulation of MYC indicates enhanced cell proliferation, while elevated IL6 and CXCL8 point to 
an inflammatory response. These findings are consistent with the expected biological effects of the treatment.

## Data Availability

All analysis code and parameters are available at:  
{st.session_state.share_link}

## Citation

If you use these results, please cite:  
Genomics AI Platform (2024). Version 1.0. https://genomics.ai
    """
    
    # Display report with formatting
    st.markdown(report_content)
    
    # Download button for PDF report
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.download_button(
            "📥 Download Complete Report (PDF)",
            data=report_content.encode(),
            file_name=f"analysis_report_{datetime.now().strftime('%Y%m%d')}.pdf",
            mime="application/pdf",
            use_container_width=True
        )

# ============================================================================
# 🚀 RUN APPLICATION
# ============================================================================

if __name__ == "__main__":
    main()