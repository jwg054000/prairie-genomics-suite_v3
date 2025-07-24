#!/usr/bin/env python3
"""
Prairie Genomics Suite - Sample Annotation Demo App

This Streamlit app demonstrates the new optional clinical data functionality,
showcasing how users can perform differential expression analysis without
uploading clinical data files.

Features demonstrated:
- Interactive sample annotation
- Pattern-based auto-detection
- Template-based experimental designs
- Flexible validation and export
- Seamless DESeq2 integration

Run with: streamlit run sample_annotation_demo.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from core.sample_annotation_manager import create_sample_annotator
from core.data_models import ExpressionData
from ui.components.sample_annotation_widget import render_sample_annotation_widget
from analysis.deseq2_engine import DESeq2Engine

# Page configuration
st.set_page_config(
    page_title="Prairie Genomics - Sample Annotation Demo",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def create_demo_data():
    """Create realistic demo expression data with various sample naming patterns"""
    # Dataset options
    datasets = {
        "Control vs Treatment (Clear Pattern)": {
            "samples": [f"Control_{i:02d}" for i in range(1, 11)] + [f"Treatment_{i:02d}" for i in range(1, 11)],
            "description": "Clear control/treatment pattern - should auto-detect perfectly"
        },
        "Case vs Normal (Medical Study)": {
            "samples": [f"Case_{i:03d}" for i in range(1, 16)] + [f"Normal_{i:03d}" for i in range(1, 16)],
            "description": "Medical case-control study with 15 samples per group"
        },
        "Time Series (4 time points)": {
            "samples": [f"T{t}_Rep{r}" for t in range(4) for r in range(1, 4)],
            "description": "Time series with 4 time points, 3 replicates each"
        },
        "Challenging Names (No Pattern)": {
            "samples": [f"Sample_{i:03d}" for i in range(1, 21)] + ["weird.name-1", "weird.name-2"],
            "description": "Difficult sample names requiring manual annotation"
        },
        "Drug Dose Response": {
            "samples": [f"Vehicle_Rep{i}" for i in range(1, 5)] + 
                       [f"Low_Dose_Rep{i}" for i in range(1, 5)] +
                       [f"Med_Dose_Rep{i}" for i in range(1, 5)] +
                       [f"High_Dose_Rep{i}" for i in range(1, 5)],
            "description": "Drug dose-response study with 4 doses, 4 replicates each"
        }
    }
    
    return datasets

def generate_expression_data(sample_names):
    """Generate realistic expression data for given sample names"""
    np.random.seed(42)  # For reproducible demo
    
    n_genes = 2000
    n_samples = len(sample_names)
    
    # Generate gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    
    # Generate count data with some realistic structure
    base_expression = np.random.negative_binomial(n=5, p=0.3, size=(n_genes, n_samples))
    
    # Add some differential expression for samples with "Treatment", "Case", "High" etc.
    de_samples = [i for i, name in enumerate(sample_names) 
                  if any(keyword in name.upper() for keyword in ["TREATMENT", "CASE", "HIGH", "T3"])]
    
    if de_samples:
        # Make some genes differentially expressed
        de_genes = np.random.choice(n_genes, size=n_genes//10, replace=False)
        for gene_idx in de_genes:
            for sample_idx in de_samples:
                # Increase expression in "treatment" samples
                base_expression[gene_idx, sample_idx] *= np.random.uniform(2, 5)
    
    # Create DataFrame
    expression_df = pd.DataFrame(
        base_expression,
        index=gene_names,
        columns=sample_names
    )
    
    return ExpressionData(
        data=expression_df,
        samples=sample_names,
        genes=gene_names,
        metadata={
            "source": "demo_data",
            "differential_genes": len(de_samples) > 0,
            "n_de_samples": len(de_samples)
        }
    )

def main():
    """Main demo application"""
    st.title("🧬 Prairie Genomics Suite")
    st.markdown("### Sample Annotation Demo - Clinical Data is Now Optional!")
    
    st.markdown("""
    **🎉 New Feature Preview:** You can now perform differential expression analysis 
    without uploading clinical data files. The system provides multiple ways to 
    annotate your samples interactively.
    """)
    
    # Sidebar for dataset selection
    st.sidebar.title("📊 Demo Configuration")
    
    datasets = create_demo_data()
    selected_dataset = st.sidebar.selectbox(
        "Choose a demo dataset:",
        list(datasets.keys()),
        help="Each dataset demonstrates different sample naming patterns"
    )
    
    dataset_info = datasets[selected_dataset]
    st.sidebar.info(f"**About this dataset:**\n{dataset_info['description']}")
    
    # Generate expression data
    with st.spinner("Generating demo expression data..."):
        expression_data = generate_expression_data(dataset_info['samples'])
    
    # Display dataset overview
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Genes", f"{expression_data.n_genes:,}")
    with col2:
        st.metric("Samples", expression_data.n_samples)
    with col3:
        st.metric("Mean Expression", f"{expression_data.data.values.mean():.1f}")
    
    # Sample names preview
    with st.expander("📋 View Sample Names"):
        sample_df = pd.DataFrame({
            "Sample ID": expression_data.samples,
            "Index": range(len(expression_data.samples))
        })
        st.dataframe(sample_df, height=200)
    
    # Main content tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🎯 Quick Annotation", "🔧 Advanced Annotation", "📊 Analysis Preview", "📈 Results"])
    
    with tab1:
        render_quick_annotation_demo(expression_data)
    
    with tab2:
        render_advanced_annotation_demo(expression_data)
    
    with tab3:
        render_analysis_preview(expression_data)
    
    with tab4:
        render_results_preview(expression_data)

def render_quick_annotation_demo(expression_data):
    """Demonstrate quick annotation features"""
    st.subheader("⚡ Quick Annotation Demo")
    st.write("Perfect for simple experimental designs. Let the system detect patterns automatically!")
    
    # Initialize annotator
    if 'demo_annotator' not in st.session_state:
        st.session_state.demo_annotator = create_sample_annotator(expression_data)
    
    annotator = st.session_state.demo_annotator
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.write("**🔍 Pattern Detection Results:**")
        
        if st.button("🧠 Analyze Sample Names", key="analyze_patterns"):
            with st.spinner("Analyzing naming patterns..."):
                suggestions = annotator.suggest_groupings()
                st.session_state.suggestions = suggestions
        
        if 'suggestions' in st.session_state:
            suggestions = st.session_state.suggestions
            
            if suggestions:
                for name, suggestion in suggestions.items():
                    confidence = suggestion.get('confidence', 0)
                    groups = suggestion.get('groups', {})
                    
                    # Color code by confidence
                    if confidence > 0.7:
                        confidence_color = "🟢"
                    elif confidence > 0.4:
                        confidence_color = "🟡"
                    else:
                        confidence_color = "🔴"
                    
                    st.markdown(f"""
                    **{confidence_color} {suggestion.get('description', name.title())}**  
                    Confidence: {confidence:.1%} | Groups: {list(groups.keys())} | Sizes: {list(groups.values())}
                    """)
                    
                    if st.button(f"Apply {name.title()}", key=f"apply_{name}"):
                        annotator.apply_suggestion(name, suggestions)
                        st.success(f"✅ Applied {name} grouping!")
                        st.rerun()
            else:
                st.info("No clear patterns detected. Try the Advanced Annotation tab for manual assignment.")
    
    with col2:
        st.write("**📊 Current Annotation Status:**")
        
        summary = annotator.get_summary()
        
        # Progress bar
        progress = summary['n_annotated'] / summary['n_samples'] if summary['n_samples'] > 0 else 0
        st.progress(progress)
        st.write(f"Progress: {summary['n_annotated']}/{summary['n_samples']} samples annotated ({progress:.1%})")
        
        # Group distribution
        if summary['n_groups'] > 0:
            group_sizes = summary['group_sizes']
            
            # Create pie chart
            fig = px.pie(
                values=list(group_sizes.values()),
                names=list(group_sizes.keys()),
                title="Group Distribution"
            )
            fig.update_layout(height=300)
            st.plotly_chart(fig, use_container_width=True)
        
        # Validation status
        is_valid, errors = annotator.validate_annotations()
        if is_valid:
            st.success("✅ Annotations are valid!")
        elif summary['n_annotated'] > 0:
            st.warning(f"⚠️ Issues to resolve: {len(errors)}")
            for error in errors[:3]:  # Show first 3 errors
                st.write(f"• {error}")

def render_advanced_annotation_demo(expression_data):
    """Demonstrate advanced annotation widget"""
    st.subheader("🔧 Advanced Sample Annotation")
    st.write("Full-featured annotation system with drag-and-drop, pattern detection, and templates.")
    
    # Use the full sample annotation widget
    clinical_data = render_sample_annotation_widget(expression_data, key="advanced_demo")
    
    if clinical_data:
        st.session_state.demo_clinical_data = clinical_data
        st.balloons()  # Celebrate successful annotation!

def render_analysis_preview(expression_data):
    """Preview analysis capabilities"""
    st.subheader("🧪 Analysis Preview")
    st.write("See how your annotated samples would be used in differential expression analysis.")
    
    if 'demo_clinical_data' not in st.session_state:
        st.info("👆 Complete sample annotation in the tabs above to see analysis preview.")
        return
    
    clinical_data = st.session_state.demo_clinical_data
    
    # Analysis setup preview
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**📋 Analysis Configuration:**")
        
        condition_column = st.selectbox(
            "Condition column:",
            clinical_data.data.columns,
            index=0 if 'Condition' in clinical_data.data.columns else 0
        )
        
        conditions = clinical_data.data[condition_column].unique()
        
        if len(conditions) >= 2:
            reference_condition = st.selectbox(
                "Reference condition:",
                conditions,
                help="This group will be used as the reference for comparisons"
            )
            
            st.write("**Planned Comparisons:**")
            for condition in conditions:
                if condition != reference_condition:
                    st.write(f"• {condition} vs {reference_condition}")
        
    with col2:
        st.write("**📊 Sample Distribution:**")
        
        # Count samples per condition
        condition_counts = clinical_data.data[condition_column].value_counts()
        
        fig = px.bar(
            x=condition_counts.index,
            y=condition_counts.values,
            title="Samples per Condition",
            labels={'x': 'Condition', 'y': 'Number of Samples'}
        )
        fig.update_layout(height=300)
        st.plotly_chart(fig, use_container_width=True)
    
    # Analysis readiness check
    st.write("**🔬 Analysis Readiness Check:**")
    
    try:
        engine = DESeq2Engine()
        count_matrix, metadata = engine.prepare_data(
            expression_data=expression_data,
            clinical_data=clinical_data,
            condition_column=condition_column
        )
        
        st.success("✅ Data is ready for DESeq2 analysis!")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Final Genes", f"{count_matrix.shape[0]:,}")
        with col2:
            st.metric("Final Samples", count_matrix.shape[1])
        with col3:
            overlap = len(set(count_matrix.columns) & set(metadata.index))
            st.metric("Sample Overlap", f"{overlap}/{count_matrix.shape[1]}")
        
        # Store prepared data for results preview
        st.session_state.prepared_data = (count_matrix, metadata, condition_column)
        
    except Exception as e:
        st.error(f"❌ Analysis preparation failed: {str(e)}")

def render_results_preview(expression_data):
    """Show what the analysis results would look like"""
    st.subheader("📈 Analysis Results Preview")
    st.write("Simulated view of what your differential expression results would look like.")
    
    if 'prepared_data' not in st.session_state:
        st.info("👆 Complete the Analysis Preview step to see results simulation.")
        return
    
    count_matrix, metadata, condition_column = st.session_state.prepared_data
    
    # Simulate DE results
    with st.spinner("Simulating differential expression analysis..."):
        simulated_results = simulate_de_results(count_matrix, metadata, condition_column)
    
    if simulated_results is None:
        st.error("Could not simulate results - check your data configuration.")
        return
    
    # Results overview
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        total_genes = len(simulated_results)
        st.metric("Total Genes", f"{total_genes:,}")
    
    with col2:
        significant = sum((simulated_results['padj'] < 0.05) & (abs(simulated_results['log2FoldChange']) > 1))
        st.metric("Significant DE", f"{significant:,}")
    
    with col3:
        upregulated = sum((simulated_results['padj'] < 0.05) & (simulated_results['log2FoldChange'] > 1))
        st.metric("Upregulated", f"{upregulated:,}")
    
    with col4:
        downregulated = sum((simulated_results['padj'] < 0.05) & (simulated_results['log2FoldChange'] < -1))
        st.metric("Downregulated", f"{downregulated:,}")
    
    # Volcano plot
    st.write("**🌋 Volcano Plot:**")
    
    # Add significance categories
    simulated_results['significance'] = 'Not Significant'
    simulated_results.loc[
        (simulated_results['log2FoldChange'] > 1) & (simulated_results['padj'] < 0.05), 
        'significance'
    ] = 'Upregulated'
    simulated_results.loc[
        (simulated_results['log2FoldChange'] < -1) & (simulated_results['padj'] < 0.05), 
        'significance'
    ] = 'Downregulated'
    
    # Handle zeros and infinities in p-values
    simulated_results['neg_log10_padj'] = -np.log10(simulated_results['padj'].replace(0, 1e-300))
    
    # Create volcano plot
    fig = px.scatter(
        simulated_results,
        x='log2FoldChange',
        y='neg_log10_padj',
        color='significance',
        title='Differential Expression Results (Simulated)',
        hover_data=['gene'],
        color_discrete_map={
            'Upregulated': '#d62728',
            'Downregulated': '#1f77b4',
            'Not Significant': '#808080'
        }
    )
    
    # Add cutoff lines
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
    fig.add_vline(x=1, line_dash="dash", line_color="gray")
    fig.add_vline(x=-1, line_dash="dash", line_color="gray")
    
    fig.update_layout(
        xaxis_title="log₂(Fold Change)",
        yaxis_title="-log₁₀(Adjusted P-value)",
        height=500
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Top DE genes table
    st.write("**🧬 Top Differentially Expressed Genes:**")
    
    top_genes = simulated_results[
        (simulated_results['padj'] < 0.05) & 
        (abs(simulated_results['log2FoldChange']) > 1)
    ].sort_values('padj').head(10)
    
    if len(top_genes) > 0:
        display_cols = ['gene', 'log2FoldChange', 'pvalue', 'padj', 'significance']
        st.dataframe(
            top_genes[display_cols].round(4),
            use_container_width=True
        )
    else:
        st.info("No significantly differentially expressed genes found in simulation.")

def simulate_de_results(count_matrix, metadata, condition_column):
    """Simulate realistic DE analysis results"""
    try:
        conditions = metadata[condition_column].unique()
        if len(conditions) < 2:
            return None
        
        n_genes = count_matrix.shape[0]
        gene_names = count_matrix.index.tolist()
        
        # Simulate realistic DE results
        np.random.seed(42)  # For reproducible demo
        
        # Most genes are not DE
        log2fc = np.random.normal(0, 0.5, n_genes)
        
        # Make some genes significantly DE
        n_de_genes = n_genes // 10  # 10% DE genes
        de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
        
        # Half upregulated, half downregulated
        log2fc[de_indices[:n_de_genes//2]] = np.random.uniform(1.5, 4, n_de_genes//2)
        log2fc[de_indices[n_de_genes//2:]] = np.random.uniform(-4, -1.5, n_de_genes - n_de_genes//2)
        
        # P-values (smaller for truly DE genes)
        pvalues = np.random.uniform(0.01, 1.0, n_genes)
        pvalues[de_indices] = np.random.uniform(0.001, 0.05, n_de_genes)
        
        # Adjusted p-values (slightly higher)
        padj = pvalues * np.random.uniform(1.0, 5.0, n_genes)
        padj = np.minimum(padj, 1.0)  # Cap at 1.0
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'gene': gene_names,
            'log2FoldChange': log2fc,
            'pvalue': pvalues,
            'padj': padj
        })
        
        return results_df
        
    except Exception as e:
        st.error(f"Error simulating results: {e}")
        return None

if __name__ == "__main__":
    main()