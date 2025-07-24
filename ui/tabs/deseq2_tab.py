"""
Prairie Genomics Suite - DESeq2 Analysis Tab

This tab provides the complete DESeq2 differential expression analysis interface,
extracted from the monolithic structure for better modularity and maintainability.

Features:
- Experimental design setup
- DESeq2 parameter configuration
- Real-time analysis execution
- Results visualization and export
- Integration with other analysis engines
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any
import logging

# Import core modules
from core.session_manager import get_session_manager
from analysis.deseq2_engine import get_deseq2_engine
from core.async_analysis_manager import get_async_manager, AnalysisStatus
from core.enhanced_cache import get_enhanced_cache
from ui.components.sample_matching_diagnostics import render_sample_matching_diagnostics
from config import get_config

logger = logging.getLogger(__name__)

class DESeq2Tab:
    """
    DESeq2 Analysis Tab Component
    
    Provides complete interface for differential expression analysis using DESeq2.
    Handles experimental design, parameter selection, and results visualization.
    """
    
    def __init__(self):
        self.session = get_session_manager()
        self.deseq2_engine = get_deseq2_engine()
        self.async_manager = get_async_manager()
        self.cache = get_enhanced_cache()
        self.config = get_config("analysis")["deseq2"]
        
        # Track active analysis task
        self.current_task = None
    
    def render(self):
        """Render the complete DESeq2 analysis tab"""
        st.header("🧬 DESeq2 Differential Expression Analysis")
        st.markdown("Perform robust differential expression analysis using DESeq2 with R integration and Python fallbacks.")
        
        # Show cache statistics
        self._render_performance_metrics()
        
        # Check data requirements
        if not self._check_data_requirements():
            return
        
        # Sample matching diagnostics (helps debug sample ID issues)
        with st.expander("🔍 Sample Matching Diagnostics", expanded=False):
            st.write("Use this tool to diagnose and fix sample ID mismatch issues.")
            render_sample_matching_diagnostics("deseq2_diagnostics")
        
        # Check for active analysis task
        if self.current_task and not self.current_task.is_complete:
            self._render_active_analysis()
            return
        
        # Main analysis interface
        col1, col2 = st.columns([2, 1])
        
        with col1:
            self._render_experimental_design()
            self._render_analysis_parameters()
            self._render_analysis_execution()
        
        with col2:
            self._render_analysis_status()
            self._render_quick_results()
        
        # Results section
        if self.session.has_analysis_results("deseq2"):
            self._render_results_section()
    
    def _check_data_requirements(self) -> bool:
        """Check if required data is available"""
        if not self.session.has_data("expression_data"):
            st.error("❌ Expression data required. Please upload data in the Data Import tab.")
            return False
        
        if not self.session.has_data("clinical_data"):
            st.error("❌ Clinical data required for group comparisons. Please upload in the Data Import tab.")
            return False
        
        # Check sample overlap
        expr_data = self.session.get_data("expression_data")
        clinical_data = self.session.get_data("clinical_data")
        
        common_samples = set(expr_data.samples) & set(clinical_data.samples)
        if len(common_samples) < 6:  # Minimum for DESeq2
            st.error(f"❌ Insufficient overlapping samples ({len(common_samples)}). Need at least 6 samples for DESeq2.")
            return False
        
        return True
    
    def _render_experimental_design(self):
        """Render experimental design setup"""
        st.subheader("🎯 Experimental Design")
        
        clinical_data = self.session.get_data("clinical_data")
        
        with st.expander("Design Setup", expanded=True):
            # Condition variable selection
            categorical_vars = clinical_data.data.select_dtypes(include=['object', 'category']).columns.tolist()
            
            if not categorical_vars:
                st.error("No categorical variables found in clinical data for group comparisons.")
                return
            
            condition_column = st.selectbox(
                "Select condition variable",
                categorical_vars,
                help="Primary variable for differential expression comparison"
            )
            
            if condition_column:
                # Show condition distribution
                condition_counts = clinical_data.data[condition_column].value_counts()
                
                col1, col2 = st.columns([2, 1])
                with col1:
                    st.write("**Condition Distribution:**")
                    for condition, count in condition_counts.items():
                        st.write(f"• {condition}: {count} samples")
                
                with col2:
                    st.bar_chart(condition_counts)
                
                # Group selection for comparison
                conditions = condition_counts.index.tolist()
                if len(conditions) >= 2:
                    reference_group = st.selectbox(
                        "Reference group (control)",
                        conditions,
                        help="Reference condition for fold change calculation"
                    )
                    
                    treatment_group = st.selectbox(
                        "Treatment group",
                        [c for c in conditions if c != reference_group],
                        help="Treatment condition to compare against reference"
                    )
                    
                    # Store design in session
                    design_info = {
                        "condition_column": condition_column,
                        "reference_group": reference_group,
                        "treatment_group": treatment_group,
                        "conditions": conditions
                    }
                    self.session.store_analysis_params("deseq2_design", design_info)
                    
                    # Sample size check
                    ref_samples = (clinical_data.data[condition_column] == reference_group).sum()
                    treat_samples = (clinical_data.data[condition_column] == treatment_group).sum()
                    
                    if ref_samples < 3 or treat_samples < 3:
                        st.warning(f"⚠️ Small sample sizes: {reference_group}={ref_samples}, {treatment_group}={treat_samples}. Consider using at least 3 samples per group.")
                
            # Optional covariates
            st.write("**Optional Covariates:**")
            available_covariates = [col for col in clinical_data.data.columns if col != condition_column]
            
            if available_covariates:
                selected_covariates = st.multiselect(
                    "Select covariates to include in model",
                    available_covariates,
                    help="Additional variables to control for (e.g., batch, age, sex)"
                )
                
                if selected_covariates:
                    design_info["covariates"] = selected_covariates
                    self.session.store_analysis_params("deseq2_design", design_info)
    
    def _render_analysis_parameters(self):
        """Render DESeq2 analysis parameters"""
        st.subheader("⚙️ Analysis Parameters")
        
        with st.expander("DESeq2 Parameters"):
            col1, col2 = st.columns(2)
            
            with col1:
                # Statistical thresholds
                fc_cutoff = st.number_input(
                    "Log2 Fold Change Cutoff",
                    min_value=0.0,
                    max_value=5.0,
                    value=self.config.get("fc_cutoff", 1.0),
                    step=0.1,
                    help="Minimum absolute log2 fold change"
                )
                
                p_cutoff = st.number_input(
                    "Adjusted P-value Cutoff",
                    min_value=0.001,
                    max_value=0.2,
                    value=self.config.get("p_cutoff", 0.05),
                    step=0.01,
                    format="%.3f",
                    help="FDR-adjusted p-value threshold"
                )
                
                # Analysis method
                method = st.selectbox(
                    "Analysis Method",
                    ["Auto (R preferred)", "R (DESeq2)", "Python (fallback)"],
                    help="Choose analysis backend"
                )
            
            with col2:
                # Gene filtering
                min_count = st.number_input(
                    "Minimum Count Threshold",
                    min_value=1,
                    max_value=100,
                    value=self.config.get("min_count", 10),
                    help="Minimum read count per gene"
                )
                
                min_samples = st.number_input(
                    "Minimum Samples",
                    min_value=1,
                    max_value=10,
                    value=3,
                    help="Minimum samples with counts above threshold"
                )
                
                # Normalization
                normalize_method = st.selectbox(
                    "Normalization Method",
                    ["median_of_ratios", "TMM", "quantile"],
                    index=0,
                    help="Method for count normalization"
                )
            
            # Store parameters
            analysis_params = {
                "fc_cutoff": fc_cutoff,
                "p_cutoff": p_cutoff,
                "method": method,
                "min_count": min_count,
                "min_samples": min_samples,
                "normalize_method": normalize_method
            }
            self.session.store_analysis_params("deseq2_params", analysis_params)
    
    def _render_analysis_execution(self):
        """Render analysis execution section"""
        st.subheader("🚀 Run Analysis")
        
        # Check if design is complete
        if not self.session.has_analysis_params("deseq2_design"):
            st.info("Complete experimental design setup above to enable analysis.")
            return
        
        design_info = self.session.get_analysis_params("deseq2_design")
        
        # Analysis summary
        with st.expander("Analysis Summary", expanded=True):
            st.write(f"**Comparison:** {design_info['treatment_group']} vs {design_info['reference_group']}")
            st.write(f"**Condition Variable:** {design_info['condition_column']}")
            
            if "covariates" in design_info:
                st.write(f"**Covariates:** {', '.join(design_info['covariates'])}")
            
            # Engine status
            engine_status = self.deseq2_engine.get_status()
            method = "R (DESeq2)" if engine_status["r_available"] else "Python (fallback)"
            st.write(f"**Analysis Method:** {method}")
        
        # Run analysis button
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col2:
            if st.button("🚀 Run DESeq2 Analysis", type="primary", use_container_width=True):
                self._execute_deseq2_analysis_async()
    
    def _execute_deseq2_analysis_async(self):
        """Execute DESeq2 analysis asynchronously"""
        try:
            # Get data and parameters
            expression_data = self.session.get_data("expression_data")
            clinical_data = self.session.get_data("clinical_data")
            design_info = self.session.get_analysis_params("deseq2_design")
            analysis_params = self.session.get_analysis_params("deseq2_params")
            
            # Start async analysis
            self.current_task = self.async_manager.run_deseq2_async(
                expression_data=expression_data,
                clinical_data=clinical_data,
                condition_column=design_info["condition_column"],
                **analysis_params
            )
            
            st.success("🚀 Analysis started! Progress will be displayed below.")
            st.rerun()  # Refresh to show progress interface
            
        except Exception as e:
            st.error(f"❌ Failed to start analysis: {e}")
            logger.error(f"DESeq2 analysis startup failed: {e}")
    
    def _render_performance_metrics(self):
        """Render performance and cache metrics"""
        with st.expander("⚡ Performance Metrics", expanded=False):
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # Cache statistics
                cache_stats = self.cache.get_stats()
                st.metric("Cache Entries", cache_stats['memory_entries'])
                st.metric("Cache Size", f"{cache_stats['memory_size_mb']:.1f} MB")
            
            with col2:
                # Analysis statistics
                active_tasks = len(self.async_manager.active_tasks)
                st.metric("Active Tasks", active_tasks)
                st.metric("Engine Status", "R Available" if self.deseq2_engine.is_available() else "Python Fallback")
            
            with col3:
                # Quick actions
                if st.button("🧹 Clear Cache", use_container_width=True):
                    self.cache.clear_all()
                    st.success("Cache cleared!")
                    st.rerun()
                
                if st.button("📊 Show Task Manager", use_container_width=True):
                    with st.expander("Task Manager", expanded=True):
                        self.async_manager.render_task_manager_ui()
    
    def _render_active_analysis(self):
        """Render active analysis progress"""
        st.subheader("🔄 Analysis in Progress")
        
        # Update task status
        current_task = self.async_manager.get_task_status(self.current_task.task_id)
        if current_task:
            self.current_task = current_task
        
        # Show progress
        col1, col2 = st.columns([3, 1])
        
        with col1:
            self.async_manager.render_task_progress(self.current_task)
        
        with col2:
            if st.button("❌ Cancel Analysis", use_container_width=True):
                if self.async_manager.cancel_task(self.current_task.task_id):
                    st.warning("Analysis cancelled")
                    self.current_task = None
                    st.rerun()
        
        # Check if analysis completed
        if self.current_task.status == AnalysisStatus.SUCCESS:
            st.success("🎉 Analysis completed successfully!")
            
            # Store results in session
            if self.current_task.result:
                self.session.store_analysis_results("deseq2", self.current_task.result)
                st.info("Results saved to session. Refresh to view results.")
            
            self.current_task = None
            
            # Auto-refresh after short delay
            import time
            time.sleep(2)
            st.rerun()
        
        elif self.current_task.status == AnalysisStatus.FAILED:
            st.error(f"❌ Analysis failed: {self.current_task.error}")
            self.current_task = None
        
        # Auto-refresh every 2 seconds while running
        if not self.current_task or not self.current_task.is_complete:
            import time
            time.sleep(2)
            st.rerun()
    
    def _render_analysis_status(self):
        """Render analysis status sidebar"""
        st.subheader("📊 Analysis Status")
        
        # Engine status
        engine_status = self.deseq2_engine.get_status()
        
        r_status = "✅" if engine_status["r_available"] else "❌"
        st.write(f"{r_status} **R Integration**")
        if engine_status["r_available"]:
            st.write("   DESeq2 R package available")
        else:
            st.write("   Using Python fallback")
        
        # Data status
        expr_ready = self.session.has_data("expression_data")
        clinical_ready = self.session.has_data("clinical_data")
        design_ready = self.session.has_analysis_params("deseq2_design")
        
        expr_icon = "✅" if expr_ready else "❌"
        clinical_icon = "✅" if clinical_ready else "❌"
        design_icon = "✅" if design_ready else "❌"
        
        st.write(f"{expr_icon} **Expression Data**")
        st.write(f"{clinical_icon} **Clinical Data**")
        st.write(f"{design_icon} **Design Setup**")
        
        # Analysis readiness
        st.divider()
        if all([expr_ready, clinical_ready, design_ready]):
            st.success("🎉 Ready to run!")
        else:
            st.info("Complete setup above")
    
    def _render_quick_results(self):
        """Render quick results preview"""
        if not self.session.has_analysis_results("deseq2"):
            return
        
        st.subheader("📈 Quick Results")
        
        results = self.session.get_analysis_results("deseq2")
        
        # Key metrics
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Total Genes", results.n_genes)
            st.metric("Significant", results.n_significant)
        
        with col2:
            st.metric("Upregulated", results.n_upregulated)
            st.metric("Downregulated", results.n_downregulated)
        
        # Top genes preview
        if hasattr(results, 'results') and not results.results.empty:
            st.write("**Top 5 Significant Genes:**")
            top_genes = results.results.head(5)[['log2FoldChange', 'padj']]
            st.dataframe(top_genes, use_container_width=True)
    
    def _render_results_section(self):
        """Render comprehensive results section"""
        st.subheader("📊 Analysis Results")
        
        results = self.session.get_analysis_results("deseq2")
        
        # Results tabs
        result_tabs = st.tabs([
            "📋 Results Table", 
            "🌋 Volcano Plot", 
            "🔥 Heatmap", 
            "📁 Export"
        ])
        
        with result_tabs[0]:
            self._render_results_table(results)
        
        with result_tabs[1]:
            self._render_volcano_plot(results)
        
        with result_tabs[2]:
            self._render_heatmap(results)
        
        with result_tabs[3]:
            self._render_export_options(results)
    
    def _render_results_table(self, results):
        """Render results table with filtering"""
        st.write("**DESeq2 Results Table**")
        
        if not hasattr(results, 'results') or results.results.empty:
            st.info("No results to display")
            return
        
        # Filtering options
        col1, col2, col3 = st.columns(3)
        
        with col1:
            show_significant_only = st.checkbox("Show significant only", value=True)
        
        with col2:
            min_abs_fc = st.number_input("Min |log2FC|", min_value=0.0, value=1.0, step=0.1)
        
        with col3:
            max_padj = st.number_input("Max adj. p-value", min_value=0.001, max_value=1.0, value=0.05, step=0.01, format="%.3f")
        
        # Filter results
        df = results.results.copy()
        
        if show_significant_only:
            df = df[df['padj'] <= max_padj]
        
        df = df[abs(df['log2FoldChange']) >= min_abs_fc]
        
        # Display results
        st.write(f"Showing {len(df)} of {len(results.results)} genes")
        
        # Format for display
        display_cols = ['log2FoldChange', 'pvalue', 'padj', 'baseMean']
        if all(col in df.columns for col in display_cols):
            df_display = df[display_cols].copy()
            df_display['log2FoldChange'] = df_display['log2FoldChange'].round(3)
            df_display['pvalue'] = df_display['pvalue'].apply(lambda x: f"{x:.2e}" if x < 0.001 else f"{x:.4f}")
            df_display['padj'] = df_display['padj'].apply(lambda x: f"{x:.2e}" if x < 0.001 else f"{x:.4f}")
            df_display['baseMean'] = df_display['baseMean'].round(1)
            
            st.dataframe(df_display, use_container_width=True)
        else:
            st.dataframe(df, use_container_width=True)
    
    def _render_volcano_plot(self, results):
        """Render volcano plot"""
        st.write("**Volcano Plot**")
        
        if hasattr(results, 'plots') and 'volcano' in results.plots:
            # Use engine-generated plot
            st.plotly_chart(results.plots['volcano'], use_container_width=True)
        else:
            # Create simple volcano plot
            if hasattr(results, 'results') and not results.results.empty:
                df = results.results.copy()
                
                # Prepare data
                df['neg_log10_padj'] = -np.log10(df['padj'].replace(0, 1e-300))
                df['significant'] = (df['padj'] <= 0.05) & (abs(df['log2FoldChange']) >= 1.0)
                
                # Create scatter plot
                import plotly.express as px
                
                fig = px.scatter(
                    df.reset_index(),
                    x='log2FoldChange',
                    y='neg_log10_padj',
                    color='significant',
                    hover_data=['index'],
                    labels={
                        'log2FoldChange': 'Log2 Fold Change',
                        'neg_log10_padj': '-Log10(Adjusted P-value)',
                        'index': 'Gene'
                    },
                    title="Volcano Plot"
                )
                
                # Add threshold lines
                fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="red", opacity=0.5)
                fig.add_vline(x=1, line_dash="dash", line_color="red", opacity=0.5)
                fig.add_vline(x=-1, line_dash="dash", line_color="red", opacity=0.5)
                
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No results available for volcano plot")
    
    def _render_heatmap(self, results):
        """Render heatmap of top genes"""
        st.write("**Expression Heatmap (Top Significant Genes)**")
        
        if not hasattr(results, 'results') or results.results.empty:
            st.info("No results available for heatmap")
            return
        
        # Get top significant genes
        top_genes = results.get_significant_genes(max_genes=50)
        
        if len(top_genes) == 0:
            st.info("No significant genes found for heatmap")
            return
        
        # Get expression data
        expression_data = self.session.get_data("expression_data")
        clinical_data = self.session.get_data("clinical_data")
        design_info = self.session.get_analysis_params("deseq2_design")
        
        # Filter to common samples and top genes
        common_samples = list(set(expression_data.samples) & set(clinical_data.samples))
        heatmap_data = expression_data.data.loc[top_genes, common_samples]
        
        # Create annotation for conditions
        condition_column = design_info["condition_column"]
        sample_conditions = clinical_data.data.loc[common_samples, condition_column]
        
        # Simple heatmap with plotly
        import plotly.express as px
        
        # Z-score normalize
        heatmap_data_norm = heatmap_data.T  # Transpose for plotly
        heatmap_data_norm = (heatmap_data_norm - heatmap_data_norm.mean()) / heatmap_data_norm.std()
        
        fig = px.imshow(
            heatmap_data_norm.values,
            x=heatmap_data_norm.columns,
            y=heatmap_data_norm.index,
            aspect="auto",
            color_continuous_scale="RdBu_r",
            title=f"Expression Heatmap (Top {len(top_genes)} Genes)"
        )
        
        fig.update_layout(
            xaxis_title="Genes",
            yaxis_title="Samples",
            height=max(400, len(common_samples) * 15)
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    def _render_export_options(self, results):
        """Render export options"""
        st.write("**Export Results**")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Export full results
            if st.button("📁 Export Full Results (CSV)", use_container_width=True):
                if hasattr(results, 'results'):
                    csv = results.results.to_csv()
                    st.download_button(
                        label="Download CSV",
                        data=csv,
                        file_name="deseq2_results.csv",
                        mime="text/csv"
                    )
        
        with col2:
            # Export significant genes only
            if st.button("⭐ Export Significant Genes", use_container_width=True):
                significant_genes = results.get_significant_genes()
                if len(significant_genes) > 0:
                    sig_results = results.results.loc[significant_genes]
                    csv = sig_results.to_csv()
                    st.download_button(
                        label="Download Significant Genes",
                        data=csv,
                        file_name="deseq2_significant_genes.csv",
                        mime="text/csv"
                    )
                else:
                    st.info("No significant genes to export")
        
        # Gene lists for pathway analysis
        st.write("**Gene Lists for Downstream Analysis:**")
        
        if hasattr(results, 'results') and not results.results.empty:
            all_genes = results.results.index.tolist()
            significant_genes = results.get_significant_genes()
            upregulated = results.get_upregulated_genes()
            downregulated = results.get_downregulated_genes()
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.write(f"• **All tested genes:** {len(all_genes)}")
                st.write(f"• **Significant genes:** {len(significant_genes)}")
            
            with col2:
                st.write(f"• **Upregulated:** {len(upregulated)}")
                st.write(f"• **Downregulated:** {len(downregulated)}")
            
            # Store gene lists for other analyses
            gene_lists = {
                "all_genes": all_genes,
                "significant": significant_genes,
                "upregulated": upregulated,
                "downregulated": downregulated
            }
            self.session.store_analysis_results("deseq2_gene_lists", gene_lists)
            
            st.success("✅ Gene lists saved for pathway and literature analysis!")


def render_deseq2_tab():
    """Render the DESeq2 analysis tab (main entry point)"""
    tab = DESeq2Tab()
    tab.render()