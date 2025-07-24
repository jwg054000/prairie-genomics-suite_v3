"""
Prairie Genomics Suite - Data Import Tab

This tab handles all data import functionality including expression data,
clinical data, and sample annotations. Extracted from the monolithic structure
to provide a clean, modular interface.

Features:
- File upload with validation
- Data preview and quality checks  
- Sample matching and annotation
- Multiple file format support
- Session state management
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import logging

# Import core modules
from core.data_models import ExpressionData, ClinicalData
from core.session_manager import get_session_manager
from core.utils import validate_dataframe, safe_convert_numeric, clean_clinical_data
from core.optimized_data_loader import get_optimized_loader
from core.async_analysis_manager import get_async_manager
from core.sample_annotation_manager import create_sample_annotator
from ui.components.sample_annotation_widget import render_sample_annotation_widget
from config import get_config

logger = logging.getLogger(__name__)

class DataImportTab:
    """
    Data Import Tab Component
    
    Handles all data import operations with validation and preview capabilities.
    Replaces the data import section from the monolithic application.
    """
    
    def __init__(self):
        self.session = get_session_manager()
        self.optimized_loader = get_optimized_loader()
        self.async_manager = get_async_manager()
        self.config = get_config("validation")
    
    def render(self):
        """Render the complete data import tab"""
        st.header("🔬 Data Import & Preparation")
        st.markdown("Import your genomics data files for analysis. All file formats are validated automatically.")
        
        # Create columns for better layout
        col1, col2 = st.columns([2, 1])
        
        with col1:
            self._render_file_upload_section()
        
        with col2:
            self._render_data_summary()
        
        # Data preview section
        if (self.session.has_data("expression_data") or 
            self.session.has_data("clinical_data")):
            self._render_data_preview()
            self._render_sample_matching()
        
        # Sample annotation section (new feature)
        if self.session.has_data("expression_data"):
            self._render_sample_annotation_section()
    
    def _render_file_upload_section(self):
        """Render file upload interface"""
        st.subheader("📁 File Upload")
        
        # Expression data upload
        with st.expander("Expression Data", expanded=True):
            expr_file = st.file_uploader(
                "Upload Expression Data",
                type=['csv', 'tsv', 'txt', 'xlsx'],
                help="Upload gene expression matrix (genes × samples)",
                key="expression_file"
            )
            
            if expr_file:
                self._process_expression_file(expr_file)
        
        # Clinical data upload (now optional)
        with st.expander("Clinical Data (Optional)"):
            st.info("💡 **New!** Clinical data is now optional. You can use the Sample Annotation tool below instead.")
            
            clinical_file = st.file_uploader(
                "Upload Clinical Data",
                type=['csv', 'tsv', 'txt', 'xlsx'],
                help="Upload clinical annotations (samples × variables) - Optional if using Sample Annotation",
                key="clinical_file"
            )
            
            if clinical_file:
                self._process_clinical_file(clinical_file)
        
        # Sample sheet upload (optional)
        with st.expander("Sample Annotations (Optional)"):
            sample_file = st.file_uploader(
                "Upload Sample Sheet",
                type=['csv', 'tsv', 'txt', 'xlsx'],
                help="Upload sample annotations and metadata",
                key="sample_file"
            )
            
            if sample_file:
                self._process_sample_file(sample_file)
    
    def _process_expression_file(self, file):
        """Process uploaded expression data file"""
        try:
            # Read file based on extension
            file_ext = Path(file.name).suffix.lower()
            
            if file_ext == '.xlsx':
                df = pd.read_excel(file, index_col=0)
            else:
                # Try to detect separator
                delimiter = '\t' if file_ext in ['.tsv', '.txt'] else ','
                df = pd.read_csv(file, index_col=0, delimiter=delimiter)
            
            # Validate expression data
            valid, errors = validate_dataframe(
                df, 
                min_rows=100,  # Reasonable gene count
                min_cols=3,    # Minimum samples
                name="Expression data"
            )
            
            if not valid:
                st.error(f"Expression data validation failed: {errors}")
                return
            
            # Convert to numeric
            numeric_df = df.apply(safe_convert_numeric)
            
            # Create ExpressionData object
            expression_data = ExpressionData(
                data=numeric_df,
                samples=numeric_df.columns.tolist(),
                genes=numeric_df.index.tolist(),
                metadata={
                    "filename": file.name,
                    "upload_time": pd.Timestamp.now(),
                    "shape": numeric_df.shape
                }
            )
            
            # Store in session
            self.session.store_data("expression_data", expression_data)
            
            st.success(f"✅ Expression data loaded: {expression_data.n_genes} genes × {expression_data.n_samples} samples")
            
            # Show basic statistics
            with st.expander("Expression Data Summary"):
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Genes", expression_data.n_genes)
                with col2:
                    st.metric("Samples", expression_data.n_samples)
                with col3:
                    mean_expr = numeric_df.values[numeric_df.values > 0].mean()
                    st.metric("Mean Expression", f"{mean_expr:.2f}")
        
        except Exception as e:
            st.error(f"Error processing expression file: {e}")
            logger.error(f"Expression file processing failed: {e}")
    
    def _process_clinical_file(self, file):
        """Process uploaded clinical data file"""
        try:
            file_ext = Path(file.name).suffix.lower()
            
            if file_ext == '.xlsx':
                df = pd.read_excel(file, index_col=0)
            else:
                delimiter = '\t' if file_ext in ['.tsv', '.txt'] else ','
                df = pd.read_csv(file, index_col=0, delimiter=delimiter)
            
            # Clean clinical data first to handle NaN columns
            original_shape = df.shape
            cleaned_df = clean_clinical_data(
                df,
                remove_nan_columns=True,
                min_valid_fraction=0.1
            )
            
            # Use flexible validation for clinical data
            valid, errors = validate_dataframe(
                cleaned_df,
                min_rows=3,
                min_cols=1,
                name="Clinical data",
                allow_nan_columns=True,
                strict_validation=False
            )
            
            if not valid:
                st.warning(f"Clinical data validation issues: {errors}")
                st.info("💡 You can still use this data - issues will be handled automatically during analysis.")
            
            # Show cleaning summary if any changes were made
            if cleaned_df.shape != original_shape:
                st.info(f"📋 Data cleaning applied: {original_shape[1]} → {cleaned_df.shape[1]} columns")
            
            df = cleaned_df
            
            # Create ClinicalData object
            clinical_data = ClinicalData(
                data=df,
                samples=df.index.tolist(),
                variables=df.columns.tolist(),
                metadata={
                    "filename": file.name,
                    "upload_time": pd.Timestamp.now(),
                    "shape": df.shape
                }
            )
            
            # Store in session
            self.session.store_data("clinical_data", clinical_data)
            
            st.success(f"✅ Clinical data loaded: {clinical_data.n_samples} samples × {clinical_data.n_variables} variables")
            
            # Show variable types
            with st.expander("Clinical Data Summary"):
                numeric_vars = df.select_dtypes(include=[np.number]).columns.tolist()
                categorical_vars = df.select_dtypes(include=['object', 'category']).columns.tolist()
                
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Numeric Variables:**")
                    if numeric_vars:
                        for var in numeric_vars[:10]:  # Show first 10
                            st.write(f"• {var}")
                        if len(numeric_vars) > 10:
                            st.write(f"... and {len(numeric_vars) - 10} more")
                
                with col2:
                    st.write("**Categorical Variables:**")
                    if categorical_vars:
                        for var in categorical_vars[:10]:
                            st.write(f"• {var}")
                        if len(categorical_vars) > 10:
                            st.write(f"... and {len(categorical_vars) - 10} more")
        
        except Exception as e:
            st.error(f"Error processing clinical file: {e}")
            logger.error(f"Clinical file processing failed: {e}")
    
    def _process_sample_file(self, file):
        """Process uploaded sample annotation file"""
        try:
            file_ext = Path(file.name).suffix.lower()
            
            if file_ext == '.xlsx':
                df = pd.read_excel(file, index_col=0)
            else:
                delimiter = '\t' if file_ext in ['.tsv', '.txt'] else ','
                df = pd.read_csv(file, index_col=0, delimiter=delimiter)
            
            # Store sample annotations
            self.session.store_data("sample_annotations", df)
            
            st.success(f"✅ Sample annotations loaded: {len(df)} samples")
            
            with st.expander("Sample Annotation Summary"):
                st.write(f"**Annotation columns:** {', '.join(df.columns.tolist())}")
                
        except Exception as e:
            st.error(f"Error processing sample file: {e}")
            logger.error(f"Sample file processing failed: {e}")
    
    def _render_data_summary(self):
        """Render data summary sidebar"""
        st.subheader("📊 Data Status")
        
        # Check what data is loaded
        has_expression = self.session.has_data("expression_data")
        has_clinical = self.session.has_data("clinical_data")
        has_samples = self.session.has_data("sample_annotations")
        
        # Status indicators
        expr_status = "✅" if has_expression else "❌"
        clinical_status = "✅" if has_clinical else "❌"
        sample_status = "✅" if has_samples else "⚪"
        
        st.write(f"{expr_status} **Expression Data**")
        if has_expression:
            expr_data = self.session.get_data("expression_data")
            st.write(f"   {expr_data.n_genes} genes × {expr_data.n_samples} samples")
        
        st.write(f"{clinical_status} **Clinical Data**")
        if has_clinical:
            clinical_data = self.session.get_data("clinical_data")
            st.write(f"   {clinical_data.n_samples} samples × {clinical_data.n_variables} variables")
        
        st.write(f"{sample_status} **Sample Annotations**")
        if has_samples:
            sample_data = self.session.get_data("sample_annotations")
            st.write(f"   {len(sample_data)} annotated samples")
        
        # Analysis readiness
        st.divider()
        if has_expression and has_clinical:
            st.success("🎉 Ready for analysis!")
        elif has_expression:
            st.warning("⚠️ Upload clinical data to enable full analysis")
        else:
            st.info("📁 Upload expression data to begin")
    
    def _render_data_preview(self):
        """Render data preview section"""
        st.subheader("🔍 Data Preview")
        
        # Tab selection for different data types
        preview_tabs = st.tabs(["Expression Data", "Clinical Data", "Sample Matching"])
        
        with preview_tabs[0]:
            if self.session.has_data("expression_data"):
                expr_data = self.session.get_data("expression_data")
                
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.write("**Expression Data Preview (first 10 genes × 10 samples):**")
                    preview_df = expr_data.data.iloc[:10, :10]
                    st.dataframe(preview_df, use_container_width=True)
                
                with col2:
                    st.write("**Expression Distribution:**")
                    # Simple histogram of mean expression
                    mean_expr = expr_data.data.mean(axis=1)
                    st.bar_chart(mean_expr.value_counts().sort_index())
        
        with preview_tabs[1]:
            if self.session.has_data("clinical_data"):
                clinical_data = self.session.get_data("clinical_data")
                
                st.write("**Clinical Data Preview:**")
                st.dataframe(clinical_data.data.head(10), use_container_width=True)
                
                # Show missing data summary
                missing_data = clinical_data.data.isnull().sum()
                if missing_data.sum() > 0:
                    st.write("**Missing Data Summary:**")
                    missing_df = pd.DataFrame({
                        'Variable': missing_data.index,
                        'Missing Count': missing_data.values,
                        'Missing %': (missing_data.values / len(clinical_data.data) * 100).round(1)
                    })
                    missing_df = missing_df[missing_df['Missing Count'] > 0]
                    st.dataframe(missing_df, use_container_width=True)
        
        with preview_tabs[2]:
            self._render_sample_matching_preview()
    
    def _render_sample_matching_preview(self):
        """Render sample matching preview"""
        if not (self.session.has_data("expression_data") and self.session.has_data("clinical_data")):
            st.info("Upload both expression and clinical data to see sample matching")
            return
        
        expr_data = self.session.get_data("expression_data")
        clinical_data = self.session.get_data("clinical_data")
        
        # Sample matching analysis
        expr_samples = set(expr_data.samples)
        clinical_samples = set(clinical_data.samples)
        
        common_samples = expr_samples & clinical_samples
        expr_only = expr_samples - clinical_samples
        clinical_only = clinical_samples - expr_samples
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Common Samples", len(common_samples))
            
        with col2:
            st.metric("Expression Only", len(expr_only))
            
        with col3:
            st.metric("Clinical Only", len(clinical_only))
        
        # Detailed sample lists
        if len(common_samples) > 0:
            with st.expander(f"Common Samples ({len(common_samples)})"):
                st.write(sorted(list(common_samples))[:20])  # Show first 20
                if len(common_samples) > 20:
                    st.write(f"... and {len(common_samples) - 20} more")
        
        if len(expr_only) > 0:
            with st.expander(f"Expression Only ({len(expr_only)})"):
                st.write(sorted(list(expr_only))[:10])
                
        if len(clinical_only) > 0:
            with st.expander(f"Clinical Only ({len(clinical_only)})"):
                st.write(sorted(list(clinical_only))[:10])
    
    def _render_sample_matching(self):
        """Render sample matching and filtering options"""
        if not (self.session.has_data("expression_data") and self.session.has_data("clinical_data")):
            return
        
        st.subheader("🔗 Sample Matching & Filtering")
        
        expr_data = self.session.get_data("expression_data")
        clinical_data = self.session.get_data("clinical_data")
        
        # Get common samples
        common_samples = list(set(expr_data.samples) & set(clinical_data.samples))
        
        if len(common_samples) == 0:
            st.error("No common samples found between expression and clinical data!")
            return
        
        # Sample filtering options
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Sample Filtering Options:**")
            
            use_common_only = st.checkbox(
                f"Use only common samples ({len(common_samples)} samples)",
                value=True,
                help="Filter to samples present in both datasets"
            )
            
            if use_common_only:
                # Additional filtering options
                min_expression = st.number_input(
                    "Minimum mean expression threshold",
                    min_value=0.0,
                    max_value=10.0,
                    value=0.1,
                    step=0.1,
                    help="Filter out genes with very low expression"
                )
                
                min_samples_expressed = st.slider(
                    "Minimum samples with expression > 0",
                    min_value=1,
                    max_value=len(common_samples),
                    value=max(1, len(common_samples) // 10),
                    help="Filter genes expressed in at least N samples"
                )
        
        with col2:
            if use_common_only:
                # Apply filtering and show results
                filtered_expr = expr_data.data[common_samples]
                
                # Gene filtering
                gene_means = filtered_expr.mean(axis=1)
                gene_nonzero = (filtered_expr > 0).sum(axis=1)
                
                genes_pass_mean = gene_means >= min_expression
                genes_pass_samples = gene_nonzero >= min_samples_expressed
                genes_keep = genes_pass_mean & genes_pass_samples
                
                st.write("**Filtering Results:**")
                st.write(f"• Original genes: {len(expr_data.genes)}")
                st.write(f"• After mean filter: {genes_pass_mean.sum()}")
                st.write(f"• After sample filter: {genes_keep.sum()}")
                st.write(f"• Final samples: {len(common_samples)}")
                
                if st.button("Apply Filtering", type="primary"):
                    # Apply filtering to data
                    filtered_expr_data = ExpressionData(
                        data=filtered_expr.loc[genes_keep],
                        samples=common_samples,
                        genes=filtered_expr.loc[genes_keep].index.tolist(),
                        metadata={
                            **expr_data.metadata,
                            "filtered": True,
                            "filter_params": {
                                "min_expression": min_expression,
                                "min_samples_expressed": min_samples_expressed
                            }
                        }
                    )
                    
                    filtered_clinical_data = ClinicalData(
                        data=clinical_data.data.loc[common_samples],
                        samples=common_samples,
                        variables=clinical_data.variables,
                        metadata={
                            **clinical_data.metadata,
                            "filtered": True
                        }
                    )
                    
                    # Update session data
                    self.session.store_data("expression_data", filtered_expr_data)
                    self.session.store_data("clinical_data", filtered_clinical_data)
                    
                    st.success("✅ Data filtering applied successfully!")
                    st.rerun()
    
    def _render_sample_annotation_section(self):
        """Render the sample annotation section for flexible clinical data creation"""
        st.header("🏷️ Sample Annotation")
        st.markdown("**New Feature!** Create experimental groups without uploading clinical data files.")
        
        # Check if we already have clinical data
        has_clinical = self.session.has_data("clinical_data")
        expression_data = self.session.get_data("expression_data")
        
        if has_clinical:
            st.info("💡 You already have clinical data uploaded. You can still use sample annotation to create additional groupings or override existing ones.")
        
        # Create tabs for different annotation workflows
        tab1, tab2 = st.tabs(["🎯 Quick Annotation", "🔧 Advanced Annotation"])
        
        with tab1:
            self._render_quick_annotation(expression_data)
        
        with tab2:
            self._render_advanced_annotation(expression_data)
    
    def _render_quick_annotation(self, expression_data):
        """Render quick annotation interface for simple experimental designs"""
        st.subheader("⚡ Quick Setup")
        st.write("Perfect for simple experimental designs like case-control or treatment-control studies.")
        
        # Quick setup options
        col1, col2 = st.columns(2)
        
        with col1:
            design_type = st.selectbox(
                "Experimental design:",
                ["Case vs Control", "Treatment vs Control", "Before vs After", "Custom Groups"],
                key="quick_design_type"
            )
            
            if design_type == "Custom Groups":
                group1_name = st.text_input("Group 1 name:", "Group_A", key="group1_name")
                group2_name = st.text_input("Group 2 name:", "Group_B", key="group2_name")
            else:
                group_names = {
                    "Case vs Control": ("Case", "Control"),
                    "Treatment vs Control": ("Treatment", "Control"), 
                    "Before vs After": ("Before", "After")
                }
                group1_name, group2_name = group_names[design_type]
        
        with col2:
            assignment_method = st.selectbox(
                "Assignment method:",
                ["Alternating (A,B,A,B...)", "First half vs Second half", "Auto-detect from names"],
                key="assignment_method"
            )
            
            # Show preview
            n_samples = len(expression_data.samples)
            if assignment_method == "Alternating (A,B,A,B...)":
                group1_count = (n_samples + 1) // 2
                group2_count = n_samples // 2
            elif assignment_method == "First half vs Second half":
                group1_count = n_samples // 2
                group2_count = n_samples - group1_count
            else:  # Auto-detect
                group1_count = "?"
                group2_count = "?"
            
            st.info(f"📊 **Preview:**\n- {group1_name}: {group1_count} samples\n- {group2_name}: {group2_count} samples")
        
        # Apply quick annotation
        if st.button("✨ Create Groups", type="primary", key="create_quick_groups"):
            try:
                # Create sample annotator
                annotator = create_sample_annotator(expression_data.data)
                
                # Apply the selected method
                if assignment_method == "Auto-detect from names":
                    suggestions = annotator.suggest_groupings()
                    if suggestions:
                        # Use the best suggestion
                        best_suggestion = max(suggestions.values(), key=lambda x: x.get('confidence', 0))
                        annotator.assign_samples(best_suggestion['assignments'])
                        annotator.annotation_method = "auto_detected"
                    else:
                        # Fall back to alternating
                        assignments = {}
                        for i, sample in enumerate(annotator.get_samples()):
                            group = group1_name if i % 2 == 0 else group2_name
                            assignments[sample] = group
                        annotator.assign_samples(assignments)
                        annotator.annotation_method = "alternating_fallback"
                
                elif assignment_method == "Alternating (A,B,A,B...)":
                    assignments = {}
                    for i, sample in enumerate(annotator.get_samples()):
                        group = group1_name if i % 2 == 0 else group2_name
                        assignments[sample] = group
                    annotator.assign_samples(assignments)
                    annotator.annotation_method = "alternating"
                
                else:  # First half vs Second half
                    mid_point = n_samples // 2
                    assignments = {}
                    for i, sample in enumerate(annotator.get_samples()):
                        group = group1_name if i < mid_point else group2_name
                        assignments[sample] = group
                    annotator.assign_samples(assignments)
                    annotator.annotation_method = "split_half"
                
                # Convert to clinical data and store
                clinical_data = annotator.export_to_clinical_data("Condition")
                self.session.store_data("clinical_data", clinical_data)
                
                st.success(f"✅ Created sample groups successfully!")
                st.success(f"📊 Groups: {annotator.get_group_sizes()}")
                
                # Show preview
                with st.expander("Preview Created Clinical Data"):
                    st.dataframe(clinical_data.data.head(10))
                
                st.rerun()
                
            except Exception as e:
                st.error(f"Error creating groups: {str(e)}")
                logger.error(f"Quick annotation failed: {e}")
    
    def _render_advanced_annotation(self, expression_data):
        """Render advanced annotation interface using the full widget"""
        st.subheader("🔧 Advanced Sample Annotation")
        st.write("Full-featured annotation system with pattern detection, templates, and manual assignment.")
        
        # Use the full sample annotation widget
        try:
            clinical_data = render_sample_annotation_widget(
                expression_data.data, 
                key="advanced_annotation"
            )
            
            if clinical_data:
                # Store the created clinical data
                self.session.store_data("clinical_data", clinical_data)
                st.success("🎉 Sample annotation completed and ready for analysis!")
                
                # Auto-advance to next step
                if st.button("➡️ Continue to DE Analysis", type="primary"):
                    st.switch_page("pages/deseq2_analysis.py")
                
        except Exception as e:
            st.error(f"Error in advanced annotation: {str(e)}")
            logger.error(f"Advanced annotation failed: {e}")


def render_data_import_tab():
    """Render the data import tab (main entry point)"""
    tab = DataImportTab()
    tab.render()