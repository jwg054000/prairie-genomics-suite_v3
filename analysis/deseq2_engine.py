"""
Prairie Genomics Suite - DESeq2 Analysis Engine

This module provides comprehensive DESeq2 differential expression analysis
using R integration (rpy2) with Python fallbacks. It's designed to be
independent and testable, similar to how you'd organize functions in an
R package for genomics analysis.

Features:
- Full DESeq2 workflow with R integration
- Python fallbacks for when R is unavailable
- Publication-quality visualizations
- Batch correction support
- Multiple contrast analysis
- Comprehensive result export

Author: Prairie Genomics Team
"""

import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import warnings
import logging

# Import core modules
from core.data_models import ExpressionData, ClinicalData, DESeq2Results
from core.utils import validate_dataframe, filter_expression_data, handle_exception, clean_clinical_data, create_minimal_clinical_data
from config import get_config, is_feature_enabled

# Set up logging
logger = logging.getLogger(__name__)

# Lazy import for optional dependencies
def get_plotting_libs():
    """Lazy import plotting libraries"""
    try:
        import plotly.express as px
        import plotly.graph_objects as go
        return px, go
    except ImportError:
        return None, None


class DESeq2Engine:
    """
    DESeq2 Analysis Engine with R integration and Python fallbacks
    
    This class provides a clean interface to DESeq2 analysis, handling
    both R integration (preferred) and Python fallbacks. Think of it
    as your main DESeq2 analysis function from R, but with robust
    error handling and fallbacks.
    """
    
    def __init__(self):
        """Initialize the DESeq2 engine"""
        self.r_available = False
        self.r = None
        self.ro = None
        self.config = get_config("analysis")["deseq2"]
        
        # Initialize R integration if available
        if is_feature_enabled("enable_r_integration"):
            self._initialize_r()
    
    def _initialize_r(self):
        """Initialize R environment and load required packages"""
        try:
            import rpy2
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri, numpy2ri
            
            # Activate automatic conversion
            pandas2ri.activate()
            numpy2ri.activate()
            
            self.r = ro.r
            self.ro = ro
            self.pandas2ri = pandas2ri
            
            # Load R script templates
            script_path = Path(__file__).parent.parent / "deseq2_templates.R"
            if script_path.exists():
                self.r.source(str(script_path))
                self.r_available = True
                logger.info("DESeq2 R integration initialized successfully")
            else:
                logger.warning("R script templates not found. Some features may be limited.")
                
        except ImportError as e:
            logger.warning(f"rpy2 not available: {e}. DESeq2 analysis will use Python fallback.")
        except Exception as e:
            logger.error(f"R integration failed: {e}")
    
    def is_available(self) -> bool:
        """Check if R integration is available"""
        return self.r_available
    
    def get_status(self) -> Dict[str, Any]:
        """Get status information about the engine"""
        return {
            "r_available": self.r_available,
            "method": "DESeq2 (R)" if self.r_available else "Python fallback",
            "config": self.config
        }
    
    @handle_exception
    def prepare_data(self, expression_data: ExpressionData, 
                    clinical_data: ClinicalData,
                    condition_column: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepare expression and clinical data for DESeq2 analysis with enhanced debugging
        
        Args:
            expression_data: Expression data container
            clinical_data: Clinical data container
            condition_column: Column name for experimental conditions
            
        Returns:
            Tuple of (prepared_expression, prepared_clinical)
        """
        logger.info("🔍 Starting enhanced data preparation with comprehensive debugging")
        
        # Debug input data
        logger.info(f"📊 Input expression data: {expression_data.data.shape} ({type(expression_data)})")
        logger.info(f"📊 Input clinical data: {clinical_data.data.shape} ({type(clinical_data)})")
        logger.info(f"🎯 Condition column: '{condition_column}'")
        
        # Extract and log sample information
        expr_samples = expression_data.samples if hasattr(expression_data, 'samples') else list(expression_data.data.columns)
        clinical_samples = clinical_data.samples if hasattr(clinical_data, 'samples') else list(clinical_data.data.index)
        
        logger.info(f"🧬 Expression samples ({len(expr_samples)}): {expr_samples[:5]}{'...' if len(expr_samples) > 5 else ''}")
        logger.info(f"🏥 Clinical samples ({len(clinical_samples)}): {clinical_samples[:5]}{'...' if len(clinical_samples) > 5 else ''}")
        logger.info(f"🔍 Expression sample types: {[type(s).__name__ for s in expr_samples[:3]]}")
        logger.info(f"🔍 Clinical sample types: {[type(s).__name__ for s in clinical_samples[:3]]}")
        
        # Validate inputs
        expr_valid, expr_errors = validate_dataframe(
            expression_data.data, 
            min_rows=100, 
            min_cols=3,
            name="Expression data"
        )
        if not expr_valid:
            logger.error(f"❌ Expression data validation failed: {expr_errors}")
            raise ValueError(f"Expression data validation failed: {expr_errors}")
        
        # Clean clinical data first to handle NaN columns
        original_clinical_shape = clinical_data.data.shape
        cleaned_clinical = clean_clinical_data(
            clinical_data.data,
            remove_nan_columns=True,
            min_valid_fraction=0.1
        )
        
        if cleaned_clinical.shape != original_clinical_shape:
            logger.info(f"🧹 Clinical data cleaned: {original_clinical_shape} → {cleaned_clinical.shape}")
        
        # Use flexible validation for clinical data
        clinical_valid, clinical_errors = validate_dataframe(
            cleaned_clinical,
            min_rows=3,
            required_columns=[condition_column],
            name="Clinical data",
            allow_nan_columns=True,
            strict_validation=False
        )
        
        # If clinical data validation fails, try to create minimal clinical data
        if not clinical_valid:
            logger.warning(f"⚠️ Clinical data validation failed: {clinical_errors}")
            logger.info("🔧 Attempting to create minimal clinical data from expression samples")
            try:
                cleaned_clinical = create_minimal_clinical_data(
                    expression_data.data,
                    condition_column=condition_column,
                    default_conditions=["Control", "Treatment"]
                )
                logger.info("✅ Successfully created minimal clinical data")
            except Exception as e:
                logger.error(f"❌ Failed to create minimal clinical data: {str(e)}")
                raise ValueError(f"Clinical data validation failed and could not create minimal data: {clinical_errors}. Error: {str(e)}")
        
        # Prepare count matrix
        logger.info("🔢 Preparing count matrix...")
        count_matrix = self._prepare_count_matrix(expression_data.data)
        logger.info(f"📊 Count matrix prepared: {count_matrix.shape}")
        
        # Prepare metadata using cleaned clinical data
        logger.info("📋 Preparing metadata...")
        metadata = self._prepare_metadata(cleaned_clinical, condition_column)
        logger.info(f"📊 Metadata prepared: {metadata.shape}")
        
        # Enhanced sample matching with detailed debugging
        logger.info("🔗 Starting enhanced sample matching...")
        count_matrix_samples = list(count_matrix.columns)
        metadata_samples = list(metadata.index)
        
        logger.info(f"🧬 Count matrix samples ({len(count_matrix_samples)}): {count_matrix_samples[:5]}{'...' if len(count_matrix_samples) > 5 else ''}")
        logger.info(f"📋 Metadata samples ({len(metadata_samples)}): {metadata_samples[:5]}{'...' if len(metadata_samples) > 5 else ''}")
        
        # Convert to sets for overlap analysis
        count_set = set(count_matrix_samples)
        metadata_set = set(metadata_samples)
        
        # Find overlaps and mismatches
        common_samples = list(count_set & metadata_set)
        expr_only = list(count_set - metadata_set)
        clinical_only = list(metadata_set - count_set)
        
        logger.info(f"✅ Common samples ({len(common_samples)}): {common_samples[:5]}{'...' if len(common_samples) > 5 else ''}")
        
        if expr_only:
            logger.warning(f"⚠️ Expression-only samples ({len(expr_only)}): {expr_only[:5]}{'...' if len(expr_only) > 5 else ''}")
        
        if clinical_only:
            logger.warning(f"⚠️ Clinical-only samples ({len(clinical_only)}): {clinical_only[:5]}{'...' if len(clinical_only) > 5 else ''}")
        
        # Enhanced error handling for insufficient overlap
        if len(common_samples) < 3:
            logger.error("❌ Insufficient overlapping samples for analysis")
            
            # Provide detailed diagnostic information
            logger.error("🔍 DETAILED SAMPLE MISMATCH ANALYSIS:")
            logger.error(f"   Expression samples: {count_matrix_samples}")
            logger.error(f"   Clinical samples: {metadata_samples}")
            
            # Check for common issues and apply fixes
            if len(count_matrix_samples) == len(metadata_samples):
                logger.error("   ℹ️ Same number of samples - likely a naming or ordering issue")
                
                # Check for exact matches after stripping whitespace
                stripped_expr = [str(s).strip() for s in count_matrix_samples]
                stripped_clinical = [str(s).strip() for s in metadata_samples]
                
                if set(stripped_expr) == set(stripped_clinical):
                    logger.error("   🔧 SOLUTION: Whitespace issue detected - applying fix")
                    
                    # Create mapping to fix whitespace issues
                    mapping = dict(zip(metadata_samples, stripped_clinical))
                    metadata.index = [mapping.get(idx, idx) for idx in metadata.index]
                    metadata_samples = list(metadata.index)
                    common_samples = list(set(count_matrix_samples) & set(metadata_samples))
                    
                    logger.info(f"   ✅ Fixed whitespace issue - new overlap: {len(common_samples)} samples")
                
                # Check for case sensitivity issues
                elif set(s.lower() for s in count_matrix_samples) == set(s.lower() for s in metadata_samples):
                    logger.error("   🔧 SOLUTION: Case sensitivity issue detected - applying fix")
                    
                    # Create case-insensitive mapping
                    expr_lower_to_orig = {s.lower(): s for s in count_matrix_samples}
                    new_index = []
                    for clinical_sample in metadata_samples:
                        matching_expr = expr_lower_to_orig.get(clinical_sample.lower())
                        if matching_expr:
                            new_index.append(matching_expr)
                        else:
                            new_index.append(clinical_sample)
                    
                    metadata.index = new_index
                    common_samples = list(set(count_matrix_samples) & set(new_index))
                    
                    logger.info(f"   ✅ Fixed case sensitivity issue - new overlap: {len(common_samples)} samples")
            
            # If still insufficient samples after attempted fixes
            if len(common_samples) < 3:
                error_msg = (
                    f"Insufficient overlapping samples: {len(common_samples)} found, need at least 3.\n"
                    f"Expression samples ({len(count_matrix_samples)}): {count_matrix_samples}\n"
                    f"Clinical samples ({len(metadata_samples)}): {metadata_samples}\n"
                    f"Common samples: {common_samples}\n"
                    f"Expression-only: {expr_only}\n"
                    f"Clinical-only: {clinical_only}\n\n"
                    f"SUGGESTED SOLUTIONS:\n"
                    f"1. Check sample naming - ensure identical names in both datasets\n"
                    f"2. Remove extra whitespace or fix capitalization\n"
                    f"3. Use the Sample Annotation tool to create matching clinical data\n"
                    f"4. Verify that both datasets contain the same samples"
                )
                raise ValueError(error_msg)
        
        # Subset to common samples
        logger.info(f"✂️ Subsetting data to {len(common_samples)} common samples...")
        count_matrix = count_matrix[common_samples]
        metadata = metadata.loc[common_samples]
        
        # Validate final data
        logger.info("✅ Final validation...")
        logger.info(f"📊 Final count matrix: {count_matrix.shape}")
        logger.info(f"📊 Final metadata: {metadata.shape}")
        logger.info(f"🎯 Condition distribution: {metadata[condition_column].value_counts().to_dict()}")
        
        # Check for minimum samples per condition
        condition_counts = metadata[condition_column].value_counts()
        min_samples_per_condition = condition_counts.min()
        
        if min_samples_per_condition < 3:
            logger.warning(f"⚠️ Some conditions have fewer than 3 samples: {condition_counts.to_dict()}")
            logger.warning("   This may affect statistical power but analysis will proceed")
        
        logger.info(f"🎉 Data preparation completed successfully: {count_matrix.shape[0]} genes, {len(common_samples)} samples")
        return count_matrix, metadata
    
    def _prepare_count_matrix(self, expression_data: pd.DataFrame) -> pd.DataFrame:
        """Prepare count matrix for DESeq2 (ensure integer counts)"""
        # Round to integers (DESeq2 requires count data)
        count_matrix = expression_data.round().astype(int)
        
        # Remove genes with all zeros
        count_matrix = count_matrix.loc[count_matrix.sum(axis=1) > 0]
        
        # Filter low-expression genes
        min_counts = self.config.get("min_counts", 10)
        min_samples = self.config.get("min_samples", 7)
        
        count_matrix = filter_expression_data(
            count_matrix,
            min_expression=min_counts,
            min_samples=min(min_samples, count_matrix.shape[1] // 2)
        )
        
        return count_matrix
    
    def _prepare_metadata(self, clinical_data: pd.DataFrame, condition_column: str) -> pd.DataFrame:
        """Prepare metadata for DESeq2 analysis"""
        metadata = clinical_data.copy()
        
        # Ensure condition column exists
        if condition_column not in metadata.columns:
            raise ValueError(f"Condition column '{condition_column}' not found in metadata")
        
        # Ensure factor levels are valid R names
        metadata[condition_column] = metadata[condition_column].astype(str)
        metadata[condition_column] = metadata[condition_column].str.replace(
            "[^A-Za-z0-9_]", "_", regex=True
        )
        
        # Remove samples with missing conditions
        metadata = metadata.dropna(subset=[condition_column])
        
        return metadata
    
    @handle_exception
    def run_analysis(self, expression_data: ExpressionData,
                    clinical_data: ClinicalData,
                    condition_column: str = "Condition",
                    reference_condition: Optional[str] = None,
                    batch_correction: bool = False,
                    batch_column: str = "Batch",
                    fc_cutoff: float = None,
                    p_cutoff: float = None) -> DESeq2Results:
        """
        Run complete DESeq2 differential expression analysis
        
        Args:
            expression_data: Expression data container
            clinical_data: Clinical data container  
            condition_column: Column defining experimental conditions
            reference_condition: Reference condition for comparisons
            batch_correction: Whether to apply batch correction
            batch_column: Column defining batches
            fc_cutoff: Fold change cutoff for significance
            p_cutoff: P-value cutoff for significance
            
        Returns:
            DESeq2Results container with all analysis results
        """
        # Use config defaults if not specified
        if fc_cutoff is None:
            fc_cutoff = self.config.get("fc_cutoff", 1.0)
        if p_cutoff is None:
            p_cutoff = self.config.get("p_cutoff", 0.05)
        
        # Prepare data
        count_matrix, metadata = self.prepare_data(expression_data, clinical_data, condition_column)
        
        # Run analysis (R or Python fallback)
        if self.is_available():
            results = self._run_r_analysis(
                count_matrix, metadata, condition_column, 
                reference_condition, batch_correction, batch_column
            )
        else:
            results = self._run_python_analysis(
                count_matrix, metadata, condition_column
            )
        
        # Create results container
        deseq2_results = DESeq2Results(
            results=results["results"],
            contrasts=list(results["results"].keys()),
            normalized_counts=results.get("normalized_counts"),
            metadata={
                "method": results["method"],
                "condition_column": condition_column,
                "reference_condition": reference_condition,
                "batch_correction": batch_correction,
                "batch_column": batch_column,
                "fc_cutoff": fc_cutoff,
                "p_cutoff": p_cutoff,
                "n_genes": count_matrix.shape[0],
                "n_samples": count_matrix.shape[1]
            }
        )
        
        return deseq2_results
    
    def _run_r_analysis(self, count_matrix: pd.DataFrame, 
                       metadata: pd.DataFrame,
                       condition_column: str,
                       reference_condition: Optional[str],
                       batch_correction: bool,
                       batch_column: str) -> Dict[str, Any]:
        """Run DESeq2 analysis using R"""
        try:
            # Convert to R objects
            r_count_matrix = self.pandas2ri.py2rpy(count_matrix)
            r_metadata = self.pandas2ri.py2rpy(metadata)
            
            # Set up R environment
            self.ro.globalenv['count_matrix'] = r_count_matrix
            self.ro.globalenv['metadata'] = r_metadata
            self.ro.globalenv['condition_column'] = condition_column
            self.ro.globalenv['reference_condition'] = reference_condition
            self.ro.globalenv['batch_correction'] = batch_correction
            self.ro.globalenv['batch_column'] = batch_column
            
            # Run DESeq2 workflow
            r_script = f"""
            # Run complete DESeq2 workflow
            design_formula <- if (batch_correction && '{batch_column}' %in% colnames(metadata)) {{
                paste("~", '{batch_column}', "+", '{condition_column}')
            }} else {{
                paste("~", '{condition_column}')
            }}
            
            workflow_results <- run_complete_deseq2_workflow(
                count_matrix = count_matrix,
                metadata = metadata,
                reference_condition = reference_condition,
                apply_batch_correction = batch_correction,
                batch_column = '{batch_column}',
                create_plots = TRUE
            )
            """
            
            self.r(r_script)
            
            # Extract results
            results = {}
            
            # Get DESeq2 results for each contrast
            r_results = self.r('workflow_results$results')
            for contrast_name in r_results.names:
                contrast_df = self.pandas2ri.rpy2py(r_results.rx2(contrast_name))
                # Ensure gene column is properly named
                if contrast_df.index.name != 'gene':
                    contrast_df = contrast_df.reset_index()
                    if 'index' in contrast_df.columns:
                        contrast_df = contrast_df.rename(columns={'index': 'gene'})
                results[contrast_name] = contrast_df
            
            # Get normalized counts
            try:
                normalized_counts = self.pandas2ri.rpy2py(self.r('workflow_results$dds'))
                normalized_counts = self.pandas2ri.rpy2py(self.r('counts(workflow_results$dds, normalized=TRUE)'))
            except:
                normalized_counts = None
            
            return {
                'results': results,
                'normalized_counts': normalized_counts,
                'method': 'DESeq2 (R)'
            }
            
        except Exception as e:
            logger.error(f"DESeq2 R analysis failed: {e}")
            raise
    
    def _run_python_analysis(self, count_matrix: pd.DataFrame,
                           metadata: pd.DataFrame,
                           condition_column: str) -> Dict[str, Any]:
        """Fallback analysis using Python statistical tests"""
        logger.info("Using Python fallback for differential expression analysis")
        
        # Import stats analyzer
        from analysis.stats_analyzer import BasicStatsAnalyzer
        analyzer = BasicStatsAnalyzer()
        
        # Get unique conditions
        conditions = metadata[condition_column].unique()
        results = {}
        
        # Perform pairwise comparisons
        for i, cond1 in enumerate(conditions):
            for cond2 in conditions[i+1:]:
                contrast_name = f"{cond2}_vs_{cond1}"
                
                # Get sample groups
                group1_samples = metadata[metadata[condition_column] == cond1].index.tolist()
                group2_samples = metadata[metadata[condition_column] == cond2].index.tolist()
                
                # Filter to available samples
                group1_samples = [s for s in group1_samples if s in count_matrix.columns]
                group2_samples = [s for s in group2_samples if s in count_matrix.columns]
                
                if len(group1_samples) >= 2 and len(group2_samples) >= 2:
                    # Run statistical analysis
                    contrast_results = analyzer.ttest_analysis(
                        count_matrix, group1_samples, group2_samples
                    )
                    results[contrast_name] = contrast_results
        
        return {
            'results': results,
            'normalized_counts': None,
            'method': 'T-test (Python fallback)'
        }
    
    @handle_exception
    def create_volcano_plot(self, results_df: pd.DataFrame, 
                          title: str = "Volcano Plot",
                          fc_cutoff: float = None,
                          p_cutoff: float = None) -> Any:
        """
        Create volcano plot from DESeq2 results
        
        Args:
            results_df: DESeq2 results DataFrame
            title: Plot title
            fc_cutoff: Fold change cutoff for significance
            p_cutoff: P-value cutoff for significance
            
        Returns:
            Plotly figure or R ggplot object
        """
        # Use config defaults if not specified
        if fc_cutoff is None:
            fc_cutoff = self.config.get("fc_cutoff", 1.0)
        if p_cutoff is None:
            p_cutoff = self.config.get("p_cutoff", 0.05)
        
        # Try R first, then Python fallback
        if self.is_available():
            return self._create_r_volcano_plot(results_df, title, fc_cutoff, p_cutoff)
        else:
            return self._create_python_volcano_plot(results_df, title, fc_cutoff, p_cutoff)
    
    def _create_r_volcano_plot(self, results_df: pd.DataFrame, title: str,
                              fc_cutoff: float, p_cutoff: float) -> Any:
        """Create volcano plot using R ggplot2"""
        try:
            # Convert to R and create plot
            r_df = self.pandas2ri.py2rpy(results_df)
            self.ro.globalenv['results_df'] = r_df
            self.ro.globalenv['plot_title'] = title
            self.ro.globalenv['fc_cutoff'] = fc_cutoff
            self.ro.globalenv['p_cutoff'] = p_cutoff
            
            volcano_plot = self.r('''
                create_enhanced_volcano_plot(
                    results_df, 
                    title = plot_title,
                    fc_cutoff = fc_cutoff,
                    p_cutoff = p_cutoff
                )
            ''')
            
            return volcano_plot
            
        except Exception as e:
            logger.warning(f"R volcano plot failed, using Python fallback: {e}")
            return self._create_python_volcano_plot(results_df, title, fc_cutoff, p_cutoff)
    
    def _create_python_volcano_plot(self, results_df: pd.DataFrame, title: str,
                                  fc_cutoff: float, p_cutoff: float):
        """Create volcano plot using Python/Plotly"""
        px, go = get_plotting_libs()
        if px is None:
            logger.error("Plotly not available for volcano plot")
            return None
        
        df = results_df.copy()
        
        # Ensure required columns exist
        if 'log2FoldChange' not in df.columns or 'padj' not in df.columns:
            logger.error("Required columns (log2FoldChange, padj) not found in results")
            return None
        
        # Add significance categories
        df['significance'] = 'Not Significant'
        df.loc[(df['log2FoldChange'] > fc_cutoff) & (df['padj'] < p_cutoff), 'significance'] = 'Up-regulated'
        df.loc[(df['log2FoldChange'] < -fc_cutoff) & (df['padj'] < p_cutoff), 'significance'] = 'Down-regulated'
        
        # Handle zeros and infinities in p-values
        df['neg_log10_padj'] = -np.log10(df['padj'].replace(0, 1e-300))
        
        # Create plot
        fig = px.scatter(
            df, 
            x='log2FoldChange', 
            y='neg_log10_padj',
            color='significance',
            title=title,
            hover_data=['gene'] if 'gene' in df.columns else None,
            color_discrete_map={
                'Up-regulated': '#d62728',
                'Down-regulated': '#1f77b4', 
                'Not Significant': '#808080'
            }
        )
        
        # Add cutoff lines
        fig.add_hline(y=-np.log10(p_cutoff), line_dash="dash", line_color="gray")
        fig.add_vline(x=fc_cutoff, line_dash="dash", line_color="gray") 
        fig.add_vline(x=-fc_cutoff, line_dash="dash", line_color="gray")
        
        fig.update_layout(
            xaxis_title="log₂(Fold Change)",
            yaxis_title="-log₁₀(Adjusted P-value)",
            showlegend=True
        )
        
        return fig
    
    def get_significant_genes(self, results: DESeq2Results,
                            contrast: str,
                            fc_cutoff: float = None,
                            p_cutoff: float = None) -> List[str]:
        """
        Get significantly differentially expressed genes from results
        
        Args:
            results: DESeq2Results container
            contrast: Name of the contrast
            fc_cutoff: Fold change cutoff
            p_cutoff: P-value cutoff
            
        Returns:
            List of significant gene names
        """
        # Use config defaults if not specified
        if fc_cutoff is None:
            fc_cutoff = self.config.get("fc_cutoff", 1.0)
        if p_cutoff is None:
            p_cutoff = self.config.get("p_cutoff", 0.05)
        
        return results.get_significant_genes(contrast, fc_cutoff, p_cutoff)
    
    def export_results(self, results: DESeq2Results, 
                      output_dir: str = ".") -> Dict[str, str]:
        """
        Export DESeq2 results to files
        
        Args:
            results: DESeq2Results container
            output_dir: Output directory
            
        Returns:
            Dictionary of exported file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        exported_files = {}
        
        # Export each contrast
        for contrast_name, contrast_results in results.results.items():
            filename = output_path / f"{contrast_name}_DESeq2_results.csv"
            contrast_results.to_csv(filename, index=False)
            exported_files[contrast_name] = str(filename)
        
        # Export normalized counts if available
        if results.normalized_counts is not None:
            norm_filename = output_path / "normalized_counts.csv"
            results.normalized_counts.to_csv(norm_filename)
            exported_files["normalized_counts"] = str(norm_filename)
        
        # Export metadata
        metadata_filename = output_path / "analysis_metadata.json"
        import json
        with open(metadata_filename, 'w') as f:
            json.dump(results.metadata, f, indent=2, default=str)
        exported_files["metadata"] = str(metadata_filename)
        
        logger.info(f"Exported {len(exported_files)} files to {output_dir}")
        return exported_files


# Create a singleton instance for easy import
_deseq2_engine = None

def get_deseq2_engine() -> DESeq2Engine:
    """Get the global DESeq2 engine instance"""
    global _deseq2_engine
    if _deseq2_engine is None:
        _deseq2_engine = DESeq2Engine()
    return _deseq2_engine