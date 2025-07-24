"""
Prairie Genomics Suite - Survival Analysis Engine

This module provides comprehensive survival analysis capabilities with both
R integration (survminer, survival packages) and Python fallbacks (lifelines).
It's designed for molecular survival analysis with gene expression integration.

Features:
- Kaplan-Meier survival curves with R survminer
- Cox proportional hazards modeling
- Gene expression-based stratification (median, tertiles, quartiles)
- Risk tables and confidence intervals
- Python fallbacks with lifelines
- Publication-quality visualizations

Author: Prairie Genomics Team
"""

import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
import warnings
import logging

# Import core modules
from core.data_models import SurvivalData, ExpressionData, ClinicalData
from core.utils import validate_dataframe, safe_convert_numeric, handle_exception
from config import get_config, is_feature_enabled

# Set up logging
logger = logging.getLogger(__name__)

# Lazy import for optional dependencies
def get_plotting_libs():
    """Lazy import plotting libraries"""
    try:
        import plotly.express as px
        import plotly.graph_objects as go
        import matplotlib.pyplot as plt
        return px, go, plt
    except ImportError:
        return None, None, None


class SurvivalEngine:
    """
    Comprehensive Survival Analysis Engine
    
    Provides both R-based (preferred) and Python-based survival analysis
    with molecular integration capabilities. Similar to survival package
    in R but with enhanced genomics features.
    """
    
    def __init__(self):
        """Initialize the survival analysis engine"""
        self.r_available = False
        self.lifelines_available = False
        self.r = None
        self.ro = None
        self.config = get_config("analysis")["survival"]
        
        # Initialize R integration if available
        if is_feature_enabled("enable_r_integration"):
            self._initialize_r()
        
        # Check for lifelines (Python survival analysis)
        self._check_lifelines()
    
    def _initialize_r(self):
        """Initialize R environment and load survival packages"""
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
            
            # Load R script templates (includes survival functions)
            script_path = Path(__file__).parent.parent / "deseq2_templates.R"
            if script_path.exists():
                self.r.source(str(script_path))
                # Load visualization packages which include survival packages
                try:
                    self.r('load_visualization_packages()')
                    self.r_available = True
                    logger.info("Survival R integration initialized successfully")
                except Exception as e:
                    logger.warning(f"Some R survival packages may not be available: {e}")
            else:
                logger.warning("R script templates not found. R survival analysis unavailable.")
                
        except ImportError as e:
            logger.warning(f"rpy2 not available: {e}. Survival analysis will use Python fallback.")
        except Exception as e:
            logger.error(f"R integration failed: {e}")
    
    def _check_lifelines(self):
        """Check if lifelines is available for Python survival analysis"""
        try:
            import lifelines
            self.lifelines_available = True
            logger.info("Lifelines available for Python survival analysis")
        except ImportError:
            logger.warning("Lifelines not available. Install with: pip install lifelines")
    
    def is_r_available(self) -> bool:
        """Check if R survival analysis is available"""
        return self.r_available
    
    def is_python_available(self) -> bool:
        """Check if Python survival analysis is available"""
        return self.lifelines_available
    
    def get_status(self) -> Dict[str, Any]:
        """Get status information about the engine"""
        return {
            "r_available": self.r_available,
            "lifelines_available": self.lifelines_available,
            "preferred_method": "R (survminer)" if self.r_available else "Python (lifelines)" if self.lifelines_available else "None",
            "config": self.config
        }
    
    @handle_exception
    def prepare_survival_data(self, clinical_data: ClinicalData,
                            time_column: str,
                            event_column: str,
                            expression_data: Optional[ExpressionData] = None,
                            gene_of_interest: Optional[str] = None) -> SurvivalData:
        """
        Prepare survival data from clinical data with optional gene integration
        
        Args:
            clinical_data: Clinical data container
            time_column: Column name for survival times
            event_column: Column name for event indicators (0/1)
            expression_data: Optional expression data for stratification
            gene_of_interest: Optional gene for expression-based stratification
            
        Returns:
            SurvivalData container
        """
        # Validate clinical data
        clinical_valid, clinical_errors = validate_dataframe(
            clinical_data.data,
            min_rows=10,  # Need reasonable sample size for survival
            required_columns=[time_column, event_column],
            name="Clinical data"
        )
        if not clinical_valid:
            raise ValueError(f"Clinical data validation failed: {clinical_errors}")
        
        # Extract survival data
        survival_df = clinical_data.data[[time_column, event_column]].copy()
        
        # Convert and validate survival times and events
        survival_df[time_column] = safe_convert_numeric(survival_df[time_column])
        survival_df[event_column] = safe_convert_numeric(survival_df[event_column])
        
        # Remove invalid entries
        survival_df = survival_df.dropna()
        survival_df = survival_df[survival_df[time_column] > 0]  # Positive times only
        survival_df = survival_df[survival_df[event_column].isin([0, 1])]  # Valid events only
        
        if len(survival_df) < 10:
            raise ValueError(f"Insufficient valid survival data: {len(survival_df)} samples")
        
        # Prepare covariates
        covariates = None
        if len(clinical_data.data.columns) > 2:
            # Include other clinical variables as potential covariates
            covariate_columns = [col for col in clinical_data.data.columns 
                               if col not in [time_column, event_column]]
            if covariate_columns:
                covariates = clinical_data.data[covariate_columns].loc[survival_df.index]
        
        # Create SurvivalData container
        survival_data = SurvivalData(
            time=survival_df[time_column].values,
            event=survival_df[event_column].values,
            samples=survival_df.index.tolist(),
            covariates=covariates
        )
        
        logger.info(f"Prepared survival data: {len(survival_data.samples)} samples, "
                   f"{survival_data.event.sum()} events")
        
        return survival_data
    
    @handle_exception
    def kaplan_meier_analysis(self, survival_data: SurvivalData,
                            group_column: Optional[str] = None,
                            expression_data: Optional[ExpressionData] = None,
                            gene_of_interest: Optional[str] = None,
                            split_method: str = "median",
                            confidence_interval: float = None) -> Dict[str, Any]:
        """
        Perform Kaplan-Meier survival analysis
        
        Args:
            survival_data: Survival data container
            group_column: Column for grouping (in covariates)
            expression_data: Optional expression data for gene-based grouping
            gene_of_interest: Gene name for expression-based stratification
            split_method: Method for splitting by expression ('median', 'tertiles', 'quartiles')
            confidence_interval: Confidence level (uses config default if None)
            
        Returns:
            Dictionary with analysis results and plots
        """
        if confidence_interval is None:
            confidence_interval = self.config.get("confidence_interval", 0.95)
        
        # Try R analysis first, then Python fallback
        if self.is_r_available():
            return self._kaplan_meier_r(
                survival_data, group_column, expression_data, 
                gene_of_interest, split_method, confidence_interval
            )
        elif self.is_python_available():
            return self._kaplan_meier_python(
                survival_data, group_column, expression_data,
                gene_of_interest, split_method, confidence_interval
            )
        else:
            raise RuntimeError("No survival analysis backend available. Install lifelines or rpy2.")
    
    def _kaplan_meier_r(self, survival_data: SurvivalData,
                       group_column: Optional[str],
                       expression_data: Optional[ExpressionData],
                       gene_of_interest: Optional[str],
                       split_method: str,
                       confidence_interval: float) -> Dict[str, Any]:
        """Kaplan-Meier analysis using R survminer"""
        try:
            # Prepare survival DataFrame for R
            survival_df = pd.DataFrame({
                'time': survival_data.time,
                'event': survival_data.event
            }, index=survival_data.samples)
            
            # Add grouping information
            if gene_of_interest and expression_data:
                # Gene expression-based grouping
                if gene_of_interest in expression_data.genes:
                    gene_expr = expression_data.data.loc[gene_of_interest, survival_data.samples]
                    
                    if split_method == "median":
                        threshold = gene_expr.median()
                        survival_df['group'] = (gene_expr > threshold).map({True: 'High', False: 'Low'})
                    elif split_method == "tertiles":
                        tertiles = gene_expr.quantile([0.33, 0.67])
                        survival_df['group'] = pd.cut(
                            gene_expr, 
                            bins=[-np.inf, tertiles.iloc[0], tertiles.iloc[1], np.inf],
                            labels=['Low', 'Medium', 'High']
                        )
                    elif split_method == "quartiles":
                        quartiles = gene_expr.quantile([0.25, 0.5, 0.75])
                        survival_df['group'] = pd.cut(
                            gene_expr,
                            bins=[-np.inf] + quartiles.tolist() + [np.inf],
                            labels=['Q1', 'Q2', 'Q3', 'Q4']
                        )
                    
                    group_column = 'group'
                else:
                    logger.warning(f"Gene {gene_of_interest} not found in expression data")
            
            elif group_column and survival_data.covariates is not None:
                # Use existing clinical grouping
                if group_column in survival_data.covariates.columns:
                    survival_df['group'] = survival_data.covariates[group_column]
                    group_column = 'group'
                else:
                    logger.warning(f"Group column {group_column} not found")
                    group_column = None
            
            # Convert to R and run analysis
            r_survival_df = self.pandas2ri.py2rpy(survival_df)
            self.ro.globalenv['survival_data'] = r_survival_df
            self.ro.globalenv['time_column'] = 'time'
            self.ro.globalenv['event_column'] = 'event'
            self.ro.globalenv['group_column'] = group_column or 'group'
            self.ro.globalenv['split_method'] = split_method
            
            # Add expression data if available for advanced analysis
            if expression_data:
                # Match samples
                common_samples = list(set(survival_data.samples) & set(expression_data.samples))
                expr_subset = expression_data.data[common_samples]
                r_expression = self.pandas2ri.py2rpy(expr_subset)
                self.ro.globalenv['expression_data'] = r_expression
                self.ro.globalenv['gene_of_interest'] = gene_of_interest
            else:
                self.ro.globalenv['expression_data'] = self.ro.NULL
                self.ro.globalenv['gene_of_interest'] = self.ro.NULL
            
            # Run R survival analysis
            survival_result = self.r('''
                create_survival_analysis(
                    survival_data = survival_data,
                    expression_data = expression_data,
                    gene_of_interest = gene_of_interest,
                    time_column = time_column,
                    event_column = event_column,
                    group_column = group_column,
                    split_method = split_method
                )
            ''')
            
            return {
                'success': True,
                'method': 'R (survminer)',
                'plot': survival_result.rx2('plot'),
                'survival_fit': survival_result.rx2('survival_fit'),
                'cox_model': survival_result.rx2('cox_model'),
                'summary': survival_result.rx2('summary'),
                'p_value': None,  # Extract from survival fit if needed
                'median_survival': None,  # Extract from survival fit if needed
                'confidence_interval': confidence_interval
            }
            
        except Exception as e:
            logger.error(f"R survival analysis failed: {e}")
            if self.is_python_available():
                logger.info("Falling back to Python survival analysis")
                return self._kaplan_meier_python(
                    survival_data, group_column, expression_data,
                    gene_of_interest, split_method, confidence_interval
                )
            else:
                raise
    
    def _kaplan_meier_python(self, survival_data: SurvivalData,
                           group_column: Optional[str],
                           expression_data: Optional[ExpressionData],
                           gene_of_interest: Optional[str],
                           split_method: str,
                           confidence_interval: float) -> Dict[str, Any]:
        """Kaplan-Meier analysis using Python lifelines"""
        try:
            from lifelines import KaplanMeierFitter
            from lifelines.statistics import logrank_test
            px, go, plt = get_plotting_libs()
            
            # Prepare survival DataFrame
            survival_df = pd.DataFrame({
                'time': survival_data.time,
                'event': survival_data.event
            }, index=survival_data.samples)
            
            # Add grouping information (same logic as R version)
            groups = None
            if gene_of_interest and expression_data:
                if gene_of_interest in expression_data.genes:
                    gene_expr = expression_data.data.loc[gene_of_interest, survival_data.samples]
                    
                    if split_method == "median":
                        threshold = gene_expr.median()
                        survival_df['group'] = (gene_expr > threshold).map({True: 'High', False: 'Low'})
                    elif split_method == "tertiles":
                        tertiles = gene_expr.quantile([0.33, 0.67])
                        survival_df['group'] = pd.cut(
                            gene_expr,
                            bins=[-np.inf, tertiles.iloc[0], tertiles.iloc[1], np.inf],
                            labels=['Low', 'Medium', 'High']
                        )
                    
                    groups = survival_df['group'].unique()
            
            elif group_column and survival_data.covariates is not None:
                if group_column in survival_data.covariates.columns:
                    survival_df['group'] = survival_data.covariates[group_column]
                    groups = survival_df['group'].unique()
            
            # Perform analysis
            kmf = KaplanMeierFitter(alpha=1-confidence_interval)
            
            if groups is not None and len(groups) > 1:
                # Group-based analysis
                if plt:
                    fig, ax = plt.subplots(figsize=(10, 6))
                    
                    group_fits = {}
                    for group in groups:
                        mask = survival_df['group'] == group
                        group_data = survival_df[mask]
                        
                        if len(group_data) >= 5:  # Minimum sample size
                            kmf.fit(
                                durations=group_data['time'],
                                event_observed=group_data['event'],
                                label=f'{group} (n={len(group_data)})'
                            )
                            group_fits[group] = kmf.copy()
                            
                            if plt:
                                kmf.plot_survival_function(ax=ax, ci_show=True)
                    
                    if plt:
                        ax.set_title(f'Survival Analysis: {gene_of_interest or "Groups"}')
                        ax.set_xlabel(f'Time ({self.config.get("time_unit", "days")})')
                        ax.set_ylabel('Survival Probability')
                        ax.grid(True, alpha=0.3)
                        ax.legend()
                    
                    # Log-rank test
                    if len(groups) == 2:
                        group1_data = survival_df[survival_df['group'] == groups[0]]
                        group2_data = survival_df[survival_df['group'] == groups[1]]
                        
                        logrank_result = logrank_test(
                            group1_data['time'], group2_data['time'],
                            group1_data['event'], group2_data['event']
                        )
                        p_value = logrank_result.p_value
                    else:
                        p_value = None
                    
                    return {
                        'success': True,
                        'method': 'Python (lifelines)',
                        'plot': fig if plt else None,
                        'group_fits': group_fits,
                        'p_value': p_value,
                        'logrank_test': logrank_result if len(groups) == 2 else None,
                        'confidence_interval': confidence_interval
                    }
            
            else:
                # Overall survival
                kmf.fit(durations=survival_df['time'], event_observed=survival_df['event'])
                
                if plt:
                    fig, ax = plt.subplots(figsize=(10, 6))
                    kmf.plot_survival_function(ax=ax, ci_show=True)
                    ax.set_title('Overall Survival Curve')
                    ax.set_xlabel(f'Time ({self.config.get("time_unit", "days")})')
                    ax.set_ylabel('Survival Probability')
                    ax.grid(True, alpha=0.3)
                
                return {
                    'success': True,
                    'method': 'Python (lifelines)',
                    'plot': fig if plt else None,
                    'survival_fit': kmf,
                    'median_survival': kmf.median_survival_time_,
                    'confidence_interval': confidence_interval
                }
                
        except ImportError:
            raise RuntimeError("Lifelines not available. Install with: pip install lifelines")
        except Exception as e:
            logger.error(f"Python survival analysis failed: {e}")
            raise
    
    @handle_exception  
    def cox_regression(self, survival_data: SurvivalData,
                      covariates: List[str],
                      expression_data: Optional[ExpressionData] = None,
                      genes_of_interest: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Perform Cox proportional hazards regression
        
        Args:
            survival_data: Survival data container
            covariates: List of covariate column names
            expression_data: Optional expression data for molecular covariates
            genes_of_interest: Optional genes to include as covariates
            
        Returns:
            Dictionary with Cox regression results
        """
        try:
            from lifelines import CoxPHFitter
            
            # Prepare regression DataFrame
            cox_df = pd.DataFrame({
                'time': survival_data.time,
                'event': survival_data.event
            }, index=survival_data.samples)
            
            # Add clinical covariates
            if survival_data.covariates is not None:
                available_covariates = [col for col in covariates 
                                      if col in survival_data.covariates.columns]
                for covar in available_covariates:
                    cox_df[covar] = survival_data.covariates[covar]
            
            # Add gene expression covariates
            if expression_data and genes_of_interest:
                for gene in genes_of_interest:
                    if gene in expression_data.genes:
                        gene_expr = expression_data.data.loc[gene, survival_data.samples]
                        cox_df[f'{gene}_expr'] = gene_expr
            
            # Remove rows with missing data
            cox_df = cox_df.dropna()
            
            if len(cox_df) < 10:
                raise ValueError("Insufficient data for Cox regression after removing missing values")
            
            # Fit Cox model
            cph = CoxPHFitter()
            cph.fit(cox_df, duration_col='time', event_col='event')
            
            return {
                'success': True,
                'method': 'Python (lifelines)',
                'model': cph,
                'summary': cph.summary,
                'concordance': cph.concordance_index_,
                'hazard_ratios': cph.hazard_ratios_,
                'confidence_intervals': cph.confidence_intervals_,
                'p_values': cph.summary['p']
            }
            
        except ImportError:
            raise RuntimeError("Lifelines not available for Cox regression. Install with: pip install lifelines")
        except Exception as e:
            logger.error(f"Cox regression failed: {e}")
            raise
    
    def create_risk_table(self, survival_fits: Dict[str, Any], 
                         time_points: Optional[List[float]] = None) -> pd.DataFrame:
        """
        Create a risk table for survival analysis
        
        Args:
            survival_fits: Dictionary of survival fits by group
            time_points: Time points for risk table (auto-generated if None)
            
        Returns:
            Risk table DataFrame
        """
        if not self.is_python_available():
            logger.warning("Risk table creation requires lifelines")
            return pd.DataFrame()
        
        try:
            # Auto-generate time points if not provided
            if time_points is None:
                max_time = max(fit.timeline.max() for fit in survival_fits.values())
                time_points = np.linspace(0, max_time, 11)[1:]  # Exclude 0
            
            risk_table = []
            
            for group_name, fit in survival_fits.items():
                at_risk = []
                for t in time_points:
                    # Number at risk at time t
                    n_at_risk = (fit.durations >= t).sum()
                    at_risk.append(n_at_risk)
                
                risk_table.append({
                    'Group': group_name,
                    **{f'Time_{int(t)}': n for t, n in zip(time_points, at_risk)}
                })
            
            return pd.DataFrame(risk_table)
            
        except Exception as e:
            logger.error(f"Risk table creation failed: {e}")
            return pd.DataFrame()
    
    def export_results(self, results: Dict[str, Any], 
                      output_dir: str = ".") -> Dict[str, str]:
        """
        Export survival analysis results
        
        Args:
            results: Survival analysis results dictionary
            output_dir: Output directory
            
        Returns:
            Dictionary of exported file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        exported_files = {}
        
        try:
            # Export summary statistics
            if 'summary' in results:
                summary_file = output_path / "survival_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write(str(results['summary']))
                exported_files['summary'] = str(summary_file)
            
            # Export plot if available
            if 'plot' in results and results['plot'] is not None:
                plot_file = output_path / "survival_plot.png"
                results['plot'].savefig(plot_file, dpi=300, bbox_inches='tight')
                exported_files['plot'] = str(plot_file)
            
            # Export Cox model results if available
            if 'model' in results:
                model_file = output_path / "cox_model_summary.csv"
                results['model'].summary.to_csv(model_file)
                exported_files['cox_model'] = str(model_file)
            
            logger.info(f"Exported {len(exported_files)} survival analysis files")
            return exported_files
            
        except Exception as e:
            logger.error(f"Export failed: {e}")
            return {}


# Create singleton instance
_survival_engine = None

def get_survival_engine() -> SurvivalEngine:
    """Get the global survival engine instance"""
    global _survival_engine
    if _survival_engine is None:
        _survival_engine = SurvivalEngine()
    return _survival_engine