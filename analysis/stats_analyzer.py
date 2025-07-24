"""
Prairie Genomics Suite - Statistical Analysis Engine

This module provides basic statistical analysis functions for genomics data
when R/DESeq2 is not available. These functions serve as fallbacks and
provide essential statistical testing capabilities.

Similar to base R statistical functions, these provide core analysis
methods that can be used independently or as fallbacks.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Tuple, Optional
import warnings
from scipy import stats
import logging

from core.utils import validate_dataframe, safe_convert_numeric
from config import get_config

logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)


class BasicStatsAnalyzer:
    """
    Basic statistical analysis for genomics data
    
    Provides essential statistical tests and analyses when more advanced
    tools (like DESeq2) are not available. Similar to base R stats functions.
    """
    
    def __init__(self):
        """Initialize the statistical analyzer"""
        self.config = get_config("analysis")
    
    def ttest_analysis(self, data: pd.DataFrame, 
                      group1_samples: List[str], 
                      group2_samples: List[str],
                      alternative: str = 'two-sided',
                      equal_var: bool = False) -> pd.DataFrame:
        """
        Perform t-test analysis between two groups
        
        Args:
            data: Expression data (genes x samples)
            group1_samples: Sample names for group 1
            group2_samples: Sample names for group 2
            alternative: Type of alternative hypothesis
            equal_var: Whether to assume equal variances
            
        Returns:
            DataFrame with statistical results
        """
        logger.info(f"Running t-test: {len(group1_samples)} vs {len(group2_samples)} samples")
        
        # Validate inputs
        valid, errors = validate_dataframe(data, min_rows=1, min_cols=2, name="Expression data")
        if not valid:
            raise ValueError(f"Data validation failed: {errors}")
        
        # Filter to available samples
        available_samples1 = [s for s in group1_samples if s in data.columns]
        available_samples2 = [s for s in group2_samples if s in data.columns]
        
        if len(available_samples1) < 2 or len(available_samples2) < 2:
            raise ValueError("Need at least 2 samples per group for t-test")
        
        # Extract data for each group
        group1_data = data[available_samples1]
        group2_data = data[available_samples2]
        
        # Calculate statistics for each gene
        results = []
        
        for gene in data.index:
            try:
                # Get expression values
                values1 = group1_data.loc[gene].values
                values2 = group2_data.loc[gene].values
                
                # Remove NaN values
                values1 = values1[~np.isnan(values1)]
                values2 = values2[~np.isnan(values2)]
                
                if len(values1) < 2 or len(values2) < 2:
                    # Not enough data for this gene
                    results.append({
                        'gene': gene,
                        'log2FoldChange': np.nan,
                        'pvalue': np.nan,
                        'padj': np.nan,
                        'baseMean': np.nan,
                        'stat': np.nan
                    })
                    continue
                
                # Calculate means
                mean1 = np.mean(values1)
                mean2 = np.mean(values2)
                baseMean = (mean1 + mean2) / 2
                
                # Calculate fold change (avoid log of zero)
                if mean1 > 0 and mean2 > 0:
                    log2FoldChange = np.log2(mean2 / mean1)
                elif mean2 > 0:
                    log2FoldChange = np.log2(mean2 / 0.01)  # Pseudo-count
                elif mean1 > 0:
                    log2FoldChange = -np.log2(mean1 / 0.01)  # Pseudo-count
                else:
                    log2FoldChange = 0
                
                # Perform t-test
                statistic, pvalue = stats.ttest_ind(
                    values1, values2, 
                    alternative=alternative,
                    equal_var=equal_var
                )
                
                results.append({
                    'gene': gene,
                    'log2FoldChange': log2FoldChange,
                    'pvalue': pvalue,
                    'padj': pvalue,  # Will be corrected later
                    'baseMean': baseMean,
                    'stat': statistic
                })
                
            except Exception as e:
                # Handle problematic genes
                logger.warning(f"Error analyzing gene {gene}: {e}")
                results.append({
                    'gene': gene,
                    'log2FoldChange': np.nan,
                    'pvalue': np.nan,
                    'padj': np.nan,
                    'baseMean': np.nan,
                    'stat': np.nan
                })
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction (Benjamini-Hochberg)
        valid_pvalues = results_df['pvalue'].dropna()
        if len(valid_pvalues) > 0:
            from statsmodels.stats.multitest import multipletests
            _, padj_values, _, _ = multipletests(valid_pvalues, method='fdr_bh')
            
            # Map back to full DataFrame
            padj_map = dict(zip(valid_pvalues.index, padj_values))
            for idx in results_df.index:
                if idx in padj_map:
                    results_df.loc[idx, 'padj'] = padj_map[idx]
        
        # Sort by p-value
        results_df = results_df.sort_values('pvalue')
        
        logger.info(f"T-test completed: {len(results_df)} genes analyzed")
        return results_df
    
    def wilcoxon_analysis(self, data: pd.DataFrame,
                         group1_samples: List[str],
                         group2_samples: List[str]) -> pd.DataFrame:
        """
        Perform Wilcoxon rank-sum test (non-parametric alternative to t-test)
        
        Args:
            data: Expression data (genes x samples)
            group1_samples: Sample names for group 1
            group2_samples: Sample names for group 2
            
        Returns:
            DataFrame with statistical results
        """
        logger.info(f"Running Wilcoxon test: {len(group1_samples)} vs {len(group2_samples)} samples")
        
        # Filter to available samples
        available_samples1 = [s for s in group1_samples if s in data.columns]
        available_samples2 = [s for s in group2_samples if s in data.columns]
        
        if len(available_samples1) < 3 or len(available_samples2) < 3:
            raise ValueError("Need at least 3 samples per group for Wilcoxon test")
        
        # Extract data for each group
        group1_data = data[available_samples1]
        group2_data = data[available_samples2]
        
        results = []
        
        for gene in data.index:
            try:
                # Get expression values
                values1 = group1_data.loc[gene].values
                values2 = group2_data.loc[gene].values
                
                # Remove NaN values
                values1 = values1[~np.isnan(values1)]
                values2 = values2[~np.isnan(values2)]
                
                if len(values1) < 3 or len(values2) < 3:
                    results.append({
                        'gene': gene,
                        'log2FoldChange': np.nan,
                        'pvalue': np.nan,
                        'padj': np.nan,
                        'baseMean': np.nan,
                        'stat': np.nan
                    })
                    continue
                
                # Calculate means and fold change
                mean1 = np.median(values1)  # Use median for non-parametric
                mean2 = np.median(values2)
                baseMean = (mean1 + mean2) / 2
                
                if mean1 > 0 and mean2 > 0:
                    log2FoldChange = np.log2(mean2 / mean1)
                elif mean2 > 0:
                    log2FoldChange = np.log2(mean2 / 0.01)
                elif mean1 > 0:
                    log2FoldChange = -np.log2(mean1 / 0.01)
                else:
                    log2FoldChange = 0
                
                # Perform Wilcoxon rank-sum test
                statistic, pvalue = stats.ranksums(values1, values2)
                
                results.append({
                    'gene': gene,
                    'log2FoldChange': log2FoldChange,
                    'pvalue': pvalue,
                    'padj': pvalue,  # Will be corrected later
                    'baseMean': baseMean,
                    'stat': statistic
                })
                
            except Exception as e:
                logger.warning(f"Error in Wilcoxon test for gene {gene}: {e}")
                results.append({
                    'gene': gene,
                    'log2FoldChange': np.nan,
                    'pvalue': np.nan,
                    'padj': np.nan,
                    'baseMean': np.nan,
                    'stat': np.nan
                })
        
        # Convert to DataFrame and apply multiple testing correction
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        valid_pvalues = results_df['pvalue'].dropna()
        if len(valid_pvalues) > 0:
            from statsmodels.stats.multitest import multipletests
            _, padj_values, _, _ = multipletests(valid_pvalues, method='fdr_bh')
            
            padj_map = dict(zip(valid_pvalues.index, padj_values))
            for idx in results_df.index:
                if idx in padj_map:
                    results_df.loc[idx, 'padj'] = padj_map[idx]
        
        results_df = results_df.sort_values('pvalue')
        logger.info(f"Wilcoxon test completed: {len(results_df)} genes analyzed")
        return results_df
    
    def correlation_analysis(self, data: pd.DataFrame,
                           method: str = 'pearson') -> pd.DataFrame:
        """
        Perform correlation analysis between samples
        
        Args:
            data: Expression data (genes x samples)
            method: Correlation method ('pearson', 'spearman', 'kendall')
            
        Returns:
            Correlation matrix
        """
        logger.info(f"Computing {method} correlation matrix")
        
        valid_methods = ['pearson', 'spearman', 'kendall']
        if method not in valid_methods:
            raise ValueError(f"Method must be one of {valid_methods}")
        
        # Calculate correlation matrix
        corr_matrix = data.corr(method=method)
        
        return corr_matrix
    
    def pca_analysis(self, data: pd.DataFrame, 
                    n_components: int = 3,
                    scale: bool = True) -> Dict[str, Any]:
        """
        Perform Principal Component Analysis
        
        Args:
            data: Expression data (genes x samples)
            n_components: Number of components to compute
            scale: Whether to scale the data
            
        Returns:
            Dictionary with PCA results
        """
        logger.info(f"Running PCA with {n_components} components")
        
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
        except ImportError:
            logger.error("scikit-learn not available for PCA")
            return None
        
        # Prepare data (transpose so samples are rows)
        pca_data = data.T
        
        # Remove genes with no variance
        pca_data = pca_data.loc[:, pca_data.var() > 0]
        
        # Scale data if requested
        if scale:
            scaler = StandardScaler()
            pca_data_scaled = scaler.fit_transform(pca_data)
        else:
            pca_data_scaled = pca_data.values
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(pca_data_scaled)
        
        # Create results DataFrame
        pca_df = pd.DataFrame(
            pca_result,
            index=pca_data.index,
            columns=[f'PC{i+1}' for i in range(n_components)]
        )
        
        return {
            'pca_data': pca_df,
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'explained_variance': pca.explained_variance_,
            'components': pca.components_,
            'feature_names': pca_data.columns.tolist()
        }
    
    def summary_statistics(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate summary statistics for expression data
        
        Args:
            data: Expression data (genes x samples)
            
        Returns:
            DataFrame with summary statistics
        """
        logger.info("Calculating summary statistics")
        
        summary_stats = []
        
        for gene in data.index:
            gene_data = data.loc[gene]
            
            stats_dict = {
                'gene': gene,
                'mean': gene_data.mean(),
                'median': gene_data.median(),
                'std': gene_data.std(),
                'min': gene_data.min(),
                'max': gene_data.max(),
                'q25': gene_data.quantile(0.25),
                'q75': gene_data.quantile(0.75),
                'n_samples': len(gene_data),
                'n_nonzero': (gene_data > 0).sum(),
                'fraction_nonzero': (gene_data > 0).mean()
            }
            
            summary_stats.append(stats_dict)
        
        return pd.DataFrame(summary_stats)
    
    def differential_variability_analysis(self, data: pd.DataFrame,
                                        group1_samples: List[str],
                                        group2_samples: List[str]) -> pd.DataFrame:
        """
        Test for differential variability between groups
        
        Args:
            data: Expression data (genes x samples)
            group1_samples: Sample names for group 1
            group2_samples: Sample names for group 2
            
        Returns:
            DataFrame with variability test results
        """
        logger.info("Testing for differential variability")
        
        # Filter to available samples
        available_samples1 = [s for s in group1_samples if s in data.columns]
        available_samples2 = [s for s in group2_samples if s in data.columns]
        
        group1_data = data[available_samples1]
        group2_data = data[available_samples2]
        
        results = []
        
        for gene in data.index:
            try:
                values1 = group1_data.loc[gene].dropna()
                values2 = group2_data.loc[gene].dropna()
                
                if len(values1) < 3 or len(values2) < 3:
                    results.append({
                        'gene': gene,
                        'var1': np.nan,
                        'var2': np.nan,
                        'var_ratio': np.nan,
                        'pvalue': np.nan
                    })
                    continue
                
                var1 = np.var(values1, ddof=1)
                var2 = np.var(values2, ddof=1)
                
                # F-test for equal variances
                if var1 > 0 and var2 > 0:
                    f_stat = var1 / var2 if var1 > var2 else var2 / var1
                    df1 = len(values1) - 1
                    df2 = len(values2) - 1
                    pvalue = 2 * (1 - stats.f.cdf(f_stat, df1, df2))
                    var_ratio = var2 / var1
                else:
                    pvalue = np.nan
                    var_ratio = np.nan
                
                results.append({
                    'gene': gene,
                    'var1': var1,
                    'var2': var2,
                    'var_ratio': var_ratio,
                    'pvalue': pvalue
                })
                
            except Exception as e:
                logger.warning(f"Error in variability test for gene {gene}: {e}")
                results.append({
                    'gene': gene,
                    'var1': np.nan,
                    'var2': np.nan,
                    'var_ratio': np.nan,
                    'pvalue': np.nan
                })
        
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        valid_pvalues = results_df['pvalue'].dropna()
        if len(valid_pvalues) > 0:
            from statsmodels.stats.multitest import multipletests
            _, padj_values, _, _ = multipletests(valid_pvalues, method='fdr_bh')
            
            padj_map = dict(zip(valid_pvalues.index, padj_values))
            results_df['padj'] = np.nan
            for idx in results_df.index:
                if idx in padj_map:
                    results_df.loc[idx, 'padj'] = padj_map[idx]
        
        return results_df.sort_values('pvalue')


# Create singleton instance
_stats_analyzer = None

def get_stats_analyzer() -> BasicStatsAnalyzer:
    """Get the global stats analyzer instance"""
    global _stats_analyzer
    if _stats_analyzer is None:
        _stats_analyzer = BasicStatsAnalyzer()
    return _stats_analyzer