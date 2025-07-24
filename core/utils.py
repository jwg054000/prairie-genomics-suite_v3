"""
Prairie Genomics Suite - Core Utilities

This module contains utility functions that are shared across the entire
application. These are like the helper functions you'd put in a utils.R
file in an R package.
"""

import pandas as pd
import numpy as np
import streamlit as st
import io
import warnings
from typing import Dict, List, Optional, Union, Any, Tuple
from pathlib import Path
import re
import hashlib
from datetime import datetime
import logging

from config import get_config

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_dataframe(df: pd.DataFrame, 
                      min_rows: int = 1, 
                      min_cols: int = 1,
                      required_columns: List[str] = None,
                      name: str = "DataFrame",
                      allow_nan_columns: bool = False,
                      strict_validation: bool = True) -> Tuple[bool, List[str]]:
    """
    Validate a DataFrame meets minimum requirements
    
    Args:
        df: DataFrame to validate
        min_rows: Minimum number of rows required
        min_cols: Minimum number of columns required
        required_columns: List of required column names
        name: Name for error messages
        allow_nan_columns: Whether to allow columns with all NaN values
        strict_validation: Whether to apply strict validation rules
    
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    warnings = []
    
    if df is None:
        errors.append(f"{name} is None")
        return False, errors
    
    if not isinstance(df, pd.DataFrame):
        errors.append(f"{name} is not a pandas DataFrame")
        return False, errors
    
    if len(df) < min_rows:
        errors.append(f"{name} has {len(df)} rows, minimum {min_rows} required")
    
    if len(df.columns) < min_cols:
        errors.append(f"{name} has {len(df.columns)} columns, minimum {min_cols} required")
    
    if required_columns:
        missing_cols = set(required_columns) - set(df.columns)
        if missing_cols:
            errors.append(f"{name} missing required columns: {list(missing_cols)}")
    
    # Check for all NaN columns (more flexible for clinical data)
    all_nan_cols = df.columns[df.isna().all()].tolist()
    if all_nan_cols:
        if allow_nan_columns:
            warnings.append(f"{name} has columns with all NaN values that will be removed: {all_nan_cols[:5]}{'...' if len(all_nan_cols) > 5 else ''}")
        elif strict_validation:
            errors.append(f"{name} has columns with all NaN values: {all_nan_cols[:5]}{'...' if len(all_nan_cols) > 5 else ''}")
        else:
            warnings.append(f"{name} has columns with all NaN values that will be ignored: {all_nan_cols[:5]}{'...' if len(all_nan_cols) > 5 else ''}")
    
    # Log warnings if any
    if warnings:
        logger.warning(f"{name} validation warnings: {'; '.join(warnings)}")
    
    return len(errors) == 0, errors


def clean_clinical_data(df: pd.DataFrame, 
                       remove_nan_columns: bool = True,
                       min_valid_fraction: float = 0.1) -> pd.DataFrame:
    """
    Clean clinical data by removing problematic columns and handling missing values
    
    Args:
        df: Clinical data DataFrame
        remove_nan_columns: Whether to remove columns with all NaN values
        min_valid_fraction: Minimum fraction of non-NaN values required to keep column
        
    Returns:
        Cleaned DataFrame
    """
    if df is None or df.empty:
        return df
    
    df_clean = df.copy()
    original_cols = len(df_clean.columns)
    
    if remove_nan_columns:
        # Remove columns with all NaN values
        all_nan_cols = df_clean.columns[df_clean.isna().all()].tolist()
        if all_nan_cols:
            df_clean = df_clean.drop(columns=all_nan_cols)
            logger.info(f"Removed {len(all_nan_cols)} columns with all NaN values")
    
    # Remove columns with too few valid values
    if min_valid_fraction > 0:
        valid_fraction = df_clean.notna().mean()
        cols_to_remove = valid_fraction[valid_fraction < min_valid_fraction].index.tolist()
        if cols_to_remove:
            df_clean = df_clean.drop(columns=cols_to_remove)
            logger.info(f"Removed {len(cols_to_remove)} columns with <{min_valid_fraction*100}% valid values")
    
    # Log summary
    final_cols = len(df_clean.columns)
    if final_cols != original_cols:
        logger.info(f"Clinical data cleaning: {original_cols} → {final_cols} columns")
    
    return df_clean


def create_minimal_clinical_data(expression_data, 
                                condition_column: str = "Condition",
                                default_conditions: List[str] = None) -> pd.DataFrame:
    """
    Create minimal clinical data from sample names when no clinical data is provided
    
    Args:
        expression_data: Expression data with sample names
        condition_column: Name for the condition column
        default_conditions: List of condition names to cycle through
        
    Returns:
        Minimal clinical DataFrame
    """
    if default_conditions is None:
        default_conditions = ["Control", "Treatment"]
    
    # Get sample names
    if hasattr(expression_data, 'samples'):
        sample_names = expression_data.samples
    elif hasattr(expression_data, 'columns'):
        sample_names = list(expression_data.columns)
    else:
        raise ValueError("Cannot determine sample names from expression data")
    
    # Create alternating conditions
    n_samples = len(sample_names)
    conditions = []
    for i in range(n_samples):
        condition_idx = i % len(default_conditions)
        conditions.append(default_conditions[condition_idx])
    
    # Create minimal clinical data
    clinical_data = pd.DataFrame({
        condition_column: conditions,
        'sample_id': sample_names,
        'created_by': 'auto_generated'
    }, index=sample_names)
    
    logger.info(f"Created minimal clinical data: {n_samples} samples, {len(default_conditions)} conditions")
    return clinical_data


def clean_gene_names(genes: List[str]) -> List[str]:
    """
    Clean gene names by removing common prefixes/suffixes and standardizing
    
    Args:
        genes: List of gene names to clean
        
    Returns:
        List of cleaned gene names
    """
    cleaned = []
    
    for gene in genes:
        if pd.isna(gene) or gene == "":
            continue
            
        # Convert to string and strip whitespace
        gene = str(gene).strip()
        
        # Remove common prefixes (ENSG, ENST, etc.)
        gene = re.sub(r'^(ENSG|ENST|ENSMUSG|ENSMUSP)\d+\|?', '', gene)
        
        # Remove version numbers (e.g., GENE.1 -> GENE)
        gene = re.sub(r'\.\d+$', '', gene)
        
        # Remove extra quotes or brackets
        gene = re.sub(r'^["\'\[\]]+|["\'\[\]]+$', '', gene)
        
        # Convert to uppercase for consistency
        gene = gene.upper()
        
        if gene:  # Only add non-empty genes
            cleaned.append(gene)
    
    return list(set(cleaned))  # Remove duplicates


def detect_file_format(file_path: Union[str, Path]) -> str:
    """
    Detect file format from file extension
    
    Args:
        file_path: Path to the file
        
    Returns:
        File format string (csv, tsv, xlsx, etc.)
    """
    path = Path(file_path)
    suffix = path.suffix.lower()
    
    format_map = {
        '.csv': 'csv',
        '.tsv': 'tsv', 
        '.txt': 'tsv',  # Assume tab-separated for .txt
        '.xlsx': 'xlsx',
        '.xls': 'xls'
    }
    
    return format_map.get(suffix, 'csv')  # Default to CSV


def read_data_file(file_content: Union[io.StringIO, io.BytesIO], 
                  file_name: str,
                  sep: str = None) -> pd.DataFrame:
    """
    Read data file with automatic format detection and error handling
    
    Args:
        file_content: File content from Streamlit file_uploader
        file_name: Name of the uploaded file
        sep: Separator character (auto-detected if None)
        
    Returns:
        pandas DataFrame
    """
    file_format = detect_file_format(file_name)
    
    try:
        if file_format in ['xlsx', 'xls']:
            df = pd.read_excel(file_content, index_col=0)
        else:
            # For text files, try to detect separator
            if sep is None:
                # Read first few lines to detect separator
                if isinstance(file_content, io.BytesIO):
                    file_content.seek(0)
                    first_lines = file_content.read(1024).decode('utf-8')
                    file_content.seek(0)
                else:
                    file_content.seek(0)
                    first_lines = file_content.read(1024)
                    file_content.seek(0)
                
                # Count separators in first line
                first_line = first_lines.split('\n')[0]
                comma_count = first_line.count(',')
                tab_count = first_line.count('\t')
                
                sep = ',' if comma_count > tab_count else '\t'
            
            df = pd.read_csv(file_content, sep=sep, index_col=0)
        
        return df
        
    except Exception as e:
        st.error(f"Error reading file {file_name}: {str(e)}")
        raise


def filter_expression_data(expr_data: pd.DataFrame, 
                         min_expression: float = 1.0,
                         min_samples: int = 3,
                         min_variance: float = 0.0) -> pd.DataFrame:
    """
    Filter expression data to remove low-expression and low-variance genes
    
    Args:
        expr_data: Expression data (genes x samples)
        min_expression: Minimum expression level
        min_samples: Minimum number of samples with expression > min_expression
        min_variance: Minimum variance across samples
        
    Returns:
        Filtered expression data
    """
    logger.info(f"Filtering expression data: {expr_data.shape[0]} genes, {expr_data.shape[1]} samples")
    
    # Filter by expression level
    expressed_mask = (expr_data > min_expression).sum(axis=1) >= min_samples
    logger.info(f"Genes passing expression filter: {expressed_mask.sum()}")
    
    # Filter by variance
    if min_variance > 0:
        variance_mask = expr_data.var(axis=1) >= min_variance
        logger.info(f"Genes passing variance filter: {variance_mask.sum()}")
        final_mask = expressed_mask & variance_mask
    else:
        final_mask = expressed_mask
    
    filtered_data = expr_data.loc[final_mask]
    logger.info(f"Final filtered data: {filtered_data.shape[0]} genes, {filtered_data.shape[1]} samples")
    
    return filtered_data


def normalize_sample_names(sample_names: List[str]) -> List[str]:
    """
    Normalize sample names by removing invalid characters and ensuring uniqueness
    
    Args:
        sample_names: List of sample names
        
    Returns:
        List of normalized sample names
    """
    normalized = []
    
    for i, name in enumerate(sample_names):
        if pd.isna(name) or name == "":
            name = f"Sample_{i+1}"
        
        # Convert to string and clean
        name = str(name).strip()
        
        # Replace invalid characters with underscores
        name = re.sub(r'[^\w\-\.]', '_', name)
        
        # Ensure it doesn't start with a number
        if name[0].isdigit():
            name = f"S_{name}"
        
        normalized.append(name)
    
    # Ensure uniqueness
    seen = set()
    unique_names = []
    
    for name in normalized:
        original_name = name
        counter = 1
        
        while name in seen:
            name = f"{original_name}_{counter}"
            counter += 1
        
        seen.add(name)
        unique_names.append(name)
    
    return unique_names


def match_samples(expr_samples: List[str], 
                 clinical_samples: List[str]) -> Tuple[List[str], List[str], List[str]]:
    """
    Match samples between expression and clinical data
    
    Args:
        expr_samples: Sample names from expression data
        clinical_samples: Sample names from clinical data
        
    Returns:
        Tuple of (common_samples, expr_only, clinical_only)
    """
    expr_set = set(expr_samples)
    clinical_set = set(clinical_samples)
    
    common = list(expr_set & clinical_set)
    expr_only = list(expr_set - clinical_set)
    clinical_only = list(clinical_set - expr_set)
    
    return common, expr_only, clinical_only


def create_summary_stats(df: pd.DataFrame, name: str = "Dataset") -> Dict[str, Any]:
    """
    Create summary statistics for a dataset
    
    Args:
        df: DataFrame to summarize
        name: Name of the dataset
        
    Returns:
        Dictionary of summary statistics
    """
    summary = {
        'name': name,
        'shape': df.shape,
        'n_rows': len(df),
        'n_cols': len(df.columns),
        'memory_usage_mb': df.memory_usage(deep=True).sum() / 1024 / 1024,
        'missing_values': df.isna().sum().sum(),
        'missing_percent': (df.isna().sum().sum() / df.size) * 100
    }
    
    # Numeric columns summary
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        summary['numeric_summary'] = {
            'n_numeric_cols': len(numeric_cols),
            'mean': df[numeric_cols].mean().mean(),
            'median': df[numeric_cols].median().median(),
            'std': df[numeric_cols].std().mean()
        }
    
    # Categorical columns summary
    cat_cols = df.select_dtypes(include=['object', 'category']).columns
    if len(cat_cols) > 0:
        summary['categorical_summary'] = {
            'n_categorical_cols': len(cat_cols),
            'unique_values': {col: df[col].nunique() for col in cat_cols}
        }
    
    return summary


def safe_convert_numeric(series: pd.Series, 
                        errors: str = 'coerce') -> pd.Series:
    """
    Safely convert a series to numeric, handling common issues
    
    Args:
        series: Pandas series to convert
        errors: How to handle conversion errors ('coerce', 'ignore', 'raise')
        
    Returns:
        Converted series
    """
    # Handle common string representations
    if series.dtype == 'object':
        # Replace common non-numeric strings
        series = series.astype(str)
        series = series.str.replace(',', '')  # Remove thousands separators
        series = series.str.replace('$', '')  # Remove currency symbols
        series = series.str.replace('%', '')  # Remove percent symbols
        series = series.replace(['', 'NA', 'N/A', 'null', 'NULL'], np.nan)
    
    return pd.to_numeric(series, errors=errors)


def calculate_md5_hash(data: Union[pd.DataFrame, str, bytes]) -> str:
    """
    Calculate MD5 hash of data for caching/comparison purposes
    
    Args:
        data: Data to hash
        
    Returns:
        MD5 hash string
    """
    if isinstance(data, pd.DataFrame):
        data_string = data.to_string()
    elif isinstance(data, str):
        data_string = data
    elif isinstance(data, bytes):
        return hashlib.md5(data).hexdigest()
    else:
        data_string = str(data)
    
    return hashlib.md5(data_string.encode()).hexdigest()


def format_number(num: Union[int, float], 
                 precision: int = 2,
                 use_thousands_sep: bool = True) -> str:
    """
    Format numbers for display
    
    Args:
        num: Number to format
        precision: Decimal places
        use_thousands_sep: Whether to use thousands separator
        
    Returns:
        Formatted number string
    """
    if pd.isna(num):
        return "N/A"
    
    if isinstance(num, (int, np.integer)):
        if use_thousands_sep:
            return f"{num:,}"
        else:
            return str(num)
    
    if isinstance(num, (float, np.floating)):
        if use_thousands_sep:
            return f"{num:,.{precision}f}"
        else:
            return f"{num:.{precision}f}"
    
    return str(num)


def create_download_link(data: Union[pd.DataFrame, str, bytes],
                        filename: str,
                        link_text: str = "Download") -> str:
    """
    Create a download link for data
    
    Args:
        data: Data to download
        filename: Name for downloaded file
        link_text: Text for the download link
        
    Returns:
        HTML download link
    """
    if isinstance(data, pd.DataFrame):
        if filename.endswith('.csv'):
            data_bytes = data.to_csv().encode()
            mime_type = 'text/csv'
        elif filename.endswith('.xlsx'):
            output = io.BytesIO()
            data.to_excel(output, index=True)
            data_bytes = output.getvalue()
            mime_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        else:
            data_bytes = data.to_string().encode()
            mime_type = 'text/plain'
    elif isinstance(data, str):
        data_bytes = data.encode()
        mime_type = 'text/plain'
    else:
        data_bytes = data
        mime_type = 'application/octet-stream'
    
    # Use Streamlit's download_button
    return st.download_button(
        label=link_text,
        data=data_bytes,
        file_name=filename,
        mime=mime_type
    )


def handle_exception(func):
    """
    Decorator for consistent exception handling
    
    Usage:
        @handle_exception
        def my_function():
            # function code
    """
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Error in {func.__name__}: {str(e)}")
            st.error(f"An error occurred in {func.__name__}: {str(e)}")
            return None
    
    return wrapper


def progress_bar(iterable, desc: str = "Processing"):
    """
    Simple progress bar for iterables
    
    Args:
        iterable: Iterable to wrap with progress bar
        desc: Description text
        
    Yields:
        Items from the iterable
    """
    total = len(iterable) if hasattr(iterable, '__len__') else None
    
    if total:
        progress_text = st.empty()
        progress_bar = st.progress(0)
        
        for i, item in enumerate(iterable):
            progress = (i + 1) / total
            progress_bar.progress(progress)
            progress_text.text(f"{desc}: {i + 1}/{total}")
            yield item
        
        progress_bar.empty()
        progress_text.empty()
    else:
        # If we can't determine length, just show spinning indicator
        with st.spinner(desc):
            for item in iterable:
                yield item