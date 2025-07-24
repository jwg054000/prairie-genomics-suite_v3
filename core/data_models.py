"""
Prairie Genomics Suite - Core Data Models

This module defines the core data structures and types used throughout
the application. These models provide type safety and data validation.

Think of these as similar to S4 classes in R - they define the structure
and validation rules for our genomics data.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any, Tuple
import pandas as pd
import numpy as np
from datetime import datetime
from enum import Enum


class AnalysisType(Enum):
    """Enumeration of supported analysis types"""
    DESEQ2 = "deseq2"
    SURVIVAL = "survival"
    PATHWAY = "pathway"
    IMMUNE = "immune"
    LITERATURE = "literature"
    VISUALIZATION = "visualization"


class FileFormat(Enum):
    """Supported file formats for data import"""
    CSV = "csv"
    TSV = "tsv"
    TXT = "txt"
    XLSX = "xlsx"
    XLS = "xls"


@dataclass
class ExpressionData:
    """
    Expression data container (similar to ExpressionSet in R/Bioconductor)
    
    Attributes:
        data: Gene expression matrix (genes x samples)
        genes: Gene identifiers 
        samples: Sample identifiers
        metadata: Additional metadata about the dataset
    """
    data: pd.DataFrame
    genes: List[str] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate and initialize data after creation"""
        if self.data is not None:
            self.genes = list(self.data.index) if not self.genes else self.genes
            self.samples = list(self.data.columns) if not self.samples else self.samples
            
    @property
    def n_genes(self) -> int:
        """Number of genes in the dataset"""
        return len(self.genes)
        
    @property 
    def n_samples(self) -> int:
        """Number of samples in the dataset"""
        return len(self.samples)
        
    def filter_genes(self, min_expression: float = 1.0, min_samples: int = 3) -> 'ExpressionData':
        """Filter low-expression genes"""
        if self.data is None:
            return self
            
        # Filter genes with sufficient expression
        keep_genes = (self.data > min_expression).sum(axis=1) >= min_samples
        filtered_data = self.data.loc[keep_genes]
        
        return ExpressionData(
            data=filtered_data,
            genes=list(filtered_data.index),
            samples=self.samples,
            metadata=self.metadata
        )


@dataclass
class ClinicalData:
    """
    Clinical/metadata container for samples
    
    Attributes:
        data: Clinical data DataFrame (samples x variables)
        samples: Sample identifiers
        variables: Clinical variable names
        metadata: Additional metadata
    """
    data: pd.DataFrame
    samples: List[str] = field(default_factory=list)
    variables: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate and initialize data after creation"""
        if self.data is not None:
            self.samples = list(self.data.index) if not self.samples else self.samples
            self.variables = list(self.data.columns) if not self.variables else self.variables
            
    @property
    def n_samples(self) -> int:
        """Number of samples"""
        return len(self.samples)
        
    @property
    def n_variables(self) -> int:
        """Number of clinical variables"""
        return len(self.variables)


@dataclass
class SurvivalData:
    """
    Survival analysis data container
    
    Attributes:
        time: Survival times
        event: Event indicators (0=censored, 1=event)
        samples: Sample identifiers
        covariates: Additional covariates for Cox models
    """
    time: Union[pd.Series, np.ndarray, List[float]]
    event: Union[pd.Series, np.ndarray, List[int]]
    samples: List[str] = field(default_factory=list)
    covariates: Optional[pd.DataFrame] = None
    
    def __post_init__(self):
        """Validate survival data"""
        if len(self.time) != len(self.event):
            raise ValueError("Time and event vectors must have the same length")
            
        if not self.samples:
            self.samples = [f"Sample_{i}" for i in range(len(self.time))]


@dataclass 
class DESeq2Results:
    """
    DESeq2 analysis results container
    
    Attributes:
        results: Main results DataFrame with log2FC, pvalue, padj
        contrasts: List of contrasts analyzed
        normalized_counts: Normalized count matrix
        metadata: Analysis metadata (parameters, etc.)
    """
    results: Dict[str, pd.DataFrame]  # One DataFrame per contrast
    contrasts: List[str] = field(default_factory=list)
    normalized_counts: Optional[pd.DataFrame] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_significant_genes(self, contrast: str, 
                            fc_cutoff: float = 1.0, 
                            p_cutoff: float = 0.05) -> List[str]:
        """Get significantly differentially expressed genes"""
        if contrast not in self.results:
            return []
            
        result_df = self.results[contrast]
        significant = result_df[
            (abs(result_df['log2FoldChange']) > fc_cutoff) & 
            (result_df['padj'] < p_cutoff)
        ]
        return list(significant.index)


@dataclass
class PathwayResults:
    """
    Pathway analysis results container
    
    Attributes:
        results: Pathway enrichment results
        gene_sets: Gene sets used for analysis
        analysis_method: Method used (ORA, GSEA, etc.)
        database: Database used (KEGG, GO, etc.)
    """
    results: pd.DataFrame
    gene_sets: Dict[str, List[str]] = field(default_factory=dict)
    analysis_method: str = "ORA"
    database: str = "KEGG"
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_significant_pathways(self, p_cutoff: float = 0.05) -> pd.DataFrame:
        """Get significantly enriched pathways"""
        if 'padj' in self.results.columns:
            return self.results[self.results['padj'] < p_cutoff]
        elif 'pvalue' in self.results.columns:
            return self.results[self.results['pvalue'] < p_cutoff]
        else:
            return self.results


@dataclass
class ImmuneResults:
    """
    Immune infiltration analysis results
    
    Attributes:
        cell_fractions: Estimated cell type fractions per sample
        p_values: P-values for deconvolution quality
        correlations: Correlations between cell types
        method: Deconvolution method used
    """
    cell_fractions: pd.DataFrame  # samples x cell_types
    p_values: Optional[pd.Series] = None
    correlations: Optional[pd.DataFrame] = None
    method: str = "CIBERSORT"
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class LiteratureSearchResult:
    """
    Individual literature search result
    
    Attributes:
        title: Paper title
        authors: Author list
        journal: Journal name
        year: Publication year
        pmid: PubMed ID
        abstract: Abstract text
        relevance_score: Computed relevance score
        keywords: Keywords/MeSH terms
    """
    title: str
    authors: str
    journal: str
    year: str
    pmid: str
    abstract: str
    relevance_score: float = 0.0
    keywords: List[str] = field(default_factory=list)
    url: Optional[str] = None


@dataclass
class LiteratureResults:
    """
    Literature search results container
    
    Attributes:
        results: Dictionary of gene/term -> list of papers
        summary: Summary statistics
        search_terms: Terms that were searched
        search_parameters: Search parameters used
    """
    results: Dict[str, List[LiteratureSearchResult]]
    summary: Dict[str, Any] = field(default_factory=dict)
    search_terms: List[str] = field(default_factory=list)
    search_parameters: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def total_papers(self) -> int:
        """Total number of papers found"""
        return sum(len(papers) for papers in self.results.values())


@dataclass
class AnalysisSession:
    """
    Complete analysis session container
    
    This is like the main workspace in R - contains all data and results
    for a complete genomics analysis session.
    """
    session_id: str
    created_at: datetime = field(default_factory=datetime.now)
    
    # Input data
    expression_data: Optional[ExpressionData] = None
    clinical_data: Optional[ClinicalData] = None
    survival_data: Optional[SurvivalData] = None
    
    # Analysis results
    deseq2_results: Optional[DESeq2Results] = None
    pathway_results: Optional[PathwayResults] = None
    immune_results: Optional[ImmuneResults] = None
    literature_results: Optional[LiteratureResults] = None
    
    # Analysis metadata
    analysis_history: List[Dict[str, Any]] = field(default_factory=list)
    settings: Dict[str, Any] = field(default_factory=dict)
    
    def add_analysis(self, analysis_type: AnalysisType, results: Any, 
                    parameters: Dict[str, Any] = None):
        """Add analysis results to the session"""
        timestamp = datetime.now()
        
        # Store results
        if analysis_type == AnalysisType.DESEQ2:
            self.deseq2_results = results
        elif analysis_type == AnalysisType.PATHWAY:
            self.pathway_results = results
        elif analysis_type == AnalysisType.IMMUNE:
            self.immune_results = results
        elif analysis_type == AnalysisType.LITERATURE:
            self.literature_results = results
            
        # Add to history
        self.analysis_history.append({
            "timestamp": timestamp,
            "analysis_type": analysis_type.value,
            "parameters": parameters or {},
            "status": "completed"
        })
    
    def get_analysis_summary(self) -> Dict[str, Any]:
        """Get summary of all analyses in this session"""
        summary = {
            "session_id": self.session_id,
            "created_at": self.created_at,
            "data_loaded": {
                "expression": self.expression_data is not None,
                "clinical": self.clinical_data is not None,
                "survival": self.survival_data is not None
            },
            "analyses_completed": {
                "deseq2": self.deseq2_results is not None,
                "pathway": self.pathway_results is not None,
                "immune": self.immune_results is not None,
                "literature": self.literature_results is not None
            },
            "total_analyses": len(self.analysis_history)
        }
        
        if self.expression_data:
            summary["data_summary"] = {
                "n_genes": self.expression_data.n_genes,
                "n_samples": self.expression_data.n_samples
            }
            
        return summary


# Type aliases for commonly used types
GenesType = List[str]
SamplesType = List[str]
ExpressionMatrix = pd.DataFrame
ResultsDict = Dict[str, Any]