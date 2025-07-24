"""
Prairie Genomics Suite - Analysis Module

This module contains all the analysis engines for genomics data processing.
Each engine is designed to be independent and testable, following the
modular design principles similar to R package organization.

Components:
- deseq2_engine: Differential expression analysis with R integration
- survival_engine: Survival analysis and Kaplan-Meier curves
- pathway_engine: Pathway enrichment and network analysis  
- immune_engine: Immune cell infiltration analysis
- literature_engine: LLM-powered literature search and synthesis
- visualization_engine: Publication-quality plotting and charts
"""

__version__ = "3.0.0"
__author__ = "Prairie Genomics Team"