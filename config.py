"""
Prairie Genomics Suite - Configuration Management

Centralized configuration for the entire application.
This file manages all settings, constants, and environment variables.
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
import streamlit as st

# Version and Application Info
APP_VERSION = "4.0.0"
APP_NAME = "Prairie Genomics Suite - Enhanced Multiomics"
APP_DESCRIPTION = "Advanced Genomics Analysis Platform with Immune Cell Infiltration & Gene Binning"
APP_AUTHOR = "Prairie Genomics Team"

# File paths and directories
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
TEMP_DIR = BASE_DIR / "temp"
CACHE_DIR = BASE_DIR / ".cache"
R_TEMPLATES_PATH = BASE_DIR / "deseq2_templates.R"

# Ensure directories exist
for directory in [DATA_DIR, TEMP_DIR, CACHE_DIR]:
    directory.mkdir(exist_ok=True)

# Streamlit Configuration
STREAMLIT_CONFIG = {
    "page_title": APP_NAME,
    "page_icon": "🧬", 
    "layout": "wide",
    "initial_sidebar_state": "expanded",
    "menu_items": {
        'Get Help': 'https://github.com/prairie-genomics/suite',
        'Report a bug': 'https://github.com/prairie-genomics/suite/issues',
        'About': f"{APP_DESCRIPTION}\n\nVersion: {APP_VERSION}\nAuthor: {APP_AUTHOR}"
    }
}

# Analysis Settings
ANALYSIS_DEFAULTS = {
    "deseq2": {
        "fc_cutoff": 1.0,
        "p_cutoff": 0.05,
        "min_counts": 10,
        "min_samples": 7,
        "normalization_method": "DESeq2"
    },
    "survival": {
        "time_unit": "months",
        "confidence_interval": 0.95,
        "split_methods": ["median", "tertiles", "quartiles"]
    },
    "pathway": {
        "databases": ["KEGG", "Reactome", "GO", "MSigDB"],
        "analysis_methods": ["ORA", "GSEA"],
        "p_cutoff": 0.05,
        "q_cutoff": 0.1
    },
    "immune": {
        "reference_set": "LM22",
        "permutations": 1000,
        "cell_types": [
            "B cells naive", "B cells memory", "Plasma cells",
            "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
            "T cells CD4 memory activated", "T cells follicular helper",
            "T cells regulatory (Tregs)", "T cells gamma delta",
            "NK cells resting", "NK cells activated",
            "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
            "Dendritic cells resting", "Dendritic cells activated",
            "Mast cells resting", "Mast cells activated",
            "Eosinophils", "Neutrophils"
        ]
    }
}

# Literature Search Configuration
LITERATURE_CONFIG = {
    "apis": {
        "pubmed": {
            "tool": "PrairieGenomicsSuite",
            "email": "genomics@prairie.edu",
            "max_results_per_query": 100,
            "rate_limit_delay": 0.5
        },
        "semantic_scholar": {
            "base_url": "https://api.semanticscholar.org/graph/v1",
            "fields": ["title", "abstract", "authors", "venue", "year", "citationCount"],
            "max_results_per_query": 50
        },
        "arxiv": {
            "base_url": "http://export.arxiv.org/api/query",
            "max_results_per_query": 30
        }
    },
    "llm": {
        "openai": {
            "model": "gpt-4",
            "max_tokens": 2000,
            "temperature": 0.1
        },
        "anthropic": {
            "model": "claude-3-sonnet-20240229",
            "max_tokens": 2000,
            "temperature": 0.1
        }
    },
    "genomics_context": {
        "gene_prefixes": ["BRCA", "TP53", "EGFR", "KRAS", "PIK3CA", "APC", "PTEN"],
        "pathways": ["apoptosis", "cell cycle", "DNA repair", "immune response"],
        "conditions": ["cancer", "tumor", "oncology", "genomics", "mutation"]
    }
}

# Visualization Settings
VISUALIZATION_CONFIG = {
    "color_palettes": {
        "default": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"],
        "conditions": ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"],
        "heatmap": "RdYlBu",
        "volcano": {"up": "#d62728", "down": "#1f77b4", "ns": "#808080"}
    },
    "figure_settings": {
        "dpi": 300,
        "format": "png",
        "width": 10,
        "height": 8,
        "font_size": 12
    },
    "plotly_theme": "plotly_white"
}

# Data Validation Rules
DATA_VALIDATION = {
    "expression_data": {
        "min_genes": 100,
        "min_samples": 3,
        "max_missing_rate": 0.5,
        "required_formats": [".csv", ".tsv", ".txt", ".xlsx"]
    },
    "clinical_data": {
        "min_samples": 3,
        "required_columns": [],  # Will be determined by analysis type
        "max_missing_rate": 0.8
    },
    "survival_data": {
        "required_columns": ["time", "event"],
        "time_min": 0,
        "event_values": [0, 1]
    }
}

# Performance Settings
PERFORMANCE_CONFIG = {
    "caching": {
        "enabled": True,
        "ttl": 3600,  # 1 hour
        "max_entries": 100
    },
    "parallel_processing": {
        "enabled": True,
        "max_workers": 4
    },
    "memory": {
        "max_file_size_mb": 500,
        "chunk_size": 10000
    }
}

# Environment Variables
def get_env_var(key: str, default: Optional[str] = None) -> Optional[str]:
    """Get environment variable with optional default"""
    return os.getenv(key, default)

# API Keys (from environment variables)
API_KEYS = {
    "openai": get_env_var("OPENAI_API_KEY"),
    "anthropic": get_env_var("ANTHROPIC_API_KEY"),
    "semantic_scholar": get_env_var("SEMANTIC_SCHOLAR_API_KEY")  # Optional
}

# Debug and Development Settings
DEBUG_MODE = get_env_var("DEBUG", "false").lower() == "true"
LOG_LEVEL = get_env_var("LOG_LEVEL", "INFO")

# Feature Flags
FEATURE_FLAGS = {
    "enable_r_integration": True,
    "enable_llm_literature_search": True,
    "enable_advanced_visualizations": True,
    "enable_parallel_processing": True,
    "enable_caching": True
}

def get_config(section: str) -> Dict[str, Any]:
    """Get configuration for a specific section"""
    config_map = {
        "app": {
            "name": APP_NAME,
            "version": APP_VERSION,
            "description": APP_DESCRIPTION,
            "author": APP_AUTHOR
        },
        "streamlit": STREAMLIT_CONFIG,
        "analysis": ANALYSIS_DEFAULTS,
        "literature": LITERATURE_CONFIG,
        "visualization": VISUALIZATION_CONFIG,
        "validation": DATA_VALIDATION,
        "performance": PERFORMANCE_CONFIG,
        "api_keys": API_KEYS,
        "features": FEATURE_FLAGS
    }
    
    return config_map.get(section, {})

def is_feature_enabled(feature: str) -> bool:
    """Check if a feature flag is enabled"""
    return FEATURE_FLAGS.get(feature, False)

def validate_api_keys() -> Dict[str, bool]:
    """Validate that required API keys are available"""
    return {
        "openai": API_KEYS["openai"] is not None,
        "anthropic": API_KEYS["anthropic"] is not None,
        "semantic_scholar": True  # Optional, so always return True
    }