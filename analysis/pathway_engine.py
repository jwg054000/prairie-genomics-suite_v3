"""
Prairie Genomics Suite - Pathway Analysis Engine

This module provides comprehensive pathway analysis capabilities including
Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA)
with support for multiple databases (KEGG, Reactome, GO, MSigDB, WikiPathways).

Features:
- Multiple pathway databases (KEGG, Reactome, GO, MSigDB, WikiPathways)
- Both ORA and GSEA analysis methods
- R integration for advanced pathway networks
- Python fallbacks with GSEApy and basic enrichment
- Interactive pathway network visualization
- Export functionality for results

Author: Prairie Genomics Team
"""

import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
import warnings
import logging
from scipy import stats

# Import core modules
from core.data_models import PathwayResults, ExpressionData, DESeq2Results
from core.utils import validate_dataframe, handle_exception
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


class PathwayEngine:
    """
    Comprehensive Pathway Analysis Engine
    
    Provides pathway enrichment analysis using multiple methods and databases.
    Similar to clusterProfiler in R but with additional database support
    and Python fallbacks.
    """
    
    def __init__(self):
        """Initialize the pathway analysis engine"""
        self.r_available = False
        self.gseapy_available = False
        self.r = None
        self.ro = None
        self.config = get_config("analysis")["pathway"]
        
        # Initialize R integration if available
        if is_feature_enabled("enable_r_integration"):
            self._initialize_r()
        
        # Check for GSEApy
        self._check_gseapy()
    
    def _initialize_r(self):
        """Initialize R environment for pathway analysis"""
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
                logger.info("Pathway R integration initialized successfully")
            else:
                logger.warning("R script templates not found. R pathway analysis unavailable.")
                
        except ImportError as e:
            logger.warning(f"rpy2 not available: {e}. Pathway analysis will use Python fallback.")
        except Exception as e:
            logger.error(f"R integration failed: {e}")
    
    def _check_gseapy(self):
        """Check if GSEApy is available for comprehensive pathway analysis"""
        try:
            import gseapy
            self.gseapy_available = True
            logger.info("GSEApy available for comprehensive pathway analysis")
        except ImportError:
            logger.warning("GSEApy not available. Install with: pip install gseapy")
    
    def is_r_available(self) -> bool:
        """Check if R pathway analysis is available"""
        return self.r_available
    
    def is_gseapy_available(self) -> bool:
        """Check if GSEApy is available"""
        return self.gseapy_available
    
    def get_status(self) -> Dict[str, Any]:
        """Get status information about the engine"""
        return {
            "r_available": self.r_available,
            "gseapy_available": self.gseapy_available,
            "preferred_method": "GSEApy" if self.gseapy_available else "Basic enrichment",
            "supported_databases": self.get_supported_databases(),
            "config": self.config
        }
    
    def get_supported_databases(self) -> List[str]:
        """Get list of supported pathway databases"""
        basic_databases = ["Basic Pathways (Built-in)"]
        
        if self.gseapy_available:
            gseapy_databases = [
                "KEGG Pathways",
                "Reactome", 
                "GO Biological Process",
                "GO Molecular Function",
                "GO Cellular Component",
                "MSigDB Hallmarks",
                "MSigDB C2 (Curated)",
                "WikiPathways"
            ]
            return basic_databases + gseapy_databases
        
        return basic_databases
    
    def get_basic_pathways(self) -> Dict[str, List[str]]:
        """Return basic pathway gene sets including immune pathways"""
        return {
            "Cell Cycle": [
                "CDK1", "CDK2", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", 
                "CCND1", "CCNE1", "RB1", "E2F1"
            ],
            "Apoptosis": [
                "TP53", "BAX", "BCL2", "CASP3", "CASP8", "CASP9", "PARP1", 
                "APAF1", "CYCS", "BAK1"
            ],
            "DNA Repair": [
                "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", 
                "XRCC1", "ERCC1", "MLH1"
            ],
            "PI3K/AKT Pathway": [
                "PIK3CA", "PIK3CB", "AKT1", "AKT2", "PTEN", "MTOR", "GSK3B", 
                "FOXO1", "PDK1"
            ],
            "MAPK Pathway": [
                "KRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAF1", 
                "EGFR", "GRB2"
            ],
            "Immune Response": [
                "CD3D", "CD4", "CD8A", "IFNG", "IL2", "TNF", "CTLA4", "PDCD1", 
                "CD274", "LAG3"
            ],
            "Metabolism": [
                "GLUT1", "HK2", "PKM", "LDHA", "PDK1", "ACLY", "FASN", "CPT1A", 
                "PPARA", "SREBF1"
            ],
            "Immune Checkpoint": [
                "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA"
            ],
            "T Cell Activation": [
                "CD3D", "CD3E", "CD3G", "CD28", "ICOS", "CD40LG", "IL2"
            ],
            "Interferon Signaling": [
                "IFNG", "IFNA1", "IFNB1", "STAT1", "STAT2", "IRF1", "IRF7"
            ],
            "Complement System": [
                "C1QA", "C1QB", "C3", "C5", "CFH", "CFI", "CD55"
            ],
            "Antigen Presentation": [
                "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "TAP1"
            ]
        }
    
    @handle_exception
    def basic_enrichment_analysis(self, gene_list: List[str], 
                                background_size: int = 20000) -> PathwayResults:
        """
        Basic pathway enrichment analysis using built-in pathways
        
        Args:
            gene_list: List of genes to analyze
            background_size: Size of background gene universe
            
        Returns:
            PathwayResults container with enrichment results
        """
        logger.info(f"Running basic enrichment analysis for {len(gene_list)} genes")
        
        pathways = self.get_basic_pathways()
        results = []
        
        for pathway_name, pathway_genes in pathways.items():
            # Calculate overlap
            overlap = set(gene_list) & set(pathway_genes)
            overlap_count = len(overlap)
            
            if overlap_count == 0:
                continue
                
            # Hypergeometric test
            total_genes = len(gene_list)
            pathway_size = len(pathway_genes)
            
            # P-value calculation (hypergeometric)
            p_value = stats.hypergeom.sf(overlap_count - 1, background_size, pathway_size, total_genes)
            
            # Enrichment score
            expected = (total_genes * pathway_size) / background_size
            enrichment_score = overlap_count / expected if expected > 0 else 0
            
            results.append({
                'ID': pathway_name.replace(" ", "_").lower(),
                'Description': pathway_name,
                'GeneRatio': f"{overlap_count}/{total_genes}",
                'BgRatio': f"{pathway_size}/{background_size}",
                'pvalue': p_value,
                'p.adjust': p_value,  # Will be corrected later
                'qvalue': p_value,
                'geneID': '/'.join(sorted(overlap)),
                'Count': overlap_count,
                'Enrichment_Score': enrichment_score
            })
        
        df = pd.DataFrame(results)
        if not df.empty:
            # Multiple testing correction
            from statsmodels.stats.multitest import multipletests
            _, padj_values, _, _ = multipletests(df['pvalue'], method='fdr_bh')
            df['p.adjust'] = padj_values
            df['qvalue'] = padj_values
            
            df = df.sort_values('pvalue')
        
        return PathwayResults(
            results=df,
            gene_sets=pathways,
            analysis_method="ORA (Basic)",
            database="Built-in Pathways",
            metadata={
                "n_genes": len(gene_list),
                "n_pathways": len(pathways),
                "significant_pathways": len(df[df['p.adjust'] < 0.05]) if not df.empty else 0
            }
        )
    
    @handle_exception
    def comprehensive_analysis(self, gene_list: List[str],
                             expression_data: Optional[ExpressionData] = None,
                             databases: List[str] = None,
                             analysis_method: str = "ORA",
                             organism: str = "Human (Homo sapiens)",
                             min_gene_set: int = 5,
                             max_gene_set: int = 1000,
                             p_cutoff: float = None) -> PathwayResults:
        """
        Comprehensive pathway analysis using GSEApy with multiple databases
        
        Args:
            gene_list: List of genes to analyze
            expression_data: Optional expression data for GSEA
            databases: List of databases to query
            analysis_method: 'ORA', 'GSEA', or 'Both'
            organism: Organism for analysis
            min_gene_set: Minimum gene set size
            max_gene_set: Maximum gene set size
            p_cutoff: P-value cutoff for significance
            
        Returns:
            PathwayResults container with comprehensive results
        """
        if p_cutoff is None:
            p_cutoff = self.config.get("p_cutoff", 0.05)
        
        if databases is None:
            databases = ["KEGG Pathways", "GO Biological Process"]
        
        # Use GSEApy if available, otherwise fallback to basic analysis
        if self.gseapy_available:
            return self._gseapy_analysis(
                gene_list, expression_data, databases, analysis_method,
                organism, min_gene_set, max_gene_set, p_cutoff
            )
        else:
            logger.info("GSEApy not available. Using basic pathway analysis.")
            return self.basic_enrichment_analysis(gene_list)
    
    def _gseapy_analysis(self, gene_list: List[str],
                        expression_data: Optional[ExpressionData],
                        databases: List[str],
                        analysis_method: str,
                        organism: str,
                        min_gene_set: int,
                        max_gene_set: int,
                        p_cutoff: float) -> PathwayResults:
        """Run comprehensive pathway analysis with GSEApy"""
        try:
            import gseapy as gp
            
            all_results = []
            
            # Map organism names
            organism_map = {
                "Human (Homo sapiens)": "Human",
                "Mouse (Mus musculus)": "Mouse", 
                "Rat (Rattus norvegicus)": "Rat"
            }
            org = organism_map.get(organism, "Human")
            
            # Map database names to GSEApy gene sets
            db_map = {
                "KEGG Pathways": f"KEGG_2021_{org}",
                "Reactome": "Reactome_2022",
                "GO Biological Process": "GO_Biological_Process_2023",
                "GO Molecular Function": "GO_Molecular_Function_2023", 
                "GO Cellular Component": "GO_Cellular_Component_2023",
                "MSigDB Hallmarks": "MSigDB_Hallmark_2020",
                "MSigDB C2 (Curated)": "MSigDB_Curated_2020",
                "WikiPathways": f"WikiPathways_2024_{org}"
            }
            
            logger.info(f"Running GSEApy analysis for {len(databases)} databases")
            
            for db_name in databases:
                if db_name not in db_map:
                    logger.warning(f"Database {db_name} not supported for GSEApy analysis")
                    continue
                
                try:
                    gmt = db_map[db_name]
                    
                    if analysis_method in ["ORA", "Over-Representation Analysis (ORA)", "Both"]:
                        # Run ORA analysis
                        logger.info(f"Running ORA for {db_name}")
                        ora_results = gp.enrichr(
                            gene_list=gene_list,
                            gene_sets=gmt,
                            organism=org,
                            cutoff=p_cutoff,
                            no_plot=True
                        )
                        
                        if not ora_results.results.empty:
                            ora_df = ora_results.results.copy()
                            ora_df['Database'] = db_name
                            ora_df['Method'] = 'ORA'
                            all_results.append(ora_df)
                    
                    if (analysis_method in ["GSEA", "Gene Set Enrichment Analysis (GSEA)", "Both"] 
                        and expression_data is not None):
                        # Run GSEA analysis
                        logger.info(f"Running GSEA for {db_name}")
                        
                        # Prepare expression data for GSEA
                        expr_df = expression_data.data.copy()
                        
                        gsea_results = gp.gsea(
                            data=expr_df,
                            gene_sets=gmt,
                            permutation_num=100,
                            min_size=min_gene_set,
                            max_size=max_gene_set,
                            no_plot=True
                        )
                        
                        if not gsea_results.res2d.empty:
                            gsea_df = gsea_results.res2d.copy()
                            gsea_df['Database'] = db_name
                            gsea_df['Method'] = 'GSEA'
                            all_results.append(gsea_df)
                            
                except Exception as e:
                    logger.warning(f"Could not analyze {db_name}: {e}")
                    continue
            
            # Combine all results
            if all_results:
                combined_results = pd.concat(all_results, ignore_index=True, sort=False)
                combined_results = combined_results.sort_values('Adjusted P-value' if 'Adjusted P-value' in combined_results.columns else 'P-value')
                
                return PathwayResults(
                    results=combined_results,
                    gene_sets={},  # GSEApy handles gene sets internally
                    analysis_method=analysis_method,
                    database=', '.join(databases),
                    metadata={
                        "n_genes": len(gene_list),
                        "databases": databases,
                        "organism": organism,
                        "method": "GSEApy",
                        "significant_pathways": len(combined_results[combined_results.get('Adjusted P-value', combined_results.get('P-value', pd.Series())) < p_cutoff])
                    }
                )
            else:
                logger.warning("No results from GSEApy analysis")
                return self.basic_enrichment_analysis(gene_list)
                
        except ImportError:
            logger.error("GSEApy not available. Install with: pip install gseapy")
            return self.basic_enrichment_analysis(gene_list)
        except Exception as e:
            logger.error(f"GSEApy analysis failed: {e}")
            return self.basic_enrichment_analysis(gene_list)
    
    @handle_exception
    def create_pathway_network(self, pathway_results: PathwayResults,
                              top_pathways: int = 20,
                              color_by: str = "pvalue",
                              layout: str = "spring") -> Dict[str, Any]:
        """
        Create pathway network visualization
        
        Args:
            pathway_results: PathwayResults container
            top_pathways: Number of top pathways to include
            color_by: Column to color nodes by
            layout: Network layout algorithm
            
        Returns:
            Dictionary with network visualization
        """
        # Try R version first, then Python fallback
        if self.is_r_available():
            return self._create_r_pathway_network(pathway_results, top_pathways, color_by)
        else:
            return self._create_python_pathway_network(pathway_results, top_pathways, color_by)
    
    def _create_r_pathway_network(self, pathway_results: PathwayResults,
                                 top_pathways: int, color_by: str) -> Dict[str, Any]:
        """Create pathway network using R"""
        try:
            # Convert results to R
            r_results = self.pandas2ri.py2rpy(pathway_results.results.head(top_pathways))
            self.ro.globalenv['pathway_results'] = r_results
            self.ro.globalenv['top_pathways'] = top_pathways
            self.ro.globalenv['color_by'] = color_by
            
            # Create network using R
            network_result = self.r('''
                create_pathway_network(
                    pathway_results = pathway_results,
                    top_pathways = top_pathways,
                    color_by = color_by
                )
            ''')
            
            return {
                'success': True,
                'method': 'R (igraph)',
                'plot': network_result.rx2('plot'),
                'network': network_result.rx2('network'),
                'nodes': network_result.rx2('nodes'),
                'edges': network_result.rx2('edges')
            }
            
        except Exception as e:
            logger.warning(f"R pathway network failed: {e}")
            return self._create_python_pathway_network(pathway_results, top_pathways, color_by)
    
    def _create_python_pathway_network(self, pathway_results: PathwayResults,
                                     top_pathways: int, color_by: str) -> Dict[str, Any]:
        """Create pathway network using Python/Plotly"""
        px, go = get_plotting_libs()
        if px is None:
            logger.error("Plotly not available for pathway network")
            return {'success': False, 'error': 'Plotting libraries not available'}
        
        try:
            # Get top pathways
            df = pathway_results.results.head(top_pathways)
            
            if df.empty:
                return {'success': False, 'error': 'No pathway results to visualize'}
            
            # Create simple network visualization
            # For now, create a scatter plot with pathway significance
            y_col = color_by if color_by in df.columns else 'pvalue'
            
            if 'pvalue' in df.columns:
                df['neg_log_pvalue'] = -np.log10(df['pvalue'].replace(0, 1e-300))
                y_axis = 'neg_log_pvalue'
                y_title = '-log10(P-value)'
            else:
                y_axis = y_col
                y_title = y_col
            
            fig = px.scatter(
                df.reset_index(),
                x='index',
                y=y_axis,
                size='Count' if 'Count' in df.columns else None,
                color=y_col,
                hover_data=['Description'] if 'Description' in df.columns else None,
                title=f"Pathway Analysis Results (Top {len(df)})",
                labels={'index': 'Pathway Rank', y_axis: y_title}
            )
            
            fig.update_layout(
                xaxis_title="Pathway Rank",
                yaxis_title=y_title,
                showlegend=True
            )
            
            return {
                'success': True,
                'method': 'Python (Plotly)',
                'plot': fig,
                'network': None,
                'nodes': df,
                'edges': None
            }
            
        except Exception as e:
            logger.error(f"Python pathway network creation failed: {e}")
            return {'success': False, 'error': str(e)}
    
    @handle_exception
    def create_dot_plot(self, pathway_results: PathwayResults,
                       top_pathways: int = 20) -> Any:
        """
        Create dot plot for pathway results
        
        Args:
            pathway_results: PathwayResults container
            top_pathways: Number of top pathways to show
            
        Returns:
            Plotly figure
        """
        px, go = get_plotting_libs()
        if px is None:
            logger.error("Plotly not available for dot plot")
            return None
        
        try:
            df = pathway_results.results.head(top_pathways)
            
            if df.empty:
                logger.warning("No pathway results to plot")
                return None
            
            # Prepare data for dot plot
            if 'Description' in df.columns:
                pathway_names = df['Description'].tolist()
            elif 'Term' in df.columns:
                pathway_names = df['Term'].tolist()
            else:
                pathway_names = df.index.tolist()
            
            # Get significance values
            if 'p.adjust' in df.columns:
                p_values = df['p.adjust']
                p_title = 'Adjusted P-value'
            elif 'Adjusted P-value' in df.columns:
                p_values = df['Adjusted P-value']
                p_title = 'Adjusted P-value'
            elif 'pvalue' in df.columns:
                p_values = df['pvalue']
                p_title = 'P-value'
            else:
                p_values = None
                p_title = 'Significance'
            
            # Get gene counts
            if 'Count' in df.columns:
                gene_counts = df['Count']
            elif 'Overlap' in df.columns:
                gene_counts = df['Overlap']
            else:
                gene_counts = None
            
            # Create dot plot
            fig = go.Figure()
            
            if p_values is not None:
                neg_log_p = -np.log10(p_values.replace(0, 1e-300))
                
                fig.add_trace(go.Scatter(
                    x=neg_log_p,
                    y=pathway_names,
                    mode='markers',
                    marker=dict(
                        size=gene_counts if gene_counts is not None else 10,
                        color=neg_log_p,
                        colorscale='Viridis',
                        showscale=True,
                        colorbar=dict(title=f'-log10({p_title})')
                    ),
                    text=[f"Genes: {count}" for count in gene_counts] if gene_counts is not None else None,
                    hovertemplate='%{y}<br>-log10(p): %{x:.3f}<br>%{text}<extra></extra>'
                ))
            
            fig.update_layout(
                title=f"Pathway Enrichment Results (Top {len(df)})",
                xaxis_title=f'-log10({p_title})',
                yaxis_title='Pathways',
                height=max(400, len(pathway_names) * 25),
                margin=dict(l=200)  # Left margin for pathway names
            )
            
            return fig
            
        except Exception as e:
            logger.error(f"Dot plot creation failed: {e}")
            return None
    
    def export_results(self, pathway_results: PathwayResults,
                      output_dir: str = ".") -> Dict[str, str]:
        """
        Export pathway analysis results to files
        
        Args:
            pathway_results: PathwayResults container
            output_dir: Output directory
            
        Returns:
            Dictionary of exported file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        exported_files = {}
        
        try:
            # Export main results
            results_filename = output_path / "pathway_analysis_results.csv"
            pathway_results.results.to_csv(results_filename, index=False)
            exported_files["results"] = str(results_filename)
            
            # Export significant pathways only
            if not pathway_results.results.empty:
                p_col = 'p.adjust' if 'p.adjust' in pathway_results.results.columns else 'pvalue'
                significant = pathway_results.get_significant_pathways()
                if not significant.empty:
                    sig_filename = output_path / "significant_pathways.csv"
                    significant.to_csv(sig_filename, index=False)
                    exported_files["significant"] = str(sig_filename)
            
            # Export gene sets if available
            if pathway_results.gene_sets:
                import json
                geneset_filename = output_path / "gene_sets.json"
                with open(geneset_filename, 'w') as f:
                    json.dump(pathway_results.gene_sets, f, indent=2)
                exported_files["gene_sets"] = str(geneset_filename)
            
            # Export metadata
            metadata_filename = output_path / "pathway_analysis_metadata.json"
            import json
            with open(metadata_filename, 'w') as f:
                json.dump(pathway_results.metadata, f, indent=2, default=str)
            exported_files["metadata"] = str(metadata_filename)
            
            logger.info(f"Exported {len(exported_files)} pathway analysis files")
            return exported_files
            
        except Exception as e:
            logger.error(f"Export failed: {e}")
            return {}


# Create singleton instance
_pathway_engine = None

def get_pathway_engine() -> PathwayEngine:
    """Get the global pathway engine instance"""
    global _pathway_engine
    if _pathway_engine is None:
        _pathway_engine = PathwayEngine()
    return _pathway_engine