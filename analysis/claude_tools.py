"""
Prairie Genomics Suite - Claude Tool Calling Interface

This module provides Claude tool calling integration for conversational literature search
and genomics analysis. Based on the pattern from Perplexity article for elegant AI-driven
research discovery.

Features:
- Natural language literature search
- Context-aware genomics queries  
- Automatic synthesis and analysis
- Integration with existing analysis engines
- Conversational research workflow
"""

import json
import logging
from typing import Dict, List, Any, Optional
from pathlib import Path

# Import our existing engines
from .simple_literature_search import get_simple_literature_search
from .literature_engine import get_literature_engine
from core.session_manager import get_session_manager
from config import get_config

logger = logging.getLogger(__name__)

# =============================================================================
# CLAUDE TOOL SCHEMAS
# =============================================================================

LITERATURE_SEARCH_TOOL = {
    "name": "pubmed_literature_search",
    "description": "Searches PubMed and other scientific databases for peer-reviewed biomedical literature. Returns relevant papers with titles, abstracts, and links. Ideal for finding research on specific genes, diseases, or biological processes.",
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string", 
                "description": "Plain-language search terms (e.g., 'BRCA1 mutations breast cancer', 'immune checkpoint inhibitors', 'TP53 pathway')"
            },
            "condition_context": {
                "type": "string",
                "description": "Disease or biological context to focus the search (e.g., 'cancer', 'immunotherapy', 'aging')",
                "default": ""
            },
            "max_results": {
                "type": "integer",
                "description": "Number of papers to return (1-20)",
                "default": 5,
                "minimum": 1,
                "maximum": 20
            },
            "include_recent_only": {
                "type": "boolean",
                "description": "Focus on papers from last 3 years",
                "default": False
            }
        },
        "required": ["query"]
    }
}

GENOMICS_CONTEXT_TOOL = {
    "name": "get_genomics_analysis_context",
    "description": "Retrieves current genomics analysis results (DESeq2, pathway analysis, etc.) to provide context for literature searches. Use this to understand what genes or pathways the user is working with.",
    "input_schema": {
        "type": "object",
        "properties": {
            "analysis_type": {
                "type": "string",
                "enum": ["deseq2", "pathway", "survival", "all"],
                "description": "Type of analysis results to retrieve",
                "default": "all"
            },
            "summary_only": {
                "type": "boolean", 
                "description": "Return only summary statistics, not full results",
                "default": True
            }
        },
        "required": []
    }
}

LITERATURE_SYNTHESIS_TOOL = {
    "name": "synthesize_literature_with_genomics",
    "description": "Connects literature findings with current genomics analysis results. Identifies relevant papers for significant genes, pathways, or findings from DESeq2/pathway analysis.",
    "input_schema": {
        "type": "object",
        "properties": {
            "analysis_context": {
                "type": "string",
                "description": "Description of the genomics analysis context (e.g., 'upregulated genes in cancer vs normal', 'immune pathways enriched')"
            },
            "focus_areas": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Specific areas to focus literature search (e.g., ['therapeutic targets', 'mechanisms', 'clinical outcomes'])",
                "default": ["mechanisms", "clinical relevance"]
            },
            "gene_limit": {
                "type": "integer",
                "description": "Maximum number of top genes to search literature for",
                "default": 10,
                "minimum": 1,
                "maximum": 50
            }
        },
        "required": ["analysis_context"]
    }
}

# All available tools
CLAUDE_TOOLS = [
    LITERATURE_SEARCH_TOOL,
    GENOMICS_CONTEXT_TOOL,
    LITERATURE_SYNTHESIS_TOOL
]

# =============================================================================
# TOOL IMPLEMENTATION FUNCTIONS
# =============================================================================

def pubmed_literature_search(query: str, condition_context: str = "", 
                            max_results: int = 5, include_recent_only: bool = False) -> Dict[str, Any]:
    """
    Execute literature search using our existing engines
    
    Returns structured data that Claude can interpret and synthesize
    """
    try:
        logger.info(f"Claude tool: Literature search for '{query}'")
        
        # Use our existing literature engine with fallback
        lit_engine = get_literature_engine()
        
        # Prepare search parameters
        search_terms = [term.strip() for term in query.split(',')][:10]  # Limit terms
        
        # Add publication year filter if recent only
        publication_years = (2021, 2024) if include_recent_only else None
        
        # Execute search
        results = lit_engine.intelligent_search(
            query_terms=search_terms,
            condition_context=condition_context,
            max_results_per_term=max_results,
            databases=['pubmed', 'semantic_scholar'],
            publication_years=publication_years,
            include_reviews=True,
            semantic_expansion=True
        )
        
        # Format results for Claude
        formatted_results = {
            "search_query": query,
            "total_papers": results.total_papers,
            "papers_found": [],
            "summary": results.summary,
            "search_status": "success"
        }
        
        # Extract paper details
        for term, papers in results.results.items():
            for paper in papers[:max_results]:  # Limit per term
                paper_info = {
                    "title": paper.title,
                    "authors": paper.authors,
                    "journal": paper.journal,
                    "year": paper.year,
                    "abstract": paper.abstract[:300] + "..." if len(paper.abstract) > 300 else paper.abstract,
                    "pmid": paper.pmid,
                    "url": paper.url,
                    "relevance_score": paper.relevance_score,
                    "search_term": term
                }
                formatted_results["papers_found"].append(paper_info)
        
        # Sort by relevance
        formatted_results["papers_found"].sort(key=lambda x: x["relevance_score"], reverse=True)
        
        return formatted_results
        
    except Exception as e:
        logger.error(f"Literature search tool error: {e}")
        return {
            "search_query": query,
            "error": f"Literature search failed: {str(e)}",
            "search_status": "error",
            "suggestion": "Try a simpler query or check if literature dependencies are installed"
        }

def get_genomics_analysis_context(analysis_type: str = "all", summary_only: bool = True) -> Dict[str, Any]:
    """
    Retrieve current genomics analysis context for Claude
    
    This helps Claude understand what the user is working with
    """
    try:
        logger.info(f"Claude tool: Getting genomics context for {analysis_type}")
        
        session = get_session_manager()
        context = {
            "session_status": "active",
            "available_data": {},
            "analysis_results": {},
            "recommendations": []
        }
        
        # Check available data
        if session.has_data("expression_data"):
            expr_data = session.get_data("expression_data")
            context["available_data"]["expression"] = {
                "n_genes": expr_data.n_genes,
                "n_samples": expr_data.n_samples,
                "sample_names": expr_data.samples[:5] + ["..."] if len(expr_data.samples) > 5 else expr_data.samples
            }
        
        if session.has_data("clinical_data"):
            clinical_data = session.get_data("clinical_data")
            context["available_data"]["clinical"] = {
                "n_samples": clinical_data.n_samples,
                "variables": clinical_data.variables[:10]  # First 10 variables
            }
        
        # Check analysis results
        if analysis_type in ["deseq2", "all"] and session.has_analysis_results("deseq2"):
            deseq2_results = session.get_analysis_results("deseq2")
            
            if summary_only:
                context["analysis_results"]["deseq2"] = {
                    "total_genes": deseq2_results.n_genes,
                    "significant_genes": deseq2_results.n_significant,
                    "upregulated": deseq2_results.n_upregulated,
                    "downregulated": deseq2_results.n_downregulated,
                    "top_significant_genes": [],
                    "analysis_complete": True
                }
                
                # Get top genes for literature search
                if hasattr(deseq2_results, 'results') and not deseq2_results.results.empty:
                    top_genes = deseq2_results.results.head(10).index.tolist()
                    context["analysis_results"]["deseq2"]["top_significant_genes"] = top_genes
                    
                    context["recommendations"].append(
                        f"Consider literature search for top significant genes: {', '.join(top_genes[:5])}"
                    )
            else:
                # Full results (limited for Claude)
                context["analysis_results"]["deseq2"] = {
                    "summary": "Full DESeq2 results available",
                    "significant_genes": deseq2_results.n_significant,
                    "note": "Use summary_only=True for detailed gene lists"
                }
        
        if analysis_type in ["pathway", "all"] and session.has_analysis_results("pathway"):
            pathway_results = session.get_analysis_results("pathway")
            context["analysis_results"]["pathway"] = {
                "pathways_found": len(pathway_results.results) if hasattr(pathway_results, 'results') else 0,
                "analysis_complete": True
            }
            
            context["recommendations"].append(
                "Consider literature search for enriched pathways and their key genes"
            )
        
        if analysis_type in ["literature", "all"] and session.has_analysis_results("literature"):
            context["analysis_results"]["literature"] = {
                "previous_search_complete": True,
                "papers_found": "Previous literature search results available"
            }
        
        # Add recommendations
        if not context["analysis_results"]:
            context["recommendations"].append(
                "No analysis results available yet. Run DESeq2 or pathway analysis first."
            )
        else:
            context["recommendations"].append(
                "Current analysis results are available for literature synthesis"
            )
        
        return context
        
    except Exception as e:
        logger.error(f"Genomics context tool error: {e}")
        return {
            "error": f"Failed to get genomics context: {str(e)}",
            "session_status": "error"
        }

def synthesize_literature_with_genomics(analysis_context: str, 
                                       focus_areas: List[str] = None,
                                       gene_limit: int = 10) -> Dict[str, Any]:
    """
    Connect literature findings with genomics analysis results
    
    This tool helps Claude find relevant literature for the user's specific analysis
    """
    try:
        logger.info(f"Claude tool: Synthesizing literature with genomics context")
        
        if focus_areas is None:
            focus_areas = ["mechanisms", "clinical relevance"]
        
        session = get_session_manager()
        synthesis_results = {
            "analysis_context": analysis_context,
            "focus_areas": focus_areas,
            "literature_searches_performed": [],
            "key_findings": [],
            "recommendations": [],
            "synthesis_status": "success"
        }
        
        # Get genomics context
        genomics_context = get_genomics_analysis_context("all", summary_only=True)
        
        # Extract genes to search
        genes_to_search = []
        
        if "deseq2" in genomics_context.get("analysis_results", {}):
            deseq2_info = genomics_context["analysis_results"]["deseq2"]
            top_genes = deseq2_info.get("top_significant_genes", [])
            genes_to_search.extend(top_genes[:gene_limit])
        
        if not genes_to_search:
            return {
                "analysis_context": analysis_context,
                "error": "No significant genes found in current analysis",
                "recommendation": "Run DESeq2 analysis first to identify genes for literature search",
                "synthesis_status": "no_data"
            }
        
        # Perform literature searches for top genes
        search_queries = []
        
        for focus_area in focus_areas:
            if focus_area.lower() in ["mechanisms", "mechanism"]:
                query = f"{', '.join(genes_to_search[:5])} mechanism pathway function"
            elif focus_area.lower() in ["clinical", "clinical relevance", "therapeutic"]:
                query = f"{', '.join(genes_to_search[:5])} therapeutic target clinical outcome"
            elif focus_area.lower() in ["disease", "pathology"]:
                query = f"{', '.join(genes_to_search[:5])} disease pathology biomarker"
            else:
                query = f"{', '.join(genes_to_search[:5])} {focus_area}"
            
            search_queries.append(query)
        
        # Execute searches
        for query in search_queries:
            search_result = pubmed_literature_search(
                query=query,
                condition_context=analysis_context,
                max_results=3,  # Limited for synthesis
                include_recent_only=True
            )
            
            if search_result.get("search_status") == "success":
                synthesis_results["literature_searches_performed"].append({
                    "query": query,
                    "papers_found": len(search_result.get("papers_found", [])),
                    "top_papers": search_result.get("papers_found", [])[:2]  # Top 2 per search
                })
        
        # Generate key findings
        total_papers = sum(len(search["top_papers"]) for search in synthesis_results["literature_searches_performed"])
        
        synthesis_results["key_findings"] = [
            f"Found {total_papers} relevant papers for {len(genes_to_search)} significant genes",
            f"Literature search focused on: {', '.join(focus_areas)}",
            f"Top genes searched: {', '.join(genes_to_search[:5])}"
        ]
        
        # Generate recommendations
        synthesis_results["recommendations"] = [
            "Review papers for novel therapeutic targets among your significant genes",
            "Look for mechanistic insights that could explain your differential expression results",
            "Consider follow-up experiments suggested in recent publications",
            "Check if any of your significant genes are in current clinical trials"
        ]
        
        return synthesis_results
        
    except Exception as e:
        logger.error(f"Literature synthesis tool error: {e}")
        return {
            "analysis_context": analysis_context,
            "error": f"Literature synthesis failed: {str(e)}",
            "synthesis_status": "error"
        }

# =============================================================================
# TOOL EXECUTION DISPATCHER
# =============================================================================

def execute_claude_tool(tool_name: str, tool_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Execute a Claude tool by name with input parameters
    
    This is the main dispatcher function that Claude will call
    """
    try:
        if tool_name == "pubmed_literature_search":
            return pubmed_literature_search(**tool_input)
        
        elif tool_name == "get_genomics_analysis_context":
            return get_genomics_analysis_context(**tool_input)
        
        elif tool_name == "synthesize_literature_with_genomics":
            return synthesize_literature_with_genomics(**tool_input)
        
        else:
            return {
                "error": f"Unknown tool: {tool_name}",
                "available_tools": [tool["name"] for tool in CLAUDE_TOOLS]
            }
            
    except Exception as e:
        logger.error(f"Tool execution error for {tool_name}: {e}")
        return {
            "error": f"Tool execution failed: {str(e)}",
            "tool_name": tool_name,
            "input_received": tool_input
        }

# =============================================================================
# TOOL REGISTRATION AND UTILITIES
# =============================================================================

def get_available_tools() -> List[Dict[str, Any]]:
    """Get list of all available Claude tools"""
    return CLAUDE_TOOLS

def get_tool_status() -> Dict[str, Any]:
    """Get status of tool dependencies"""
    status = {
        "tools_available": len(CLAUDE_TOOLS),
        "dependencies": {
            "literature_engine": False,
            "session_manager": False,
            "simple_search": False
        },
        "ready": False
    }
    
    try:
        # Test literature engine
        lit_engine = get_literature_engine()
        status["dependencies"]["literature_engine"] = True
    except:
        pass
    
    try:
        # Test session manager
        session = get_session_manager()
        status["dependencies"]["session_manager"] = True
    except:
        pass
    
    try:
        # Test simple search
        simple_search = get_simple_literature_search()
        status["dependencies"]["simple_search"] = True
    except:
        pass
    
    status["ready"] = all(status["dependencies"].values())
    
    return status