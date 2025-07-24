"""
Prairie Genomics Suite - Simple Literature Search

Basic literature search functionality that works without API keys or complex dependencies.
This serves as a fallback when the advanced LiteratureIntelligenceEngine is not available.

Features:
- Simple PubMed search using public E-utilities API (no key required)
- Basic result parsing and display
- Works with minimal dependencies
- Clear error handling and user guidance
"""

import requests
import pandas as pd
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Any
import logging
import time
from urllib.parse import quote

logger = logging.getLogger(__name__)

class SimpleLiteratureSearch:
    """
    Simple literature search using PubMed E-utilities API
    
    This is a fallback search engine that works without API keys
    and provides basic literature search functionality.
    """
    
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.rate_limit_delay = 0.5  # Respect NCBI rate limits
    
    def search_pubmed_simple(self, query_terms: List[str], 
                           condition_context: str = "",
                           max_results_per_term: int = 5) -> Dict[str, List[Dict]]:
        """
        Simple PubMed search for multiple terms
        
        Args:
            query_terms: List of genes/terms to search
            condition_context: Disease/condition context
            max_results_per_term: Maximum results per term
            
        Returns:
            Dictionary with search results for each term
        """
        logger.info(f"Starting simple PubMed search for {len(query_terms)} terms")
        
        results = {}
        
        for term in query_terms[:10]:  # Limit to 10 terms to avoid rate limits
            try:
                papers = self._search_single_term(term, condition_context, max_results_per_term)
                results[term] = papers
                
                # Rate limiting
                time.sleep(self.rate_limit_delay)
                
            except Exception as e:
                logger.warning(f"Search failed for term {term}: {e}")
                results[term] = []
        
        return results
    
    def _search_single_term(self, term: str, condition: str, max_results: int) -> List[Dict]:
        """Search PubMed for a single term"""
        try:
            # Construct search query
            search_terms = [f'"{term}"[Title/Abstract]']
            if condition and condition.strip():
                search_terms.append(f'"{condition}"[Title/Abstract]')
            
            query = " AND ".join(search_terms)
            
            # Step 1: Search for PMIDs
            search_url = f"{self.base_url}esearch.fcgi"
            search_params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "retmode": "json"
            }
            
            response = requests.get(search_url, params=search_params, timeout=10)
            response.raise_for_status()
            
            search_data = response.json()
            pmids = search_data.get("esearchresult", {}).get("idlist", [])
            
            if not pmids:
                return []
            
            # Step 2: Fetch article details
            return self._fetch_article_details(pmids)
            
        except Exception as e:
            logger.error(f"PubMed search failed for {term}: {e}")
            return []
    
    def _fetch_article_details(self, pmids: List[str]) -> List[Dict]:
        """Fetch detailed information for PubMed IDs"""
        try:
            # Fetch article summaries
            summary_url = f"{self.base_url}esummary.fcgi"
            summary_params = {
                "db": "pubmed",
                "id": ",".join(pmids),
                "retmode": "json"
            }
            
            response = requests.get(summary_url, params=summary_params, timeout=10)
            response.raise_for_status()
            
            summary_data = response.json()
            results = []
            
            for pmid in pmids:
                try:
                    article_data = summary_data.get("result", {}).get(pmid, {})
                    
                    if article_data:
                        # Extract basic information
                        authors = self._extract_authors(article_data.get("authors", []))
                        
                        paper = {
                            "title": article_data.get("title", "No title available"),
                            "authors": authors,
                            "journal": article_data.get("fulljournalname", "Unknown journal"),
                            "year": article_data.get("pubdate", "Unknown year")[:4] if article_data.get("pubdate") else "Unknown year",
                            "pmid": pmid,
                            "abstract": "Abstract not available in summary",  # Would need separate API call
                            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            "relevance_score": 1.0,  # Simple fixed score
                            "database": "PubMed (Simple)"
                        }
                        
                        results.append(paper)
                        
                except Exception as e:
                    logger.warning(f"Error processing article {pmid}: {e}")
                    continue
            
            return results
            
        except Exception as e:
            logger.error(f"Failed to fetch article details: {e}")
            return []
    
    def _extract_authors(self, authors_list: List[Dict]) -> str:
        """Extract and format author names"""
        try:
            if not authors_list:
                return "Authors not available"
            
            author_names = []
            for author in authors_list[:3]:  # Limit to first 3 authors
                name = author.get("name", "")
                if name:
                    author_names.append(name)
            
            if len(authors_list) > 3:
                author_names.append("et al.")
            
            return ", ".join(author_names) if author_names else "Authors not available"
            
        except Exception:
            return "Authors not available"
    
    def create_example_results(self, query_terms: List[str]) -> Dict[str, List[Dict]]:
        """
        Create example literature results for demonstration
        
        This provides working examples when real search is not available
        """
        example_papers = {
            "TP53": [
                {
                    "title": "p53 mutations in human cancers: functional consequences and therapeutic opportunities",
                    "authors": "Smith, J., Johnson, A., et al.",
                    "journal": "Nature Reviews Cancer",
                    "year": "2023",
                    "pmid": "example_001",
                    "abstract": "p53 is one of the most frequently mutated genes in human cancer...",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/example_001/",
                    "relevance_score": 2.8,
                    "database": "Example Data"
                },
                {
                    "title": "Therapeutic targeting of p53 in cancer treatment",
                    "authors": "Wilson, M., Davis, K.",
                    "journal": "Cell",
                    "year": "2023",
                    "pmid": "example_002",
                    "abstract": "Recent advances in understanding p53 biology have opened new therapeutic avenues...",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/example_002/",
                    "relevance_score": 2.5,
                    "database": "Example Data"
                }
            ],
            "BRCA1": [
                {
                    "title": "BRCA1 mutations and breast cancer susceptibility",
                    "authors": "Brown, L., Taylor, R., et al.",
                    "journal": "New England Journal of Medicine",
                    "year": "2023",
                    "pmid": "example_003",
                    "abstract": "BRCA1 germline mutations confer high risk for breast and ovarian cancers...",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/example_003/",
                    "relevance_score": 2.9,
                    "database": "Example Data"
                }
            ],
            "MYC": [
                {
                    "title": "MYC oncogene in cancer progression and therapy",
                    "authors": "Anderson, P., Miller, S.",
                    "journal": "Science",
                    "year": "2023",
                    "pmid": "example_004",
                    "abstract": "The MYC oncogene plays central roles in cancer development and progression...",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/example_004/",
                    "relevance_score": 2.7,
                    "database": "Example Data"
                }
            ]
        }
        
        results = {}
        for term in query_terms:
            term_upper = term.upper()
            if term_upper in example_papers:
                results[term] = example_papers[term_upper]
            else:
                # Create a generic example for any gene
                results[term] = [{
                    "title": f"Functional analysis of {term} in human disease",
                    "authors": "Example, A., Research, B.",
                    "journal": "Journal of Molecular Biology",
                    "year": "2023",
                    "pmid": f"example_{term.lower()}",
                    "abstract": f"This study investigates the role of {term} in disease mechanisms...",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/example_{term.lower()}/",
                    "relevance_score": 2.0,
                    "database": "Example Data"
                }]
        
        return results
    
    def get_status(self) -> Dict[str, Any]:
        """Get status of the simple literature search"""
        try:
            # Test connection to PubMed
            test_url = f"{self.base_url}einfo.fcgi"
            response = requests.get(test_url, timeout=5)
            pubmed_available = response.status_code == 200
        except:
            pubmed_available = False
        
        return {
            "simple_search_available": True,
            "pubmed_api_available": pubmed_available,
            "requires_api_keys": False,
            "features": [
                "Basic PubMed search",
                "Example literature results",
                "No API keys required"
            ]
        }


# Global instance
_simple_search = None

def get_simple_literature_search() -> SimpleLiteratureSearch:
    """Get the global simple literature search instance"""
    global _simple_search
    if _simple_search is None:
        _simple_search = SimpleLiteratureSearch()
    return _simple_search