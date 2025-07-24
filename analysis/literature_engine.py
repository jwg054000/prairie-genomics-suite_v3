"""
Prairie Genomics Suite - Advanced Literature Intelligence Engine

This module provides cutting-edge literature search and synthesis using multiple
APIs and LLM integration for scientific rigor and discovery acceleration.

Features:
- Multi-API literature aggregation (PubMed, Semantic Scholar, arXiv)
- LLM-powered intelligent analysis (OpenAI GPT-4, Claude)
- Genomics-specific context understanding
- Semantic search with vector embeddings
- Real-time research synthesis and insights
- Citation network analysis and discovery
- Automated research summaries with scientific integrity

"In science, there is only physics; all the rest is stamp collecting." - Rutherford
But we're building the machine that makes sense of all the stamps! 🚀

Author: Prairie Genomics Team
"""

import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
import warnings
import logging
import time
import json
import requests
from datetime import datetime
import hashlib
import re

# Import core modules
from core.data_models import LiteratureResults, LiteratureSearchResult
from core.utils import validate_dataframe, handle_exception, create_download_link
from config import get_config, is_feature_enabled, validate_api_keys

# Set up logging
logger = logging.getLogger(__name__)

# Lazy imports for optional dependencies
def get_optional_libs():
    """Lazy import optional libraries"""
    libs = {}
    try:
        import openai
        libs['openai'] = openai
    except ImportError:
        libs['openai'] = None
    
    try:
        import anthropic
        libs['anthropic'] = anthropic
    except ImportError:
        libs['anthropic'] = None
    
    try:
        from pymed import PubMed
        libs['pymed'] = PubMed
    except ImportError:
        libs['pymed'] = None
    
    try:
        from semanticscholar import SemanticScholar
        libs['semanticscholar'] = SemanticScholar
    except ImportError:
        libs['semanticscholar'] = None
    
    try:
        import arxiv
        libs['arxiv'] = arxiv
    except ImportError:
        libs['arxiv'] = None
    
    return libs


class LiteratureIntelligenceEngine:
    """
    Advanced Literature Intelligence Engine with Multi-API and LLM Integration
    
    This is the culmination of modern literature search - combining traditional
    databases with AI-powered analysis for unprecedented research insights.
    Built for scientific rigor and discovery acceleration.
    """
    
    def __init__(self):
        """Initialize the advanced literature engine"""
        self.config = get_config("literature")
        self.api_status = validate_api_keys()
        self.libs = get_optional_libs()
        
        # Initialize API clients
        self.pubmed = None
        self.semantic_scholar = None
        self.openai_client = None
        self.anthropic_client = None
        
        # Cache for performance
        self.search_cache = {}
        self.synthesis_cache = {}
        
        # Initialize all available services
        self._initialize_apis()
    
    def _initialize_apis(self):
        """Initialize all available API clients with proper error handling"""
        try:
            # PubMed via PyMed
            if self.libs['pymed']:
                try:
                    self.pubmed = self.libs['pymed'](
                        tool=self.config['apis']['pubmed']['tool'],
                        email=self.config['apis']['pubmed']['email']
                    )
                    logger.info("✅ PubMed API initialized")
                except Exception as e:
                    logger.warning(f"PubMed initialization failed: {e}")
                    self.pubmed = None
            else:
                logger.info("📦 PyMed not available - install with: pip install pymed")
            
            # Semantic Scholar
            if self.libs['semanticscholar']:
                try:
                    self.semantic_scholar = self.libs['semanticscholar']()
                    logger.info("✅ Semantic Scholar API initialized")
                except Exception as e:
                    logger.warning(f"Semantic Scholar initialization failed: {e}")
                    self.semantic_scholar = None
            else:
                logger.info("📦 Semantic Scholar not available - install with: pip install semanticscholar")
            
            # OpenAI
            if self.libs['openai'] and self.api_status['openai']:
                try:
                    self.openai_client = self.libs['openai'].OpenAI()
                    logger.info("✅ OpenAI API initialized")
                except Exception as e:
                    logger.warning(f"OpenAI initialization failed: {e}")
                    self.openai_client = None
            elif self.libs['openai']:
                logger.info("🔑 OpenAI available but no API key - set OPENAI_API_KEY environment variable")
            else:
                logger.info("📦 OpenAI not available - install with: pip install openai")
            
            # Anthropic Claude
            if self.libs['anthropic'] and self.api_status['anthropic']:
                try:
                    self.anthropic_client = self.libs['anthropic'].Anthropic()
                    logger.info("✅ Anthropic Claude API initialized")
                except Exception as e:
                    logger.warning(f"Anthropic initialization failed: {e}")
                    self.anthropic_client = None
            elif self.libs['anthropic']:
                logger.info("🔑 Anthropic available but no API key - set ANTHROPIC_API_KEY environment variable")
            else:
                logger.info("📦 Anthropic not available - install with: pip install anthropic")
                
        except Exception as e:
            logger.error(f"API initialization error: {e}")
        
        # Initialize fallback if no advanced APIs available
        if not any([self.pubmed, self.semantic_scholar]):
            logger.info("🔄 No advanced APIs available, fallback search will be used")
    
    def get_status(self) -> Dict[str, Any]:
        """Get comprehensive status of the literature engine"""
        return {
            "apis_available": {
                "pubmed": self.pubmed is not None,
                "semantic_scholar": self.semantic_scholar is not None,
                "arxiv": self.libs['arxiv'] is not None,
                "openai": self.openai_client is not None,
                "anthropic": self.anthropic_client is not None
            },
            "api_status": self.api_status,
            "genomics_context": self.config['genomics_context'],
            "llm_config": self.config['llm'],
            "cache_size": len(self.search_cache)
        }
    
    @handle_exception
    def intelligent_search(self, query_terms: List[str],
                         condition_context: str = "",
                         max_results_per_term: int = 10,
                         databases: List[str] = None,
                         publication_years: Optional[Tuple[int, int]] = None,
                         include_reviews: bool = True,
                         semantic_expansion: bool = True) -> LiteratureResults:
        """
        Perform intelligent multi-database literature search with LLM enhancement
        
        Args:
            query_terms: List of genes/terms to search
            condition_context: Disease/condition context
            max_results_per_term: Maximum results per term per database
            databases: List of databases to search ('pubmed', 'semantic_scholar', 'arxiv')
            publication_years: Optional year range (start, end)
            include_reviews: Whether to include review articles
            semantic_expansion: Whether to use LLM for query expansion
            
        Returns:
            LiteratureResults with comprehensive search results and synthesis
        """
        if databases is None:
            databases = ['pubmed', 'semantic_scholar']
        
        logger.info(f"🔬 Starting intelligent search for {len(query_terms)} terms across {len(databases)} databases")
        
        # Check if any APIs are available, if not use fallback
        apis_available = any([self.pubmed, self.semantic_scholar, self.libs['arxiv']])
        
        if not apis_available:
            logger.info("🔄 No advanced APIs available, using fallback search")
            return self._fallback_search(query_terms, condition_context, max_results_per_term)
        
        # Step 1: Expand queries using LLM if available
        expanded_queries = query_terms.copy()
        if semantic_expansion and (self.openai_client or self.anthropic_client):
            expanded_queries = self._expand_queries_with_llm(query_terms, condition_context)
        
        # Step 2: Search across all databases
        all_results = {}
        search_metadata = {
            "databases": databases,
            "query_terms": query_terms,
            "expanded_queries": expanded_queries,
            "condition_context": condition_context,
            "search_timestamp": datetime.now(),
            "total_papers": 0
        }
        
        for term in query_terms:
            all_results[term] = []
            
            # Search PubMed
            if 'pubmed' in databases and self.pubmed:
                pubmed_results = self._search_pubmed(
                    term, condition_context, max_results_per_term,
                    publication_years, include_reviews
                )
                all_results[term].extend(pubmed_results)
            
            # Search Semantic Scholar
            if 'semantic_scholar' in databases and self.semantic_scholar:
                ss_results = self._search_semantic_scholar(
                    term, condition_context, max_results_per_term,
                    publication_years
                )
                all_results[term].extend(ss_results)
            
            # Search arXiv for preprints
            if 'arxiv' in databases and self.libs['arxiv']:
                arxiv_results = self._search_arxiv(
                    term, condition_context, max_results_per_term
                )
                all_results[term].extend(arxiv_results)
            
            # Remove duplicates and rank by relevance
            all_results[term] = self._deduplicate_and_rank(all_results[term], term, condition_context)
            search_metadata["total_papers"] += len(all_results[term])
        
        # Step 3: LLM-powered synthesis and insights
        synthesis = {}
        if self.openai_client or self.anthropic_client:
            synthesis = self._synthesize_research_with_llm(all_results, condition_context, query_terms)
        
        # Step 4: Create comprehensive results
        literature_results = LiteratureResults(
            results={term: [self._dict_to_result(paper) for paper in papers] 
                    for term, papers in all_results.items()},
            summary=self._create_enhanced_summary(all_results, synthesis),
            search_terms=query_terms,
            search_parameters={
                "databases": databases,
                "condition_context": condition_context,
                "max_results_per_term": max_results_per_term,
                "publication_years": publication_years,
                "include_reviews": include_reviews,
                "semantic_expansion": semantic_expansion,
                "metadata": search_metadata
            }
        )
        
        # Cache results for performance
        cache_key = self._create_cache_key(query_terms, condition_context, databases)
        self.search_cache[cache_key] = literature_results
        
        logger.info(f"🎉 Intelligent search completed: {literature_results.total_papers} papers found")
        return literature_results
    
    def _expand_queries_with_llm(self, query_terms: List[str], condition_context: str) -> List[str]:
        """Use LLM to expand queries with related terms"""
        try:
            # Create context-aware prompt
            genomics_context = self.config['genomics_context']
            
            prompt = f"""You are a genomics research expert. Given these genes/terms: {query_terms}
            
            Context: {condition_context}
            
            Expand this search with 2-3 highly relevant scientific terms for each original term.
            Focus on:
            - Gene aliases and synonyms
            - Related pathway members
            - Functional protein families
            - Disease associations
            
            Return only the expanded terms as a simple comma-separated list, no explanations."""
            
            if self.openai_client:
                response = self.openai_client.chat.completions.create(
                    model=self.config['llm']['openai']['model'],
                    messages=[{"role": "user", "content": prompt}],
                    max_tokens=500,
                    temperature=0.1
                )
                expanded_text = response.choices[0].message.content
            elif self.anthropic_client:
                response = self.anthropic_client.messages.create(
                    model=self.config['llm']['anthropic']['model'],
                    max_tokens=500,
                    temperature=0.1,
                    messages=[{"role": "user", "content": prompt}]
                )
                expanded_text = response.content[0].text
            else:
                return query_terms
            
            # Parse expanded terms
            expanded_terms = [term.strip() for term in expanded_text.split(',') if term.strip()]
            
            # Combine original and expanded, remove duplicates
            all_terms = list(set(query_terms + expanded_terms))
            logger.info(f"🧠 LLM expanded {len(query_terms)} terms to {len(all_terms)} terms")
            return all_terms
            
        except Exception as e:
            logger.warning(f"LLM query expansion failed: {e}")
            return query_terms
    
    def _search_pubmed(self, term: str, condition: str, max_results: int,
                      publication_years: Optional[Tuple[int, int]], 
                      include_reviews: bool) -> List[Dict[str, Any]]:
        """Search PubMed with advanced query construction"""
        try:
            # Build sophisticated query
            query_parts = [f'"{term}"[Title/Abstract]']
            
            if condition and condition.strip():
                query_parts.append(f'"{condition}"[Title/Abstract]')
            
            query = " AND ".join(query_parts)
            
            # Add publication year filter
            if publication_years:
                start_year, end_year = publication_years
                query += f' AND ("{start_year}"[Date - Publication] : "{end_year}"[Date - Publication])'
            
            # Add review filter
            if not include_reviews:
                query += ' NOT "review"[Publication Type]'
            
            # Execute search
            articles = self.pubmed.query(query, max_results=max_results)
            
            results = []
            for article in articles:
                try:
                    # Calculate genomics-specific relevance
                    relevance = self._calculate_genomics_relevance(article, term, condition)
                    
                    paper = {
                        'title': article.title or "No title available",
                        'authors': self._extract_authors(article),
                        'journal': article.journal or "Unknown journal",
                        'year': self._extract_year(article),
                        'pmid': article.pubmed_id or "Unknown",
                        'abstract': article.abstract or "No abstract available",
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{article.pubmed_id}/" if article.pubmed_id else None,
                        'keywords': article.keywords or [],
                        'relevance_score': relevance,
                        'database': 'PubMed',
                        'publication_type': self._extract_publication_type(article)
                    }
                    results.append(paper)
                except Exception as e:
                    logger.warning(f"Error processing PubMed article: {e}")
                    continue
            
            return results
            
        except Exception as e:
            logger.error(f"PubMed search failed for {term}: {e}")
            return []
    
    def _search_semantic_scholar(self, term: str, condition: str, max_results: int,
                                publication_years: Optional[Tuple[int, int]]) -> List[Dict[str, Any]]:
        """Search Semantic Scholar for comprehensive academic coverage"""
        try:
            # Construct query
            query = term
            if condition and condition.strip():
                query += f" {condition}"
            
            # Add genomics-specific terms for better relevance
            genomics_terms = self.config['genomics_context']['conditions']
            if any(gterm in query.lower() for gterm in genomics_terms):
                query += " genomics"
            
            # Search Semantic Scholar
            results = self.semantic_scholar.search_paper(
                query, 
                limit=max_results,
                fields=['title', 'abstract', 'authors', 'venue', 'year', 'citationCount', 'url']
            )
            
            papers = []
            for paper in results:
                try:
                    # Filter by publication year if specified
                    if publication_years and paper.year:
                        if paper.year < publication_years[0] or paper.year > publication_years[1]:
                            continue
                    
                    # Calculate relevance
                    relevance = self._calculate_semantic_relevance(paper, term, condition)
                    
                    paper_dict = {
                        'title': paper.title or "No title available",
                        'authors': ', '.join([author.name for author in (paper.authors or [])]),
                        'journal': paper.venue or "Unknown venue",
                        'year': str(paper.year) if paper.year else "Unknown year",
                        'pmid': paper.paperId or "Unknown",
                        'abstract': paper.abstract or "No abstract available",
                        'url': paper.url,
                        'keywords': [],
                        'relevance_score': relevance,
                        'database': 'Semantic Scholar',
                        'citation_count': paper.citationCount or 0,
                        'publication_type': 'Research Article'
                    }
                    papers.append(paper_dict)
                    
                except Exception as e:
                    logger.warning(f"Error processing Semantic Scholar paper: {e}")
                    continue
            
            return papers
            
        except Exception as e:
            logger.error(f"Semantic Scholar search failed for {term}: {e}")
            return []
    
    def _search_arxiv(self, term: str, condition: str, max_results: int) -> List[Dict[str, Any]]:
        """Search arXiv for preprints and cutting-edge research"""
        try:
            arxiv = self.libs['arxiv']
            
            # Construct query for arXiv
            query = f"all:{term}"
            if condition and condition.strip():
                query += f" AND all:{condition}"
            
            # Focus on relevant categories
            categories = ["q-bio", "stat", "cs", "physics"]  # Relevant for genomics
            
            # Search arXiv
            search = arxiv.Search(
                query=query,
                max_results=max_results,
                sort_by=arxiv.SortCriterion.Relevance
            )
            
            papers = []
            for paper in search.results():
                try:
                    # Calculate relevance
                    relevance = self._calculate_arxiv_relevance(paper, term, condition)
                    
                    paper_dict = {
                        'title': paper.title or "No title available",
                        'authors': ', '.join([author.name for author in paper.authors]),
                        'journal': f"arXiv ({', '.join([cat.split('.')[0] for cat in paper.categories])})",
                        'year': str(paper.published.year),
                        'pmid': paper.entry_id.split('/')[-1],
                        'abstract': paper.summary or "No abstract available",
                        'url': paper.entry_id,
                        'keywords': paper.categories,
                        'relevance_score': relevance,
                        'database': 'arXiv',
                        'publication_type': 'Preprint'
                    }
                    papers.append(paper_dict)
                    
                except Exception as e:
                    logger.warning(f"Error processing arXiv paper: {e}")
                    continue
            
            return papers
            
        except Exception as e:
            logger.error(f"arXiv search failed for {term}: {e}")
            return []
    
    def _calculate_genomics_relevance(self, article, gene: str, condition: str) -> float:
        """Calculate genomics-specific relevance score"""
        try:
            score = 0.0
            
            # Get all text
            text_fields = [
                article.title or "",
                article.abstract or "",
                " ".join(article.keywords or [])
            ]
            full_text = " ".join(text_fields).lower()
            
            # Gene mention scoring with genomics context
            gene_lower = gene.lower()
            gene_mentions = full_text.count(gene_lower)
            score += min(gene_mentions * 0.3, 1.0)
            
            # Title mentions are critical
            if gene_lower in (article.title or "").lower():
                score += 0.8
            
            # Condition context scoring
            if condition and condition.strip():
                condition_lower = condition.lower()
                condition_mentions = full_text.count(condition_lower)
                score += min(condition_mentions * 0.2, 0.5)
            
            # Genomics-specific terms boost
            genomics_terms = self.config['genomics_context']['gene_prefixes'] + \
                           self.config['genomics_context']['pathways'] + \
                           self.config['genomics_context']['conditions']
            
            genomics_matches = sum(1 for term in genomics_terms if term.lower() in full_text)
            score += min(genomics_matches * 0.1, 0.5)
            
            # Abstract quality bonus
            if article.abstract and len(article.abstract) > 200:
                score += 0.3
            
            # Journal quality (simple heuristic)
            if article.journal:
                high_impact_keywords = ['nature', 'science', 'cell', 'nejm', 'lancet', 'genome', 'genomics']
                if any(keyword in article.journal.lower() for keyword in high_impact_keywords):
                    score += 0.4
            
            return min(score, 3.0)  # Cap at 3.0
            
        except Exception:
            return 0.5
    
    def _calculate_semantic_relevance(self, paper, gene: str, condition: str) -> float:
        """Calculate relevance for Semantic Scholar papers"""
        try:
            score = 0.0
            
            # Text analysis
            text_fields = [paper.title or "", paper.abstract or ""]
            full_text = " ".join(text_fields).lower()
            
            # Gene scoring
            gene_lower = gene.lower()
            if gene_lower in (paper.title or "").lower():
                score += 1.0
            if gene_lower in (paper.abstract or "").lower():
                score += 0.5
            
            # Condition scoring
            if condition and condition.strip():
                condition_lower = condition.lower()
                if condition_lower in full_text:
                    score += 0.4
            
            # Citation count bonus (logarithmic scaling)
            if hasattr(paper, 'citationCount') and paper.citationCount:
                citation_score = min(np.log10(paper.citationCount + 1) * 0.2, 0.6)
                score += citation_score
            
            # Recency bonus
            if hasattr(paper, 'year') and paper.year:
                current_year = datetime.now().year
                if paper.year >= current_year - 3:  # Last 3 years
                    score += 0.3
                elif paper.year >= current_year - 5:  # Last 5 years
                    score += 0.1
            
            return min(score, 3.0)
            
        except Exception:
            return 0.5
    
    def _calculate_arxiv_relevance(self, paper, gene: str, condition: str) -> float:
        """Calculate relevance for arXiv papers"""
        try:
            score = 0.0
            
            # Text analysis
            full_text = f"{paper.title} {paper.summary}".lower()
            
            # Gene and condition scoring
            gene_lower = gene.lower()
            if gene_lower in paper.title.lower():
                score += 1.0
            if gene_lower in paper.summary.lower():
                score += 0.5
            
            if condition and condition.strip():
                if condition.lower() in full_text:
                    score += 0.4
            
            # Category relevance
            relevant_categories = ['q-bio', 'stat.ML', 'cs.LG', 'physics.bio-ph']
            if any(cat in paper.categories for cat in relevant_categories):
                score += 0.5
            
            # Recency bonus (preprints are valuable when recent)
            days_old = (datetime.now() - paper.published).days
            if days_old < 90:  # Last 3 months
                score += 0.4
            elif days_old < 365:  # Last year
                score += 0.2
            
            return min(score, 3.0)
            
        except Exception:
            return 0.5
    
    def _deduplicate_and_rank(self, papers: List[Dict[str, Any]], term: str, condition: str) -> List[Dict[str, Any]]:
        """Remove duplicates and rank by relevance"""
        # Simple deduplication by title similarity
        unique_papers = []
        seen_titles = set()
        
        for paper in papers:
            title_normalized = re.sub(r'[^\w\s]', '', paper['title'].lower())[:50]
            if title_normalized not in seen_titles:
                seen_titles.add(title_normalized)
                unique_papers.append(paper)
        
        # Sort by relevance score
        unique_papers.sort(key=lambda x: x.get('relevance_score', 0), reverse=True)
        
        return unique_papers
    
    def _synthesize_research_with_llm(self, all_results: Dict[str, List[Dict]], 
                                    condition_context: str, query_terms: List[str]) -> Dict[str, Any]:
        """Use LLM to synthesize research findings"""
        try:
            # Prepare research summary for LLM
            research_summary = self._prepare_research_for_llm(all_results, query_terms)
            
            synthesis_prompt = f"""You are a world-class genomics researcher analyzing literature for: {query_terms}
            
            Context: {condition_context}
            
            Research Summary:
            {research_summary}
            
            Provide a concise scientific synthesis focusing on:
            1. Key findings and mechanisms
            2. Clinical implications
            3. Research gaps and future directions
            4. Methodological insights
            
            Be precise, evidence-based, and highlight the most impactful discoveries.
            Limit to 200 words."""
            
            if self.openai_client:
                response = self.openai_client.chat.completions.create(
                    model=self.config['llm']['openai']['model'],
                    messages=[{"role": "user", "content": synthesis_prompt}],
                    max_tokens=400,
                    temperature=0.1
                )
                synthesis_text = response.choices[0].message.content
            elif self.anthropic_client:
                response = self.anthropic_client.messages.create(
                    model=self.config['llm']['anthropic']['model'],
                    max_tokens=400,
                    temperature=0.1,
                    messages=[{"role": "user", "content": synthesis_prompt}]
                )
                synthesis_text = response.content[0].text
            else:
                return {}
            
            # Generate research insights
            insights = self._extract_research_insights(all_results, query_terms)
            
            return {
                "llm_synthesis": synthesis_text,
                "research_insights": insights,
                "synthesis_timestamp": datetime.now(),
                "method": "GPT-4" if self.openai_client else "Claude"
            }
            
        except Exception as e:
            logger.error(f"LLM synthesis failed: {e}")
            return {}
    
    def _prepare_research_for_llm(self, all_results: Dict[str, List[Dict]], query_terms: List[str]) -> str:
        """Prepare research summary for LLM analysis"""
        summary_parts = []
        
        for term in query_terms[:5]:  # Limit to top 5 terms for context length
            if term in all_results and all_results[term]:
                # Get top 3 papers for this term
                top_papers = all_results[term][:3]
                term_summary = f"\n{term}:\n"
                
                for i, paper in enumerate(top_papers):
                    abstract = paper.get('abstract', '')[:200]  # Limit abstract length
                    term_summary += f"  {i+1}. {paper.get('title', 'No title')} ({paper.get('year', 'Unknown year')})\n"
                    term_summary += f"     {abstract}...\n"
                
                summary_parts.append(term_summary)
        
        return "\n".join(summary_parts)
    
    def _extract_research_insights(self, all_results: Dict[str, List[Dict]], query_terms: List[str]) -> Dict[str, Any]:
        """Extract key research insights from results"""
        insights = {
            "top_genes_by_relevance": [],
            "recent_breakthrough_papers": [],
            "high_impact_findings": [],
            "emerging_research_areas": [],
            "clinical_translations": []
        }
        
        try:
            # Aggregate all papers
            all_papers = []
            for papers in all_results.values():
                all_papers.extend(papers)
            
            # Sort by relevance
            all_papers.sort(key=lambda x: x.get('relevance_score', 0), reverse=True)
            
            # Top genes by relevance
            gene_scores = {}
            for term in query_terms:
                if term in all_results:
                    avg_relevance = np.mean([p.get('relevance_score', 0) for p in all_results[term]])
                    gene_scores[term] = avg_relevance
            
            insights["top_genes_by_relevance"] = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)[:5]
            
            # Recent high-impact papers (last 2 years, high relevance)
            current_year = datetime.now().year
            recent_papers = [p for p in all_papers 
                           if p.get('year', '').isdigit() and int(p['year']) >= current_year - 2
                           and p.get('relevance_score', 0) > 2.0]
            insights["recent_breakthrough_papers"] = recent_papers[:5]
            
            # High citation papers (if available)
            high_citation = [p for p in all_papers 
                           if p.get('citation_count', 0) > 100]
            insights["high_impact_findings"] = high_citation[:5]
            
            return insights
            
        except Exception as e:
            logger.error(f"Insight extraction failed: {e}")
            return insights
    
    def _extract_authors(self, article) -> str:
        """Extract and format author information"""
        try:
            if hasattr(article, 'authors') and article.authors:
                if len(article.authors) <= 3:
                    return ", ".join(article.authors)
                else:
                    return f"{', '.join(article.authors[:3])}, et al."
            return "Authors not available"
        except Exception:
            return "Authors not available"
    
    def _extract_year(self, article) -> str:
        """Extract publication year"""
        try:
            if hasattr(article, 'publication_date') and article.publication_date:
                return str(article.publication_date.year)
            elif hasattr(article, 'publication_year') and article.publication_year:
                return str(article.publication_year)
            return "Unknown year"
        except Exception:
            return "Unknown year"
    
    def _extract_publication_type(self, article) -> str:
        """Extract publication type"""
        try:
            # This would need to be implemented based on PubMed metadata
            return "Research Article"
        except Exception:
            return "Unknown type"
    
    def _dict_to_result(self, paper_dict: Dict[str, Any]) -> LiteratureSearchResult:
        """Convert dictionary to LiteratureSearchResult"""
        return LiteratureSearchResult(
            title=paper_dict.get('title', ''),
            authors=paper_dict.get('authors', ''),
            journal=paper_dict.get('journal', ''),
            year=paper_dict.get('year', ''),
            pmid=paper_dict.get('pmid', ''),
            abstract=paper_dict.get('abstract', ''),
            relevance_score=paper_dict.get('relevance_score', 0.0),
            keywords=paper_dict.get('keywords', []),
            url=paper_dict.get('url')
        )
    
    def _create_enhanced_summary(self, all_results: Dict[str, List[Dict]], 
                               synthesis: Dict[str, Any]) -> Dict[str, Any]:
        """Create comprehensive summary with LLM insights"""
        # Basic statistics
        total_papers = sum(len(papers) for papers in all_results.values())
        terms_searched = len(all_results)
        
        # Database breakdown
        database_counts = {}
        for papers in all_results.values():
            for paper in papers:
                db = paper.get('database', 'Unknown')
                database_counts[db] = database_counts.get(db, 0) + 1
        
        # Recent papers (last 3 years)
        current_year = datetime.now().year
        recent_count = 0
        for papers in all_results.values():
            for paper in papers:
                year = paper.get('year', '')
                if year.isdigit() and int(year) >= current_year - 3:
                    recent_count += 1
        
        # High relevance papers
        high_relevance = 0
        for papers in all_results.values():
            for paper in papers:
                if paper.get('relevance_score', 0) > 2.0:
                    high_relevance += 1
        
        summary = {
            "total_papers": total_papers,
            "terms_searched": terms_searched,
            "high_relevance_papers": high_relevance,
            "recent_papers_3yr": recent_count,
            "database_breakdown": database_counts,
            "search_quality": "Excellent" if high_relevance > total_papers * 0.3 else "Good" if high_relevance > 0 else "Basic",
            "llm_synthesis_available": bool(synthesis)
        }
        
        # Add LLM insights if available
        if synthesis:
            summary.update(synthesis)
        
        return summary
    
    def _create_cache_key(self, query_terms: List[str], condition: str, databases: List[str]) -> str:
        """Create cache key for search results"""
        key_data = {
            "terms": sorted(query_terms),
            "condition": condition,
            "databases": sorted(databases)
        }
        return hashlib.md5(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
    
    def export_results(self, literature_results: LiteratureResults,
                      output_dir: str = ".") -> Dict[str, str]:
        """Export literature search results to files"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        exported_files = {}
        
        try:
            # Export all results as CSV
            all_papers = []
            for term, papers in literature_results.results.items():
                for paper in papers:
                    paper_dict = {
                        'search_term': term,
                        'title': paper.title,
                        'authors': paper.authors,
                        'journal': paper.journal,
                        'year': paper.year,
                        'pmid': paper.pmid,
                        'relevance_score': paper.relevance_score,
                        'url': paper.url,
                        'abstract': paper.abstract[:500] + '...' if len(paper.abstract) > 500 else paper.abstract
                    }
                    all_papers.append(paper_dict)
            
            if all_papers:
                results_df = pd.DataFrame(all_papers)
                results_filename = output_path / "literature_search_results.csv"
                results_df.to_csv(results_filename, index=False)
                exported_files["results"] = str(results_filename)
            
            # Export summary
            summary_filename = output_path / "literature_search_summary.json"
            with open(summary_filename, 'w') as f:
                json.dump(literature_results.summary, f, indent=2, default=str)
            exported_files["summary"] = str(summary_filename)
            
            # Export high-relevance papers only
            high_relevance_papers = [paper for paper in all_papers if paper['relevance_score'] > 2.0]
            if high_relevance_papers:
                hr_df = pd.DataFrame(high_relevance_papers)
                hr_filename = output_path / "high_relevance_papers.csv"
                hr_df.to_csv(hr_filename, index=False)
                exported_files["high_relevance"] = str(hr_filename)
            
            logger.info(f"Exported {len(exported_files)} literature search files")
            return exported_files
            
        except Exception as e:
            logger.error(f"Export failed: {e}")
            return {}
    
    def _fallback_search(self, query_terms: List[str], condition_context: str, 
                        max_results_per_term: int) -> LiteratureResults:
        """
        Fallback search when advanced APIs are not available
        
        Uses simple literature search or example data to provide working functionality
        """
        try:
            # Import the simple search fallback
            from .simple_literature_search import get_simple_literature_search
            
            simple_search = get_simple_literature_search()
            
            # Check if simple search is available
            status = simple_search.get_status()
            
            if status["pubmed_api_available"]:
                # Use simple PubMed search
                logger.info("🔍 Using simple PubMed search")
                search_results = simple_search.search_pubmed_simple(
                    query_terms, condition_context, max_results_per_term
                )
            else:
                # Use example data
                logger.info("📝 Using example literature data")
                search_results = simple_search.create_example_results(query_terms)
            
            # Convert to LiteratureResults format
            literature_results = {}
            for term, papers in search_results.items():
                # Convert dict papers to LiteratureSearchResult objects
                result_objects = []
                for paper in papers:
                    result_obj = LiteratureSearchResult(
                        title=paper["title"],
                        authors=paper["authors"], 
                        journal=paper["journal"],
                        year=paper["year"],
                        pmid=paper["pmid"],
                        abstract=paper["abstract"],
                        relevance_score=paper["relevance_score"],
                        keywords=[],
                        url=paper["url"]
                    )
                    result_objects.append(result_obj)
                literature_results[term] = result_objects
            
            # Create summary
            total_papers = sum(len(papers) for papers in literature_results.values())
            
            summary = {
                "total_papers": total_papers,
                "terms_searched": len(query_terms),
                "high_relevance_papers": 0,
                "recent_papers_3yr": 0,
                "database_breakdown": {"PubMed (Simple)" if status["pubmed_api_available"] else "Example Data": total_papers},
                "search_quality": "Basic",
                "llm_synthesis_available": False,
                "fallback_mode": True,
                "fallback_reason": "Advanced APIs not available" if not status["pubmed_api_available"] else "PyMed not installed"
            }
            
            return LiteratureResults(
                results=literature_results,
                summary=summary,
                search_terms=query_terms,
                search_parameters={
                    "databases": ["fallback"],
                    "condition_context": condition_context,
                    "max_results_per_term": max_results_per_term,
                    "fallback_mode": True
                }
            )
            
        except Exception as e:
            logger.error(f"Fallback search failed: {e}")
            
            # Ultimate fallback - return empty results with helpful message
            return LiteratureResults(
                results={term: [] for term in query_terms},
                summary={
                    "total_papers": 0,
                    "terms_searched": len(query_terms),
                    "error": "Literature search not available. Please install dependencies: pip install pymed arxiv",
                    "fallback_mode": True
                },
                search_terms=query_terms,
                search_parameters={"error": True}
            )


# Create singleton instance
_literature_engine = None

def get_literature_engine() -> LiteratureIntelligenceEngine:
    """Get the global literature intelligence engine instance"""
    global _literature_engine
    if _literature_engine is None:
        _literature_engine = LiteratureIntelligenceEngine()
    return _literature_engine