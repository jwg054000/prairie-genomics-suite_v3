"""
Prairie Genomics Suite - Literature Search Tab

This tab provides advanced literature search capabilities with LLM integration,
multi-database search, and intelligent synthesis. Extracted from the monolithic
structure to provide a focused, powerful research interface.

Features:
- Multi-database literature search (PubMed, Semantic Scholar, arXiv)
- LLM-powered query expansion and synthesis
- Gene-context literature discovery
- Research insights and trend analysis
- Export and citation management
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any
import logging

# Import core modules
from core.session_manager import get_session_manager
from analysis.literature_engine import get_literature_engine
from config import get_config

logger = logging.getLogger(__name__)

class LiteratureTab:
    """
    Literature Search Tab Component
    
    Provides comprehensive literature search and analysis interface with
    LLM integration for intelligent research discovery and synthesis.
    """
    
    def __init__(self):
        self.session = get_session_manager()
        self.literature_engine = get_literature_engine()
        self.config = get_config("literature")
    
    def render(self):
        """Render the complete literature search tab"""
        st.header("📚 Intelligent Literature Search")
        st.markdown("Advanced literature discovery with AI-powered analysis across multiple scientific databases.")
        
        # Engine status check
        self._render_engine_status()
        
        # Main search interface
        col1, col2 = st.columns([2, 1])
        
        with col1:
            self._render_search_interface()
        
        with col2:
            self._render_search_sidebar()
        
        # Results section
        if self.session.has_analysis_results("literature"):
            self._render_results_section()
    
    def _render_engine_status(self):
        """Render literature engine status"""
        engine_status = self.literature_engine.get_status()
        
        # Quick status indicators
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            pubmed_icon = "✅" if engine_status["apis_available"]["pubmed"] else "🔄"
            status_text = "PubMed" if engine_status["apis_available"]["pubmed"] else "PubMed (Fallback)"
            st.write(f"{pubmed_icon} **{status_text}**")
        
        with col2:
            semantic_icon = "✅" if engine_status["apis_available"]["semantic_scholar"] else "❌"
            st.write(f"{semantic_icon} **Semantic Scholar**")
        
        with col3:
            arxiv_icon = "✅" if engine_status["apis_available"]["arxiv"] else "❌"
            st.write(f"{arxiv_icon} **arXiv**")
        
        with col4:
            llm_icon = "✅" if (engine_status["apis_available"]["openai"] or engine_status["apis_available"]["anthropic"]) else "❌"
            st.write(f"{llm_icon} **AI Analysis**")
        
        # Show helpful status information
        unavailable_apis = [api for api, available in engine_status["apis_available"].items() if not available]
        
        if len(unavailable_apis) == len(engine_status["apis_available"]):
            # All APIs unavailable - show fallback info
            st.info("🔄 **Fallback Mode Active** - Using basic literature search with example data. For full functionality, install: `pip install pymed arxiv openai anthropic`")
        elif unavailable_apis:
            # Some APIs available
            available_count = len(engine_status["apis_available"]) - len(unavailable_apis)
            st.success(f"✅ **{available_count} of {len(engine_status['apis_available'])} APIs active** - Literature search ready!")
            
            # Show what's missing
            missing_packages = []
            if not engine_status["apis_available"]["pubmed"]:
                missing_packages.append("pymed")
            if not engine_status["apis_available"]["semantic_scholar"]:
                missing_packages.append("semanticscholar")
            if not engine_status["apis_available"]["arxiv"]:
                missing_packages.append("arxiv")
            if not engine_status["apis_available"]["openai"]:
                missing_packages.append("openai + API key")
            if not engine_status["apis_available"]["anthropic"]:
                missing_packages.append("anthropic + API key")
            
            if missing_packages:
                st.info(f"💡 **Optional enhancements available:** Install {', '.join(missing_packages)} for additional features")
        else:
            # All APIs available
            st.success("🎉 **All literature APIs active** - Full functionality available!")
        
        # Add session reset button
        if st.button("🔄 Reset Session", help="Clear cached data and restart"):
            self.session.clear_session()
            st.success("Session reset! Please refresh the page.")
            st.rerun()
    
    def _render_search_interface(self):
        """Render main search interface"""
        st.subheader("🔍 Search Configuration")
        
        # Search source selection
        search_tabs = st.tabs([
            "🧬 Gene-Based Search", 
            "📝 Manual Search", 
            "🎯 Quick Examples"
        ])
        
        with search_tabs[0]:
            self._render_gene_based_search()
        
        with search_tabs[1]:
            self._render_manual_search()
        
        with search_tabs[2]:
            self._render_quick_examples()
    
    def _render_gene_based_search(self):
        """Render gene-based search from analysis results"""
        st.write("**Search Literature for Genes from Analysis Results**")
        
        # Check for available gene lists
        gene_sources = []
        
        if self.session.has_analysis_results("deseq2_gene_lists"):
            gene_lists = self.session.get_analysis_results("deseq2_gene_lists")
            gene_sources.append(("DESeq2 Results", gene_lists))
        
        if self.session.has_data("expression_data"):
            expr_data = self.session.get_data("expression_data")
            top_variable_genes = expr_data.data.var(axis=1).nlargest(100).index.tolist()
            gene_sources.append(("Top Variable Genes", {"top_variable": top_variable_genes}))
        
        if not gene_sources:
            st.info("No analysis results available. Run DESeq2 analysis first or use manual search.")
            return
        
        # Gene list selection
        source_options = [source[0] for source in gene_sources]
        selected_source = st.selectbox("Select gene source", source_options)
        
        if selected_source:
            # Get selected gene lists
            source_data = next(data for name, data in gene_sources if name == selected_source)
            
            if selected_source == "DESeq2 Results":
                gene_list_type = st.selectbox(
                    "Select gene list type",
                    ["significant", "upregulated", "downregulated", "all_genes"]
                )
                
                selected_genes = source_data.get(gene_list_type, [])
                
                # Limit number of genes for literature search
                max_genes = st.slider(
                    "Maximum genes to search",
                    min_value=5,
                    max_value=min(50, len(selected_genes)),
                    value=min(20, len(selected_genes)),
                    help="Limit genes to avoid API limits and focus on top results"
                )
                
                selected_genes = selected_genes[:max_genes]
                
            elif selected_source == "Top Variable Genes":
                selected_genes = source_data["top_variable"]
                max_genes = st.slider("Maximum genes to search", 5, 50, 20)
                selected_genes = selected_genes[:max_genes]
            
            # Context information
            condition_context = ""
            if self.session.has_analysis_params("deseq2_design"):
                design_info = self.session.get_analysis_params("deseq2_design")
                condition_context = f"{design_info['treatment_group']} vs {design_info['reference_group']}"
            
            context_input = st.text_input(
                "Research context (optional)",
                value=condition_context,
                help="Provide disease/condition context to improve search relevance"
            )
            
            # Display selected genes
            if selected_genes:
                st.write(f"**Selected genes ({len(selected_genes)}):**")
                genes_display = ", ".join(selected_genes[:10])
                if len(selected_genes) > 10:
                    genes_display += f"... and {len(selected_genes) - 10} more"
                st.write(genes_display)
                
                # Store search parameters
                self.session.store_analysis_params("literature_search", {
                    "query_terms": selected_genes,
                    "condition_context": context_input,
                    "search_type": "gene_based"
                })
    
    def _render_manual_search(self):
        """Render manual search interface"""
        st.write("**Manual Literature Search**")
        
        # Manual gene/term input
        manual_terms = st.text_area(
            "Enter genes or search terms (one per line)",
            height=100,
            placeholder="TP53\nBRCA1\napoptosis\ncell cycle",
            help="Enter gene names, pathways, or biological processes"
        )
        
        # Context
        condition_context = st.text_input(
            "Research context",
            placeholder="cancer, immunotherapy, aging, etc.",
            help="Provide biological or clinical context"
        )
        
        if manual_terms:
            # Parse terms
            terms = [term.strip() for term in manual_terms.split('\n') if term.strip()]
            
            if len(terms) > 50:
                st.warning(f"⚠️ Too many terms ({len(terms)}). Consider limiting to 50 or fewer.")
                terms = terms[:50]
            
            st.write(f"**Search terms ({len(terms)}):** {', '.join(terms[:5])}{'...' if len(terms) > 5 else ''}")
            
            # Store search parameters
            self.session.store_analysis_params("literature_search", {
                "query_terms": terms,
                "condition_context": condition_context,
                "search_type": "manual"
            })
    
    def _render_quick_examples(self):
        """Render quick example searches"""
        st.write("**Quick Example Searches**")
        st.markdown("Click any example to run a demonstration search:")
        
        examples = [
            {
                "name": "🧬 Oncogenes in Cancer",
                "genes": ["TP53", "MYC", "RAS", "EGFR", "PI3K"],
                "context": "cancer oncogenes"
            },
            {
                "name": "🛡️ Immune Checkpoint Genes",
                "genes": ["PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT"],
                "context": "immune checkpoint immunotherapy"
            },
            {
                "name": "🔄 Cell Cycle Regulation",
                "genes": ["CDK1", "CDK2", "CCND1", "RB1", "E2F1"],
                "context": "cell cycle regulation"
            },
            {
                "name": "💀 Apoptosis Pathway",
                "genes": ["BAX", "BCL2", "CASP3", "TP53", "APAF1"],
                "context": "apoptosis programmed cell death"
            }
        ]
        
        for example in examples:
            if st.button(example["name"], use_container_width=True):
                # Store example search
                self.session.store_analysis_params("literature_search", {
                    "query_terms": example["genes"],
                    "condition_context": example["context"],
                    "search_type": "example"
                })
                st.success(f"✅ Loaded example: {example['name']}")
                st.rerun()
    
    def _render_search_sidebar(self):
        """Render search configuration sidebar"""
        st.subheader("⚙️ Search Settings")
        
        # Database selection
        available_databases = ["pubmed", "semantic_scholar", "arxiv"]
        
        selected_databases = st.multiselect(
            "Select databases",
            available_databases,
            default=["pubmed", "semantic_scholar"],
            help="Choose which databases to search"
        )
        
        # Search parameters
        max_results_per_term = st.slider(
            "Results per gene/term",
            min_value=5,
            max_value=50,
            value=10,
            help="Maximum papers to retrieve per search term"
        )
        
        # Publication year filter
        use_year_filter = st.checkbox("Filter by publication year")
        
        publication_years = None
        if use_year_filter:
            year_range = st.slider(
                "Publication year range",
                min_value=2000,
                max_value=2024,
                value=(2020, 2024),
                help="Limit search to recent publications"
            )
            publication_years = year_range
        
        # Advanced options
        with st.expander("Advanced Options"):
            include_reviews = st.checkbox("Include review articles", value=True)
            semantic_expansion = st.checkbox("AI query expansion", value=True, help="Use LLM to expand search terms")
        
        # Execute search
        st.divider()
        
        if self.session.has_analysis_params("literature_search"):
            search_params = self.session.get_analysis_params("literature_search")
            
            st.write("**Ready to Search:**")
            st.write(f"• Terms: {len(search_params['query_terms'])}")
            st.write(f"• Databases: {len(selected_databases)}")
            st.write(f"• Context: {search_params['condition_context'][:30]}{'...' if len(search_params['condition_context']) > 30 else ''}")
            
            if st.button("🚀 Start Literature Search", type="primary", use_container_width=True):
                self._execute_literature_search(
                    search_params,
                    selected_databases,
                    max_results_per_term,
                    publication_years,
                    include_reviews,
                    semantic_expansion
                )
        else:
            st.info("Configure search terms above")
    
    def _execute_literature_search(self, search_params, databases, max_results_per_term,
                                 publication_years, include_reviews, semantic_expansion):
        """Execute literature search"""
        try:
            # Show progress
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            status_text.text("Initializing literature search...")
            progress_bar.progress(0.1)
            
            # Prepare search parameters
            query_terms = search_params["query_terms"]
            condition_context = search_params["condition_context"]
            
            status_text.text(f"Searching {len(databases)} databases for {len(query_terms)} terms...")
            progress_bar.progress(0.3)
            
            # Execute intelligent search
            results = self.literature_engine.intelligent_search(
                query_terms=query_terms,
                condition_context=condition_context,
                max_results_per_term=max_results_per_term,
                databases=databases,
                publication_years=publication_years,
                include_reviews=include_reviews,
                semantic_expansion=semantic_expansion
            )
            
            progress_bar.progress(0.8)
            status_text.text("Processing and synthesizing results...")
            
            # Store results
            self.session.store_analysis_results("literature", results)
            
            progress_bar.progress(1.0)
            status_text.text("Literature search complete!")
            
            # Show appropriate success message based on mode
            if hasattr(results, 'search_parameters') and results.search_parameters.get('fallback_mode'):
                if 'error' in results.search_parameters:
                    st.error("❌ Literature search failed. Please install required dependencies.")
                else:
                    fallback_reason = results.summary.get('fallback_reason', 'using basic search')
                    st.info(f"🔄 Literature search completed in fallback mode ({fallback_reason}). Found {results.total_papers} papers.")
                    if results.total_papers > 0:
                        st.info("💡 These are example results. Install `pip install pymed arxiv` for real literature search.")
            else:
                st.success(f"✅ Literature search completed! Found {results.total_papers} papers across {len(query_terms)} terms.")
            
            # Auto-refresh to show results
            st.rerun()
            
        except Exception as e:
            st.error(f"❌ Literature search failed: {e}")
            logger.error(f"Literature search failed: {e}")
    
    def _render_results_section(self):
        """Render comprehensive results section"""
        st.subheader("📊 Literature Search Results")
        
        results = self.session.get_analysis_results("literature")
        
        # Results overview
        self._render_results_overview(results)
        
        # Results tabs
        result_tabs = st.tabs([
            "📋 Papers by Gene",
            "🧠 AI Synthesis", 
            "📈 Research Insights",
            "📁 Export & Citations"
        ])
        
        with result_tabs[0]:
            self._render_papers_by_gene(results)
        
        with result_tabs[1]:
            self._render_ai_synthesis(results)
        
        with result_tabs[2]:
            self._render_research_insights(results)
        
        with result_tabs[3]:
            self._render_export_options(results)
    
    def _render_results_overview(self, results):
        """Render results overview metrics"""
        # Check if in fallback mode
        is_fallback = results.search_parameters.get("fallback_mode", False)
        
        if is_fallback:
            fallback_reason = results.summary.get("fallback_reason", "basic search mode")
            st.info(f"🔄 **Fallback Mode**: Results from {fallback_reason}")
        
        with st.container():
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                papers_label = "Example Papers" if is_fallback and results.total_papers > 0 else "Total Papers"
                st.metric(papers_label, results.total_papers)
            
            with col2:
                st.metric("Search Terms", len(results.search_terms))
            
            with col3:
                high_relevance = sum(
                    len([p for p in papers if p.relevance_score > 2.0])
                    for papers in results.results.values()
                )
                st.metric("High Relevance", high_relevance)
            
            with col4:
                databases = results.search_parameters.get("databases", [])
                db_count = len(databases) if not is_fallback else 1
                db_label = "Fallback" if is_fallback else "Databases"
                st.metric(db_label, db_count)
    
    def _render_papers_by_gene(self, results):
        """Render papers organized by gene/term"""
        st.write("**Literature Results by Search Term**")
        
        # Term selection
        term_options = list(results.results.keys())
        selected_term = st.selectbox("Select term to view papers", term_options)
        
        if selected_term and selected_term in results.results:
            papers = results.results[selected_term]
            
            if not papers:
                st.info(f"No papers found for {selected_term}")
                return
            
            st.write(f"**{len(papers)} papers found for {selected_term}**")
            
            # Sort by relevance
            papers_sorted = sorted(papers, key=lambda x: x.relevance_score, reverse=True)
            
            # Display papers
            for i, paper in enumerate(papers_sorted[:20]):  # Show top 20
                with st.expander(f"[{i+1}] {paper.title} (Relevance: {paper.relevance_score:.2f})"):
                    col1, col2 = st.columns([2, 1])
                    
                    with col1:
                        st.write(f"**Authors:** {paper.authors}")
                        st.write(f"**Journal:** {paper.journal} ({paper.year})")
                        
                        if paper.abstract:
                            st.write("**Abstract:**")
                            st.write(paper.abstract[:500] + "..." if len(paper.abstract) > 500 else paper.abstract)
                    
                    with col2:
                        if paper.url:
                            st.link_button("📖 Read Paper", paper.url)
                        
                        st.write(f"**PMID:** {paper.pmid}")
                        st.write(f"**Relevance:** {paper.relevance_score:.2f}")
                        
                        if paper.keywords:
                            st.write(f"**Keywords:** {', '.join(paper.keywords[:5])}")
            
            if len(papers_sorted) > 20:
                st.info(f"Showing top 20 of {len(papers_sorted)} papers. Use export to get complete results.")
    
    def _render_ai_synthesis(self, results):
        """Render AI-powered synthesis"""
        st.write("**AI-Powered Research Synthesis**")
        
        if "llm_synthesis" in results.summary:
            synthesis_text = results.summary["llm_synthesis"]
            
            st.markdown("### 🧠 Research Summary")
            st.write(synthesis_text)
            
            # Synthesis metadata
            if "synthesis_timestamp" in results.summary:
                st.caption(f"Generated: {results.summary['synthesis_timestamp']}")
            
            if "method" in results.summary:
                st.caption(f"AI Model: {results.summary['method']}")
        else:
            st.info("AI synthesis not available. This may be due to:")
            st.write("• No LLM APIs configured")
            st.write("• Insufficient search results")
            st.write("• API rate limits")
    
    def _render_research_insights(self, results):
        """Render research insights and trends"""
        st.write("**Research Insights & Analytics**")
        
        if "research_insights" in results.summary:
            insights = results.summary["research_insights"]
            
            # Top genes by relevance
            if "top_genes_by_relevance" in insights:
                st.subheader("🎯 Top Genes by Research Relevance")
                
                top_genes = insights["top_genes_by_relevance"]
                if top_genes:
                    genes_df = pd.DataFrame(top_genes, columns=["Gene", "Avg Relevance"])
                    genes_df["Avg Relevance"] = genes_df["Avg Relevance"].round(3)
                    st.dataframe(genes_df, use_container_width=True)
            
            # Recent breakthrough papers
            if "recent_breakthrough_papers" in insights:
                st.subheader("🚀 Recent High-Impact Papers")
                
                recent_papers = insights["recent_breakthrough_papers"]
                if recent_papers:
                    for paper in recent_papers[:5]:
                        with st.expander(f"{paper['title']} ({paper['year']})"):
                            st.write(f"**Authors:** {paper['authors']}")
                            st.write(f"**Journal:** {paper['journal']}")
                            st.write(f"**Relevance Score:** {paper['relevance_score']:.2f}")
                            if paper.get('url'):
                                st.link_button("📖 Read Paper", paper['url'])
            
            # High-impact findings
            if "high_impact_findings" in insights:
                st.subheader("📈 High-Citation Papers")
                
                high_impact = insights["high_impact_findings"]
                if high_impact:
                    st.write(f"Found {len(high_impact)} highly-cited papers in your search results")
                    
                    for paper in high_impact[:3]:
                        citation_count = paper.get('citation_count', 'Unknown')
                        st.write(f"• **{paper['title']}** ({citation_count} citations)")
        
        # Database distribution
        st.subheader("📊 Database Coverage")
        
        database_counts = {}
        for papers in results.results.values():
            for paper in papers:
                # Extract database from paper metadata if available
                database = getattr(paper, 'database', 'Unknown')
                database_counts[database] = database_counts.get(database, 0) + 1
        
        if database_counts:
            db_df = pd.DataFrame(list(database_counts.items()), columns=["Database", "Papers"])
            st.bar_chart(db_df.set_index("Database"))
    
    def _render_export_options(self, results):
        """Render export and citation options"""
        st.write("**Export & Citation Management**")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📁 Export Results")
            
            # Export all results
            if st.button("📊 Export All Results (CSV)", use_container_width=True):
                # Flatten results for CSV export
                all_papers = []
                for term, papers in results.results.items():
                    for paper in papers:
                        paper_dict = {
                            "search_term": term,
                            "title": paper.title,
                            "authors": paper.authors,
                            "journal": paper.journal,
                            "year": paper.year,
                            "pmid": paper.pmid,
                            "relevance_score": paper.relevance_score,
                            "url": paper.url,
                            "abstract": paper.abstract
                        }
                        all_papers.append(paper_dict)
                
                if all_papers:
                    csv_data = pd.DataFrame(all_papers).to_csv(index=False)
                    st.download_button(
                        label="Download CSV",
                        data=csv_data,
                        file_name="literature_search_results.csv",
                        mime="text/csv"
                    )
            
            # Export high-relevance only
            if st.button("⭐ Export High-Relevance Papers", use_container_width=True):
                high_rel_papers = []
                for term, papers in results.results.items():
                    for paper in papers:
                        if paper.relevance_score > 2.0:
                            paper_dict = {
                                "search_term": term,
                                "title": paper.title,
                                "authors": paper.authors,
                                "journal": paper.journal,
                                "year": paper.year,
                                "pmid": paper.pmid,
                                "relevance_score": paper.relevance_score,
                                "url": paper.url
                            }
                            high_rel_papers.append(paper_dict)
                
                if high_rel_papers:
                    csv_data = pd.DataFrame(high_rel_papers).to_csv(index=False)
                    st.download_button(
                        label="Download High-Relevance Papers",
                        data=csv_data,
                        file_name="high_relevance_papers.csv",
                        mime="text/csv"
                    )
                else:
                    st.info("No high-relevance papers found")
        
        with col2:
            st.subheader("📚 Citation Export")
            
            # Generate citation text
            if st.button("📋 Generate Citation List", use_container_width=True):
                citations = []
                paper_count = 0
                
                for term, papers in results.results.items():
                    for paper in papers[:5]:  # Top 5 per term
                        if paper.relevance_score > 1.5:  # Only relevant papers
                            citation = f"{paper.authors}. {paper.title}. {paper.journal}. {paper.year}."
                            if paper.pmid and paper.pmid != "Unknown":
                                citation += f" PMID: {paper.pmid}."
                            citations.append(citation)
                            paper_count += 1
                
                if citations:
                    citation_text = "\n\n".join(citations)
                    st.download_button(
                        label="Download Citations",
                        data=citation_text,
                        file_name="literature_citations.txt",
                        mime="text/plain"
                    )
                    st.success(f"Generated {paper_count} citations")
        
        # Save for other analyses
        st.subheader("🔗 Integration with Other Analyses")
        
        # Extract genes for pathway analysis
        all_genes = []
        for term in results.search_terms:
            # Simple heuristic: if term looks like a gene symbol, include it
            if term.isupper() and len(term) <= 10:
                all_genes.append(term)
        
        if all_genes:
            literature_gene_lists = {
                "literature_genes": all_genes,
                "high_relevance_genes": all_genes[:20]  # Top genes for pathway analysis
            }
            self.session.store_analysis_results("literature_gene_lists", literature_gene_lists)
            
            st.success(f"✅ Saved {len(all_genes)} genes for pathway analysis!")
            st.info("💡 Use these gene lists in the Pathway Analysis tab for functional enrichment.")


def render_literature_tab():
    """Render the literature search tab (main entry point)"""
    tab = LiteratureTab()
    tab.render()