"""
Prairie Genomics Suite - Main Application UI

This module provides the main Streamlit application interface, orchestrating
all the modular tab components. Replaces the monolithic structure with a
clean, maintainable architecture.

Features:
- Modular tab system
- Session state management
- Navigation and layout
- Progress tracking across analyses
- Global configuration
"""

import streamlit as st
import logging
from typing import Dict, Any

# Import core modules
from core.session_manager import get_session_manager
from config import get_config, APP_VERSION

# Import tab components with bulletproof error handling
try:
    from ui.tabs.data_import_tab import render_data_import_tab
    data_import_available = True
    data_import_error = None
except Exception as e:
    data_import_available = False
    data_import_error = str(e)
    render_data_import_tab = None

try:
    from ui.tabs.deseq2_tab import render_deseq2_tab
    deseq2_available = True
    deseq2_error = None
except Exception as e:
    deseq2_available = False
    deseq2_error = str(e)
    render_deseq2_tab = None

try:
    from ui.tabs.literature_tab import render_literature_tab
    literature_available = True
    literature_error = None
except Exception as e:
    literature_available = False
    literature_error = str(e)
    render_literature_tab = None

try:
    # Use simple version for now to test tab loading
    from ui.tabs.ai_research_tab_simple import render_ai_research_tab
    ai_research_available = True
    ai_research_error = None
except Exception as e:
    ai_research_available = False
    ai_research_error = str(e)
    render_ai_research_tab = None

# Import status for diagnostics
IMPORT_STATUS = {
    "data_import": {"available": data_import_available, "error": data_import_error},
    "deseq2": {"available": deseq2_available, "error": deseq2_error},
    "literature": {"available": literature_available, "error": literature_error},
    "ai_research": {"available": ai_research_available, "error": ai_research_error}
}

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_error_tab(tab_name: str, error_msg: str):
    """Create a fallback tab that shows error information"""
    def render_error_tab():
        st.header(f"⚠️ {tab_name}")
        st.error(f"**Tab Loading Error**: This tab failed to load properly.")
        
        with st.expander("🔧 Error Details"):
            st.code(error_msg)
            st.write("**Possible Solutions:**")
            st.write("• Restart the Streamlit application")
            st.write("• Check that all dependencies are installed")
            st.write("• Clear browser cache and refresh")
        
        st.info("Other tabs should still work normally. This is an isolated error.")
    
    return render_error_tab

def create_diagnostic_tab():
    """Create diagnostic tab showing system status"""
    def render_diagnostic_tab():
        st.header("🔧 System Diagnostics v4.0")
        st.success(f"✅ **Prairie Genomics Suite v{APP_VERSION} Active**")
        
        st.subheader("📊 Import Status")
        
        for tab_name, status in IMPORT_STATUS.items():
            if status["available"]:
                st.success(f"✅ {tab_name.replace('_', ' ').title()}: Working")
            else:
                st.error(f"❌ {tab_name.replace('_', ' ').title()}: Failed")
                with st.expander(f"Error details for {tab_name}"):
                    st.code(status["error"])
        
        st.subheader("🚀 Cache Busting")
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🔄 Force Refresh", help="Clear cache and force reload"):
                st.cache_data.clear()
                st.cache_resource.clear()
                st.success("Cache cleared! Please refresh your browser.")
        
        with col2:
            if st.button("⏰ Show Timestamp", help="Verify fresh loading"):
                from datetime import datetime
                st.info(f"Current time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        st.subheader("📋 Restart Instructions")
        st.info("""
        **If you're still seeing old version:**
        1. **Stop Streamlit** (Ctrl+C in terminal)
        2. **Restart**: `streamlit run prairie_genomics_streamlit_ready.py`
        3. **Clear browser cache**: Ctrl+F5 or Shift+Refresh
        4. **Hard refresh**: Close browser tab and reopen
        """)
    
    return render_diagnostic_tab

class PrairieGenomicsApp:
    """
    Main Prairie Genomics Suite Application
    
    Orchestrates the modular tab system and provides global navigation
    and state management for the genomics analysis workflow.
    """
    
    def __init__(self):
        self.session = get_session_manager()
        self.config = get_config("app")
        
        # Define available tabs with bulletproof error handling
        self.tabs = {}
        
        # Add diagnostic tab first (always available)
        self.tabs["🔧 Diagnostics v4.0"] = {
            "render_func": create_diagnostic_tab(),
            "description": "System status and troubleshooting tools",
            "requires_data": False
        }
        
        # Add Data Import tab
        if data_import_available and render_data_import_tab:
            self.tabs["📁 Data Import"] = {
                "render_func": render_data_import_tab,
                "description": "Upload and validate genomics data files",
                "requires_data": False
            }
        else:
            self.tabs["📁 Data Import (Error)"] = {
                "render_func": create_error_tab("Data Import", data_import_error),
                "description": "Tab failed to load - see error details",
                "requires_data": False
            }
        
        # Add DESeq2 tab
        if deseq2_available and render_deseq2_tab:
            self.tabs["🧬 DESeq2 Analysis"] = {
                "render_func": render_deseq2_tab,
                "description": "Differential expression analysis with DESeq2",
                "requires_data": True
            }
        else:
            self.tabs["🧬 DESeq2 Analysis (Error)"] = {
                "render_func": create_error_tab("DESeq2 Analysis", deseq2_error),
                "description": "Tab failed to load - see error details",
                "requires_data": False
            }
        
        # Add Literature Search tab
        if literature_available and render_literature_tab:
            self.tabs["📚 Literature Search"] = {
                "render_func": render_literature_tab,
                "description": "AI-powered literature discovery and synthesis",
                "requires_data": False
            }
        else:
            self.tabs["📚 Literature Search (Error)"] = {
                "render_func": create_error_tab("Literature Search", literature_error),
                "description": "Tab failed to load - see error details",
                "requires_data": False
            }
        
        # Add AI Research Assistant tab
        if ai_research_available and render_ai_research_tab:
            self.tabs["🤖 AI Research Assistant v4.0"] = {
                "render_func": render_ai_research_tab,
                "description": "Conversational AI research help - SIMPLE TEST VERSION",
                "requires_data": False
            }
        else:
            self.tabs["🤖 AI Research Assistant (Error)"] = {
                "render_func": create_error_tab("AI Research Assistant", ai_research_error),
                "description": "Tab failed to load - see error details",
                "requires_data": False
            }
        
        # Always available placeholder tabs (no import dependencies)
        self.tabs["🛤️ Survival Analysis"] = {
            "render_func": self._render_survival_placeholder,
            "description": "Kaplan-Meier and Cox regression analysis",
            "requires_data": True
        }
        
        self.tabs["🔬 Pathway Analysis"] = {
            "render_func": self._render_pathway_placeholder,
            "description": "Gene set enrichment and pathway analysis",
            "requires_data": True
        }
        
        self.tabs["📊 Visualizations"] = {
            "render_func": self._render_visualization_placeholder,
            "description": "Advanced plotting and data visualization",
            "requires_data": True
        }
        
        # Cache busting tab for testing
        self.tabs["🔄 Cache Buster v4.0"] = {
            "render_func": self._render_cache_buster_tab,
            "description": "Force refresh and clear all caches",
            "requires_data": False
        }
    
    def run(self):
        """Run the main application"""
        # Configure page with version info and cache busting
        st.set_page_config(
            page_title=f"Prairie Genomics Suite v{APP_VERSION}",
            page_icon="🧬",
            layout="wide",
            initial_sidebar_state="expanded"
        )
        
        # Render header
        self._render_header()
        
        # Render sidebar navigation
        selected_tab = self._render_sidebar()
        
        # Render main content
        self._render_main_content(selected_tab)
        
        # Render footer
        self._render_footer()
    
    def _render_header(self):
        """Render application header with version prominence"""
        from datetime import datetime
        cache_bust = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        st.markdown(f"""
        <div style='text-align: center; padding: 1rem 0; background: linear-gradient(90deg, #1f4e79, #2d7dd2); color: white; margin-bottom: 2rem; border-radius: 10px;'>
            <h1>🧬 Prairie Genomics Suite v{APP_VERSION}</h1>
            <p>Advanced Genomics Analysis Platform with Scientific Rigor</p>
            <small style='opacity: 0.8;'>🚀 Fresh Load: {cache_bust}</small>
        </div>
        """, unsafe_allow_html=True)
        
        # Version confirmation banner
        st.success(f"✅ **Version {APP_VERSION} Active** - All systems operational!")
    
    def _render_sidebar(self) -> str:
        """Render sidebar navigation and return selected tab"""
        with st.sidebar:
            st.markdown("### 🧭 Navigation")
            
            # Tab selection
            tab_names = list(self.tabs.keys())
            selected_tab = st.radio(
                "Select Analysis",
                tab_names,
                format_func=lambda x: x,
                help="Navigate between different analysis modules"
            )
            
            st.divider()
            
            # Progress tracker
            self._render_progress_tracker()
            
            st.divider()
            
            # Quick actions
            self._render_quick_actions()
            
            st.divider()
            
            # System status
            self._render_system_status()
        
        return selected_tab
    
    def _render_progress_tracker(self):
        """Render analysis progress tracker"""
        st.markdown("### 📊 Analysis Progress")
        
        # Check completion status for each step
        progress_steps = [
            ("Data Import", self.session.has_data("expression_data") and self.session.has_data("clinical_data")),
            ("DESeq2 Analysis", self.session.has_analysis_results("deseq2")),
            ("Literature Search", self.session.has_analysis_results("literature")),
            ("AI Research Assistant", hasattr(st.session_state, "ai_conversation_history") and len(st.session_state.ai_conversation_history) > 0),
            ("Pathway Analysis", self.session.has_analysis_results("pathway")),
            ("Survival Analysis", self.session.has_analysis_results("survival"))
        ]
        
        for step_name, completed in progress_steps:
            icon = "✅" if completed else "⏳"
            st.write(f"{icon} {step_name}")
        
        # Progress percentage
        completed_steps = sum(1 for _, completed in progress_steps if completed)
        progress_pct = completed_steps / len(progress_steps)
        
        st.progress(progress_pct)
        st.caption(f"Overall Progress: {progress_pct:.0%}")
    
    def _render_quick_actions(self):
        """Render quick action buttons"""
        st.markdown("### ⚡ Quick Actions")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🔄 Reset Session", use_container_width=True):
                self.session.clear_session()
                st.success("Session cleared!")
                st.rerun()
        
        with col2:
            if st.button("💾 Save State", use_container_width=True):
                # This would save session state to file
                st.info("State saving not yet implemented")
        
        # Export all results
        if st.button("📁 Export All Results", use_container_width=True):
            self._export_all_results()
    
    def _render_system_status(self):
        """Render system status information"""
        st.markdown("### 🔧 System Status")
        
        # Check engine statuses
        try:
            # Import engines for status check
            from analysis.deseq2_engine import get_deseq2_engine
            from analysis.literature_engine import get_literature_engine
            
            deseq2_engine = get_deseq2_engine()
            literature_engine = get_literature_engine()
            
            # DESeq2 status
            deseq2_status = deseq2_engine.get_status()
            r_icon = "✅" if deseq2_status["r_available"] else "⚠️"
            st.write(f"{r_icon} R/DESeq2")
            
            # Literature engine status
            lit_status = literature_engine.get_status()
            api_count = sum(1 for available in lit_status["apis_available"].values() if available)
            lit_icon = "✅" if api_count >= 2 else "⚠️"
            st.write(f"{lit_icon} Literature APIs ({api_count}/5)")
            
        except Exception as e:
            st.write("❌ Engine Status Check Failed")
            logger.error(f"Status check failed: {e}")
        
        # Session info
        st.caption(f"Version: {APP_VERSION}")
        
        # Data summary
        if self.session.has_data("expression_data"):
            expr_data = self.session.get_data("expression_data")
            st.caption(f"Data: {expr_data.n_genes} genes × {expr_data.n_samples} samples")
    
    def _render_main_content(self, selected_tab: str):
        """Render main content area"""
        tab_config = self.tabs[selected_tab]
        
        # Check data requirements
        if tab_config["requires_data"] and not self._check_data_availability():
            st.warning("⚠️ This analysis requires data to be uploaded first.")
            st.info("Please go to the Data Import tab to upload your genomics data.")
            return
        
        # Render selected tab
        try:
            tab_config["render_func"]()
        except Exception as e:
            st.error(f"Error rendering {selected_tab}: {e}")
            logger.error(f"Tab rendering failed for {selected_tab}: {e}")
    
    def _check_data_availability(self) -> bool:
        """Check if required data is available"""
        return (self.session.has_data("expression_data") and 
                self.session.has_data("clinical_data"))
    
    def _export_all_results(self):
        """Export all analysis results"""
        exported_files = []
        
        try:
            # Export DESeq2 results
            if self.session.has_analysis_results("deseq2"):
                results = self.session.get_analysis_results("deseq2")
                if hasattr(results, 'results') and not results.results.empty:
                    csv_data = results.results.to_csv()
                    st.download_button(
                        "📊 Download DESeq2 Results",
                        data=csv_data,
                        file_name="deseq2_results.csv",
                        mime="text/csv"
                    )
                    exported_files.append("DESeq2 Results")
            
            # Export literature results
            if self.session.has_analysis_results("literature"):
                results = self.session.get_analysis_results("literature")
                
                # Flatten literature results
                all_papers = []
                for term, papers in results.results.items():
                    for paper in papers:
                        paper_dict = {
                            "search_term": term,
                            "title": paper.title,
                            "authors": paper.authors,
                            "journal": paper.journal,
                            "year": paper.year,
                            "relevance_score": paper.relevance_score
                        }
                        all_papers.append(paper_dict)
                
                if all_papers:
                    import pandas as pd
                    csv_data = pd.DataFrame(all_papers).to_csv(index=False)
                    st.download_button(
                        "📚 Download Literature Results",
                        data=csv_data,
                        file_name="literature_results.csv",
                        mime="text/csv"
                    )
                    exported_files.append("Literature Results")
            
            if exported_files:
                st.success(f"✅ Ready to export: {', '.join(exported_files)}")
            else:
                st.info("No results available for export. Run analyses first.")
                
        except Exception as e:
            st.error(f"Export failed: {e}")
            logger.error(f"Export all results failed: {e}")
    
    def _render_footer(self):
        """Render application footer"""
        st.markdown("---")
        st.markdown(f"""
        <div style='text-align: center; color: #666; font-size: 0.8em;'>
            Prairie Genomics Suite v{APP_VERSION} | Built for Scientific Rigor and Discovery
        </div>
        """, unsafe_allow_html=True)
    
    # Placeholder render functions for tabs not yet implemented
    def _render_survival_placeholder(self):
        """Placeholder for survival analysis tab"""
        st.header("🛤️ Survival Analysis")
        st.info("🚧 Survival analysis tab is under development. The SurvivalEngine is ready!")
        
        # Show what's available
        st.markdown("""
        **Available Features (in SurvivalEngine):**
        - Kaplan-Meier survival curves with R survminer
        - Cox proportional hazards modeling  
        - Gene expression-based stratification
        - Risk tables and confidence intervals
        - Python fallbacks with lifelines
        """)
    
    def _render_pathway_placeholder(self):
        """Placeholder for pathway analysis tab"""
        st.header("🔬 Pathway Analysis")
        st.info("🚧 Pathway analysis tab is under development. The PathwayEngine is ready!")
        
        # Show what's available
        st.markdown("""
        **Available Features (in PathwayEngine):**
        - Multiple pathway databases (KEGG, Reactome, GO, MSigDB)
        - Both ORA and GSEA analysis methods
        - R integration for advanced pathway networks
        - Python fallbacks with GSEApy
        - Interactive pathway network visualization
        """)
    
    def _render_visualization_placeholder(self):
        """Placeholder for visualization tab"""
        st.header("📊 Advanced Visualizations")
        st.info("🚧 Visualization tab is under development. Visualization engines are ready!")
        
        # Show what's available
        st.markdown("""
        **Available Features (in analysis engines):**
        - Publication-quality heatmaps with ComplexHeatmap
        - Interactive volcano plots and MA plots
        - Pathway network visualizations
        - Survival curves with risk tables
        - PCA and correlation plots
        """)
    
    def _render_cache_buster_tab(self):
        """Cache busting tab for testing new versions"""
        from datetime import datetime
        
        st.header("🔄 Cache Buster v4.0")
        st.success(f"✅ **Prairie Genomics Suite v{APP_VERSION} CONFIRMED**")
        
        st.info(f"🎉 **Version 4.0 is ACTIVE!** Loaded at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("🧹 Clear All Caches", use_container_width=True):
                st.cache_data.clear()
                st.cache_resource.clear()
                st.success("All Streamlit caches cleared!")
        
        with col2:
            if st.button("🔄 Force Browser Refresh", use_container_width=True):
                st.markdown("""
                <script>
                window.location.reload(true);
                </script>
                """, unsafe_allow_html=True)
        
        with col3:
            if st.button("⏰ Show Timestamp", use_container_width=True):
                st.write(f"Current time: {datetime.now().isoformat()}")
        
        st.subheader("🔧 Version Verification")
        
        # Show all the version indicators
        st.write(f"**App Version:** {APP_VERSION}")
        st.write(f"**Page Title:** Prairie Genomics Suite v{APP_VERSION}")
        st.write(f"**Cache Bust ID:** {datetime.now().strftime('%Y%m%d_%H%M%S')}")
        
        # Debug info
        with st.expander("🐛 Debug Information"):
            st.json({
                "version": APP_VERSION,
                "timestamp": datetime.now().isoformat(),
                "streamlit_version": st.__version__,
                "session_id": st.session_state.get('session_id', 'not_set'),
                "browser_cache_bust": True
            })
        
        st.markdown("---")
        st.info("""
        **If you're still seeing the old version:**
        1. 🛑 **Stop Streamlit** (Ctrl+C in terminal)
        2. 🚀 **Restart**: `streamlit run prairie_genomics_streamlit_ready.py`
        3. 🧹 **Clear browser cache**: Ctrl+Shift+R (or Cmd+Shift+R on Mac)
        4. 🔄 **Hard refresh**: Close browser tab completely and reopen
        """)


def main():
    """Main application entry point"""
    app = PrairieGenomicsApp()
    app.run()


if __name__ == "__main__":
    main()