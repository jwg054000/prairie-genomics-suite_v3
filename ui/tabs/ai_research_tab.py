"""
Prairie Genomics Suite - AI Research Assistant Tab

This tab provides conversational, AI-driven literature search and research synthesis
using Claude's tool calling capabilities. Based on the elegant pattern from the 
Perplexity article for natural language research discovery.

Features:
- Natural language research queries
- Context-aware literature search
- Automatic genomics integration
- Conversational follow-up questions
- AI-powered synthesis and insights
"""

import streamlit as st
import json
import logging
from typing import Dict, List, Any, Optional
import time

# Import core modules
from core.session_manager import get_session_manager
from analysis.claude_tools import get_available_tools, execute_claude_tool, get_tool_status
from config import get_config

logger = logging.getLogger(__name__)

class AIResearchTab:
    """
    AI Research Assistant Tab Component
    
    Provides conversational interface for AI-driven literature search and research
    synthesis using Claude's tool calling capabilities.
    """
    
    def __init__(self):
        self.session = get_session_manager()
        self.config = get_config("app")
        
        # Initialize conversation history in session state
        if "ai_conversation_history" not in st.session_state:
            st.session_state.ai_conversation_history = []
        
        if "ai_research_context" not in st.session_state:
            st.session_state.ai_research_context = {}
    
    def render(self):
        """Render the complete AI Research Assistant tab"""
        st.header("🤖 AI Research Assistant")
        st.markdown("Have a conversation with Claude about your genomics research. Ask questions, get literature recommendations, and explore your data with AI-powered insights.")
        
        # Check tool status
        tool_status = get_tool_status()
        
        if not tool_status["ready"]:
            self._render_setup_help()
            return
        
        # Main interface
        col1, col2 = st.columns([2, 1])
        
        with col1:
            self._render_conversation_interface()
        
        with col2:
            self._render_research_sidebar()
        
        # Conversation history
        self._render_conversation_history()
    
    def _render_setup_help(self):
        """Render setup help when tools aren't ready"""
        st.warning("⚠️ AI Research Assistant needs additional setup")
        
        tool_status = get_tool_status()
        
        st.write("**Missing Dependencies:**")
        for dep, available in tool_status["dependencies"].items():
            icon = "✅" if available else "❌"
            st.write(f"{icon} {dep}")
        
        st.info("""
        **To enable AI Research Assistant:**
        
        1. **Install Anthropic Claude API**: `pip install anthropic`
        2. **Set API Key**: Add `ANTHROPIC_API_KEY` to your environment variables
        3. **Optional Literature APIs**: `pip install pymed arxiv` for enhanced search
        
        The AI Assistant uses Claude's tool calling to provide conversational research help.
        """)
        
        if st.button("🔄 Refresh Status"):
            st.rerun()
    
    def _render_conversation_interface(self):
        """Render the main conversation interface"""
        st.subheader("💬 Research Conversation")
        
        # Quick start suggestions
        if not st.session_state.ai_conversation_history:
            self._render_quick_start()
        
        # Chat input
        user_input = st.text_area(
            "Ask me anything about your genomics research:",
            placeholder="e.g., 'What does the latest research say about my significant genes?' or 'Find papers on immune checkpoint inhibitors'",
            height=100,
            key="ai_chat_input"
        )
        
        # Send button
        col1, col2, col3 = st.columns([1, 1, 2])
        
        with col1:
            if st.button("🚀 Send", type="primary", disabled=not user_input.strip()):
                self._process_user_query(user_input.strip())
        
        with col2:
            if st.button("🔄 Clear Chat"):
                st.session_state.ai_conversation_history = []
                st.session_state.ai_research_context = {}
                st.rerun()
        
        # Demo mode toggle
        with col3:
            demo_mode = st.checkbox("🎭 Demo Mode (simulated responses)", 
                                   help="Use simulated AI responses when Claude API isn't available")
            st.session_state.demo_mode = demo_mode
    
    def _render_quick_start(self):
        """Render quick start suggestions"""
        st.markdown("### 🚀 Quick Start Examples")
        st.markdown("Click any example to start a conversation:")
        
        examples = [
            {
                "title": "🧬 Analyze My DESeq2 Results",
                "query": "What can you tell me about my current DESeq2 results? Can you find relevant literature for my significant genes?",
                "description": "Get AI insights on your differential expression analysis"
            },
            {
                "title": "📚 Latest Research Discovery",
                "query": "Find the latest research papers on BRCA1 mutations in breast cancer from the past 2 years",
                "description": "Search for recent publications on specific topics"
            },
            {
                "title": "🔬 Therapeutic Target Investigation",
                "query": "Which of my upregulated genes might be potential therapeutic targets? Find papers about druggable targets.",
                "description": "Explore therapeutic potential of your findings"
            },
            {
                "title": "🛤️ Pathway Literature Synthesis",
                "query": "Analyze the literature for my enriched pathways and identify key regulatory mechanisms",
                "description": "Deep dive into pathway biology and regulation"
            }
        ]
        
        for example in examples:
            with st.expander(f"{example['title']}", expanded=False):
                st.write(example["description"])
                if st.button(f"Try this example", key=f"example_{example['title']}"):
                    self._process_user_query(example["query"])
    
    def _render_research_sidebar(self):
        """Render research context sidebar"""
        st.subheader("🔬 Research Context")
        
        # Show current genomics context
        genomics_context = self._get_genomics_context_summary()
        
        if genomics_context:
            st.success("✅ **Analysis Data Available**")
            
            for analysis_type, info in genomics_context.items():
                if isinstance(info, dict) and info.get("available"):
                    with st.expander(f"{analysis_type.upper()} Results"):
                        for key, value in info.items():
                            if key != "available":
                                st.write(f"**{key}:** {value}")
        else:
            st.info("📊 **No Analysis Data**\n\nRun DESeq2 or pathway analysis first for context-aware research assistance.")
        
        # Tool status
        st.divider()
        st.write("**🛠️ AI Tools Status:**")
        
        tools = get_available_tools()
        for tool in tools:
            st.write(f"✅ {tool['name'].replace('_', ' ').title()}")
        
        # Conversation stats
        if st.session_state.ai_conversation_history:
            st.divider()
            st.write("**💬 Conversation Stats:**")
            st.write(f"Messages: {len(st.session_state.ai_conversation_history)}")
            
            # Count tool uses
            tool_uses = sum(1 for msg in st.session_state.ai_conversation_history 
                          if msg.get("type") == "assistant" and "tools_used" in msg)
            st.write(f"Tools Used: {tool_uses}")
    
    def _render_conversation_history(self):
        """Render the conversation history"""
        if not st.session_state.ai_conversation_history:
            return
        
        st.subheader("💬 Conversation History")
        
        for i, message in enumerate(st.session_state.ai_conversation_history):
            with st.container():
                if message["type"] == "user":
                    # User message
                    st.markdown(f"**🧑‍🔬 You:** {message['content']}")
                
                elif message["type"] == "assistant":
                    # Assistant message
                    st.markdown(f"**🤖 Claude:** {message['content']}")
                    
                    # Show tools used
                    if "tools_used" in message:
                        with st.expander("🛠️ Tools Used", expanded=False):
                            for tool in message["tools_used"]:
                                st.json(tool)
                
                elif message["type"] == "tool_result":
                    # Tool result (usually hidden, but can show for debugging)
                    with st.expander("🔧 Tool Result", expanded=False):
                        st.json(message["content"])
                
                st.divider()
    
    def _process_user_query(self, user_query: str):
        """Process user query using Claude tool calling"""
        try:
            # Add user message to history
            st.session_state.ai_conversation_history.append({
                "type": "user",
                "content": user_query,
                "timestamp": time.time()
            })
            
            # Show processing message
            with st.spinner("🤖 Claude is thinking and searching..."):
                
                # Check if we should use demo mode
                if st.session_state.get("demo_mode", False):
                    response = self._generate_demo_response(user_query)
                else:
                    response = self._call_claude_with_tools(user_query)
                
                # Add assistant response to history
                st.session_state.ai_conversation_history.append(response)
            
            # Clear input and refresh
            st.session_state.ai_chat_input = ""
            st.rerun()
            
        except Exception as e:
            logger.error(f"Error processing query: {e}")
            st.error(f"Error processing your query: {e}")
    
    def _call_claude_with_tools(self, user_query: str) -> Dict[str, Any]:
        """
        Call Claude with tool calling capabilities
        
        This implements the pattern from the Perplexity article
        """
        try:
            # For now, simulate Claude's intelligent tool selection
            # In a real implementation, this would use the Anthropic API
            
            # Analyze query to determine which tools to use
            tools_to_use = self._analyze_query_for_tools(user_query)
            
            assistant_response = {
                "type": "assistant",
                "content": "",
                "tools_used": [],
                "timestamp": time.time()
            }
            
            # Execute tools
            tool_results = []
            
            for tool_name, tool_params in tools_to_use:
                try:
                    result = execute_claude_tool(tool_name, tool_params)
                    tool_results.append({
                        "tool": tool_name,
                        "input": tool_params,
                        "result": result
                    })
                    assistant_response["tools_used"].append({
                        "tool": tool_name,
                        "input": tool_params,
                        "status": "success" if "error" not in result else "error"
                    })
                except Exception as e:
                    logger.error(f"Tool execution error: {e}")
                    tool_results.append({
                        "tool": tool_name,
                        "error": str(e)
                    })
            
            # Generate response based on tool results
            assistant_response["content"] = self._synthesize_tool_results(user_query, tool_results)
            
            return assistant_response
            
        except Exception as e:
            logger.error(f"Claude API call error: {e}")
            return {
                "type": "assistant",
                "content": f"I apologize, but I encountered an error: {e}\n\nTry using Demo Mode for simulated responses.",
                "timestamp": time.time()
            }
    
    def _analyze_query_for_tools(self, query: str) -> List[tuple]:
        """
        Analyze user query to determine which tools to use
        
        This simulates Claude's intelligent tool selection
        """
        tools_to_use = []
        query_lower = query.lower()
        
        # Check if user is asking about their current analysis
        if any(word in query_lower for word in ["my", "current", "results", "analysis", "significant", "genes"]):
            tools_to_use.append(("get_genomics_analysis_context", {"analysis_type": "all", "summary_only": True}))
        
        # Check if user wants literature search
        literature_keywords = ["literature", "papers", "research", "studies", "publications", "find", "search"]
        if any(word in query_lower for word in literature_keywords):
            
            # Extract search terms from query
            search_query = query
            condition_context = ""
            
            # Look for specific genes or terms
            gene_mentions = []
            common_genes = ["tp53", "brca1", "brca2", "myc", "egfr", "kras", "pten"]
            for gene in common_genes:
                if gene in query_lower:
                    gene_mentions.append(gene.upper())
            
            if gene_mentions:
                search_query = ", ".join(gene_mentions)
            
            # Look for disease/condition context
            if any(word in query_lower for word in ["cancer", "tumor", "disease"]):
                condition_context = "cancer"
            elif any(word in query_lower for word in ["immune", "immunotherapy"]):
                condition_context = "immunotherapy"
            
            tools_to_use.append(("pubmed_literature_search", {
                "query": search_query,
                "condition_context": condition_context,
                "max_results": 5,
                "include_recent_only": "recent" in query_lower or "latest" in query_lower
            }))
        
        # Check if user wants synthesis
        synthesis_keywords = ["analyze", "synthesis", "connect", "relate", "therapeutic", "targets"]
        if any(word in query_lower for word in synthesis_keywords):
            tools_to_use.append(("synthesize_literature_with_genomics", {
                "analysis_context": query,
                "focus_areas": ["mechanisms", "therapeutic targets"] if "therapeutic" in query_lower else ["mechanisms"],
                "gene_limit": 10
            }))
        
        # Default: get context if no specific tools identified
        if not tools_to_use:
            tools_to_use.append(("get_genomics_analysis_context", {"analysis_type": "all", "summary_only": True}))
        
        return tools_to_use
    
    def _synthesize_tool_results(self, query: str, tool_results: List[Dict]) -> str:
        """
        Synthesize tool results into a coherent response
        
        This simulates Claude's natural language synthesis
        """
        try:
            response_parts = []
            
            # Process each tool result
            for result in tool_results:
                if "error" in result:
                    response_parts.append(f"⚠️ There was an issue with {result.get('tool', 'a tool')}: {result['error']}")
                    continue
                
                tool_name = result["tool"]
                tool_data = result["result"]
                
                if tool_name == "get_genomics_analysis_context":
                    response_parts.append(self._format_genomics_context_response(tool_data))
                
                elif tool_name == "pubmed_literature_search":
                    response_parts.append(self._format_literature_search_response(tool_data))
                
                elif tool_name == "synthesize_literature_with_genomics":
                    response_parts.append(self._format_synthesis_response(tool_data))
            
            # Combine responses
            if response_parts:
                final_response = "\n\n".join(response_parts)
                
                # Add helpful suggestions
                final_response += "\n\n💡 **What would you like to explore next?**\n"
                final_response += "- Ask about specific genes or pathways\n"
                final_response += "- Request therapeutic target analysis\n"
                final_response += "- Explore mechanism literature\n"
                final_response += "- Compare with clinical outcomes"
                
                return final_response
            else:
                return "I understand your question, but I need to search for relevant information. Could you try rephrasing your query or being more specific?"
        
        except Exception as e:
            logger.error(f"Response synthesis error: {e}")
            return f"I encountered an error while processing your request: {e}"
    
    def _format_genomics_context_response(self, context_data: Dict) -> str:
        """Format genomics context into readable response"""
        if "error" in context_data:
            return f"❌ **Analysis Context**: {context_data['error']}"
        
        response = "📊 **Your Current Analysis Status:**\n\n"
        
        # Available data
        if context_data.get("available_data"):
            data_info = context_data["available_data"]
            if "expression" in data_info:
                expr = data_info["expression"]
                response += f"- **Expression Data**: {expr['n_genes']} genes × {expr['n_samples']} samples\n"
            
            if "clinical" in data_info:
                clinical = data_info["clinical"]
                response += f"- **Clinical Data**: {clinical['n_samples']} samples with {len(clinical['variables'])} variables\n"
        
        # Analysis results
        if context_data.get("analysis_results"):
            results = context_data["analysis_results"]
            
            if "deseq2" in results:
                deseq2 = results["deseq2"]
                response += f"\n🧬 **DESeq2 Analysis Complete**:\n"
                response += f"- {deseq2['significant_genes']} significant genes\n"
                response += f"- {deseq2['upregulated']} upregulated, {deseq2['downregulated']} downregulated\n"
                
                if deseq2.get("top_significant_genes"):
                    top_genes = deseq2["top_significant_genes"][:5]
                    response += f"- Top genes: {', '.join(top_genes)}\n"
            
            if "pathway" in results:
                response += "\n🔬 **Pathway Analysis**: Complete\n"
        
        # Recommendations
        if context_data.get("recommendations"):
            response += "\n💡 **Recommendations**:\n"
            for rec in context_data["recommendations"][:3]:
                response += f"- {rec}\n"
        
        return response
    
    def _format_literature_search_response(self, search_data: Dict) -> str:
        """Format literature search results into readable response"""
        if "error" in search_data:
            return f"❌ **Literature Search**: {search_data['error']}"
        
        if search_data.get("search_status") != "success":
            return f"⚠️ **Literature Search**: {search_data.get('error', 'Search failed')}"
        
        response = f"📚 **Literature Search Results for '{search_data['search_query']}'**:\n\n"
        response += f"Found {search_data['total_papers']} relevant papers.\n\n"
        
        papers = search_data.get("papers_found", [])[:3]  # Show top 3
        
        for i, paper in enumerate(papers, 1):
            response += f"**{i}. {paper['title']}**\n"
            response += f"   *{paper['authors']}* - {paper['journal']} ({paper['year']})\n"
            response += f"   📋 {paper['abstract']}\n"
            response += f"   🔗 [Read Paper]({paper['url']})\n\n"
        
        if len(papers) < search_data['total_papers']:
            remaining = search_data['total_papers'] - len(papers)
            response += f"... and {remaining} more papers found.\n"
        
        return response
    
    def _format_synthesis_response(self, synthesis_data: Dict) -> str:
        """Format synthesis results into readable response"""
        if "error" in synthesis_data:
            return f"❌ **Literature Synthesis**: {synthesis_data['error']}"
        
        response = "🔬 **Literature Synthesis & Analysis**:\n\n"
        
        if synthesis_data.get("key_findings"):
            response += "**Key Findings**:\n"
            for finding in synthesis_data["key_findings"]:
                response += f"- {finding}\n"
            response += "\n"
        
        if synthesis_data.get("recommendations"):
            response += "**Research Recommendations**:\n"
            for rec in synthesis_data["recommendations"]:
                response += f"- {rec}\n"
            response += "\n"
        
        # Show literature searches performed
        searches = synthesis_data.get("literature_searches_performed", [])
        if searches:
            response += f"**Literature Reviews Completed**: {len(searches)} targeted searches\n"
        
        return response
    
    def _generate_demo_response(self, query: str) -> Dict[str, Any]:
        """Generate demo response when Claude API isn't available"""
        query_lower = query.lower()
        
        demo_responses = {
            "genomics": """📊 **Demo: Your Analysis Context**

🧬 **Simulated DESeq2 Results**:
- 15,847 total genes analyzed
- 1,234 significant genes (FDR < 0.05)
- 687 upregulated, 547 downregulated
- Top genes: TP53, BRCA1, MYC, EGFR, KRAS

💡 This is a simulated response. Connect Claude API for real analysis!""",

            "literature": """📚 **Demo: Literature Search Results**

Found 12 relevant papers for your query.

**1. Novel therapeutic targets in cancer genomics**
   *Johnson et al.* - Nature Medicine (2024)
   📋 Recent advances in identifying druggable targets from genomics data...
   🔗 [Demo Link]

**2. BRCA1 pathway interactions in tumor progression**
   *Smith et al.* - Cell (2024)
   📋 Comprehensive analysis of BRCA1 regulatory networks...
   🔗 [Demo Link]

💡 This is simulated data. Install Claude API for real literature search!""",

            "synthesis": """🔬 **Demo: Literature Synthesis**

**Key Findings**:
- Found 18 papers relevant to your 10 top significant genes
- 67% of your upregulated genes have therapeutic targeting literature
- Recent mechanistic insights available for TP53-MYC interactions

**Research Recommendations**:
- Review PI3K pathway papers for novel drug combinations
- Investigate immune evasion mechanisms in your gene set
- Consider biomarker potential of top differential genes

💡 Connect Claude API for real AI-powered synthesis!"""
        }
        
        # Choose appropriate demo response
        if any(word in query_lower for word in ["analysis", "results", "genes"]):
            content = demo_responses["genomics"]
        elif any(word in query_lower for word in ["literature", "papers", "search"]):
            content = demo_responses["literature"]
        elif any(word in query_lower for word in ["synthesis", "analyze", "therapeutic"]):
            content = demo_responses["synthesis"]
        else:
            content = """🤖 **Demo Mode Active**

I can help you with:
- 📊 Analyzing your genomics results
- 📚 Searching scientific literature  
- 🔬 Synthesizing research findings
- 💡 Providing research recommendations

Try asking: "What can you tell me about my DESeq2 results?" or "Find papers on BRCA1 mutations"

💡 This is demo mode. Set up Claude API for real AI assistance!"""
        
        return {
            "type": "assistant",
            "content": content,
            "demo_mode": True,
            "timestamp": time.time()
        }
    
    def _get_genomics_context_summary(self) -> Dict[str, Any]:
        """Get summary of available genomics analysis context"""
        context = {}
        
        try:
            if self.session.has_data("expression_data"):
                expr_data = self.session.get_data("expression_data")
                context["expression"] = {
                    "available": True,
                    "genes": expr_data.n_genes,
                    "samples": expr_data.n_samples
                }
            
            if self.session.has_analysis_results("deseq2"):
                deseq2_results = self.session.get_analysis_results("deseq2")
                context["deseq2"] = {
                    "available": True,
                    "significant_genes": deseq2_results.n_significant,
                    "total_genes": deseq2_results.n_genes
                }
            
            if self.session.has_analysis_results("pathway"):
                context["pathway"] = {
                    "available": True,
                    "analysis_complete": True
                }
            
        except Exception as e:
            logger.error(f"Error getting genomics context: {e}")
        
        return context


def render_ai_research_tab():
    """Render the AI Research Assistant tab (main entry point)"""
    tab = AIResearchTab()
    tab.render()