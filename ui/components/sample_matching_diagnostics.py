"""
Sample Matching Diagnostics Component

This component provides detailed information about sample matching between
expression and clinical data, helping users debug sample ID mismatch issues.
"""

import streamlit as st
import pandas as pd
from typing import List, Dict, Tuple, Optional
import logging

from core.session_manager import get_session_manager
from core.data_models import ExpressionData, ClinicalData

logger = logging.getLogger(__name__)


def render_sample_matching_diagnostics(key: str = "sample_diagnostics") -> Optional[Dict]:
    """
    Render comprehensive sample matching diagnostics UI component.
    
    Args:
        key: Unique key for the component
        
    Returns:
        Dictionary with diagnostic information if data is available
    """
    st.subheader("🔍 Sample Matching Diagnostics")
    st.write("This tool helps identify and resolve sample ID mismatches between expression and clinical data.")
    
    # Get data from session
    session = get_session_manager()
    
    if not session.has_data("expression_data"):
        st.info("📤 Upload expression data first to see sample diagnostics.")
        return None
    
    if not session.has_data("clinical_data"):
        st.info("📤 Upload clinical data or create sample annotations to see matching diagnostics.")
        return None
    
    # Get the data
    expression_data = session.get_data("expression_data")
    clinical_data = session.get_data("clinical_data")
    
    # Extract sample information
    expr_samples = expression_data.samples if hasattr(expression_data, 'samples') else list(expression_data.data.columns)
    clinical_samples = clinical_data.samples if hasattr(clinical_data, 'samples') else list(clinical_data.data.index)
    
    # Create diagnostic information
    diagnostics = analyze_sample_matching(expr_samples, clinical_samples)
    
    # Display overview
    st.write("### 📊 Sample Overview")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("🧬 Expression Samples", len(expr_samples))
    
    with col2:
        st.metric("🏥 Clinical Samples", len(clinical_samples))
    
    with col3:
        st.metric("✅ Common Samples", diagnostics['n_common'])
    
    # Display matching status
    if diagnostics['n_common'] >= 6:
        st.success(f"🎉 **Excellent!** {diagnostics['n_common']} samples match perfectly. Ready for analysis!")
    elif diagnostics['n_common'] >= 3:
        st.warning(f"⚠️ **Acceptable** - {diagnostics['n_common']} samples match. Analysis possible but limited statistical power.")
    else:
        st.error(f"❌ **Problem!** Only {diagnostics['n_common']} samples match. Need at least 3 for analysis.")
    
    # Detailed analysis tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🔍 Sample Lists", "🔗 Matching Analysis", "🔧 Suggested Fixes", "📋 Export"])
    
    with tab1:
        render_sample_lists(expr_samples, clinical_samples, diagnostics)
    
    with tab2:
        render_matching_analysis(diagnostics)
    
    with tab3:
        render_suggested_fixes(diagnostics, expr_samples, clinical_samples)
    
    with tab4:
        render_export_options(expr_samples, clinical_samples, diagnostics)
    
    return diagnostics


def analyze_sample_matching(expr_samples: List[str], clinical_samples: List[str]) -> Dict:
    """
    Analyze sample matching between expression and clinical data.
    
    Args:
        expr_samples: List of expression sample names
        clinical_samples: List of clinical sample names
        
    Returns:
        Dictionary with diagnostic information
    """
    expr_set = set(expr_samples)
    clinical_set = set(clinical_samples)
    
    # Basic overlap analysis
    common = expr_set & clinical_set
    expr_only = expr_set - clinical_set
    clinical_only = clinical_set - expr_set
    
    # Check for common issues
    has_whitespace_issue = False
    has_case_issue = False
    has_character_issue = False
    
    if len(expr_samples) == len(clinical_samples) and len(common) == 0:
        # Check whitespace
        stripped_expr = set(str(s).strip() for s in expr_samples)
        stripped_clinical = set(str(s).strip() for s in clinical_samples)
        
        if stripped_expr == stripped_clinical:
            has_whitespace_issue = True
        
        # Check case sensitivity
        elif set(s.lower() for s in expr_samples) == set(s.lower() for s in clinical_samples):
            has_case_issue = True
        
        # Check common character substitutions
        else:
            expr_normalized = set(s.replace('-', '_').replace(' ', '_') for s in expr_samples)
            clinical_normalized = set(s.replace('-', '_').replace(' ', '_') for s in clinical_samples)
            
            if expr_normalized == clinical_normalized:
                has_character_issue = True
    
    # Calculate similarity scores
    similarities = calculate_sample_similarities(expr_samples, clinical_samples)
    
    return {
        'n_expr': len(expr_samples),
        'n_clinical': len(clinical_samples),
        'n_common': len(common),
        'n_expr_only': len(expr_only),
        'n_clinical_only': len(clinical_only),
        'common_samples': sorted(list(common)),
        'expr_only_samples': sorted(list(expr_only)),
        'clinical_only_samples': sorted(list(clinical_only)),
        'has_whitespace_issue': has_whitespace_issue,
        'has_case_issue': has_case_issue,
        'has_character_issue': has_character_issue,
        'similarities': similarities,
        'same_count': len(expr_samples) == len(clinical_samples)
    }


def calculate_sample_similarities(expr_samples: List[str], clinical_samples: List[str]) -> List[Dict]:
    """Calculate similarity scores between unmatched samples."""
    from difflib import SequenceMatcher
    
    similarities = []
    
    # Find best matches for unmatched clinical samples
    expr_set = set(expr_samples)
    clinical_set = set(clinical_samples)
    
    clinical_only = clinical_set - expr_set
    expr_only = expr_set - clinical_set
    
    for clinical_sample in clinical_only:
        best_match = None
        best_score = 0.0
        
        for expr_sample in expr_only:
            score = SequenceMatcher(None, str(clinical_sample), str(expr_sample)).ratio()
            
            if score > best_score:
                best_match = expr_sample
                best_score = score
        
        if best_match and best_score > 0.5:  # Only show reasonably similar matches
            similarities.append({
                'clinical': clinical_sample,
                'expression': best_match,
                'similarity': best_score
            })
    
    return sorted(similarities, key=lambda x: x['similarity'], reverse=True)


def render_sample_lists(expr_samples: List[str], clinical_samples: List[str], diagnostics: Dict):
    """Render detailed sample lists with type information."""
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**🧬 Expression Data Samples:**")
        
        # Show sample types
        if expr_samples:
            sample_types = [type(s).__name__ for s in expr_samples[:3]]
            st.write(f"*Types: {', '.join(set(sample_types))}*")
        
        # Display samples with status indicators
        with st.expander(f"View all {len(expr_samples)} expression samples", expanded=len(expr_samples) <= 10):
            for i, sample in enumerate(expr_samples):
                if sample in diagnostics['common_samples']:
                    st.write(f"{i+1:2d}. ✅ `{sample}`")
                else:
                    st.write(f"{i+1:2d}. ⚠️ `{sample}` (no match)")
    
    with col2:
        st.write("**🏥 Clinical Data Samples:**")
        
        # Show sample types
        if clinical_samples:
            sample_types = [type(s).__name__ for s in clinical_samples[:3]]
            st.write(f"*Types: {', '.join(set(sample_types))}*")
        
        # Display samples with status indicators
        with st.expander(f"View all {len(clinical_samples)} clinical samples", expanded=len(clinical_samples) <= 10):
            for i, sample in enumerate(clinical_samples):
                if sample in diagnostics['common_samples']:
                    st.write(f"{i+1:2d}. ✅ `{sample}`")
                else:
                    st.write(f"{i+1:2d}. ⚠️ `{sample}` (no match)")


def render_matching_analysis(diagnostics: Dict):
    """Render detailed matching analysis."""
    st.write("**🔗 Sample Matching Results:**")
    
    # Summary metrics
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("✅ Matching", diagnostics['n_common'])
        if diagnostics['common_samples']:
            with st.expander("View matching samples"):
                for sample in diagnostics['common_samples']:
                    st.write(f"• `{sample}`")
    
    with col2:
        st.metric("🧬 Expression Only", diagnostics['n_expr_only'])
        if diagnostics['expr_only_samples']:
            with st.expander("View expression-only samples"):
                for sample in diagnostics['expr_only_samples']:
                    st.write(f"• `{sample}`")
    
    with col3:
        st.metric("🏥 Clinical Only", diagnostics['n_clinical_only'])
        if diagnostics['clinical_only_samples']:
            with st.expander("View clinical-only samples"):
                for sample in diagnostics['clinical_only_samples']:
                    st.write(f"• `{sample}`")
    
    # Similarity analysis for potential matches
    if diagnostics['similarities']:
        st.write("**🔍 Potential Matches (Similar Names):**")
        
        for sim in diagnostics['similarities'][:5]:  # Show top 5 matches
            confidence = "High" if sim['similarity'] > 0.8 else "Medium" if sim['similarity'] > 0.6 else "Low"
            st.write(f"• `{sim['clinical']}` ↔ `{sim['expression']}` ({sim['similarity']:.1%} similarity, {confidence} confidence)")


def render_suggested_fixes(diagnostics: Dict, expr_samples: List[str], clinical_samples: List[str]):
    """Render suggested fixes for sample matching issues."""
    st.write("**🔧 Suggested Solutions:**")
    
    if diagnostics['n_common'] >= 6:
        st.success("✅ No fixes needed - sample matching is excellent!")
        return
    
    fixes = []
    
    # Issue-specific suggestions
    if diagnostics['has_whitespace_issue']:
        fixes.append({
            'title': '🧹 Fix Whitespace Issues',
            'description': 'Remove extra spaces from sample names in your data files',
            'severity': 'High',
            'icon': '🟢'
        })
    
    if diagnostics['has_case_issue']:
        fixes.append({
            'title': '🔤 Fix Case Sensitivity',
            'description': 'Ensure consistent capitalization in sample names',
            'severity': 'High',
            'icon': '🟢'
        })
    
    if diagnostics['has_character_issue']:
        fixes.append({
            'title': '🔤 Standardize Characters',
            'description': 'Replace hyphens with underscores or vice versa consistently',
            'severity': 'Medium',
            'icon': '🟡'
        })
    
    # General suggestions
    if diagnostics['same_count'] and diagnostics['n_common'] == 0:
        fixes.append({
            'title': '📝 Check Sample Names',
            'description': 'You have the same number of samples - likely a naming consistency issue',
            'severity': 'High',
            'icon': '🟢'
        })
    
    if diagnostics['similarities']:
        fixes.append({
            'title': '🔍 Review Similar Names',
            'description': f'Found {len(diagnostics["similarities"])} potentially matching pairs - check for typos',
            'severity': 'Medium',
            'icon': '🟡'
        })
    
    # Always suggest the annotation tool
    fixes.append({
        'title': '🏷️ Use Sample Annotation Tool',
        'description': 'Create clinical data directly from your expression data samples',
        'severity': 'Recommended',
        'icon': '⭐'
    })
    
    fixes.append({
        'title': '📁 Verify Data Files',
        'description': 'Ensure both expression and clinical files contain the same samples',
        'severity': 'Basic',
        'icon': '📋'
    })
    
    # Display fixes
    for fix in fixes:
        with st.container():
            col1, col2 = st.columns([1, 10])
            
            with col1:
                st.write(fix['icon'])
            
            with col2:
                st.write(f"**{fix['title']}** ({fix['severity']})")
                st.write(fix['description'])
            
            st.write("")  # Add spacing


def render_export_options(expr_samples: List[str], clinical_samples: List[str], diagnostics: Dict):
    """Render options to export diagnostic information."""
    st.write("**📋 Export Diagnostic Information:**")
    
    # Create diagnostic report
    report = create_diagnostic_report(expr_samples, clinical_samples, diagnostics)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.download_button(
            label="📄 Download Diagnostic Report",
            data=report,
            file_name="sample_matching_diagnostics.txt",
            mime="text/plain",
            help="Download a detailed report of sample matching analysis"
        )
    
    with col2:
        if st.button("📋 Copy to Clipboard"):
            st.write("```")
            st.code(report, language=None)
            st.write("```")
            st.info("Copy the text above to share with support or for documentation")


def create_diagnostic_report(expr_samples: List[str], clinical_samples: List[str], diagnostics: Dict) -> str:
    """Create a detailed diagnostic report as text."""
    from datetime import datetime
    
    report = f"""
SAMPLE MATCHING DIAGNOSTIC REPORT
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

=== SUMMARY ===
Expression samples: {diagnostics['n_expr']}
Clinical samples: {diagnostics['n_clinical']}
Common samples: {diagnostics['n_common']}
Expression-only: {diagnostics['n_expr_only']}
Clinical-only: {diagnostics['n_clinical_only']}

Analysis Status: {'✅ READY' if diagnostics['n_common'] >= 6 else '⚠️ LIMITED' if diagnostics['n_common'] >= 3 else '❌ INSUFFICIENT'}

=== DETECTED ISSUES ===
Whitespace issues: {'Yes' if diagnostics['has_whitespace_issue'] else 'No'}
Case sensitivity issues: {'Yes' if diagnostics['has_case_issue'] else 'No'}
Character differences: {'Yes' if diagnostics['has_character_issue'] else 'No'}
Same sample count: {'Yes' if diagnostics['same_count'] else 'No'}

=== SAMPLE LISTS ===

Expression samples:
{chr(10).join(f'  {i+1:2d}. {sample}' for i, sample in enumerate(expr_samples))}

Clinical samples:
{chr(10).join(f'  {i+1:2d}. {sample}' for i, sample in enumerate(clinical_samples))}

=== MATCHING ANALYSIS ===

Common samples ({len(diagnostics['common_samples'])}):
{chr(10).join(f'  • {sample}' for sample in diagnostics['common_samples'])}

Expression-only samples ({len(diagnostics['expr_only_samples'])}):
{chr(10).join(f'  • {sample}' for sample in diagnostics['expr_only_samples'])}

Clinical-only samples ({len(diagnostics['clinical_only_samples'])}):
{chr(10).join(f'  • {sample}' for sample in diagnostics['clinical_only_samples'])}

=== POTENTIAL MATCHES ===
{chr(10).join(f'  • {sim["clinical"]} ↔ {sim["expression"]} ({sim["similarity"]:.1%})' for sim in diagnostics['similarities'][:10])}

=== RECOMMENDATIONS ===
1. If whitespace issues detected: Remove extra spaces from sample names
2. If case issues detected: Ensure consistent capitalization
3. If character issues detected: Standardize use of hyphens, underscores, spaces
4. Consider using the Sample Annotation tool to create matching clinical data
5. Verify that both files contain the same set of samples

=== END OF REPORT ===
"""
    
    return report.strip()