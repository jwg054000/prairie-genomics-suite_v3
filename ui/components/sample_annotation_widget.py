"""
Prairie Genomics Suite - Sample Annotation Widget

This module provides Streamlit UI components for interactive sample annotation.
Users can assign samples to groups using drag-and-drop interface, apply 
automatic suggestions, and validate their annotations before analysis.

Features:
- Interactive sample assignment interface
- Automatic grouping suggestions
- Visual validation and preview
- Template-based designs
- Export capabilities

Author: Prairie Genomics Team
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List, Optional, Any, Union
from pathlib import Path
import logging

from core.sample_annotation_manager import SampleAnnotationManager, create_sample_annotator
from core.data_models import ExpressionData, ClinicalData

# Set up logging
logger = logging.getLogger(__name__)


def render_sample_annotation_widget(expression_data: Union[ExpressionData, pd.DataFrame], 
                                  key: str = "sample_annotation") -> Optional[ClinicalData]:
    """
    Render the main sample annotation widget
    
    Args:
        expression_data: Expression data container or DataFrame
        key: Unique key for the widget
        
    Returns:
        ClinicalData object if annotations are complete, None otherwise
    """
    st.subheader("📝 Sample Annotation")
    st.write("Assign your samples to experimental groups for differential expression analysis.")
    
    # Initialize annotation manager
    if f"annotator_{key}" not in st.session_state:
        st.session_state[f"annotator_{key}"] = create_sample_annotator(expression_data)
    
    annotator = st.session_state[f"annotator_{key}"]
    
    # Create tabs for different annotation methods
    tab1, tab2, tab3, tab4 = st.tabs([
        "🤖 Smart Suggestions", 
        "✏️ Manual Assignment", 
        "📋 Templates",
        "📊 Validation"
    ])
    
    with tab1:
        render_smart_suggestions(annotator, key)
    
    with tab2:
        render_manual_assignment(annotator, key)
    
    with tab3:
        render_template_assignment(annotator, key)
    
    with tab4:
        clinical_data = render_validation_and_export(annotator, key)
        return clinical_data
    
    return None


def render_smart_suggestions(annotator: SampleAnnotationManager, key: str):
    """Render smart grouping suggestions"""
    st.subheader("🧠 Smart Grouping Suggestions")
    st.write("AI-powered analysis of your sample names to suggest optimal groupings.")
    
    # Generate suggestions button
    if st.button("🔍 Analyze Sample Names", key=f"analyze_{key}"):
        with st.spinner("Analyzing sample naming patterns..."):
            suggestions = annotator.suggest_groupings()
            st.session_state[f"suggestions_{key}"] = suggestions
    
    # Display suggestions if available
    if f"suggestions_{key}" in st.session_state:
        suggestions = st.session_state[f"suggestions_{key}"]
        
        if not suggestions:
            st.info("No clear patterns detected in sample names. Try manual assignment or templates.")
            return
        
        st.write(f"Found **{len(suggestions)}** potential grouping strategies:")
        
        # Create columns for suggestions
        cols = st.columns(min(len(suggestions), 3))
        
        for i, (suggestion_name, suggestion_data) in enumerate(suggestions.items()):
            col = cols[i % len(cols)]
            
            with col:
                # Create suggestion card
                confidence = suggestion_data.get('confidence', 0)
                confidence_color = "green" if confidence > 0.7 else "orange" if confidence > 0.4 else "red"
                
                st.markdown(f"""
                <div style="border: 1px solid #ddd; padding: 10px; border-radius: 5px; margin-bottom: 10px;">
                    <h4>{suggestion_data.get('description', suggestion_name.title())}</h4>
                    <p><strong>Confidence:</strong> <span style="color: {confidence_color};">{confidence:.1%}</span></p>
                    <p><strong>Groups:</strong> {list(suggestion_data.get('groups', {}).keys())}</p>
                    <p><strong>Sizes:</strong> {list(suggestion_data.get('groups', {}).values())}</p>
                </div>
                """, unsafe_allow_html=True)
                
                # Apply suggestion button
                if st.button(f"Apply {suggestion_name.title()}", 
                           key=f"apply_{suggestion_name}_{key}",
                           help=f"Apply this grouping strategy (confidence: {confidence:.1%})"):
                    annotator.apply_suggestion(suggestion_name, suggestions)
                    st.success(f"Applied {suggestion_name} grouping!")
                    st.rerun()
        
        # Show preview of best suggestion
        if suggestions:
            best_suggestion_name = max(suggestions.keys(), 
                                     key=lambda x: suggestions[x].get('confidence', 0))
            best_suggestion = suggestions[best_suggestion_name]
            
            st.subheader("📊 Preview: Best Suggestion")
            st.write(f"**{best_suggestion.get('description')}** (Confidence: {best_suggestion.get('confidence', 0):.1%})")
            
            # Create visualization
            preview_df = pd.DataFrame([
                {'Sample': sample, 'Group': group} 
                for sample, group in best_suggestion['assignments'].items()
            ])
            
            fig = px.histogram(preview_df, x='Group', 
                             title="Sample Distribution by Group",
                             color='Group')
            fig.update_layout(showlegend=False, height=300)
            st.plotly_chart(fig, use_container_width=True)


def render_manual_assignment(annotator: SampleAnnotationManager, key: str):
    """Render manual sample assignment interface"""
    st.subheader("✏️ Manual Sample Assignment")
    st.write("Manually assign each sample to experimental groups.")
    
    samples = annotator.get_samples()
    current_annotations = annotator.get_annotations()
    
    # Group management
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.write("**Define your experimental groups:**")
        
        # Add new group
        new_group = st.text_input("Add new group:", key=f"new_group_{key}")
        if st.button("➕ Add Group", key=f"add_group_{key}") and new_group:
            if f"groups_{key}" not in st.session_state:
                st.session_state[f"groups_{key}"] = ["Control", "Treatment"]
            if new_group not in st.session_state[f"groups_{key}"]:
                st.session_state[f"groups_{key}"].append(new_group)
                st.success(f"Added group: {new_group}")
                st.rerun()
    
    with col2:
        # Initialize default groups
        if f"groups_{key}" not in st.session_state:
            st.session_state[f"groups_{key}"] = ["Control", "Treatment"]
        
        current_groups = st.session_state[f"groups_{key}"]
        st.write("**Current groups:**")
        for group in current_groups:
            count = sum(1 for g in current_annotations.values() if g == group)
            st.write(f"• {group} ({count} samples)")
    
    # Sample assignment interface
    st.write("**Assign samples to groups:**")
    
    # Create assignment form
    with st.form(key=f"assignment_form_{key}"):
        assignments = {}
        
        # Group samples into chunks for better UI
        chunk_size = 10
        sample_chunks = [samples[i:i+chunk_size] for i in range(0, len(samples), chunk_size)]
        
        for chunk_idx, sample_chunk in enumerate(sample_chunks):
            if len(sample_chunks) > 1:
                st.write(f"**Samples {chunk_idx * chunk_size + 1} - {min((chunk_idx + 1) * chunk_size, len(samples))}:**")
            
            # Create columns for this chunk
            cols = st.columns(min(len(sample_chunk), 3))
            
            for i, sample in enumerate(sample_chunk):
                col = cols[i % len(cols)]
                
                with col:
                    current_group = current_annotations.get(sample, "")
                    selected_group = st.selectbox(
                        f"{sample}:",
                        [""] + current_groups,
                        index=current_groups.index(current_group) + 1 if current_group in current_groups else 0,
                        key=f"assign_{sample}_{key}"
                    )
                    
                    if selected_group:
                        assignments[sample] = selected_group
        
        # Submit assignments
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.form_submit_button("💾 Save Assignments"):
                annotator.assign_samples(assignments)
                st.success(f"Saved assignments for {len(assignments)} samples!")
                st.rerun()
        
        with col2:
            if st.form_submit_button("🗑️ Clear All"):
                annotator.clear_annotations()
                st.success("Cleared all assignments!")
                st.rerun()
        
        with col3:
            # Auto-balance button
            if st.form_submit_button("⚖️ Auto Balance"):
                _auto_balance_groups(annotator, current_groups)
                st.success("Auto-balanced groups!")
                st.rerun()


def render_template_assignment(annotator: SampleAnnotationManager, key: str):
    """Render template-based assignment"""
    st.subheader("📋 Template Assignment")
    st.write("Apply common experimental design templates.")
    
    # Template selection
    template_options = {
        "case_control": "Case vs Control",
        "time_series": "Time Series",
        "dose_response": "Dose Response",
        "paired": "Paired Design (Before/After)"
    }
    
    selected_template = st.selectbox(
        "Select a template:",
        list(template_options.keys()),
        format_func=lambda x: template_options[x],
        key=f"template_select_{key}"
    )
    
    # Template-specific parameters
    template_params = {}
    
    if selected_template == "case_control":
        st.write("**Case-Control Design**")
        template_params['case_ratio'] = st.slider(
            "Proportion of cases:", 0.1, 0.9, 0.5, 0.1,
            key=f"case_ratio_{key}"
        )
        st.info(f"Will assign {int(annotator.n_samples * template_params['case_ratio'])} samples to 'Case' and the rest to 'Control'")
    
    elif selected_template == "time_series":
        st.write("**Time Series Design**")
        n_timepoints = st.number_input(
            "Number of time points:", 2, 8, 4,
            key=f"timepoints_{key}"
        )
        template_params['time_points'] = [f"T{i}" for i in range(n_timepoints)]
        st.info(f"Will create {n_timepoints} time points: {', '.join(template_params['time_points'])}")
    
    elif selected_template == "dose_response":
        st.write("**Dose Response Design**")
        dose_input = st.text_input(
            "Dose levels (comma-separated):", "Low,Medium,High",
            key=f"doses_{key}"
        )
        template_params['doses'] = [d.strip() for d in dose_input.split(',') if d.strip()]
        st.info(f"Will create dose groups: {', '.join(template_params['doses'])}")
    
    elif selected_template == "paired":
        st.write("**Paired Design**")
        if annotator.n_samples % 2 != 0:
            st.error("Paired design requires an even number of samples!")
        else:
            st.info("Will assign first half to 'Before' and second half to 'After'")
    
    # Apply template button
    if st.button(f"Apply {template_options[selected_template]}", key=f"apply_template_{key}"):
        try:
            assignments = annotator.create_template_annotation(selected_template, **template_params)
            annotator.assign_samples(assignments)
            annotator.annotation_method = f"template_{selected_template}"
            st.success(f"Applied {template_options[selected_template]} template!")
            st.rerun()
        except Exception as e:
            st.error(f"Error applying template: {str(e)}")


def render_validation_and_export(annotator: SampleAnnotationManager, key: str) -> Optional[ClinicalData]:
    """Render validation and export interface"""
    st.subheader("📊 Validation & Export")
    
    # Get current status
    summary = annotator.get_summary()
    annotations = annotator.get_annotations()
    
    # Status overview
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Samples", summary['n_samples'])
    
    with col2:
        st.metric("Annotated", summary['n_annotated'])
    
    with col3:
        st.metric("Groups", summary['n_groups'])
    
    with col4:
        completion = summary['n_annotated'] / summary['n_samples'] if summary['n_samples'] > 0 else 0
        st.metric("Progress", f"{completion:.1%}")
    
    # Validation results
    if annotations:
        is_valid, errors = annotator.validate_annotations()
        
        if is_valid:
            st.success("✅ Annotations are valid and ready for analysis!")
        else:
            st.error("❌ Annotations have issues that need to be fixed:")
            for error in errors:
                st.write(f"• {error}")
    
    # Visualization
    if annotations:
        st.subheader("📈 Group Distribution")
        
        # Create distribution chart
        group_sizes = annotator.get_group_sizes()
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Bar chart
            fig_bar = px.bar(
                x=list(group_sizes.keys()), 
                y=list(group_sizes.values()),
                title="Samples per Group",
                labels={'x': 'Group', 'y': 'Number of Samples'}
            )
            fig_bar.update_layout(showlegend=False, height=300)
            st.plotly_chart(fig_bar, use_container_width=True)
        
        with col2:
            # Pie chart
            fig_pie = px.pie(
                values=list(group_sizes.values()),
                names=list(group_sizes.keys()),
                title="Group Proportions"
            )
            fig_pie.update_layout(height=300)
            st.plotly_chart(fig_pie, use_container_width=True)
        
        # Sample assignments table
        if st.checkbox("Show detailed assignments", key=f"show_details_{key}"):
            assignment_df = pd.DataFrame([
                {'Sample': sample, 'Group': group}
                for sample, group in annotations.items()
            ])
            st.dataframe(assignment_df, use_container_width=True)
    
    # Export options
    if annotations and annotator.validate_annotations()[0]:
        st.subheader("📤 Export for Analysis")
        
        condition_column = st.text_input(
            "Condition column name:", "Condition",
            key=f"condition_col_{key}"
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🧬 Use for DE Analysis", key=f"export_analysis_{key}"):
                try:
                    clinical_data = annotator.export_to_clinical_data(condition_column)
                    st.success("✅ Annotations ready for differential expression analysis!")
                    
                    # Show preview
                    with st.expander("Preview clinical data"):
                        st.dataframe(clinical_data.data.head())
                    
                    return clinical_data
                
                except Exception as e:
                    st.error(f"Export error: {str(e)}")
        
        with col2:
            # Download as CSV
            clinical_df = annotator.export_clinical_format(condition_column)
            csv_data = clinical_df.to_csv()
            
            st.download_button(
                label="💾 Download CSV",
                data=csv_data,
                file_name=f"sample_annotations_{key}.csv",
                mime="text/csv",
                key=f"download_csv_{key}"
            )
    
    else:
        st.info("Complete your sample annotations to export for analysis.")
    
    return None


def _auto_balance_groups(annotator: SampleAnnotationManager, groups: List[str]):
    """Auto-balance samples across groups"""
    samples = annotator.get_samples()
    
    # Simple round-robin assignment
    assignments = {}
    for i, sample in enumerate(samples):
        group = groups[i % len(groups)]
        assignments[sample] = group
    
    annotator.assign_samples(assignments)
    annotator.annotation_method = "auto_balanced"


def render_annotation_summary_card(annotator: SampleAnnotationManager) -> None:
    """Render a summary card of current annotations"""
    summary = annotator.get_summary()
    
    st.markdown(f"""
    <div style="border: 1px solid #ddd; padding: 15px; border-radius: 10px; margin: 10px 0;">
        <h3>📝 Annotation Summary</h3>
        <ul>
            <li><strong>Samples:</strong> {summary['n_annotated']}/{summary['n_samples']} annotated</li>
            <li><strong>Groups:</strong> {summary['n_groups']} ({', '.join(summary['groups'])})</li>
            <li><strong>Method:</strong> {summary['annotation_method'] or 'Manual'}</li>
            <li><strong>Status:</strong> {'✅ Valid' if summary['is_valid'] else '❌ Needs attention'}</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)


# Example usage function for testing
def demo_sample_annotation_widget():
    """Demo function to test the sample annotation widget"""
    # Create fake expression data
    np.random.seed(42)
    n_genes, n_samples = 100, 20
    
    gene_names = [f"Gene_{i:03d}" for i in range(n_genes)]
    sample_names = [f"Sample_{i:02d}" for i in range(n_samples)]
    
    expression_data = pd.DataFrame(
        np.random.negative_binomial(10, 0.3, size=(n_genes, n_samples)),
        index=gene_names,
        columns=sample_names
    )
    
    st.title("🧬 Sample Annotation Demo")
    
    clinical_data = render_sample_annotation_widget(expression_data, "demo")
    
    if clinical_data:
        st.success("🎉 Ready for differential expression analysis!")
        st.json(clinical_data.metadata)