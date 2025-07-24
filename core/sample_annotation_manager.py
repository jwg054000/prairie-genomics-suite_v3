"""
Prairie Genomics Suite - Sample Annotation Manager

This module provides flexible sample annotation capabilities, allowing users 
to perform differential expression analysis without requiring formal clinical 
data files. It supports interactive sample grouping, pattern detection, and 
automatic annotation strategies.

Features:
- Interactive sample-to-group assignment
- Pattern-based automatic grouping
- Manual annotation workflows
- Integration with existing clinical data
- Export to standard clinical data format

Author: Prairie Genomics Team
"""

import pandas as pd
import numpy as np
import re
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
from datetime import datetime
import logging

from core.data_models import ExpressionData, ClinicalData

# Set up logging
logger = logging.getLogger(__name__)


class SampleAnnotationManager:
    """
    Manages flexible sample annotation workflows for differential expression analysis.
    
    This class provides multiple ways to assign samples to experimental groups:
    1. Interactive manual assignment
    2. Pattern-based automatic detection
    3. Template-based common designs
    4. Integration with existing clinical data
    """
    
    def __init__(self, expression_data: Union[ExpressionData, pd.DataFrame]):
        """
        Initialize with expression data
        
        Args:
            expression_data: Expression data container or DataFrame
        """
        if isinstance(expression_data, pd.DataFrame):
            self.expression_df = expression_data
            self.samples = expression_data.columns.tolist()
        else:
            self.expression_df = expression_data.data
            self.samples = expression_data.samples
            
        self.n_samples = len(self.samples)
        self.annotations = {}  # sample_id -> group mapping
        self.annotation_method = None
        self.created_at = datetime.now()
        
        logger.info(f"Initialized SampleAnnotationManager with {self.n_samples} samples")
    
    def get_samples(self) -> List[str]:
        """Get list of all sample names"""
        return self.samples.copy()
    
    def get_annotations(self) -> Dict[str, str]:
        """Get current sample annotations"""
        return self.annotations.copy()
    
    def get_groups(self) -> List[str]:
        """Get list of unique groups"""
        return list(set(self.annotations.values()))
    
    def get_group_sizes(self) -> Dict[str, int]:
        """Get sample count for each group"""
        if not self.annotations:
            return {}
        
        groups = {}
        for group in self.annotations.values():
            groups[group] = groups.get(group, 0) + 1
        return groups
    
    def clear_annotations(self):
        """Clear all current annotations"""
        self.annotations = {}
        self.annotation_method = None
        logger.info("Cleared all sample annotations")
    
    def assign_sample(self, sample_id: str, group: str):
        """
        Assign a single sample to a group
        
        Args:
            sample_id: Sample identifier
            group: Group name
        """
        if sample_id not in self.samples:
            raise ValueError(f"Sample '{sample_id}' not found in expression data")
        
        self.annotations[sample_id] = group
        logger.debug(f"Assigned {sample_id} to group '{group}'")
    
    def assign_samples(self, assignments: Dict[str, str]):
        """
        Assign multiple samples to groups
        
        Args:
            assignments: Dictionary mapping sample_id -> group
        """
        for sample_id, group in assignments.items():
            self.assign_sample(sample_id, group)
        
        logger.info(f"Assigned {len(assignments)} samples to groups")
    
    def suggest_groupings(self) -> Dict[str, Dict[str, Any]]:
        """
        Suggest potential groupings based on sample name patterns
        
        Returns:
            Dictionary of suggested grouping strategies with metadata
        """
        suggestions = {}
        
        # Pattern-based suggestions
        pattern_suggestion = self._suggest_pattern_grouping()
        if pattern_suggestion:
            suggestions['pattern'] = pattern_suggestion
        
        # Numeric-based suggestions
        numeric_suggestion = self._suggest_numeric_grouping()
        if numeric_suggestion:
            suggestions['numeric'] = numeric_suggestion
        
        # Alternating suggestions
        alternating_suggestion = self._suggest_alternating_grouping()
        suggestions['alternating'] = alternating_suggestion
        
        # Split-half suggestions
        split_suggestion = self._suggest_split_grouping()
        suggestions['split_half'] = split_suggestion
        
        logger.info(f"Generated {len(suggestions)} grouping suggestions")
        return suggestions
    
    def _suggest_pattern_grouping(self) -> Optional[Dict[str, Any]]:
        """Suggest groupings based on common naming patterns"""
        # Common patterns to look for
        patterns = [
            (r'(.+)[_\-](control|ctrl|con)([_\-].*)?$', 'Control'),
            (r'(.+)[_\-](treatment|treat|trt|exp)([_\-].*)?$', 'Treatment'),
            (r'(.+)[_\-](case|patient|disease)([_\-].*)?$', 'Case'),
            (r'(.+)[_\-](normal|healthy|norm)([_\-].*)?$', 'Normal'),
            (r'(.+)[_\-](tumor|cancer|malignant)([_\-].*)?$', 'Tumor'),
            (r'(.+)[_\-](before|pre|baseline)([_\-].*)?$', 'Before'),
            (r'(.+)[_\-](after|post|followup)([_\-].*)?$', 'After'),
        ]
        
        best_match = None
        best_count = 0
        
        for pattern, group_name in patterns:
            matches = []
            non_matches = []
            
            for sample in self.samples:
                if re.search(pattern, sample, re.IGNORECASE):
                    matches.append(sample)
                else:
                    non_matches.append(sample)
            
            # Check if this pattern captures a reasonable portion of samples
            if len(matches) >= 2 and len(non_matches) >= 2 and len(matches) > best_count:
                best_count = len(matches)
                
                # Create grouping
                assignments = {}
                for sample in matches:
                    assignments[sample] = group_name
                
                # Assign remaining samples to "Other" or infer opposite group
                other_group = self._infer_opposite_group(group_name)
                for sample in non_matches:
                    assignments[sample] = other_group
                
                confidence = min(0.95, len(matches) / self.n_samples * 2)  # Scale confidence
                
                best_match = {
                    'assignments': assignments,
                    'confidence': confidence,
                    'description': f"Pattern-based: {group_name} vs {other_group}",
                    'method': 'pattern',
                    'pattern': pattern,
                    'matched_samples': len(matches),
                    'groups': {group_name: len(matches), other_group: len(non_matches)}
                }
        
        return best_match
    
    def _suggest_numeric_grouping(self) -> Optional[Dict[str, Any]]:
        """Suggest groupings based on numeric patterns in sample names"""
        # Extract numbers from sample names
        sample_numbers = {}
        number_pattern = r'(\d+)'
        
        for sample in self.samples:
            numbers = re.findall(number_pattern, sample)
            if numbers:
                # Use the last number found (often the most meaningful)
                sample_numbers[sample] = int(numbers[-1])
        
        if len(sample_numbers) < 4:  # Need at least 4 samples with numbers
            return None
        
        # Try different splitting strategies
        numbers = list(sample_numbers.values())
        numbers.sort()
        
        # Split at median
        median_val = np.median(numbers)
        assignments = {}
        
        for sample, number in sample_numbers.items():
            if number <= median_val:
                assignments[sample] = "Group_Low"
            else:
                assignments[sample] = "Group_High"
        
        # Assign samples without numbers to a separate group or distribute evenly
        remaining_samples = [s for s in self.samples if s not in sample_numbers]
        for i, sample in enumerate(remaining_samples):
            group = "Group_Low" if i % 2 == 0 else "Group_High"
            assignments[sample] = group
        
        group_sizes = {}
        for group in assignments.values():
            group_sizes[group] = group_sizes.get(group, 0) + 1
        
        # Calculate confidence based on balance and coverage
        coverage = len(sample_numbers) / self.n_samples
        balance = min(group_sizes.values()) / max(group_sizes.values()) if group_sizes else 0
        confidence = (coverage * 0.7 + balance * 0.3) * 0.8  # Cap at 0.8 for numeric
        
        return {
            'assignments': assignments,
            'confidence': confidence,
            'description': f"Numeric split at {median_val}",
            'method': 'numeric',
            'split_value': median_val,
            'groups': group_sizes
        }
    
    def _suggest_alternating_grouping(self) -> Dict[str, Any]:
        """Suggest alternating group assignment"""
        assignments = {}
        
        for i, sample in enumerate(self.samples):
            group = "Group_A" if i % 2 == 0 else "Group_B"
            assignments[sample] = group
        
        return {
            'assignments': assignments,
            'confidence': 0.6,  # Moderate confidence for alternating
            'description': "Alternating assignment (A, B, A, B, ...)",
            'method': 'alternating',
            'groups': {"Group_A": (self.n_samples + 1) // 2, "Group_B": self.n_samples // 2}
        }
    
    def _suggest_split_grouping(self) -> Dict[str, Any]:
        """Suggest first-half vs second-half grouping"""
        mid_point = self.n_samples // 2
        assignments = {}
        
        for i, sample in enumerate(self.samples):
            group = "First_Half" if i < mid_point else "Second_Half"
            assignments[sample] = group
        
        return {
            'assignments': assignments,
            'confidence': 0.5,  # Lower confidence for arbitrary split
            'description': "Split into first half vs second half",
            'method': 'split_half',
            'groups': {"First_Half": mid_point, "Second_Half": self.n_samples - mid_point}
        }
    
    def _infer_opposite_group(self, group_name: str) -> str:
        """Infer the opposite group name for common terms"""
        opposites = {
            'control': 'treatment',
            'ctrl': 'treatment', 
            'con': 'treatment',
            'treatment': 'control',
            'treat': 'control',
            'trt': 'control',
            'case': 'control',
            'patient': 'control',
            'disease': 'normal',
            'normal': 'disease',
            'healthy': 'disease',
            'tumor': 'normal',
            'cancer': 'normal',
            'malignant': 'normal',
            'before': 'after',
            'pre': 'post',
            'baseline': 'followup',
            'after': 'before',
            'post': 'pre',
            'followup': 'baseline'
        }
        
        group_lower = group_name.lower()
        return opposites.get(group_lower, 'Other').title()
    
    def apply_suggestion(self, suggestion_key: str, suggestions: Dict[str, Any] = None):
        """
        Apply a suggested grouping
        
        Args:
            suggestion_key: Key of the suggestion to apply
            suggestions: Optional suggestions dict (will generate if not provided)
        """
        if suggestions is None:
            suggestions = self.suggest_groupings()
        
        if suggestion_key not in suggestions:
            raise ValueError(f"Suggestion '{suggestion_key}' not found")
        
        suggestion = suggestions[suggestion_key]
        self.annotations = suggestion['assignments'].copy()
        self.annotation_method = suggestion['method']
        
        logger.info(f"Applied {suggestion_key} grouping suggestion to {len(self.annotations)} samples")
    
    def create_template_annotation(self, template_name: str, **kwargs) -> Dict[str, str]:
        """
        Create annotations from predefined templates
        
        Args:
            template_name: Name of the template
            **kwargs: Template-specific parameters
            
        Returns:
            Dictionary of sample assignments
        """
        if template_name == 'case_control':
            return self._create_case_control_template(**kwargs)
        elif template_name == 'time_series':
            return self._create_time_series_template(**kwargs)
        elif template_name == 'dose_response':
            return self._create_dose_response_template(**kwargs)
        elif template_name == 'paired':
            return self._create_paired_template(**kwargs)
        else:
            raise ValueError(f"Unknown template: {template_name}")
    
    def _create_case_control_template(self, case_ratio: float = 0.5) -> Dict[str, str]:
        """Create case-control design"""
        n_cases = int(self.n_samples * case_ratio)
        assignments = {}
        
        for i, sample in enumerate(self.samples):
            group = "Case" if i < n_cases else "Control"
            assignments[sample] = group
        
        return assignments
    
    def _create_time_series_template(self, time_points: List[str] = None) -> Dict[str, str]:
        """Create time series design"""
        if time_points is None:
            n_points = min(4, self.n_samples // 2)
            time_points = [f"T{i}" for i in range(n_points)]
        
        assignments = {}
        points_per_time = self.n_samples // len(time_points)
        
        for i, sample in enumerate(self.samples):
            time_idx = i // points_per_time
            if time_idx >= len(time_points):
                time_idx = len(time_points) - 1
            assignments[sample] = time_points[time_idx]
        
        return assignments
    
    def _create_dose_response_template(self, doses: List[str] = None) -> Dict[str, str]:
        """Create dose-response design"""
        if doses is None:
            doses = ["Low", "Medium", "High"]
        
        assignments = {}
        samples_per_dose = self.n_samples // len(doses)
        
        for i, sample in enumerate(self.samples):
            dose_idx = i // samples_per_dose
            if dose_idx >= len(doses):
                dose_idx = len(doses) - 1
            assignments[sample] = doses[dose_idx]
        
        return assignments
    
    def _create_paired_template(self, **kwargs) -> Dict[str, str]:
        """Create paired design (before/after)"""
        if self.n_samples % 2 != 0:
            raise ValueError("Paired design requires even number of samples")
        
        assignments = {}
        for i, sample in enumerate(self.samples):
            group = "Before" if i < self.n_samples // 2 else "After"
            assignments[sample] = group
        
        return assignments
    
    def validate_annotations(self) -> Tuple[bool, List[str]]:
        """
        Validate current annotations
        
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        # Check if all samples are annotated
        unannotated = set(self.samples) - set(self.annotations.keys())
        if unannotated:
            errors.append(f"Unannotated samples: {list(unannotated)[:5]}")
        
        # Check for at least 2 groups
        groups = self.get_groups()
        if len(groups) < 2:
            errors.append(f"Need at least 2 groups, found {len(groups)}")
        
        # Check group sizes
        group_sizes = self.get_group_sizes()
        min_size = min(group_sizes.values()) if group_sizes else 0
        if min_size < 2:
            errors.append(f"All groups need at least 2 samples, minimum found: {min_size}")
        
        # Check for reasonable balance (optional warning)
        if group_sizes:
            max_size = max(group_sizes.values())
            if min_size < max_size * 0.3:  # Very imbalanced
                errors.append(f"Groups are very imbalanced: {group_sizes}")
        
        return len(errors) == 0, errors
    
    def export_clinical_format(self, condition_column: str = "Condition") -> pd.DataFrame:
        """
        Export annotations as clinical data format
        
        Args:
            condition_column: Name for the condition column
            
        Returns:
            Clinical data DataFrame
        """
        if not self.annotations:
            raise ValueError("No annotations to export")
        
        # Create clinical data DataFrame
        clinical_data = pd.DataFrame({
            condition_column: [self.annotations.get(sample, "Unknown") for sample in self.samples],
            'sample_id': self.samples,
            'annotation_method': self.annotation_method or 'manual',
            'created_at': self.created_at.isoformat(),
            'created_by': 'sample_annotation_manager'
        }, index=self.samples)
        
        logger.info(f"Exported annotations to clinical format: {clinical_data.shape}")
        return clinical_data
    
    def export_to_clinical_data(self, condition_column: str = "Condition") -> ClinicalData:
        """
        Export as ClinicalData object
        
        Args:
            condition_column: Name for the condition column
            
        Returns:
            ClinicalData object
        """
        clinical_df = self.export_clinical_format(condition_column)
        
        return ClinicalData(
            data=clinical_df,
            samples=clinical_df.index.tolist(),
            metadata={
                'source': 'sample_annotation_manager',
                'method': self.annotation_method,
                'created_at': self.created_at.isoformat(),
                'n_groups': len(self.get_groups()),
                'group_sizes': self.get_group_sizes()
            }
        )
    
    def save_annotations(self, filepath: Union[str, Path]):
        """Save annotations to file"""
        filepath = Path(filepath)
        
        # Export as CSV
        clinical_df = self.export_clinical_format()
        clinical_df.to_csv(filepath)
        
        logger.info(f"Saved annotations to {filepath}")
    
    def load_annotations(self, filepath: Union[str, Path], condition_column: str = "Condition"):
        """Load annotations from file"""
        filepath = Path(filepath)
        
        if not filepath.exists():
            raise FileNotFoundError(f"Annotation file not found: {filepath}")
        
        # Load CSV
        clinical_df = pd.read_csv(filepath, index_col=0)
        
        if condition_column not in clinical_df.columns:
            raise ValueError(f"Condition column '{condition_column}' not found in file")
        
        # Extract annotations
        self.annotations = clinical_df[condition_column].to_dict()
        self.annotation_method = clinical_df.get('annotation_method', {}).get(self.samples[0], 'loaded')
        
        logger.info(f"Loaded annotations from {filepath}: {len(self.annotations)} samples")
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary of current annotation state"""
        return {
            'n_samples': self.n_samples,
            'n_annotated': len(self.annotations),
            'n_groups': len(self.get_groups()),
            'groups': self.get_groups(),
            'group_sizes': self.get_group_sizes(),
            'annotation_method': self.annotation_method,
            'is_valid': self.validate_annotations()[0],
            'created_at': self.created_at.isoformat()
        }


def create_sample_annotator(expression_data: Union[ExpressionData, pd.DataFrame]) -> SampleAnnotationManager:
    """
    Convenience function to create a sample annotation manager
    
    Args:
        expression_data: Expression data container or DataFrame
        
    Returns:
        SampleAnnotationManager instance
    """
    return SampleAnnotationManager(expression_data)