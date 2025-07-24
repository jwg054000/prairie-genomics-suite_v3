"""
Prairie Genomics Suite - Session Management

This module manages Streamlit session state in a centralized, organized way.
It provides a clean interface for storing and retrieving analysis data,
similar to how you'd manage your R workspace environment.
"""

import streamlit as st
from typing import Any, Dict, List, Optional, Type, TypeVar, Tuple
from datetime import datetime
import uuid

from .data_models import (
    AnalysisSession, 
    ExpressionData, 
    ClinicalData, 
    SurvivalData,
    DESeq2Results,
    PathwayResults,
    ImmuneResults,
    LiteratureResults,
    AnalysisType
)
from config import get_config

T = TypeVar('T')


class SessionManager:
    """
    Centralized session state management for the Prairie Genomics Suite.
    
    This class provides a clean interface to Streamlit's session state,
    organizing data and results in a structured way. Think of it as
    your R workspace manager.
    """
    
    def __init__(self):
        """Initialize the session manager"""
        self._initialize_session()
    
    def _initialize_session(self):
        """Initialize session state with default values"""
        
        # Core session info
        if 'session_id' not in st.session_state:
            st.session_state.session_id = str(uuid.uuid4())
            
        if 'session_created' not in st.session_state:
            st.session_state.session_created = datetime.now()
            
        # Analysis session container
        if 'analysis_session' not in st.session_state:
            st.session_state.analysis_session = AnalysisSession(
                session_id=st.session_state.session_id,
                created_at=st.session_state.session_created
            )
        
        # UI state
        if 'current_tab' not in st.session_state:
            st.session_state.current_tab = "Data Import"
            
        if 'analysis_settings' not in st.session_state:
            st.session_state.analysis_settings = get_config("analysis")
            
        # Data loading flags
        if 'data_loaded' not in st.session_state:
            st.session_state.data_loaded = {
                'expression': False,
                'clinical': False,
                'survival': False
            }
            
        # Analysis completion flags
        if 'analysis_completed' not in st.session_state:
            st.session_state.analysis_completed = {
                'deseq2': False,
                'pathway': False,
                'immune': False,
                'literature': False
            }
            
        # Caching for performance
        if 'cached_results' not in st.session_state:
            st.session_state.cached_results = {}
            
        # Error tracking
        if 'error_log' not in st.session_state:
            st.session_state.error_log = []
    
    @property 
    def session(self) -> AnalysisSession:
        """Get the current analysis session"""
        return st.session_state.analysis_session
    
    @property
    def session_id(self) -> str:
        """Get the current session ID"""
        return st.session_state.session_id
    
    # Data Management Methods
    
    def set_expression_data(self, data: ExpressionData):
        """Set expression data and update session"""
        self.session.expression_data = data
        st.session_state.data_loaded['expression'] = True
        self._log_action("expression_data_loaded", {
            "n_genes": data.n_genes,
            "n_samples": data.n_samples
        })
    
    def get_expression_data(self) -> Optional[ExpressionData]:
        """Get expression data from session"""
        return self.session.expression_data
    
    def set_clinical_data(self, data: ClinicalData):
        """Set clinical data and update session"""
        self.session.clinical_data = data
        st.session_state.data_loaded['clinical'] = True
        self._log_action("clinical_data_loaded", {
            "n_samples": data.n_samples,
            "n_variables": data.n_variables
        })
    
    def get_clinical_data(self) -> Optional[ClinicalData]:
        """Get clinical data from session"""
        return self.session.clinical_data
    
    def set_survival_data(self, data: SurvivalData):
        """Set survival data and update session"""
        self.session.survival_data = data
        st.session_state.data_loaded['survival'] = True
        self._log_action("survival_data_loaded", {
            "n_samples": len(data.samples)
        })
    
    def get_survival_data(self) -> Optional[SurvivalData]:
        """Get survival data from session"""
        return self.session.survival_data
    
    # Analysis Results Management
    
    def set_analysis_results(self, analysis_type: AnalysisType, results: Any, 
                           parameters: Dict[str, Any] = None):
        """Set analysis results and update session"""
        self.session.add_analysis(analysis_type, results, parameters)
        st.session_state.analysis_completed[analysis_type.value] = True
        self._log_action(f"{analysis_type.value}_completed", parameters)
    
    def get_deseq2_results(self) -> Optional[DESeq2Results]:
        """Get DESeq2 results from session"""
        return self.session.deseq2_results
    
    def get_pathway_results(self) -> Optional[PathwayResults]:
        """Get pathway analysis results from session"""
        return self.session.pathway_results
    
    def get_immune_results(self) -> Optional[ImmuneResults]:
        """Get immune infiltration results from session"""
        return self.session.immune_results
    
    def get_literature_results(self) -> Optional[LiteratureResults]:
        """Get literature search results from session"""
        return self.session.literature_results
    
    # Data availability check methods
    
    def has_data(self, data_type: str) -> bool:
        """Check if specific data type is available"""
        if data_type == "expression_data":
            return self.session.expression_data is not None
        elif data_type == "clinical_data":
            return self.session.clinical_data is not None
        elif data_type == "survival_data":
            return self.session.survival_data is not None
        else:
            return False
    
    def has_analysis_results(self, analysis_type: str) -> bool:
        """Check if specific analysis results are available"""
        if analysis_type == "deseq2":
            return self.session.deseq2_results is not None
        elif analysis_type == "pathway":
            return self.session.pathway_results is not None
        elif analysis_type == "immune":
            return self.session.immune_results is not None
        elif analysis_type == "literature":
            return self.session.literature_results is not None
        else:
            return False
    
    def store_data(self, data_type: str, data: Any):
        """Store data in session (generic method)"""
        if data_type == "expression_data":
            self.set_expression_data(data)
        elif data_type == "clinical_data":
            self.set_clinical_data(data)
        else:
            # Store as generic data
            if not hasattr(st.session_state, 'custom_data'):
                st.session_state.custom_data = {}
            st.session_state.custom_data[data_type] = data
    
    def get_data(self, data_type: str) -> Any:
        """Get data from session (generic method)"""
        if data_type == "expression_data":
            return self.get_expression_data()
        elif data_type == "clinical_data":
            return self.get_clinical_data()
        else:
            # Get custom data
            custom_data = getattr(st.session_state, 'custom_data', {})
            return custom_data.get(data_type)
    
    def store_analysis_results(self, analysis_type: str, results: Any):
        """Store analysis results in session"""
        if analysis_type == "deseq2":
            self.session.deseq2_results = results
        elif analysis_type == "pathway":
            self.session.pathway_results = results
        elif analysis_type == "literature":
            self.session.literature_results = results
        else:
            # Store as custom results
            if not hasattr(st.session_state, 'custom_results'):
                st.session_state.custom_results = {}
            st.session_state.custom_results[analysis_type] = results
    
    def get_analysis_results(self, analysis_type: str) -> Any:
        """Get analysis results from session"""
        if analysis_type == "deseq2":
            return self.get_deseq2_results()
        elif analysis_type == "pathway":
            return self.get_pathway_results()
        elif analysis_type == "literature":
            return self.get_literature_results()
        else:
            # Get custom results
            custom_results = getattr(st.session_state, 'custom_results', {})
            return custom_results.get(analysis_type)
    
    def store_analysis_params(self, param_key: str, params: Dict[str, Any]):
        """Store analysis parameters"""
        if not hasattr(st.session_state, 'analysis_params'):
            st.session_state.analysis_params = {}
        st.session_state.analysis_params[param_key] = params
    
    def get_analysis_params(self, param_key: str) -> Optional[Dict[str, Any]]:
        """Get analysis parameters"""
        params = getattr(st.session_state, 'analysis_params', {})
        return params.get(param_key)
    
    def has_analysis_params(self, param_key: str) -> bool:
        """Check if analysis parameters exist"""
        params = getattr(st.session_state, 'analysis_params', {})
        return param_key in params
    
    # UI State Management
    
    def set_current_tab(self, tab_name: str):
        """Set the current active tab"""
        st.session_state.current_tab = tab_name
        self._log_action("tab_changed", {"tab": tab_name})
    
    def get_current_tab(self) -> str:
        """Get the current active tab"""
        return st.session_state.current_tab
    
    def set_analysis_setting(self, analysis_type: str, setting: str, value: Any):
        """Set an analysis setting"""
        if analysis_type not in st.session_state.analysis_settings:
            st.session_state.analysis_settings[analysis_type] = {}
        st.session_state.analysis_settings[analysis_type][setting] = value
    
    def get_analysis_setting(self, analysis_type: str, setting: str, default: Any = None) -> Any:
        """Get an analysis setting"""
        return st.session_state.analysis_settings.get(analysis_type, {}).get(setting, default)
    
    # Data Validation and Status
    
    def is_data_loaded(self, data_type: str) -> bool:
        """Check if a specific data type is loaded"""
        return st.session_state.data_loaded.get(data_type, False)
    
    def is_analysis_completed(self, analysis_type: str) -> bool:
        """Check if a specific analysis is completed"""
        return st.session_state.analysis_completed.get(analysis_type, False)
    
    def get_data_status(self) -> Dict[str, bool]:
        """Get status of all data loading"""
        return st.session_state.data_loaded.copy()
    
    def get_analysis_status(self) -> Dict[str, bool]:
        """Get status of all analyses"""
        return st.session_state.analysis_completed.copy()
    
    def can_run_analysis(self, analysis_type: str) -> Tuple[bool, List[str]]:
        """
        Check if an analysis can be run based on data availability
        
        Returns:
            Tuple of (can_run, missing_requirements)
        """
        requirements = {
            'deseq2': ['expression', 'clinical'],
            'pathway': ['expression'],  # Can use DE results if available
            'immune': ['expression'],
            'literature': [],  # Can run independently
            'survival': ['clinical']  # May also need expression for stratification
        }
        
        required_data = requirements.get(analysis_type, [])
        missing = [req for req in required_data if not self.is_data_loaded(req)]
        
        return len(missing) == 0, missing
    
    # Caching Support
    
    def cache_result(self, key: str, value: Any, ttl_seconds: int = 3600):
        """Cache a result with optional TTL"""
        cache_entry = {
            'value': value,
            'timestamp': datetime.now(),
            'ttl': ttl_seconds
        }
        st.session_state.cached_results[key] = cache_entry
    
    def get_cached_result(self, key: str) -> Optional[Any]:
        """Get a cached result if still valid"""
        if key not in st.session_state.cached_results:
            return None
            
        cache_entry = st.session_state.cached_results[key]
        
        # Check if expired
        if cache_entry['ttl'] > 0:
            elapsed = (datetime.now() - cache_entry['timestamp']).total_seconds()
            if elapsed > cache_entry['ttl']:
                del st.session_state.cached_results[key]
                return None
        
        return cache_entry['value']
    
    def clear_cache(self):
        """Clear all cached results"""
        st.session_state.cached_results = {}
    
    # Error Handling
    
    def log_error(self, error: str, context: Dict[str, Any] = None):
        """Log an error with context"""
        error_entry = {
            'timestamp': datetime.now(),
            'error': error,
            'context': context or {},
            'session_id': self.session_id
        }
        st.session_state.error_log.append(error_entry)
        
        # Keep only last 50 errors to prevent memory issues
        if len(st.session_state.error_log) > 50:
            st.session_state.error_log = st.session_state.error_log[-50:]
    
    def get_errors(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent errors"""
        return st.session_state.error_log[-limit:]
    
    def clear_errors(self):
        """Clear error log"""
        st.session_state.error_log = []
    
    # Session Management
    
    def reset_session(self):
        """Reset the entire session (useful for starting over)"""
        # Clear all session state except basic Streamlit keys
        keys_to_keep = ['session_id']
        keys_to_remove = [key for key in st.session_state.keys() if key not in keys_to_keep]
        
        for key in keys_to_remove:
            del st.session_state[key]
            
        # Reinitialize
        self._initialize_session()
    
    def clear_session(self):
        """Alias for reset_session for compatibility"""
        self.reset_session()
    
    def export_session_summary(self) -> Dict[str, Any]:
        """Export a summary of the current session"""
        return {
            'session_info': {
                'session_id': self.session_id,
                'created_at': st.session_state.session_created,
                'current_tab': self.get_current_tab()
            },
            'data_status': self.get_data_status(),
            'analysis_status': self.get_analysis_status(),
            'session_summary': self.session.get_analysis_summary(),
            'error_count': len(st.session_state.error_log)
        }
    
    # Private helper methods
    
    def _log_action(self, action: str, context: Dict[str, Any] = None):
        """Log an action for debugging/analytics"""
        # This could be extended to send to external analytics
        pass


# Global session manager instance
_session_manager = None

def get_session_manager() -> SessionManager:
    """Get the global session manager instance"""
    global _session_manager
    if _session_manager is None:
        _session_manager = SessionManager()
    return _session_manager


# Convenience functions for common operations
def get_expression_data() -> Optional[ExpressionData]:
    """Convenience function to get expression data"""
    return get_session_manager().get_expression_data()

def get_clinical_data() -> Optional[ClinicalData]:
    """Convenience function to get clinical data"""
    return get_session_manager().get_clinical_data()

def is_ready_for_analysis(analysis_type: str) -> bool:
    """Convenience function to check if ready for analysis"""
    can_run, _ = get_session_manager().can_run_analysis(analysis_type)
    return can_run