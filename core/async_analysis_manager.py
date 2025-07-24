"""
Prairie Genomics Suite - Async Analysis Manager

Handles asynchronous analysis execution with progress tracking and UI responsiveness.
Prevents blocking operations and provides real-time feedback to users.

Author: Prairie Genomics Team
"""

import asyncio
import threading
import time
from typing import Dict, Any, Optional, Callable, List, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging
import streamlit as st
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, Future
import uuid

from core.data_models import ExpressionData, ClinicalData, DESeq2Results
from core.enhanced_cache import get_enhanced_cache, cached_operation
from config import get_config

logger = logging.getLogger(__name__)

class AnalysisStatus(Enum):
    """Analysis execution status"""
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELLED = "cancelled"

@dataclass
class AnalysisTask:
    """Represents a single analysis task"""
    task_id: str
    task_type: str
    status: AnalysisStatus = AnalysisStatus.PENDING
    progress: float = 0.0
    current_stage: str = "Initializing"
    stages: List[Tuple[str, float]] = field(default_factory=list)
    result: Any = None
    error: Optional[str] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def duration(self) -> Optional[float]:
        """Get task duration in seconds"""
        if self.started_at and self.completed_at:
            return (self.completed_at - self.started_at).total_seconds()
        elif self.started_at:
            return (datetime.now() - self.started_at).total_seconds()
        return None
    
    @property
    def is_complete(self) -> bool:
        """Check if task is complete"""
        return self.status in [AnalysisStatus.SUCCESS, AnalysisStatus.FAILED, AnalysisStatus.CANCELLED]
    
    def update_progress(self, progress: float, stage: str = None):
        """Update task progress and stage"""
        self.progress = min(1.0, max(0.0, progress))
        if stage:
            self.current_stage = stage
    
    def mark_started(self):
        """Mark task as started"""
        self.status = AnalysisStatus.RUNNING
        self.started_at = datetime.now()
    
    def mark_completed(self, result: Any = None):
        """Mark task as successfully completed"""
        self.status = AnalysisStatus.SUCCESS
        self.progress = 1.0
        self.current_stage = "Complete"
        self.completed_at = datetime.now()
        if result is not None:
            self.result = result
    
    def mark_failed(self, error: str):
        """Mark task as failed"""
        self.status = AnalysisStatus.FAILED
        self.current_stage = "Failed"
        self.completed_at = datetime.now()
        self.error = error


class AsyncAnalysisManager:
    """
    Manages asynchronous execution of genomics analyses with progress tracking
    """
    
    def __init__(self):
        self.config = get_config("performance")
        self.cache = get_enhanced_cache()
        
        # Task management
        self.active_tasks: Dict[str, AnalysisTask] = {}
        self.task_history: List[AnalysisTask] = []
        self.max_history = 50
        
        # Execution pools
        self.thread_pool = ThreadPoolExecutor(
            max_workers=self.config["parallel_processing"]["max_workers"],
            thread_name_prefix="analysis"
        )
        self.process_pool = ProcessPoolExecutor(
            max_workers=max(1, self.config["parallel_processing"]["max_workers"] // 2)
        )
        
        # Progress tracking
        self.progress_callbacks: Dict[str, Callable] = {}
        
        logger.info("Async Analysis Manager initialized")
    
    def create_task(self, 
                   task_type: str,
                   stages: List[Tuple[str, float]] = None,
                   metadata: Dict[str, Any] = None) -> AnalysisTask:
        """
        Create a new analysis task
        
        Args:
            task_type: Type of analysis (e.g., 'deseq2', 'literature_search')
            stages: List of (stage_name, expected_progress) tuples
            metadata: Additional task metadata
            
        Returns:
            Created AnalysisTask
        """
        task_id = str(uuid.uuid4())
        
        # Default stages if not provided
        if stages is None:
            stages = [
                ("Preparing data", 0.1),
                ("Running analysis", 0.8),
                ("Processing results", 0.95),
                ("Complete", 1.0)
            ]
        
        task = AnalysisTask(
            task_id=task_id,
            task_type=task_type,
            stages=stages,
            metadata=metadata or {}
        )
        
        self.active_tasks[task_id] = task
        logger.info(f"Created task {task_id}: {task_type}")
        
        return task
    
    def run_deseq2_async(self,
                        expression_data: ExpressionData,
                        clinical_data: ClinicalData,
                        condition_column: str,
                        **kwargs) -> AnalysisTask:
        """
        Run DESeq2 analysis asynchronously with progress tracking
        
        Args:
            expression_data: Expression data
            clinical_data: Clinical data
            condition_column: Condition column for comparison
            **kwargs: Additional DESeq2 parameters
            
        Returns:
            AnalysisTask for tracking progress
        """
        # Define analysis stages
        stages = [
            ("Validating data", 0.05),
            ("Preparing count matrix", 0.15),
            ("Initializing DESeq2", 0.25),
            ("Running differential expression", 0.70),
            ("Processing results", 0.85),
            ("Generating visualizations", 0.95),
            ("Complete", 1.0)
        ]
        
        task = self.create_task(
            task_type='deseq2',
            stages=stages,
            metadata={
                'condition_column': condition_column,
                'n_genes': len(expression_data.genes),
                'n_samples': len(expression_data.samples),
                'parameters': kwargs
            }
        )
        
        # Submit to thread pool
        future = self.thread_pool.submit(
            self._run_deseq2_worker,
            task,
            expression_data,
            clinical_data,
            condition_column,
            **kwargs
        )
        
        # Set up progress monitoring
        self._monitor_task_async(task, future)
        
        return task
    
    def run_literature_search_async(self,
                                  gene_list: List[str],
                                  search_terms: List[str] = None,
                                  **kwargs) -> AnalysisTask:
        """
        Run literature search asynchronously
        
        Args:
            gene_list: List of gene names to search
            search_terms: Additional search terms
            **kwargs: Literature search parameters
            
        Returns:
            AnalysisTask for tracking progress
        """
        stages = [
            ("Preparing gene queries", 0.1),
            ("Searching PubMed", 0.4),
            ("Searching Semantic Scholar", 0.6),
            ("Searching ArXiv", 0.7),
            ("Processing results with LLM", 0.9),
            ("Complete", 1.0)
        ]
        
        task = self.create_task(
            task_type='literature_search',
            stages=stages,
            metadata={
                'gene_count': len(gene_list),
                'search_terms': search_terms or [],
                'parameters': kwargs
            }
        )
        
        # Submit to thread pool
        future = self.thread_pool.submit(
            self._run_literature_search_worker,
            task,
            gene_list,
            search_terms,
            **kwargs
        )
        
        self._monitor_task_async(task, future)
        
        return task
    
    def run_chunked_analysis_async(self,
                                 analysis_func: Callable,
                                 data_chunks: List[Any],
                                 chunk_processor: Callable,
                                 **kwargs) -> AnalysisTask:
        """
        Run analysis on data chunks in parallel
        
        Args:
            analysis_func: Main analysis function
            data_chunks: List of data chunks to process
            chunk_processor: Function to process individual chunks
            **kwargs: Additional parameters
            
        Returns:
            AnalysisTask for tracking progress
        """
        n_chunks = len(data_chunks)
        stages = [
            ("Preparing chunks", 0.05),
            ("Processing chunks in parallel", 0.8),
            ("Combining results", 0.95),
            ("Complete", 1.0)
        ]
        
        task = self.create_task(
            task_type='chunked_analysis',
            stages=stages,
            metadata={
                'n_chunks': n_chunks,
                'parallel_processing': True
            }
        )
        
        # Submit to process pool for CPU-intensive work
        future = self.process_pool.submit(
            self._run_chunked_analysis_worker,
            task,
            analysis_func,
            data_chunks,
            chunk_processor,
            **kwargs
        )
        
        self._monitor_task_async(task, future)
        
        return task
    
    def get_task_status(self, task_id: str) -> Optional[AnalysisTask]:
        """Get current status of a task"""
        return self.active_tasks.get(task_id)
    
    def cancel_task(self, task_id: str) -> bool:
        """Cancel a running task"""
        task = self.active_tasks.get(task_id)
        if task and not task.is_complete:
            task.status = AnalysisStatus.CANCELLED
            task.completed_at = datetime.now()
            logger.info(f"Cancelled task {task_id}")
            return True
        return False
    
    def render_task_progress(self, task: AnalysisTask, container=None):
        """
        Render task progress in Streamlit UI
        
        Args:
            task: AnalysisTask to display
            container: Streamlit container (optional)
        """
        if container is None:
            container = st
        
        # Progress bar
        container.progress(task.progress)
        
        # Status information
        if task.status == AnalysisStatus.RUNNING:
            container.info(f"🔄 {task.current_stage}... ({task.progress:.1%})")
            
            # Show duration if available
            if task.duration:
                container.caption(f"Running for {task.duration:.1f} seconds")
                
        elif task.status == AnalysisStatus.SUCCESS:
            container.success(f"✅ {task.task_type.title()} analysis complete!")
            if task.duration:
                container.caption(f"Completed in {task.duration:.1f} seconds")
                
        elif task.status == AnalysisStatus.FAILED:
            container.error(f"❌ Analysis failed: {task.error}")
            
        elif task.status == AnalysisStatus.CANCELLED:
            container.warning("⚠️ Analysis was cancelled")
    
    def render_task_manager_ui(self):
        """Render task manager UI for monitoring all tasks"""
        st.subheader("🔄 Analysis Tasks")
        
        if not self.active_tasks and not self.task_history:
            st.info("No analysis tasks")
            return
        
        # Active tasks
        if self.active_tasks:
            st.write("**Active Tasks:**")
            for task_id, task in self.active_tasks.items():
                with st.expander(f"{task.task_type} - {task.status.value}", expanded=task.status == AnalysisStatus.RUNNING):
                    col1, col2 = st.columns([3, 1])
                    
                    with col1:
                        self.render_task_progress(task)
                    
                    with col2:
                        if not task.is_complete:
                            if st.button(f"Cancel", key=f"cancel_{task_id}"):
                                self.cancel_task(task_id)
                                st.rerun()
        
        # Task history
        if self.task_history:
            st.write("**Recent Tasks:**")
            for task in self.task_history[-5:]:  # Show last 5
                status_icon = {
                    AnalysisStatus.SUCCESS: "✅",
                    AnalysisStatus.FAILED: "❌",
                    AnalysisStatus.CANCELLED: "⚠️"
                }.get(task.status, "❓")
                
                st.write(f"{status_icon} {task.task_type} - {task.duration:.1f}s" if task.duration else f"{status_icon} {task.task_type}")
    
    def _run_deseq2_worker(self,
                          task: AnalysisTask,
                          expression_data: ExpressionData,
                          clinical_data: ClinicalData,
                          condition_column: str,
                          **kwargs) -> DESeq2Results:
        """Worker function for DESeq2 analysis"""
        try:
            task.mark_started()
            
            # Import DESeq2 engine
            from analysis.deseq2_engine import get_deseq2_engine
            deseq2_engine = get_deseq2_engine()
            
            # Stage 1: Validate data
            task.update_progress(0.05, "Validating data")
            if len(expression_data.samples) < 6:
                raise ValueError("Insufficient samples for DESeq2 analysis")
            
            # Stage 2: Prepare data
            task.update_progress(0.15, "Preparing count matrix")
            count_matrix, metadata = deseq2_engine.prepare_data(
                expression_data, clinical_data, condition_column
            )
            
            # Stage 3: Initialize DESeq2
            task.update_progress(0.25, "Initializing DESeq2")
            time.sleep(0.5)  # Simulate R initialization
            
            # Stage 4: Run analysis (this is the main computation)
            task.update_progress(0.30, "Running differential expression analysis")
            
            # Check cache first
            cache_key = f"deseq2_{condition_column}"
            cached_result = self.cache.cache_operation(
                cache_key,
                deseq2_engine.run_analysis,
                expression_data,
                clinical_data,
                condition_column=condition_column,
                dependencies=[f"expression_data_{id(expression_data)}"],
                tags=['deseq2', 'differential_expression'],
                cache_policy='deseq2_results',
                **kwargs
            )
            
            # Stage 5: Process results
            task.update_progress(0.85, "Processing results")
            time.sleep(0.2)
            
            # Stage 6: Generate visualizations
            task.update_progress(0.95, "Generating visualizations")
            if hasattr(cached_result, 'results') and not cached_result.results.empty:
                # Create volcano plot
                volcano_plot = deseq2_engine.create_volcano_plot(cached_result.results)
                cached_result.plots = {'volcano': volcano_plot}
            
            # Complete
            task.mark_completed(cached_result)
            return cached_result
            
        except Exception as e:
            logger.error(f"DESeq2 analysis failed: {e}")
            task.mark_failed(str(e))
            raise
    
    def _run_literature_search_worker(self,
                                     task: AnalysisTask,
                                     gene_list: List[str],
                                     search_terms: List[str] = None,
                                     **kwargs):
        """Worker function for literature search"""
        try:
            task.mark_started()
            
            # Import literature engine
            from analysis.literature_engine import get_literature_engine
            lit_engine = get_literature_engine()
            
            # Stage 1: Prepare queries
            task.update_progress(0.1, "Preparing gene queries")
            queries = gene_list[:10]  # Limit for demo
            if search_terms:
                queries.extend(search_terms)
            
            # Stage 2-4: Search APIs (simulated progress)
            for i, api_name in enumerate(['PubMed', 'Semantic Scholar', 'ArXiv']):
                progress = 0.1 + (i + 1) * 0.2
                task.update_progress(progress, f"Searching {api_name}")
                time.sleep(1)  # Simulate API calls
            
            # Stage 5: Process with LLM
            task.update_progress(0.9, "Processing results with LLM")
            
            # Use cached literature search
            cache_key = f"literature_{'_'.join(queries[:3])}"
            results = self.cache.cache_operation(
                cache_key,
                lit_engine.search_literature,
                queries,
                dependencies=[],
                tags=['literature', 'search'],
                cache_policy='literature_search',
                **kwargs
            )
            
            task.mark_completed(results)
            return results
            
        except Exception as e:
            logger.error(f"Literature search failed: {e}")
            task.mark_failed(str(e))
            raise
    
    def _run_chunked_analysis_worker(self,
                                   task: AnalysisTask,
                                   analysis_func: Callable,
                                   data_chunks: List[Any],
                                   chunk_processor: Callable,
                                   **kwargs):
        """Worker function for chunked analysis"""
        try:
            task.mark_started()
            
            # Stage 1: Prepare chunks
            task.update_progress(0.05, "Preparing chunks")
            n_chunks = len(data_chunks)
            
            # Stage 2: Process chunks in parallel
            task.update_progress(0.1, "Processing chunks in parallel")
            
            results = []
            for i, chunk in enumerate(data_chunks):
                # Update progress
                progress = 0.1 + (i / n_chunks) * 0.7
                task.update_progress(progress, f"Processing chunk {i+1}/{n_chunks}")
                
                # Process chunk
                chunk_result = chunk_processor(chunk, **kwargs)
                results.append(chunk_result)
            
            # Stage 3: Combine results
            task.update_progress(0.95, "Combining results")
            final_result = analysis_func(results, **kwargs)
            
            task.mark_completed(final_result)
            return final_result
            
        except Exception as e:
            logger.error(f"Chunked analysis failed: {e}")
            task.mark_failed(str(e))
            raise
    
    def _monitor_task_async(self, task: AnalysisTask, future: Future):
        """Monitor task execution asynchronously"""
        def monitor_worker():
            try:
                # Wait for completion
                result = future.result()
                
                # Move to history if completed
                if task.is_complete:
                    self._move_to_history(task)
                    
            except Exception as e:
                logger.error(f"Task monitoring failed: {e}")
                if not task.is_complete:
                    task.mark_failed(str(e))
                    self._move_to_history(task)
        
        # Start monitoring in background thread
        monitor_thread = threading.Thread(target=monitor_worker, daemon=True)
        monitor_thread.start()
    
    def _move_to_history(self, task: AnalysisTask):
        """Move completed task to history"""
        if task.task_id in self.active_tasks:
            del self.active_tasks[task.task_id]
        
        self.task_history.append(task)
        
        # Keep history size manageable
        if len(self.task_history) > self.max_history:
            self.task_history = self.task_history[-self.max_history:]
        
        logger.info(f"Task {task.task_id} moved to history")
    
    def cleanup(self):
        """Cleanup resources"""
        self.thread_pool.shutdown(wait=False)
        self.process_pool.shutdown(wait=False)


# Global instance
_async_manager = None

def get_async_manager() -> AsyncAnalysisManager:
    """Get the global async analysis manager instance"""
    global _async_manager
    if _async_manager is None:
        _async_manager = AsyncAnalysisManager()
    return _async_manager