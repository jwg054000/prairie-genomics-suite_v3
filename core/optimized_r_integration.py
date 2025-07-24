"""
Prairie Genomics Suite - Optimized R Integration

Process pooling and efficient data transfer for R-based genomics analyses.
Reduces R startup overhead and optimizes data serialization for large datasets.

Author: Prairie Genomics Team
"""

import pickle
import tempfile
import time
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List, Callable, Union
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, Future
import threading
import queue
import pandas as pd
import numpy as np

from config import get_config, is_feature_enabled
from core.enhanced_cache import get_enhanced_cache, cached_operation

logger = logging.getLogger(__name__)

@dataclass
class RProcessInfo:
    """Information about an R process in the pool"""
    process_id: int
    is_busy: bool = False
    last_used: float = 0.0
    initialization_time: float = 0.0
    total_jobs: int = 0

class OptimizedRIntegration:
    """
    Optimized R integration with process pooling and efficient data transfer
    """
    
    def __init__(self, pool_size: int = 2):
        self.pool_size = pool_size
        self.cache = get_enhanced_cache()
        self.config = get_config("performance")
        
        # R availability check
        self.r_available = False
        self.r_process_pool = None
        self.process_info: Dict[int, RProcessInfo] = {}
        
        # Data transfer optimization
        self.temp_dir = Path(tempfile.mkdtemp(prefix="prairie_r_"))
        self.temp_dir.mkdir(exist_ok=True)
        
        # Initialize if R integration is enabled
        if is_feature_enabled("enable_r_integration"):
            self._initialize_r_pool()
        
        logger.info(f"Optimized R integration initialized: pool_size={pool_size}, available={self.r_available}")
    
    def _initialize_r_pool(self):
        """Initialize R process pool with warm processes"""
        try:
            # Test R availability first
            import rpy2
            import rpy2.robjects as ro
            
            # Create process pool
            self.r_process_pool = ProcessPoolExecutor(
                max_workers=self.pool_size,
                initializer=self._initialize_r_worker
            )
            
            # Submit initialization tasks to warm up processes
            futures = []
            for i in range(self.pool_size):
                future = self.r_process_pool.submit(self._test_r_process, i)
                futures.append(future)
            
            # Check initialization results
            successful_inits = 0
            for i, future in enumerate(futures):
                try:
                    process_info = future.result(timeout=30)
                    self.process_info[i] = process_info
                    successful_inits += 1
                except Exception as e:
                    logger.warning(f"R process {i} initialization failed: {e}")
            
            if successful_inits > 0:
                self.r_available = True
                logger.info(f"R process pool initialized: {successful_inits}/{self.pool_size} processes ready")
            else:
                logger.error("No R processes could be initialized")
                
        except ImportError:
            logger.warning("rpy2 not available, R integration disabled")
        except Exception as e:
            logger.error(f"R process pool initialization failed: {e}")
    
    @staticmethod
    def _initialize_r_worker():
        """Initialize R worker process (runs once per process)"""
        try:
            import rpy2
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri, numpy2ri
            
            # Activate conversions
            pandas2ri.activate()
            numpy2ri.activate()
            
            # Load required R packages
            ro.r('''
            suppressPackageStartupMessages({
                library(DESeq2)
                library(dplyr)
                library(ggplot2)
            })
            ''')
            
            logger.debug("R worker process initialized successfully")
            return True
            
        except Exception as e:
            logger.error(f"R worker initialization failed: {e}")
            return False
    
    @staticmethod 
    def _test_r_process(process_id: int) -> RProcessInfo:
        """Test R process functionality"""
        start_time = time.time()
        
        try:
            # Import in worker process
            import rpy2.robjects as ro
            
            # Test basic R functionality
            result = ro.r('1 + 1')[0]
            if result != 2:
                raise ValueError("R basic test failed")
            
            initialization_time = time.time() - start_time
            
            return RProcessInfo(
                process_id=process_id,
                initialization_time=initialization_time
            )
            
        except Exception as e:
            logger.error(f"R process {process_id} test failed: {e}")
            raise
    
    def is_available(self) -> bool:
        """Check if R integration is available"""
        return self.r_available and self.r_process_pool is not None
    
    def get_status(self) -> Dict[str, Any]:
        """Get R integration status"""
        if not self.is_available():
            return {
                'available': False,
                'pool_size': 0,
                'active_processes': 0
            }
        
        active_processes = sum(1 for info in self.process_info.values() if info.is_busy)
        
        return {
            'available': True,
            'pool_size': self.pool_size,
            'active_processes': active_processes,
            'total_jobs': sum(info.total_jobs for info in self.process_info.values()),
            'process_info': self.process_info
        }
    
    @cached_operation(cache_policy='r_session')
    def run_deseq2_analysis(self,
                          count_matrix: pd.DataFrame,
                          metadata: pd.DataFrame,
                          design_formula: str,
                          contrasts: List[str] = None,
                          **kwargs) -> Dict[str, Any]:
        """
        Run DESeq2 analysis using optimized R integration
        
        Args:
            count_matrix: Gene expression count matrix
            metadata: Sample metadata
            design_formula: R design formula
            contrasts: List of contrasts to test
            **kwargs: Additional DESeq2 parameters
            
        Returns:
            Dictionary with DESeq2 results
        """
        if not self.is_available():
            raise RuntimeError("R integration not available")
        
        # Optimize data transfer
        temp_data = self._prepare_data_for_r(count_matrix, metadata)
        
        # Submit to R process pool
        future = self.r_process_pool.submit(
            self._run_deseq2_worker,
            temp_data,
            design_formula,
            contrasts,
            **kwargs
        )
        
        try:
            results = future.result(timeout=600)  # 10 minute timeout
            return results
        except Exception as e:
            logger.error(f"DESeq2 analysis failed: {e}")
            raise
        finally:
            # Cleanup temporary data
            self._cleanup_temp_data(temp_data['temp_files'])
    
    def _prepare_data_for_r(self, count_matrix: pd.DataFrame, metadata: pd.DataFrame) -> Dict[str, Any]:
        """Prepare data for efficient transfer to R process"""
        temp_files = {}
        
        try:
            # Save count matrix to temporary file
            count_file = self.temp_dir / f"counts_{int(time.time() * 1000)}.csv"
            count_matrix.to_csv(count_file, index=True)
            temp_files['counts'] = str(count_file)
            
            # Save metadata to temporary file  
            meta_file = self.temp_dir / f"metadata_{int(time.time() * 1000)}.csv"
            metadata.to_csv(meta_file, index=True)
            temp_files['metadata'] = str(meta_file)
            
            return {
                'temp_files': temp_files,
                'count_shape': count_matrix.shape,
                'metadata_shape': metadata.shape
            }
            
        except Exception as e:
            # Cleanup on error
            self._cleanup_temp_data(temp_files)
            raise RuntimeError(f"Data preparation failed: {e}")
    
    @staticmethod
    def _run_deseq2_worker(temp_data: Dict[str, Any],
                          design_formula: str,
                          contrasts: List[str] = None,
                          **kwargs) -> Dict[str, Any]:
        """Worker function to run DESeq2 in R process"""
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            
            # Load data from temporary files
            ro.r(f'''
            count_matrix <- read.csv("{temp_data['temp_files']['counts']}", row.names=1)
            metadata <- read.csv("{temp_data['temp_files']['metadata']}", row.names=1)
            ''')
            
            # Run DESeq2 analysis
            ro.r(f'''
            # Create DESeq2 dataset
            dds <- DESeqDataSetFromMatrix(
                countData = count_matrix,
                colData = metadata,
                design = {design_formula}
            )
            
            # Filter low count genes
            dds <- dds[rowSums(counts(dds)) >= 10, ]
            
            # Run DESeq2
            dds <- DESeq(dds, quiet = TRUE)
            
            # Get results
            results_list <- list()
            ''')
            
            # Get results for each contrast
            if contrasts:
                for i, contrast in enumerate(contrasts):
                    ro.r(f'''
                    res_{i} <- results(dds, contrast = {contrast})
                    res_{i} <- as.data.frame(res_{i})
                    results_list[["{contrast}"]] <- res_{i}
                    ''')
            else:
                # Default comparison
                ro.r('''
                res_default <- results(dds)
                res_default <- as.data.frame(res_default)
                results_list[["default"]] <- res_default
                ''')
            
            # Get normalized counts
            ro.r('''
            normalized_counts <- counts(dds, normalized = TRUE)
            normalized_counts <- as.data.frame(normalized_counts)
            ''')
            
            # Convert results back to Python
            results = {}
            
            # Get results for each contrast
            r_results_list = ro.r('results_list')
            for contrast_name in r_results_list.names:
                contrast_df = pandas2ri.rpy2py(r_results_list.rx2(contrast_name))
                contrast_df = contrast_df.reset_index()
                results[contrast_name] = contrast_df
            
            # Get normalized counts
            norm_counts = pandas2ri.rpy2py(ro.r('normalized_counts'))
            
            return {
                'results': results,
                'normalized_counts': norm_counts,
                'method': 'DESeq2 (R - Optimized)',
                'n_genes_analyzed': int(ro.r('nrow(dds)')[0]),
                'n_samples': int(ro.r('ncol(dds)')[0])
            }
            
        except Exception as e:
            logger.error(f"DESeq2 R worker failed: {e}")
            raise RuntimeError(f"DESeq2 analysis failed in R: {e}")
    
    def run_r_script(self, 
                    script: str,
                    data_objects: Dict[str, Any] = None,
                    return_objects: List[str] = None) -> Dict[str, Any]:
        """
        Run arbitrary R script with optimized data transfer
        
        Args:
            script: R script to execute
            data_objects: Python objects to transfer to R
            return_objects: R object names to return to Python
            
        Returns:
            Dictionary with returned R objects
        """
        if not self.is_available():
            raise RuntimeError("R integration not available")
        
        # Prepare data
        temp_data = {}
        if data_objects:
            temp_data = self._prepare_objects_for_r(data_objects)
        
        # Submit to R process pool
        future = self.r_process_pool.submit(
            self._run_r_script_worker,
            script,
            temp_data,
            return_objects or []
        )
        
        try:
            results = future.result(timeout=300)  # 5 minute timeout
            return results
        except Exception as e:
            logger.error(f"R script execution failed: {e}")
            raise
        finally:
            # Cleanup
            if temp_data and 'temp_files' in temp_data:
                self._cleanup_temp_data(temp_data['temp_files'])
    
    def _prepare_objects_for_r(self, objects: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare Python objects for transfer to R"""
        temp_files = {}
        object_info = {}
        
        for name, obj in objects.items():
            if isinstance(obj, pd.DataFrame):
                # Save DataFrame to CSV
                temp_file = self.temp_dir / f"{name}_{int(time.time() * 1000)}.csv"
                obj.to_csv(temp_file, index=True)
                temp_files[name] = str(temp_file)
                object_info[name] = {'type': 'dataframe', 'shape': obj.shape}
                
            elif isinstance(obj, np.ndarray):
                # Save array to CSV  
                temp_file = self.temp_dir / f"{name}_{int(time.time() * 1000)}.csv"
                pd.DataFrame(obj).to_csv(temp_file, index=False, header=False)
                temp_files[name] = str(temp_file)
                object_info[name] = {'type': 'array', 'shape': obj.shape}
                
            else:
                # For other objects, use pickle
                temp_file = self.temp_dir / f"{name}_{int(time.time() * 1000)}.pkl"
                with open(temp_file, 'wb') as f:
                    pickle.dump(obj, f)
                temp_files[name] = str(temp_file)
                object_info[name] = {'type': 'pickle'}
        
        return {
            'temp_files': temp_files,
            'object_info': object_info
        }
    
    @staticmethod
    def _run_r_script_worker(script: str,
                           temp_data: Dict[str, Any],
                           return_objects: List[str]) -> Dict[str, Any]:
        """Worker function to run R script"""
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            
            # Load data objects into R
            if temp_data and 'temp_files' in temp_data:
                for name, file_path in temp_data['temp_files'].items():
                    obj_info = temp_data['object_info'][name]
                    
                    if obj_info['type'] == 'dataframe':
                        ro.r(f'{name} <- read.csv("{file_path}", row.names=1)')
                    elif obj_info['type'] == 'array':
                        ro.r(f'{name} <- as.matrix(read.csv("{file_path}", header=FALSE))')
                    elif obj_info['type'] == 'pickle':
                        # Handle pickle objects (would need Python-R bridge)
                        logger.warning(f"Pickle object {name} not supported in R worker")
            
            # Execute R script
            ro.r(script)
            
            # Return requested objects
            results = {}
            for obj_name in return_objects:
                try:
                    r_obj = ro.r(obj_name)  
                    python_obj = pandas2ri.rpy2py(r_obj)
                    results[obj_name] = python_obj
                except Exception as e:
                    logger.warning(f"Could not return R object {obj_name}: {e}")
            
            return results
            
        except Exception as e:
            logger.error(f"R script worker failed: {e}")
            raise RuntimeError(f"R script execution failed: {e}")
    
    def _cleanup_temp_data(self, temp_files: Dict[str, str]):
        """Clean up temporary data files"""
        for name, file_path in temp_files.items():
            try:
                Path(file_path).unlink(missing_ok=True)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp file {file_path}: {e}")
    
    def shutdown(self):
        """Shutdown R process pool and cleanup"""
        if self.r_process_pool:
            self.r_process_pool.shutdown(wait=False)
        
        # Cleanup temp directory
        try:
            import shutil
            shutil.rmtree(self.temp_dir, ignore_errors=True)
        except Exception as e:
            logger.warning(f"Failed to cleanup temp directory: {e}")
        
        logger.info("Optimized R integration shutdown complete")


# Global instance
_optimized_r_integration = None

def get_optimized_r_integration() -> OptimizedRIntegration:
    """Get the global optimized R integration instance"""
    global _optimized_r_integration
    if _optimized_r_integration is None:
        _optimized_r_integration = OptimizedRIntegration()
    return _optimized_r_integration