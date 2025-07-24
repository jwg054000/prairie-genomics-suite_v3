#!/usr/bin/env python3
"""
Prairie Genomics Suite v3 - Startup Script

Quick startup validation and dependency checking for deployment.
Run this before starting the Streamlit app to ensure everything is ready.
"""

import sys
import os
from pathlib import Path
import importlib
import logging

# Set up basic logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def check_python_version():
    """Check if Python version is compatible"""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        logger.error(f"Python 3.8+ required, found {version.major}.{version.minor}")
        return False
    
    logger.info(f"✅ Python {version.major}.{version.minor}.{version.micro}")
    return True

def check_core_dependencies():
    """Check if core dependencies are available"""
    core_deps = {
        'streamlit': 'Streamlit web framework',
        'pandas': 'Data manipulation',
        'numpy': 'Numerical computing', 
        'plotly': 'Interactive visualizations',
        'scipy': 'Scientific computing'
    }
    
    missing = []
    
    for dep, description in core_deps.items():
        try:
            importlib.import_module(dep)
            logger.info(f"✅ {dep} - {description}")
        except ImportError:
            logger.error(f"❌ {dep} - {description}")
            missing.append(dep)
    
    return len(missing) == 0, missing

def check_optional_dependencies():
    """Check optional dependencies and their impact"""
    optional_deps = {
        'rpy2': 'R integration for enhanced DESeq2 analysis',
        'h5py': 'HDF5 support for large datasets',
        'numba': 'Performance optimization',
        'tables': 'Advanced HDF5 features'
    }
    
    available = []
    missing = []
    
    for dep, description in optional_deps.items():
        try:
            importlib.import_module(dep)
            logger.info(f"✅ {dep} - {description}")
            available.append(dep)
        except ImportError:
            logger.warning(f"⚠️  {dep} - {description} (optional)")
            missing.append(dep)
    
    return available, missing

def check_directory_structure():
    """Verify required directories exist"""
    required_dirs = ['core', 'ui', 'analysis']
    required_files = ['config.py', 'app.py']
    
    missing_dirs = []
    missing_files = []
    
    for dir_name in required_dirs:
        if not Path(dir_name).exists():
            missing_dirs.append(dir_name)
        else:
            logger.info(f"✅ Directory: {dir_name}")
    
    for file_name in required_files:
        if not Path(file_name).exists():
            missing_files.append(file_name)
        else:
            logger.info(f"✅ File: {file_name}")
    
    return len(missing_dirs) == 0 and len(missing_files) == 0, missing_dirs, missing_files

def create_missing_directories():
    """Create any missing directories"""
    dirs_to_create = ['data', 'temp', '.streamlit']
    
    for dir_name in dirs_to_create:
        dir_path = Path(dir_name)
        if not dir_path.exists():
            dir_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"✅ Created directory: {dir_name}")

def main():
    """Main startup validation"""
    logger.info("🚀 Prairie Genomics Suite v3 - Startup Check")
    logger.info("=" * 50)
    
    success = True
    
    # Check Python version
    if not check_python_version():
        success = False
    
    # Check core dependencies
    core_ok, missing_core = check_core_dependencies()
    if not core_ok:
        logger.error(f"❌ Missing core dependencies: {missing_core}")
        logger.error("Run: pip install -r requirements.txt")
        success = False
    
    # Check optional dependencies
    available_optional, missing_optional = check_optional_dependencies()
    if missing_optional:
        logger.info(f"💡 Optional features available with: pip install {' '.join(missing_optional)}")
    
    # Check directory structure
    struct_ok, missing_dirs, missing_files = check_directory_structure()
    if not struct_ok:
        logger.error(f"❌ Missing directories: {missing_dirs}")
        logger.error(f"❌ Missing files: {missing_files}")
        success = False
    
    # Create missing directories
    create_missing_directories()
    
    # Summary
    logger.info("=" * 50)
    if success:
        logger.info("🎉 Startup check passed - Ready to launch!")
        logger.info("\nTo start the application:")
        logger.info("  Main app:    streamlit run app.py")
        logger.info("  Demo:        streamlit run sample_annotation_demo.py")
        
        # Feature availability summary
        logger.info("\n🌟 Available features:")
        logger.info("  ✅ Core differential expression analysis")
        logger.info("  ✅ Optional clinical data & sample annotation")
        logger.info("  ✅ Interactive visualizations")
        logger.info("  ✅ Performance optimizations")
        
        if 'rpy2' in available_optional:
            logger.info("  ✅ Enhanced R/DESeq2 integration")
        else:
            logger.info("  ⚠️  Python-only analysis (R integration disabled)")
            
        if 'h5py' in available_optional:
            logger.info("  ✅ Large dataset support (HDF5)")
        else:
            logger.info("  ⚠️  Standard dataset support only")
            
        if 'numba' in available_optional:
            logger.info("  ✅ Performance acceleration")
        else:
            logger.info("  ⚠️  Standard performance")
        
    else:
        logger.error("❌ Startup check failed - Please fix issues above")
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)