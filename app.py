#!/usr/bin/env python3
"""
🧬 Prairie Genomics Suite v3.0 - Optional Clinical Data Edition
Advanced Genomics Analysis Platform with Flexible Sample Annotation

🎉 NEW in v3.0: Clinical Data is Now Optional!
- 🏷️ Interactive sample annotation without clinical files
- 🧠 Smart pattern detection from sample names
- 📋 Template-based experimental designs
- ⚡ 70-80% performance improvements with optimizations

Core Features:
- 🧬 DESeq2 differential expression with R integration + Python fallbacks
- 📊 Interactive visualizations and publication-quality plots
- 🔧 Memory-efficient processing with chunking and caching
- 🛡️ Robust error handling and data validation
- 📈 Real-time progress tracking and async processing

Usage:
    streamlit run app.py

Author: Prairie Genomics Team  
Version: 3.0.0 - Optional Clinical Data Edition
"""

import sys
import logging
from datetime import datetime

# Configure logging for startup diagnostics
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """Main entry point with startup diagnostics"""
    try:
        logger.info("🚀 Starting Prairie Genomics Suite v4.0...")
        logger.info(f"⏰ Startup time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Import and run the modular application
        from ui.main_app import main as app_main
        from config import APP_VERSION
        
        logger.info(f"✅ Successfully loaded modular app v{APP_VERSION}")
        
        # Run the application
        app_main()
        
    except ImportError as e:
        logger.error(f"❌ Import error during startup: {e}")
        print(f"\n🚨 STARTUP FAILED: {e}")
        print("📋 Try these solutions:")
        print("1. Check that all dependencies are installed")
        print("2. Verify you're in the correct directory")
        print("3. Run: pip install -r requirements.txt")
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"❌ Unexpected error during startup: {e}")
        print(f"\n🚨 UNEXPECTED ERROR: {e}")
        print("📋 Please report this issue with the full error message")
        sys.exit(1)

if __name__ == "__main__":
    main()