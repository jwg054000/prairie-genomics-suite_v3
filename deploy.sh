#!/bin/bash
# Prairie Genomics Suite v3 - Deployment Script

echo "🧬 Prairie Genomics Suite v3 - Deployment"
echo "==========================================="

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 not found. Please install Python 3.8+"
    exit 1
fi

# Check Python version
python_version=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "🐍 Python version: $python_version"

# Install dependencies
echo "📦 Installing dependencies..."
pip3 install -r requirements.txt

# Run startup check
echo "🔍 Running startup validation..."
python3 startup.py

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Deployment successful!"
    echo ""
    echo "🚀 To start the application:"
    echo "   Main app:  streamlit run app.py"
    echo "   Demo:      streamlit run sample_annotation_demo.py"
    echo ""
    echo "🌟 Features available:"
    echo "   ✅ Optional clinical data & smart sample annotation"
    echo "   ✅ Differential expression analysis (DESeq2)"
    echo "   ✅ Interactive visualizations"
    echo "   ✅ Performance optimizations"
    echo "   ✅ Robust error handling"
    echo ""
    echo "💡 Access the app at: http://localhost:8501"
else
    echo "❌ Deployment failed. Please check the errors above."
    exit 1
fi