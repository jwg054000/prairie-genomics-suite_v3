# MC9 Demo Dataset

This dataset is designed to demonstrate the Prairie Genomics Suite capabilities.
It's based on the validated MC9 vs MLM comparison that achieved 100% expert agreement.

## Files:
- `MC9_demo_expression.csv`: Gene expression count matrix (1000 genes × 12 samples)
- `MC9_demo_metadata.csv`: Sample metadata (OPTIONAL - the app can detect groups automatically!)
- `gene_annotations.csv`: Gene symbol mappings

## Sample Groups:
- MC9: Mouse cancer cell line (3 replicates)
- MLM: Mouse lymphoma cell line (3 replicates)  
- M1245: Mouse cell line (3 replicates)
- M242: Mouse cell line (3 replicates)

## Expected Results:
- ~100-150 differentially expressed genes for MC9 vs MLM
- Strong upregulation of cancer genes (Il6, Myc, Cish) in MC9
- Clear separation in PCA plot
- Significant pathway enrichment for cancer-related processes

## Usage:
1. Upload `MC9_demo_expression.csv` to the app
2. Let the app auto-detect sample groups (or upload metadata)
3. Run DESeq2 analysis with default parameters
4. Explore pathway analysis results

This demonstrates the same workflow that achieved 100% expert validation!
