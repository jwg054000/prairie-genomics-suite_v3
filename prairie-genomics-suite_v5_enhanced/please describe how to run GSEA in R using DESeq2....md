Gene Set Enrichment Analysis (GSEA) is a powerful computational method used in RNA-seq analysis to interpret the biological meaning of differential gene expression results. Instead of focusing on individual differentially expressed genes (DEGs), GSEA assesses whether a predefined set of genes (e.g., genes belonging to a specific biological pathway or Gene Ontology term) is coordinately up- or down-regulated in a particular condition.2 This approach helps researchers gain deeper biological insights by identifying affected cellular processes and pathways.4

To run GSEA in R using data derived from DESeq2, you will primarily use the fgsea or clusterProfiler packages, which are robust Bioconductor tools for this purpose.

### **Prerequisites: DESeq2 Differential Expression Results**

GSEA typically requires a ranked list of genes, usually based on their log2 fold change and significance (p-value or adjusted p-value) from a differential expression analysis. DESeq2 is a widely used Bioconductor package that provides these essential outputs.8

Assuming you have already performed differential expression analysis with DESeq2 and obtained your results, typically stored in a DESeqResults object (e.g., deseq2Results), you will extract the necessary columns for GSEA.

**Basic DESeq2 Workflow to get deseq2Results (recap):**

1. **Load Libraries:** library(DESeq2) 14  
2. **Prepare Count Matrix and Metadata:** Ensure your raw count matrix (genes as rows, samples as columns) and metadata (colData) are correctly formatted, with matching sample names.8  
3. **Create DESeqDataSet Object:** dds \<- DESeqDataSetFromMatrix(countData \= count\_matrix, colData \= metadata, design \= \~ condition).8  
4. **Pre-filter Low Counts:** dds \<- dds (or similar threshold).8  
5. **Run DESeq2 Analysis:** dds \<- DESeq(dds).8  
6. **Extract Results:** deseq2Results \<- results(dds, contrast=c("condition\_column", "numerator\_level", "denominator\_level")).8

### **Method 1: Gene Set Enrichment Analysis with fgsea**

fgsea (Fast Gene Set Enrichment Analysis) is an R package known for its speed and accuracy in calculating GSEA p-values.22 It requires a preranked list of genes and a collection of gene sets.

1\. Prepare Gene-level Statistics (Ranked List) from DESeq2 Output:  
GSEA typically uses a ranked list of all genes, not just the differentially expressed ones. The ranking is often based on a statistic that reflects both fold change and significance, such as the signed p-value or a combination of log2 fold change and adjusted p-value.

R

\# Load DESeq2 results (assuming 'deseq2Results' is your DESeqResults object)  
\# deseq2Results \<- results(dds, contrast=c("condition", "treated", "control")) \# Example from previous step

\# Convert DESeqResults to a data frame  
deseq2ResDF \<- as.data.frame(deseq2Results) %\>%  
  rownames\_to\_column("GeneSymbol") \# Convert row names (Ensembl IDs or Gene Symbols) to a column

\# Handle NA values in padj or log2FoldChange (e.g., set to a small non-zero value or remove)  
\# For fgsea, it's common to remove genes with NA p-values or padj  
deseq2ResDF \<- na.omit(deseq2ResDF)

\# Create a ranked list of genes  
\# A common ranking metric is the log2FoldChange, or a signed p-value  
\# For fgsea, a named numeric vector is required, where names are gene symbols and values are ranks.  
\# Let's use log2FoldChange as the ranking statistic.  
gene\_ranks \<- deseq2ResDF$log2FoldChange  
names(gene\_ranks) \<- deseq2ResDF$GeneSymbol

\# Sort the ranks in decreasing order (most upregulated to most downregulated)  
gene\_ranks \<- sort(gene\_ranks, decreasing \= TRUE)

\# View the top few ranked genes  
head(gene\_ranks)  
tail(gene\_ranks)

2\. Obtain Gene Sets (Pathways):  
Gene sets are collections of genes that share a common biological function, pathway, or other characteristic. The Molecular Signatures Database (MSigDB) is a widely used source for these.24 You can access MSigDB gene sets in R using packages like  
msigdbr or r4msigdb.

R

\# Install and load msigdbr (if not already installed)  
if (\!requireNamespace("BiocManager", quietly \= TRUE))  
    install.packages("BiocManager")  
BiocManager::install("msigdbr")  
library(msigdbr)

\# Get human gene sets from MSigDB (e.g., KEGG pathways, Hallmark gene sets)  
\# For human: msigdbr\_species() to see available species  
\# For mouse: msigdbr\_species()  
\# Example: Human KEGG pathways  
m\_df \<- msigdbr(species \= "Homo sapiens", category \= "C2", subcategory \= "CP:KEGG") \# C2: curated gene sets, CP:KEGG: KEGG pathways

\# Convert to the list format required by fgsea  
\# The list should have pathway names as names, and vectors of gene symbols as elements  
kegg\_pathways \<- m\_df %\>%  
  dplyr::select(gs\_name, gene\_symbol) %\>%  
  group\_by(gs\_name) %\>%  
  summarise(genes \= list(gene\_symbol)) %\>%  
  deframe() \# Converts to a named list

\# View a few pathway names  
head(names(kegg\_pathways))

**3\. Run fgsea Analysis:**

R

\# Install and load fgsea (if not already installed)  
BiocManager::install("fgsea")  
library(fgsea)

\# Run fgsea  
\# pathways: your list of gene sets (e.g., kegg\_pathways)  
\# stats: your ranked gene list (e.g., gene\_ranks)  
\# minSize, maxSize: minimum and maximum size of gene sets to consider (important for filtering) \[22, 24\]  
\# nperm: number of permutations to perform for p-value calculation (higher \= more accurate, but slower) \[24\]  
fgsea\_results \<- fgsea(  
  pathways \= kegg\_pathways,  
  stats \= gene\_ranks,  
  minSize \= 15,  
  maxSize \= 500,  
  nperm \= 10000 \# Recommended for more accurate p-values \[24\]  
)

\# View the top enriched pathways, sorted by adjusted p-value (padj)  
fgsea\_results\_sorted \<- fgsea\_results %\>%  
  arrange(padj) %\>%  
  as.data.frame() \# Convert to data frame for easier viewing

print("Top fgsea results:")  
head(fgsea\_results\_sorted)

Interpreting fgsea Results:  
The fgsea\_results table contains several important columns 22:

* pathway: Name of the gene set/pathway.  
* pval: Nominal p-value for enrichment.  
* padj: Adjusted p-value (FDR) to account for multiple hypothesis testing. This is the most important metric for significance.  
* ES: Enrichment Score. Positive ES indicates enrichment in upregulated genes, negative ES indicates enrichment in downregulated genes.  
* NES: Normalized Enrichment Score. Allows comparison of enrichment results across different gene sets and experiments.  
* size: Number of genes in the gene set that are present in your ranked list.  
* leadingEdge: Genes that contribute most to the enrichment score.

**4\. Visualize fgsea Results:**

R

\# Plot a specific enriched pathway (e.g., the top one)  
\# This shows the distribution of genes from the pathway along the ranked list  
plotEnrichment(  
  pathway \= kegg\_pathways\[\[fgsea\_results\_sorted$pathway\[1\]\]\],  
  stats \= gene\_ranks  
) \+  
  labs(title \= fgsea\_results\_sorted$pathway\[1\]) \+  
  theme\_minimal()

\# Plot a table of top enriched pathways  
\# This requires the 'data.table' package for plotGseaTable  
library(data.table) \# \[22\]  
top\_pathways\_up \<- fgsea\_results %\>%  
  filter(NES \> 0) %\>%  
  arrange(padj) %\>%  
  head(10) %\>%  
  pull(pathway)

top\_pathways\_down \<- fgsea\_results %\>%  
  filter(NES \< 0) %\>%  
  arrange(padj) %\>%  
  head(10) %\>%  
  pull(pathway)

\# Combine and plot  
plotGseaTable(  
  pathways \= kegg\_pathways\[c(top\_pathways\_up, top\_pathways\_down)\],  
  stats \= gene\_ranks,  
  fgseaRes \= fgsea\_results,  
  gseaParam \= 0.5 \# Adjust for plot aesthetics  
)

### **Method 2: Gene Set Enrichment Analysis with clusterProfiler**

clusterProfiler is another comprehensive Bioconductor package for functional enrichment analysis, supporting Gene Ontology (GO) and KEGG pathway enrichment, among others.2 It can directly use the output from DESeq2.

1\. Prepare Input Genes and Background Genes from DESeq2 Output:  
clusterProfiler typically takes a vector of differentially expressed gene IDs (e.g., gene symbols) and a background set of all genes analyzed.

R

\# Load DESeq2 results (assuming 'deseq2Results' is your DESeqResults object)  
\# deseq2Results \<- results(dds, contrast=c("condition", "treated", "control"))

\# Convert DESeqResults to a data frame and add gene symbols  
deseq2ResDF\_annotated \<- as.data.frame(deseq2Results) %\>%  
  rownames\_to\_column("EnsemblID") %\>%  
  \# Assuming you have a gene symbol column, or perform conversion here  
  \# For this example, let's assume EnsemblID is already gene symbol for simplicity,  
  \# or you've already converted it as per previous instructions.  
  \# If not, you'd use mapIds or biomaRt here.  
  mutate(GeneSymbol \= EnsemblID) \# Placeholder: replace with actual gene symbol conversion

\# Filter for significant differentially expressed genes (DEGs)  
\# Adjust padj and log2FoldChange thresholds as needed  
significant\_genes \<- deseq2ResDF\_annotated %\>%  
  filter(padj \< 0.05 & abs(log2FoldChange) \> 1) %\>%  
  pull(GeneSymbol) \# Extract gene symbols of DEGs

\# Define background genes (all genes that were tested in DESeq2)  
background\_genes \<- deseq2ResDF\_annotated %\>%  
  pull(GeneSymbol)

\# View the number of significant genes and background genes  
length(significant\_genes)  
length(background\_genes)

**2\. Perform Gene Ontology (GO) Enrichment Analysis:**

R

\# Install and load clusterProfiler (if not already installed)  
BiocManager::install("clusterProfiler")  
library(clusterProfiler)  
library(org.Hs.eg.db) \# For human GO annotations (use org.Mm.eg.db for mouse)

\# Convert gene symbols to Entrez IDs (required by enrichGO/enrichKEGG)  
\# This step is crucial for clusterProfiler  
entrez\_ids \<- mapIds(  
  x \= org.Hs.eg.db, \# Use org.Mm.eg.db for mouse  
  keys \= significant\_genes,  
  column \= "ENTREZID",  
  keytype \= "SYMBOL",  
  multiVals \= "first"  
)  
entrez\_ids \<- na.omit(entrez\_ids) \# Remove genes that couldn't be mapped

\# Perform GO enrichment analysis  
\# ont: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "ALL"  
\# OrgDb: organism database (e.g., org.Hs.eg.db)  
\# keyType: type of input gene IDs  
\# readable: TRUE converts Entrez IDs back to gene symbols in results  
go\_results \<- enrichGO(  
  gene \= entrez\_ids,  
  OrgDb \= org.Hs.eg.db,  
  ont \= "BP",  
  keyType \= "ENTREZID",  
  readable \= TRUE  
)

\# View top GO terms  
print("Top GO enrichment results:")  
head(as.data.frame(go\_results))

**3\. Perform KEGG Pathway Enrichment Analysis:**

R

\# Install and load ReactomePA (often used with clusterProfiler for pathways)  
BiocManager::install("ReactomePA")  
library(ReactomePA)

\# Perform KEGG enrichment analysis  
\# organism: "hsa" for human, "mmu" for mouse  
kegg\_results \<- enrichKEGG(  
  gene \= entrez\_ids,  
  organism \= "hsa", \# "mmu" for mouse  
  keyType \= "ncbi-geneid" \# or "kegg"  
)

\# View top KEGG pathways  
print("Top KEGG enrichment results:")  
head(as.data.frame(kegg\_results))

4\. Visualize clusterProfiler Results:  
clusterProfiler offers a rich set of visualization functions.7

R

\# Bar plot of enriched GO terms  
barplot(go\_results, showCategory \= 10) \+  
  labs(title \= "Top 10 Enriched Biological Processes") \+  
  theme\_minimal()

\# Dot plot of enriched KEGG pathways  
dotplot(kegg\_results, showCategory \= 10) \+  
  labs(title \= "Top 10 Enriched KEGG Pathways") \+  
  theme\_minimal()

\# Pathway map visualization (for KEGG)  
\# This will open a browser window showing the pathway with your genes highlighted  
\# You need to specify a pathway ID from your results (e.g., 'hsa04110' for Cell cycle)  
\# browseKEGG(kegg\_results, 'hsa04110')

### **Best Practices and Tips for GSEA:**

* **Input Data Quality:** Ensure your DESeq2 results are robust, with proper normalization and batch effect correction if necessary.8 GSEA is sensitive to the quality of the input gene ranks.  
* **Gene Set Selection:** Choose appropriate gene sets relevant to your biological question (e.g., specific pathways, GO terms, or custom gene sets).2  
* **Adjusted P-values:** Always use adjusted p-values (FDR) to determine significance, as you are performing multiple tests.8  
* **Interpretation:** GSEA results should be interpreted in the context of your experimental design and biological knowledge. Visualizations like enrichment plots are crucial for understanding the enrichment patterns.22  
* **Computational Resources:** For large gene sets or many permutations, GSEA can be computationally intensive. Consider using parallelization options if available in the packages or running it on a more powerful machine.8  
* **Reproducibility:** Document all steps, including package versions and parameters used, to ensure your analysis is reproducible.