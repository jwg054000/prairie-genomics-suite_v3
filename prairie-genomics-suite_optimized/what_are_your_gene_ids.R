# IDENTIFY YOUR GENE ID FORMAT
# Quick script to help determine what type of gene IDs you have

cat("🔍 GENE ID FORMAT IDENTIFIER\n")
cat("=" , rep("=", 60), "\n")

cat("📋 Please provide some example gene IDs from your data:\n\n")

# Common gene ID examples to help user identify their format
examples <- list(
  "Human Ensembl" = c("ENSG00000139618", "ENSG00000012048", "ENSG00000111640"),
  "Mouse Ensembl" = c("ENSMUSG00000051951", "ENSMUSG00000020122", "ENSMUSG00000057666"), 
  "Human Symbols" = c("BRCA1", "BRCA2", "TP53", "GAPDH"),
  "Mouse Symbols" = c("Brca1", "Brca2", "Trp53", "Gapdh"),
  "Entrez IDs" = c("672", "675", "7157", "2597"),
  "RefSeq" = c("NM_007294", "NM_000059", "NM_000546")
)

cat("🧬 COMMON GENE ID FORMATS:\n")
for (format_name in names(examples)) {
  ids <- examples[[format_name]]
  cat(sprintf("%-15s: %s\n", format_name, paste(ids, collapse = ", ")))
}

cat("\n❓ QUESTIONS TO ASK YOURSELF:\n")
cat("1. Do your gene IDs start with 'ENSG' followed by numbers? → Human Ensembl\n")
cat("2. Do your gene IDs start with 'ENSMUSG' followed by numbers? → Mouse Ensembl\n") 
cat("3. Are your gene IDs all uppercase letters? → Human gene symbols\n")
cat("4. Are your gene IDs mixed case (first letter uppercase)? → Mouse gene symbols\n")
cat("5. Are your gene IDs just numbers? → Entrez IDs\n")
cat("6. Do your gene IDs start with 'NM_' or 'NR_'? → RefSeq IDs\n")

cat("\n🎯 QUICK FIX BASED ON YOUR GENE TYPE:\n")
cat("• If HUMAN ENSEMBL (ENSG...): Set species to 'Human' in UI\n")
cat("• If MOUSE ENSEMBL (ENSMUSG...): Set species to 'Mouse' in UI\n") 
cat("• If already gene symbols: Disable 'Convert gene IDs to symbols' checkbox\n")
cat("• If other format: May need custom conversion approach\n")

cat("\n📁 EASY TEST:\n")
cat("1. Look at the first column of your data file\n")
cat("2. Copy the first few gene IDs\n")
cat("3. Compare them to the examples above\n")
cat("4. Set the species correctly in the Prairie Genomics Suite UI\n")

cat("\n💡 TROUBLESHOOTING:\n")
cat("If conversion still fails after setting correct species:\n")
cat("• Check your internet connection (needed for biomaRt backup)\n")
cat("• Verify gene IDs are current/valid Ensembl IDs\n")
cat("• Some gene IDs may not have symbol mappings\n")

cat("\n✅ Identification guide complete\n")