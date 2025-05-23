# Step 3: Extract Differentially Expressed Genes (DEGs)

# Source the previous script to load the DESeqDataSet
source("prepare_deseq.R")

# Extract results from DESeq analysis
res <- results(dds)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj),]

# Extract significant genes (padj < 0.05)
sigGenes <- subset(resOrdered, padj < 0.05)

# Print summary of the results
cat("\nDifferential Expression Analysis Results:\n")
cat("Total number of genes analyzed:", nrow(res), "\n")
cat("Number of significant genes (padj < 0.05):", nrow(sigGenes), "\n")

# Create a summary of the results
summary(res)

# Save the results
write.csv(as.data.frame(resOrdered), "all_genes_results.csv")
write.csv(as.data.frame(sigGenes), "significant_genes.csv")

# Create a volcano plot
pdf("volcano_plot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Create a heatmap of the top 50 significant genes
# First, perform variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

# Create heatmap of top 50 significant genes
pdf("top50_genes_heatmap.pdf")
pheatmap(assay(vsd)[rownames(sigGenes)[1:50], ],
         scale="row",
         show_rownames=TRUE,
         show_colnames=TRUE,
         main="Top 50 Significant Genes")
dev.off()

cat("\nAnalysis complete. Results have been saved to:\n")
cat("- all_genes_results.csv\n")
cat("- significant_genes.csv\n")
cat("- volcano_plot.pdf\n")
cat("- top50_genes_heatmap.pdf\n") 