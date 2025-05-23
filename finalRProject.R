# Final R Project: Analysis of Smoking and Aging Impact on Mouse Brain
# This script combines all previous analyses into one comprehensive workflow

# =============================================
# Step 1: Load Required Packages
# =============================================
cat("\nStep 1: Loading required packages...\n")

# Check if BiocManager is installed, if not install it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c("DESeq2", "airway", "biomaRt", "ggplot2", "pheatmap", 
                      "TabulaMurisData", "TabulaMurisSenisData", "SummarizedExperiment",
                      "VennDiagram"))

# Load the required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(TabulaMurisData)
library(TabulaMurisSenisData)
library(SummarizedExperiment)
library(VennDiagram)

# =============================================
# Step 2: Prepare Dataset and DESeq Analysis
# =============================================
cat("\nStep 2: Preparing dataset and running DESeq analysis...\n")

# Function to prepare smoking mouse dataset
prepare_smoking_dataset <- function() {
  set.seed(123)  # For reproducibility
  
  # Create mock count data
  n_genes <- 1000
  n_samples <- 6
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 1), 
                  nrow = n_genes, 
                  ncol = n_samples)
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)
  
  # Create mock sample information
  col_data <- data.frame(
    condition = factor(rep(c("control", "smoking"), each = 3)),
    row.names = colnames(counts)
  )
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
  
  return(se)
}

# Load the dataset
smokingMouse <- prepare_smoking_dataset()

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = assay(smokingMouse),
  colData = colData(smokingMouse),
  design = ~ condition
)

# Run DESeq analysis
dds <- DESeq(dds)

# =============================================
# Step 3: Extract Differentially Expressed Genes
# =============================================
cat("\nStep 3: Extracting differentially expressed genes...\n")

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

# =============================================
# Step 4: Enhanced Visualization
# =============================================
cat("\nStep 4: Creating enhanced visualizations...\n")

# Create detailed volcano plot
pdf("detailed_volcano_plot.pdf")
plot(res$log2FoldChange, -log10(res$padj),
     pch=20, main="Volcano Plot",
     xlab="log2 Fold Change",
     ylab="-log10 adjusted p-value")
abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1,1), col="red", lty=2)
dev.off()

# Create heatmap of top 50 significant genes
vsd <- vst(dds, blind=FALSE)
pdf("top50_genes_heatmap.pdf")
pheatmap(assay(vsd)[rownames(sigGenes)[1:50], ],
         scale="row",
         show_rownames=TRUE,
         show_colnames=TRUE,
         main="Top 50 Significant Genes")
dev.off()

# =============================================
# Step 5: Single-cell Data Analysis
# =============================================
cat("\nStep 5: Loading and analyzing single-cell data...\n")

# Load Tabula Muris data
tm_brain <- TabulaMurisData(tissue = "brain", celltype = "all")
# Load Tabula Muris Senis data
tms_brain <- TabulaMurisSenisData(tissue = "brain")

# =============================================
# Step 6: Cell Type Mapping
# =============================================
cat("\nStep 6: Mapping DEGs to cell types...\n")

# Function to find cell type markers
find_cell_markers <- function(single_cell_data, degs) {
  common_genes <- intersect(rownames(degs), rownames(single_cell_data))
  
  cell_markers <- data.frame(
    gene = common_genes,
    cell_type = character(length(common_genes)),
    expression_level = numeric(length(common_genes))
  )
  
  for(i in seq_along(common_genes)) {
    gene <- common_genes[i]
    expr <- assay(single_cell_data)[gene,]
    cell_markers$cell_type[i] <- names(which.max(tapply(expr, 
                                                       colData(single_cell_data)$cell_type, 
                                                       mean)))
    cell_markers$expression_level[i] <- max(tapply(expr, 
                                                  colData(single_cell_data)$cell_type, 
                                                  mean))
  }
  
  return(cell_markers)
}

# Map DEGs to cell types
cell_markers <- find_cell_markers(tm_brain, sigGenes)
write.csv(cell_markers, "cell_type_mapping.csv")

# =============================================
# Step 7: Cell Type Proportions
# =============================================
cat("\nStep 7: Estimating cell type proportions...\n")

# Create mock cell type proportions
cell_proportions <- data.frame(
  cell_type = unique(colData(tm_brain)$cell_type),
  proportion = runif(length(unique(colData(tm_brain)$cell_type)))
)
write.csv(cell_proportions, "cell_proportions.csv")

# =============================================
# Step 8: DEG Comparison
# =============================================
cat("\nStep 8: Comparing DEG lists...\n")

# Create mock aging DEG list
aging_degs <- sample(rownames(sigGenes), size = 100)
smoking_degs <- rownames(sigGenes)[1:100]

# Create Venn diagram
pdf("deg_comparison_venn.pdf")
venn.diagram(
  x = list(
    Smoking = smoking_degs,
    Aging = aging_degs
  ),
  filename = "deg_comparison_venn.pdf",
  imagetype = "pdf",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("#440154FF", "#21908CFF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3)),
  main = "Comparison of Smoking and Aging DEGs"
)
dev.off()

# =============================================
# Step 9: Gene Expression Overlap
# =============================================
cat("\nStep 9: Analyzing gene expression overlap...\n")

# Find shared genes
shared_genes <- intersect(smoking_degs, aging_degs)

# Create comparison plot for shared genes
pdf("shared_genes_expression.pdf")
expr_matrix <- matrix(rnorm(length(shared_genes) * 4), 
                     nrow = length(shared_genes),
                     ncol = 4)
colnames(expr_matrix) <- c("Smoking_1", "Smoking_2", "Aging_1", "Aging_2")
rownames(expr_matrix) <- shared_genes

pheatmap(expr_matrix,
         scale = "row",
         main = "Expression of Shared Genes",
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()

# =============================================
# Save All Results
# =============================================
cat("\nSaving all results...\n")

# Save DESeq results
write.csv(as.data.frame(resOrdered), "all_genes_results.csv")
write.csv(as.data.frame(sigGenes), "significant_genes.csv")

# Save DESeqDataSet
saveRDS(dds, "deseq_dataset.rds")

cat("\nAnalysis complete! All results have been saved to the following files:\n")
cat("- all_genes_results.csv\n")
cat("- significant_genes.csv\n")
cat("- deseq_dataset.rds\n")
cat("- detailed_volcano_plot.pdf\n")
cat("- top50_genes_heatmap.pdf\n")
cat("- cell_type_mapping.csv\n")
cat("- cell_proportions.csv\n")
cat("- deg_comparison_venn.pdf\n")
cat("- shared_genes_expression.pdf\n") 