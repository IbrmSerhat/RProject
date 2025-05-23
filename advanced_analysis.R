# Steps 4-9: Advanced Analysis and Integration

# Source the previous script to load DEGs
source("extract_degs.R")

# Step 4: Enhanced Visualization of Results
# Create a more detailed volcano plot
pdf("detailed_volcano_plot.pdf")
plot(res$log2FoldChange, -log10(res$padj),
     pch=20, main="Volcano Plot",
     xlab="log2 Fold Change",
     ylab="-log10 adjusted p-value")
abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1,1), col="red", lty=2)
dev.off()

# Step 5: Load single-cell data
cat("\nLoading single-cell data...\n")
# Load Tabula Muris data
tm_brain <- TabulaMurisData(tissue = "brain", celltype = "all")
# Load Tabula Muris Senis data
tms_brain <- TabulaMurisSenisData(tissue = "brain")

# Step 6: Map DEGs to cell types
cat("\nMapping DEGs to cell types...\n")
# Function to find cell type markers
find_cell_markers <- function(single_cell_data, degs) {
  # Get common genes between DEGs and single-cell data
  common_genes <- intersect(rownames(degs), rownames(single_cell_data))
  
  # Create a data frame to store results
  cell_markers <- data.frame(
    gene = common_genes,
    cell_type = character(length(common_genes)),
    expression_level = numeric(length(common_genes))
  )
  
  # For demonstration, we'll use a simplified approach
  # In practice, you would use more sophisticated methods
  for(i in seq_along(common_genes)) {
    gene <- common_genes[i]
    # Get expression levels for this gene
    expr <- assay(single_cell_data)[gene,]
    # Find cell type with highest expression
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

# Save cell type mapping results
write.csv(cell_markers, "cell_type_mapping.csv")

# Step 7: Estimate cell type proportions
cat("\nEstimating cell type proportions...\n")
# Note: In practice, this would be done using CIBERSORTx
# Here we'll create a mock result for demonstration
cell_proportions <- data.frame(
  cell_type = unique(colData(tm_brain)$cell_type),
  proportion = runif(length(unique(colData(tm_brain)$cell_type)))
)
write.csv(cell_proportions, "cell_proportions.csv")

# Step 8: Compare DEG lists (Venn diagram)
cat("\nCreating Venn diagram of DEGs...\n")
# For demonstration, we'll create a mock aging DEG list
# In practice, you would load your actual aging DEG list
aging_degs <- sample(rownames(sigGenes), size = 100)
smoking_degs <- rownames(sigGenes)[1:100]

# Create Venn diagram
library(VennDiagram)
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

# Step 9: Plot gene expression overlap
cat("\nCreating gene expression overlap plots...\n")
# Find shared genes
shared_genes <- intersect(smoking_degs, aging_degs)

# Create comparison plot for shared genes
pdf("shared_genes_expression.pdf")
# Create a mock expression matrix for demonstration
expr_matrix <- matrix(rnorm(length(shared_genes) * 4), 
                     nrow = length(shared_genes),
                     ncol = 4)
colnames(expr_matrix) <- c("Smoking_1", "Smoking_2", "Aging_1", "Aging_2")
rownames(expr_matrix) <- shared_genes

# Create heatmap
pheatmap(expr_matrix,
         scale = "row",
         main = "Expression of Shared Genes",
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()

cat("\nAdvanced analysis complete. Results have been saved to:\n")
cat("- detailed_volcano_plot.pdf\n")
cat("- cell_type_mapping.csv\n")
cat("- cell_proportions.csv\n")
cat("- deg_comparison_venn.pdf\n")
cat("- shared_genes_expression.pdf\n") 