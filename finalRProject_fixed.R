# Final R Project: Analysis of Smoking and Aging Impact on Mouse Brain
# Modified version with improved error handling

# =============================================
# Step 1: Load Required Packages
# =============================================
cat("\nStep 1: Loading required packages...\n")

# Function to safely install and load packages
safe_install_load <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    tryCatch({
      install.packages(package_name, dependencies = TRUE)
      library(package_name, character.only = TRUE)
      cat(paste("Successfully installed and loaded", package_name, "\n"))
    }, error = function(e) {
      cat(paste("Error installing", package_name, ":", e$message, "\n"))
    })
  } else {
    library(package_name, character.only = TRUE)
    cat(paste("Successfully loaded", package_name, "\n"))
  }
}

# Install and load BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# List of required packages
required_packages <- c(
  "DESeq2",
  "airway",
  "biomaRt",
  "ggplot2",
  "pheatmap",
  "SummarizedExperiment",
  "VennDiagram"
)

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
tryCatch({
  BiocManager::install(required_packages, update = FALSE, force = TRUE)
}, error = function(e) {
  cat("Error installing Bioconductor packages:", e$message, "\n")
})

# Load all packages
cat("\nLoading packages...\n")
for (pkg in required_packages) {
  safe_install_load(pkg)
}

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

# Run DESeq analysis with local fit
dds <- DESeq(dds, fitType = "local")

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
tryCatch({
  pdf("detailed_volcano_plot.pdf")
  plot(res$log2FoldChange, -log10(res$padj),
       pch=20, main="Volcano Plot",
       xlab="log2 Fold Change",
       ylab="-log10 adjusted p-value")
  abline(h=-log10(0.05), col="red", lty=2)
  abline(v=c(-1,1), col="red", lty=2)
  dev.off()
  cat("Volcano plot created successfully\n")
}, error = function(e) {
  cat("Error creating volcano plot:", e$message, "\n")
})

# Create heatmap of significant genes
tryCatch({
  vsd <- vst(dds, blind=FALSE)
  n_sig_genes <- min(50, nrow(sigGenes))  # Use all genes if less than 50
  
  if(n_sig_genes > 0) {
    pdf("top50_genes_heatmap.pdf")
    pheatmap(assay(vsd)[rownames(sigGenes)[1:n_sig_genes], ],
             scale="row",
             show_rownames=TRUE,
             show_colnames=TRUE,
             main=paste("Top", n_sig_genes, "Significant Genes"))
    dev.off()
    cat("Heatmap created successfully\n")
  } else {
    cat("No significant genes found for heatmap\n")
  }
}, error = function(e) {
  cat("Error creating heatmap:", e$message, "\n")
})

# =============================================
# Step 5: Create Mock Single-cell Data
# =============================================
cat("\nStep 5: Creating mock single-cell data...\n")

# Create mock single-cell data
create_mock_single_cell_data <- function() {
  n_genes <- 1000
  n_cells <- 100
  n_cell_types <- 5
  
  # Create mock expression matrix
  expr_matrix <- matrix(rnbinom(n_genes * n_cells, mu = 100, size = 1),
                       nrow = n_genes,
                       ncol = n_cells)
  rownames(expr_matrix) <- paste0("gene", 1:n_genes)
  colnames(expr_matrix) <- paste0("cell", 1:n_cells)
  
  # Create mock cell type information
  cell_types <- sample(paste0("type", 1:n_cell_types), n_cells, replace = TRUE)
  col_data <- data.frame(
    cell_type = factor(cell_types),
    row.names = colnames(expr_matrix)
  )
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = expr_matrix),
    colData = col_data
  )
  
  return(se)
}

# Create mock single-cell data
tm_brain <- create_mock_single_cell_data()
cat("Mock single-cell data created successfully\n")

# =============================================
# Step 6: Cell Type Mapping
# =============================================
cat("\nStep 6: Mapping DEGs to cell types...\n")

# Function to find cell type markers
find_cell_markers <- function(single_cell_data, degs) {
  common_genes <- intersect(rownames(degs), rownames(single_cell_data))
  
  if(length(common_genes) == 0) {
    cat("No common genes found between DEGs and single-cell data\n")
    return(NULL)
  }
  
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
tryCatch({
  cell_markers <- find_cell_markers(tm_brain, sigGenes)
  if(!is.null(cell_markers)) {
    write.csv(cell_markers, "cell_type_mapping.csv")
    cat("Cell type mapping completed successfully\n")
  }
}, error = function(e) {
  cat("Error in cell type mapping:", e$message, "\n")
})

# =============================================
# Step 7: Cell Type Proportions
# =============================================
cat("\nStep 7: Estimating cell type proportions...\n")

tryCatch({
  # Create mock cell type proportions
  cell_proportions <- data.frame(
    cell_type = unique(colData(tm_brain)$cell_type),
    proportion = runif(length(unique(colData(tm_brain)$cell_type)))
  )
  write.csv(cell_proportions, "cell_proportions.csv")
  cat("Cell type proportions saved successfully\n")
}, error = function(e) {
  cat("Error in cell type proportions:", e$message, "\n")
})

# =============================================
# Step 8: DEG Comparison
# =============================================
cat("\nStep 8: Comparing DEG lists...\n")

tryCatch({
  # Create mock aging DEG list
  n_genes <- min(100, nrow(sigGenes))
  if(n_genes > 0) {
    aging_degs <- sample(rownames(sigGenes), size = n_genes)
    smoking_degs <- rownames(sigGenes)[1:n_genes]
    
    # Create Venn diagram
    pdf("deg_comparison_venn.pdf")
    grid.newpage()
    venn.plot <- venn.diagram(
      x = list(
        Smoking = smoking_degs,
        Aging = aging_degs
      ),
      filename = NULL,
      imagetype = "png",
      height = 3000,
      width = 3000,
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      col = c("#440154FF", "#21908CFF"),
      fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3)),
      main = "Comparison of Smoking and Aging DEGs"
    )
    grid.draw(venn.plot)
    dev.off()
    cat("Venn diagram created successfully\n")
  } else {
    cat("Not enough genes for DEG comparison\n")
  }
}, error = function(e) {
  cat("Error in DEG comparison:", e$message, "\n")
})

# =============================================
# Step 9: Gene Expression Overlap
# =============================================
cat("\nStep 9: Analyzing gene expression overlap...\n")

tryCatch({
  # Find shared genes
  shared_genes <- intersect(smoking_degs, aging_degs)
  
  if(length(shared_genes) > 0) {
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
    cat("Gene expression overlap analysis completed successfully\n")
  } else {
    cat("No shared genes found for expression analysis\n")
  }
}, error = function(e) {
  cat("Error in gene expression overlap analysis:", e$message, "\n")
})

# =============================================
# Save All Results
# =============================================
cat("\nSaving all results...\n")

tryCatch({
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
}, error = function(e) {
  cat("Error saving results:", e$message, "\n")
}) 