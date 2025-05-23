# Step 2: Prepare DESeqDataSet for smoking mouse analysis

# Source the package loading script
source("load_packages.R")

# Function to load and prepare the smoking mouse dataset
prepare_smoking_dataset <- function() {
  # Note: This is a placeholder for the actual dataset loading
  # You'll need to replace this with your actual dataset loading code
  # Example structure of what the data should look like:
  
  # Create a mock SummarizedExperiment object for demonstration
  # In practice, you would load your actual data here
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

# Print summary of the dataset
cat("\nDataset Summary:\n")
cat("Number of genes:", nrow(dds), "\n")
cat("Number of samples:", ncol(dds), "\n")
cat("Design formula:", design(dds), "\n")

# Save the DESeqDataSet object
saveRDS(dds, "deseq_dataset.rds")
cat("\nDESeqDataSet has been prepared and saved as 'deseq_dataset.rds'\n") 