# Load required packages for the smoking and aging mouse brain analysis project

# Check if BiocManager is installed, if not install it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c("DESeq2", "airway", "biomaRt", "ggplot2", "pheatmap", 
                      "TabulaMurisData", "TabulaMurisSenisData", "SummarizedExperiment"))

# Load the required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(TabulaMurisData)
library(TabulaMurisSenisData)
library(SummarizedExperiment)

# Print confirmation message
cat("All required packages have been loaded successfully.\n") 