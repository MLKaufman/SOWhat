library(Seurat)

# Create a small sample dataset
# Use pbmc_small as a base or create one from scratch
# Let's create one from scratch to ensure Seurat v5 compatibility

set.seed(42)
counts <- matrix(rpois(1000 * 500, lambda = 1), nrow = 1000, ncol = 500)
rownames(counts) <- paste0("Gene", 1:1000)
colnames(counts) <- paste0("Cell", 1:500)

pbmc <- CreateSeuratObject(counts = counts, project = "PBMC_Sample")

# Normalize and scale
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 200)
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Find Neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Generate multiple clustering resolutions
resolutions <- c(0.1, 0.5, 1.0, 1.5)
for (res in resolutions) {
  pbmc <- FindClusters(pbmc, resolution = res)
}

# Save the object
saveRDS(pbmc, "sample_seurat.rds")
message("Sample Seurat object saved to sample_seurat.rds")
