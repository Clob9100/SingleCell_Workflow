setwd("C:/Users/clob9/Desktop")

library(Seurat)
library(tidyverse)
library(ggplot2)

#####################
## 1. Load the NSCLC dataset
#####################
nsclc.sparse.m <- Read10X_h5(filename = "20k_NSCLC_DTC_3p_nextgem_intron_donor_1_count_sample_feature_bc_matrix.h5")
str(nsclc.sparse.m)
# The sparse matrix contains the gene expression data for the NSCLC dataset.
# It is stored in a sparse format to save memory and computational resources.
# The matrix has three components: data, i, and j, which represent the non-zero values, row indices, and column indices, respectively.


cts <- nsclc.sparse.m$`Gene Expression` # Select the Gene Expressoin matrix
# The Gene Expression matrix contains the raw counts of gene expression for each cell in the dataset.
# It is a sparse matrix that stores the counts of each gene in each cell, with rows representing genes and columns representing cells.
# The counts represent the number of RNA molecules detected for each gene in each cell.


# Inititate the Seurat object with the raw (non-normalized data)
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 21478 features across 1575 samples
# The Seurat object contains the raw counts of gene expression data for the NSCLC dataset.
# It has 21478 features (genes) and 1575 samples (cells).
# The object is initialized with the raw counts and basic metadata, such as the project name, minimum number of cells, and minimum number of features.
# The object can be further processed and analyzed using Seurat functions to explore the gene expression patterns and identify cell populations in the dataset.


#####################
## 2. QC (filter the low quality cells)
#####################
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)
# Calculate the percentage of mitochondrial (MT) reads for each cell in the dataset.
# This metric is commonly used as a quality control measure to identify cells with high levels of mitochondrial contamination, which may indicate low-quality cells or doublets.
# The percentage of MT reads is calculated by summing the counts of mitochondrial genes and dividing by the total counts for each cell.
# The results are stored in the meta.data slot of the Seurat object, allowing for easy access and visualization of the QC metrics.


VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
# The violin plot shows the distribution of the number of RNA features detected (nFeature_RNA), the total RNA counts (nCount_RNA), and the percentage of mitochondrial (MT) reads in the dataset.
# The scatter plot shows the relationship between the total RNA counts and the number of RNA features detected for each cell.
# Cells with high total RNA counts and a high number of RNA features are likely to be high-quality cells, while cells with low counts and features may be low-quality cells or empty droplets.
# The percentage of mitochondrial reads can be used as a quality control metric to identify cells with high levels of mitochondrial contamination, which may indicate low-quality cells or doublets.


#####################
## 3. Filtering (filter the low quality cells)
#####################
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
  percent.mt < 5)
# nFeature_RNA > 200 : retains cell with at least 200 genes detected (filters out low-quality cells (e.g., empty dropllets, broken cells)
# nFeatures_RNA < 2500 : retains cells with fewer than 2500 genes detected (filters out doublets)
# percent.mt < 5 : retains cells with less than 5% mitochondrial reads (filters out low-quality cells with high mitochondrial contamination)
# These filters help to remove low-quality cells, empty droplets, and doublets from the dataset, ensuring that only high-quality cells are included in downstream analysis.


#####################
## 4. Normalization
#####################
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Apply logarithmic normalization (log1p transformation) to the data.
# For each cell, raw counts are divided by the total counts for that cell, multiplied by a scale factor (10,000 by default), and log-transformed.
# Why 10000 ? This is a default convention in scRNA-seq analysis, but the scale factor can be adjusted based on the data and the desired range of expression values.
# It mimics "couts per million" (CPM) normalization used in bulk RNA-seq analysis.

str(nsclc.seurat.obj)
# With the @commands, we can see the different slots in the Seurat object

#####################
## 5. Identify highly variable features
#####################
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)
# nfeatures = 20000 : Select the 2000 most variable features based on the variance stabilizing transformation (VST)
# The VST is a variance-stabilizing transformation that stabilizes the variance of each gene, making it more consistent across the range of expression values.
# This is particularly useful for scRNA-seq data, where the variance of lowly expressed genes is typically higher than highly expressed genes.
# The FindVariableFeatures function identifies the most variable features based on the VST, which helps to capture the biological heterogeneity in the dataset.

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
# The dots in red are the top 10 most variable genes
# the dots in black are the other genes


#####################
## 6. Scaling (removing unwanted sources of variation)
#####################
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
# Scale the data by centering and scaling each gene to have a mean of 0 and a variance of 1.
# This removes unwanted sources of variation, such as differences in sequencing depth or cell size, and ensures that all genes are on a comparable scale.
# Scaling the data helps to improve the performance of downstream analyses, such as clustering and dimensionality reduction.


#####################
## 7. Perform Linear dimensional reduction (PCA)
#####################
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))
# Perform principal component analysis (PCA) on the scaled data using the 2000 most variable features identified earlier.
# PCA reduces the dimensionality of the data by identifying the principal components that capture the most variation in the dataset.
# The RunPCA function calculates the principal components and stores the results in the Seurat object for downstream analysis.

# Visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1:5, cells = 500, balanced = TRUE)
# Heatmap of PCA loadings for the top 5 features in the first 5 principal components
# This heatmap shows the contribution of each gene to the principal components, highlighting genes that drive the separation of cells in the PCA space.
# The cells are colored by the expression of the top 5 features in the first principal component.
# This helps to visualize how the genes are contributing to the separation of cells in the PCA space.


# Determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)
# Select the number of principal components to use for downstream analysis based on the "elbow" in the plot.
# Significant drop in the explained variance indicates the optimal number of principal components to retain.
# Helps to balance the amount of variance explained by the components with the complexity of the model.
# Up to PC15, the explained variance is still decreasing, so we can use PC15 for downstream analysis.


#####################
## 8. Clustering (identify cell clusters)
#####################

# Identify the neighbors of each cell based on the PCA space
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
# The resolution parameter controls the granularity of the clustering. Higher values lead to more clusters, while lower values lead to fewer clusters.
# The optimal resolution depends on the dataset and the biological question being addressed.
# A common approach is to explore a range of resolutions and evaluate the cluster stability and biological relevance of the resulting clusters.
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1.0 ))
View(nsclc.seurat.obj@meta.data)
# you can visualize how many cluster are formed with each resolution
# you can choose the resolution that best fits your data

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
# Visualize the clustering results in the PCA space, coloring the cells by the identified clusters.
# Clusters are labeled with different colors, and each dot represents a cell in the dataset.
# This plot helps to visualize the separation of cells into distinct clusters based on their gene expression profiles.
# Resolutions of 0.5 seems to be a good fit for this dataset to identify distinct cell populations.

# Setting identidy of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
Idents(nsclc.seurat.obj)


#####################
## 9. Non-Linear dimensial reduction (UMAP)
#####################
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)
# Perform uniform manifold approximation and projection (UMAP) on the PCA space to visualize the cell clusters in a two-dimensional plot.
# UMAP is a non-linear dimensionality reduction technique that preserves the local and global structure of the data, making it ideal for visualizing complex datasets.
# The DimPlot function visualizes the cell clusters in the UMAP space, coloring the cells by the identified clusters.
# This plot helps to visualize the separation of cells into distinct clusters and identify the relationships between different cell populations in the dataset.


#####################
## 10. Finding differentally expressed features (cluster biomarkers)
#####################

# Find all markers of cluster 2
cluster2.markers <- FindMarkers(nsclc.seurat.obj, ident.1 = 2)
head(cluster2.markers)
# This function identifies differentially expressed features (genes) for cluster 2 compared to all other cells in the dataset.
# The results include the log-fold change, p-value, and adjusted p-value for each feature, as well as the average expression in cluster 2 and all other cells.
# This information helps to identify genes that are specifically upregulated or downregulated in cluster 2 compared to other cells.

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(nsclc.seurat.obj, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers)
# This function identifies differentially expressed features (genes) for cluster 5 compared to clusters 0 and 3.
# The results include the log-fold change, p-value, and adjusted p-value for each feature, as well as the average expression in cluster 5 and clusters 0 and 3.
# This information helps to identify genes that are specifically upregulated or downregulated in cluster 5 compared to clusters 0 and 3.

# Find markers for every cluster compared to all remaining cells, report only the positive ones
cluster.markers <- FindAllMarkers(nsclc.seurat.obj, only.pos = TRUE)
cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.3)
# This function identifies differentially expressed features (genes) for each cluster compared to all other cells in the dataset.
# The results include the log-fold change, p-value, and adjusted p-value for each feature, as well as the average expression in each cluster and all other cells.
# This information helps to identify cluster-specific markers that can be used to characterize the different cell populations in the dataset.
# The filter(avg_log2FC > 1.3) line selects only markers with a log-fold change greater than 1.3, which indicates a significant difference in expression between the cluster and other cells.
# This threshold can be adjusted based on the dataset and the desired level of significance.
# The results can be further analyzed to understand the biological significance of the identified markers and their potential roles in different cell populations.

cluster0.markers <- FindMarkers(nsclc.seurat.obj, ident.1 = 0, logfc.threshold = 0.25,
                                test.use = "roc", only.pos = TRUE)
# This function identifies differentially expressed features (genes) for cluster 0 compared to all other cells in the dataset using a receiver operating characteristic (ROC) test.
# The results include the log-fold change, p-value, and adjusted p-value for each feature, as well as the area under the ROC curve (AUC) and the average expression in cluster 0 and all other cells.
# The ROC test is useful for identifying markers that are specifically upregulated in cluster 0 compared to other cells, with a focus on sensitivity and specificity.
# The logfc.threshold = 0.25 parameter sets a threshold for the log-fold change to select markers with a minimum effect size.

VlnPlot(nsclc.seurat.obj, features = c("PDCD1", "CD4", "CD3E", "BCL6", "CXCR3", "MKI67"), pt.size = 0.1)
# Visualize the expression of selected markers (PDCD1, CD4, CD3E, BCL6, CXCR3) across all cells in the dataset using violin plots.
# Violin plots show the distribution of expression levels for each marker, with the width of the plot indicating the density of cells at each expression level.
# This plot helps to visualize the expression patterns of key markers in different cell populations and identify potential biomarkers for further investigation.

VlnPlot(nsclc.seurat.obj, features = c("PDCD1", "CD4", "CD3E", "BCL6", "CXCR3", "MKI67"), slot = "counts", log = TRUE)
# Visualize the raw counts of selected markers (PDCD1, CD4, CD3E, BCL6, CXCR3) across all cells in the dataset using violin plots.
# This plot shows the distribution of raw counts for each marker, with the width of the plot indicating the density of cells at each count level.
# By visualizing the raw counts, you can assess the expression levels of key markers in different cell populations and identify potential biomarkers for further analysis.


FeaturePlot(nsclc.seurat.obj, features = c("PDCD1", "CD4", "CD3E", "BCL6", "CXCR3", "MKI67"))
# Visualize the expression of selected markers (PDCD1, CD4, CD3E, BCL6, CXCR3) across all cells in the dataset using feature plots.
# Feature plots show the expression levels of each marker on a per-cell basis, with each dot representing a cell and the color indicating the expression level.
# This plot helps to visualize the expression patterns of key markers in different cell populations and identify potential biomarkers for further investigation.

cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.3) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(nsclc.seurat.obj, features = top10$gene) + NoLegend()

# This code snippet identifies the top 5 differentially expressed genes for each cluster with an average log-fold change greater than 1.3.
# The results are stored in the top10 object, which contains the top 5 genes for each cluster that show significant differential expression.
# The DoHeatmap function generates a heatmap of the expression levels of these top genes across all cells in the dataset.
# The heatmap helps to visualize the expression patterns of these genes in different cell populations and identify cluster-specific markers that can be used to characterize the cell types.


#####################
## 11. Assigning cell types to clusters
#####################

new.cluster.ids <- c("T Cell", "B Cell", "Tfh", "Macrophage", "NK Cell", "Dendritic Cell", "Fibroblast", "Endothelial Cell", "Epithelial Cell")
names(new.cluster.ids) <- levels(nsclc.seurat.obj)
nscl <- RenameIdents(nsclc.seurat.obj, new.cluster.ids)
DimPlot(nscl, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# Assigning cell types to clusters based on marker gene expression profiles and known cell type markers.
# The new.cluster.ids vector contains the assigned cell types for each cluster, which are then renamed using the RenameIdents function.
# The DimPlot function visualizes the cell types in the UMAP space, with each cell colored by its assigned cell type.
# This plot helps to visualize the distribution of cell types in the dataset and identify distinct cell populations based on their gene expression profiles.


plot <- DimPlot(nscl, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "C:/Users/clob9/Desktop/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)