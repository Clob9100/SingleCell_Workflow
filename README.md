# **Single-Cell RNA-Seq Data Analysis Pipeline**

This repository contains an R-based pipeline for analyzing single-cell RNA-Seq (scRNA-Seq) data using the **Seurat** package. The dataset used in this analysis is the **Non-Small Cell Lung Cancer (NSCLC)** dataset from **10x Genomics**, specifically the data from **Donor 1**.

---

## **What is Single-Cell RNA-Seq?**

Single-cell RNA sequencing (scRNA-Seq) is a powerful technology that allows researchers to measure gene expression at the resolution of individual cells. Unlike bulk RNA-Seq, which averages gene expression across thousands of cells, scRNA-Seq provides insights into cellular heterogeneity, enabling the identification of distinct cell types, states, and rare populations within a tissue.

### **Key Applications**:
- **Cell Type Identification**: Classify cells into distinct types based on their gene expression profiles.
- **Differential Expression Analysis**: Identify genes that are differentially expressed between cell types or conditions.
- **Trajectory Inference**: Study developmental processes or transitions between cell states.
- **Spatial Transcriptomics**: Map gene expression to tissue architecture.

---

## **Dataset Overview**

The dataset used in this pipeline is the **20k Mixture of NSCLC DTCs from 7 Donors** provided by **10x Genomics**. Specifically, the data from **Donor 1** is analyzed here.

### **Dataset Details**:
- **Format**: `.h5` (HDF5) file.
- **Content**: Contains gene expression data (intronic + exonic reads) for **20,000 cells** from a NSCLC tumor.
- **File Name**: `20k_NSCLC_DTC_3p_nextgem_intron_donor_1_count_sample_feature_bc_matrix.h5`.
- **Data Structure**:
  - **Gene Expression Matrix**: Sparse matrix storing raw counts of RNA molecules for each gene in each cell.
  - **Barcodes**: Unique identifiers for each cell.
  - **Features**: Gene names and metadata.

---

## **Pipeline Overview**

The pipeline consists of the following steps:

### **1. Load the Dataset**
- The `.h5` file is loaded using `Read10X_h5()`.
- A **Seurat object** is created to store the raw gene expression data.

### **2. Quality Control (QC)**
- **Metrics**:
  - `nFeature_RNA`: Number of genes detected per cell.
  - `nCount_RNA`: Total RNA counts per cell.
  - `percent.mt`: Percentage of mitochondrial reads (indicator of cell quality).
- **Visualization**:
  - Violin plots and scatter plots are used to assess QC metrics.
- **Filtering**:
  - Cells with fewer than 200 genes, more than 2,500 genes, or more than 5% mitochondrial reads are removed.

### **3. Normalization**
- **Log Normalization**: Raw counts are normalized using a log-transformation (`LogNormalize`) with a scale factor of 10,000.
- **Purpose**: Corrects for differences in sequencing depth and stabilizes variance.

### **4. Identify Highly Variable Features**
- **Variable Genes**: The top 2,000 most variable genes are identified using the **variance-stabilizing transformation (VST)**.
- **Purpose**: Focuses on genes that drive biological heterogeneity.

### **5. Scaling**
- **Scale Data**: Gene expression values are scaled to have a mean of 0 and a variance of 1.
- **Purpose**: Removes technical noise and ensures genes are comparable.

### **6. Dimensionality Reduction (PCA)**
- **Principal Component Analysis (PCA)**: Reduces the dimensionality of the data to capture the most significant sources of variation.
- **Elbow Plot**: Helps determine the optimal number of principal components (PCs) for downstream analysis.

### **7. Clustering**
- **Find Neighbors**: Cells are grouped based on their similarity in PCA space.
- **Find Clusters**: Cells are clustered using a shared nearest neighbor (SNN) graph and the **Louvain algorithm**.
- **Resolution**: Adjusted to balance the granularity of clusters (e.g., 0.5).

### **8. Non-Linear Dimensionality Reduction (UMAP)**
- **UMAP**: Projects the high-dimensional data into 2D space for visualization.
- **Purpose**: Visualizes cell clusters and relationships between them.

### **9. Differential Expression Analysis**
- **Find Markers**: Identifies genes that are differentially expressed between clusters.
- **Heatmaps and Violin Plots**: Visualize marker gene expression across clusters.

### **10. Assign Cell Types**
- **Cluster Annotation**: Based on known marker genes, clusters are annotated with cell types (e.g., T cells, B cells, macrophages).

---

## **Key Outputs**

1. **QC Plots**:
   - Violin plots and scatter plots for QC metrics.
2. **Variable Feature Plot**:
   - Highlights the top 10 most variable genes.
3. **PCA and UMAP Plots**:
   - Visualizes cell clusters in reduced dimensions.
4. **Heatmaps**:
   - Shows expression of top marker genes across clusters.
5. **Cell Type Annotation**:
   - UMAP plot with clusters labeled by cell type.

---

## **How to Run the Pipeline**

1. **Install Required Packages**:
   ```R
   install.packages("Seurat")
   install.packages("tidyverse")
   install.packages("ggplot2")
   ```

2. **Set Working Directory**:
   ```R
   setwd("path/to/your/directory")
   ```

3. **Run the Script**:
   - Execute the provided R script (`pipeline.R`) step by step.

4. **Save Results**:
   - Plots and analysis results are saved to the working directory.

---

## **References**

- **Seurat**: https://satijalab.org/seurat/
- **Seurat Tutorial**: https://satijalab.org/seurat/articles/pbmc3k_tutorial
- **10x Genomics Dataset**: https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard
- **Single-Cell Best Practices**: https://www.sc-best-practices.org/


---

## **Contact**

For questions or feedback, please contact:
- **Name**: [Cloé BRENNA]
- **Email**: [clob9100@gmail.com]

---

This README provides a comprehensive overview of the pipeline and its steps. You can customize it further to include additional details or specific instructions for your project.