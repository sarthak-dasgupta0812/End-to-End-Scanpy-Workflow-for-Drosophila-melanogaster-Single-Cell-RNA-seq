# End-to-End-Scanpy-Workflow-for-Drosophila-melanogaster-Single-Cell-RNA-seq
This repository implements a Scanpy-based single-cell RNA-seq pipeline for Drosophila melanogaster, covering QC, normalization, HVG selection, PCA, neighborhood graph construction, Leiden clustering, and marker detection. The objective was of  a cell-type annotation using top three marker genes per cluster with end-to-end runtime under 20 seconds.
# Scanpy-Based Single-Cell RNA Sequencing Analysis of *Drosophila melanogaster*

## Overview

This repository presents a complete, reproducible, and computationally efficient single-cell RNA sequencing (scRNA-seq) analysis pipeline implemented using the Scanpy framework for *Drosophila melanogaster*. The workflow is designed to take raw 10X Genomics data as input and transform it into biologically interpretable cell populations through principled preprocessing, graph-based clustering, and marker gene identification. In addition to biological insights, the pipeline emphasizes algorithmic transparency, runtime efficiency, and low memory usage, making it suitable for benchmarking, teaching, and exploratory single-cell studies.

## Data Loading and AnnData Object Construction

The analysis begins with loading raw 10X Genomics output files, including the sparse gene–cell count matrix, barcodes, and gene feature metadata. These components are integrated into a single Scanpy AnnData object, which serves as the central data structure throughout the analysis. The AnnData object enables efficient storage of sparse matrices while simultaneously maintaining cell-level (obs) and gene-level (var) annotations required for downstream analysis.

At initialization, the dataset consisted of **3,000 cells (n_obs)** and **17,473 genes (n_vars)**. The raw count matrix represents unnormalized UMI counts and preserves the inherent sparsity and skewness of single-cell transcriptomic data. All subsequent filtering, normalization, dimensionality reduction, and clustering steps are applied to this same AnnData object to ensure reproducibility and traceability.

## Quality Control and Filtering

Quality control (QC) is a critical step in scRNA-seq analysis to remove technical artifacts and low-quality observations. Cell-level QC focused on gene complexity, measured as the number of detected genes per cell. Cells expressing fewer than 200 genes were removed, as these likely represent empty droplets, dying cells, or technical failures.

Gene-level QC involved removing genes expressed in fewer than three cells. Such genes contribute little to population-level structure and introduce unnecessary noise. The dataset uses FlyBase gene identifiers (FBGN IDs) rather than gene symbols; therefore, all QC steps and downstream analyses were performed using FBGN IDs to ensure consistency with *Drosophila* genomic annotations. Mitochondrial genes were identified using FlyBase annotations, while non-mitochondrial ribosomal genes such as *RPL32* were excluded from mitochondrial calculations.

After QC filtering, the dataset was reduced to 2,811 cells and 11,969 genes, preserving high-quality biological signal while removing uninformative data.

## Normalization and Log Transformation

To make gene expression values comparable across cells with different sequencing depths, total-count normalization was applied. Each cell was scaled to a target sum of **10,000 counts**, a commonly used normalization constant in scRNA-seq workflows. This step corrects for differences in sequencing depth and library size.

Following normalization, a logarithmic transformation was applied using Scanpy’s `log1p` function (log(1 + x)). Log transformation compresses large expression values, stabilizes variance, and reduces the dominance of highly expressed housekeeping and ribosomal genes. After log transformation, gene expression values ranged approximately between 0 and 5–10, indicating effective compression of the raw count distribution. The raw count matrix was preserved within the AnnData object for reference.

## Highly Variable Gene Selection

To focus downstream analyses on biologically meaningful variation, highly variable genes (HVGs) were identified. Scanpy models the relationship between mean expression and variance and selects genes whose variability exceeds expectations from technical noise alone. From the filtered gene set, the top **2,000 HVGs** were selected. Restricting analysis to HVGs improves computational efficiency and enhances biological signal for dimensionality reduction, clustering, and visualization.

## Scaling and Dimensionality Reduction

Before principal component analysis (PCA), HVGs were scaled so that each gene had a mean of zero and a variance of one. Scaling prevents PCA from being dominated by highly expressed genes. Values were clipped to a maximum absolute value of **10**, ensuring numerical stability and limiting the influence of extreme outliers.

PCA was then performed to reduce the high-dimensional gene expression space into a lower-dimensional representation while retaining maximal variance. Principal components are ordered by decreasing variance explained, with PC1 capturing the largest variance, followed by PC2, and so on. A variance ratio (elbow) plot showed a characteristic decreasing curve, and the top **30 principal components** were selected for downstream analyses.

## Neighborhood Graph Construction

Using the top 30 PCs, a k-nearest neighbor (kNN) graph was constructed with **15 neighbors per cell** based on Euclidean distance in PCA space. In this graph, cells are represented as nodes and edges connect each cell to its nearest neighbors. The resulting distance and connectivity matrices are sparse **2811 × 2811** matrices, capturing local neighborhood structure while remaining memory efficient. This graph forms the foundation for clustering and visualization.

## Leiden Clustering

Cell clustering was performed using the **Leiden algorithm**, a graph-based community detection method that identifies groups of cells with dense internal connectivity and sparse external connectivity. Clustering depends on parameters such as the number of PCs, number of neighbors, distance metric, and resolution. A resolution of **1.0** was used, yielding **15 clusters** across the 2,811 cells. These clusters represent transcriptionally distinct cell populations.

## Visualization with UMAP and t-SNE

To visualize clustering results, both **UMAP** and **t-SNE** embeddings were generated from the kNN graph. UMAP preserves both local and global structure, while t-SNE emphasizes local neighborhood relationships. Both methods revealed compact and well-separated clusters, with several clusters forming distinct islands and central clusters showing moderate overlap. The agreement between UMAP and t-SNE supports the robustness of the clustering structure.

## Marker Gene Identification and Cell-Type Annotation

Biological interpretation was achieved by identifying marker genes that are overexpressed in each cluster relative to all others. Differential expression analysis was performed using the **Wilcoxon rank-sum test**, a non-parametric method suitable for scRNA-seq data. For each of the 15 clusters, the **top three marker genes** were identified based on statistical significance and effect size. These markers predominantly included protein-coding genes and long non-coding RNAs (lncRNAs). Marker genes were annotated using FlyBase resources and mapped back onto low-dimensional embeddings for visualization and interpretation.

## Runtime and Performance

A core objective of this project was to evaluate and demonstrate the computational efficiency of Scanpy when applied to a real-world scRNA-seq dataset. To achieve fine-grained runtime profiling, each major pipeline step was explicitly timed using Python’s time module, allowing precise measurement of individual computational costs rather than relying solely on total runtime.

Step-wise Runtime Breakdown

The observed runtimes for each stage of the pipeline are summarized below:

Data loading and AnnData creation: 1.21 seconds
QC metrics computation: 0.05 seconds
QC filtering: 0.20 seconds
Library size normalization: 0.06 seconds
Log1p transformation: 0.06 seconds
Highly Variable Gene (HVG) selection: 0.10 seconds
HVG subsetting: ~0.00 seconds
Scaling: 0.07 seconds
Principal Component Analysis (PCA): 2.86 seconds
kNN neighbor graph construction: 0.09 seconds
Leiden clustering: 0.14 seconds
UMAP embedding: 6.12 seconds
Marker gene detection (Wilcoxon test): 8.55 seconds

The total end-to-end pipeline runtime, computed as the sum of all individual steps, was under 20 seconds, with the majority of computation concentrated in biologically intensive steps such as dimensionality reduction, visualization, and differential expression analysis.


##Why Scanpy Achieves High Runtime Efficiency
Scanpy’s efficiency in this project arises from several design principles. First, it operates primarily on sparse matrix representations, which is crucial for scRNA-seq data where the majority of entries are zeros. This drastically reduces memory overhead and avoids unnecessary computations on absent values.
Second, Scanpy tightly integrates optimized numerical libraries such as NumPy, SciPy, and scikit-learn, which delegate heavy linear algebra operations (e.g., PCA and nearest-neighbor search) to highly optimized, low-level C and Fortran routines. This is particularly evident in PCA, where nearly 3 seconds are spent computing principal components across 2,000 HVGs—an operation that would be prohibitively slow in pure Python.
Third, Scanpy enforces early dimensionality reduction through HVG selection. By restricting downstream steps (PCA, kNN, UMAP, clustering) to the top 2,000 informative genes instead of the full 11,969-gene matrix, the computational complexity is dramatically reduced while preserving biological signal. This design choice explains the negligible runtime of HVG subsetting and the fast execution of graph-based methods.
Finally, Scanpy’s graph-based algorithms (kNN, Leiden, UMAP) operate on cell–cell graphs rather than gene-level matrices, shifting computation to a lower-dimensional space (cells × 
cells). This is why neighbor graph construction and Leiden clustering complete in well under a second despite operating on thousands of cells.

Memory (RAM) Usage

Live memory usage was monitored using the psutil library during execution. The observed resident set size (RSS) was approximately 1509.89 MB. While this may appear high relative to dataset size, it reflects the cumulative memory footprint of Python, Scanpy, scientific libraries, sparse matrices, cached PCA results, and neighbor graphs loaded simultaneously in memory. Importantly, this RAM usage remained stable throughout execution and did not scale explosively, demonstrating Scanpy’s controlled memory behavior.


## Conclusion

This repository demonstrates a fast, reproducible, and biologically interpretable scRNA-seq analysis pipeline for *Drosophila melanogaster* using Scanpy. By integrating robust preprocessing, state-of-the-art algorithms, and performance benchmarking, the workflow provides a solid reference for single-cell analysis, benchmarking studies, and educational use.
