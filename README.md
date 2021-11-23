# CD8_clustering

# Background and Rationale

Single cell RNA-seq (scRNA-seq) experiments are increasing in prevalence within biological experimentation. An advantage of scRNA-seq over bulk sequencing experiments is the ability to look at gene expression of individual cells. This level of granularity allows for analysis that would otherwise not be possible with bulk sequencing experiments.

A classic example of single cell specific analysis is cell clustering. In a bulk sequencing experiment, the gene expression data is a mixture of various cell types that cannot be separated. In scRNA-seq experiments, individual cells within a sample can be separated into specific sub-type clusters based on their gene expression profiles using classic clustering methods.

In this CD8_clustering workflow, we look specifically at CD8+ T cells, also known as "killer T cells", which have a cytotoxic function within adaptive immunity. CD8+ T cells are typically categorized into specific subtypes - naive, memory (stem cell, central, and effector), effector, and exhausted. The subtypes of CD8+ T cell have different gene expression profiles as well as different functions in the immune system.

By performing clustering on scRNA-seq data of CD8+ T cells, one can infer the subtypes of individual cells based on the characteristic gene expressions of each cluster and gain valuable information on the cellular composition of given sample(s).

# Clustering Method

This workflow makes use of the [K Nearest Neighbours clustering algorithm](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) (KNN) where the number of nearest neighbours to use can be user specified. Since scRNA-seq data is highly dimensional and sparse, dimensionality reduction is necessary in order to performing clustering. The most common dimensionality reduction technique is [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis), but other methods are also available in this workflow ([t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) and [UMAP](https://arxiv.org/abs/1802.03426)).

# Workflow

![alt text](https://github.com/pattiey/CD8_clustering/blob/main/dag.svg)

The main steps of this workflow are:

1. FTP download of raw gene expression matrix files for all relevant samples from the NCBI GEO (tar format).
2. Extract individual sample specific data files.
3. Aggregation of individual sample data into one singular gene expression matrix.
4. Filter features and cells based on user determined thresholds.
5. Clustering of cells using KNN on dimension reduced gene expression data.

# Usage

First, clone the git repository to your local machine and enter the project directory:

```
git clone git@github.com:pattiey/CD8_clustering.git
cd CD8_clustering
```

The create an environment with the required packages:

```
conda create --name cd8_clustering --file env.txt
```

Activate the environment:

```
conda activate cd8_clustering
```

Run the Snakemake workflow with relevant parameters. Be sure to set the appropriate number of cores. Ensure that the `config.yaml` file is updated with the relevant fields before running.

```
snakemake --snakefile /path/to/CD8_clustering/Snakemake/Snakefile --configfile /path/to/CD8_clustering/Snakemake/config.yaml --cores 4
```

# Input

The input for this workflow is controlled by the `config.yaml` file. Here is a description of the fields of the `config.yaml` file.

| Field | Description |
| ----- | ----------- |
| `FTP_URL` | FTP download file from NCBI GEO of raw gene expression data |
| `DATA_DIR` | Directory where data is to be stored |
| `SCRIPTS_DIR` | Directory where project scripts are stored |
| `OUTPUT_DIR` | Directory where output files are to be stored |
| `PROJECT` | Name of experiment/project |
| `SAMPLES` | Sample names from GEO |

Other user specified parameters are also available to adjust in the Snakefile.

| Snakemake rule | Parameter | Description |
| -------------- | --------- | ----------- |
| `filter_cells` | `mito` | maximum percentage threshold of mitochondrial expression to filter |
| `filter_cells` | `ribo` | maximum percentage threshold of ribosomal expression to filter |
| `filter_cells` | `nFeature_lo` | minimum number of features present in cells to keep |
| `filter_cells` |`nFeature_hi` | maximum number of features present in cells to keep |
| `filter_cells` | `nCount_lo` | minimum number of counts required to keep a cell |
| `filter_cells` | `nCount_hi` | maximum number of counts required to keep a cell|
| `cluster_cells` | `reduction` | method of dimensionality reduction for clustering, pca, tsne, or umap |
| `cluster_cells` | `k` | Number of neighbours to use for K nearest neighbours clustering |
| `cluster_cells` | `num_features` | Number of features to use for SCTransform |
| `cluster_cells` | `resolution` | Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. |
| `cluster_cells` | `dims` | Number of dimensions to use for clustering. |


The example config file in this repository contains the samples and URL data from the NCBI Gene Expression Omnibus (GEO) [GSE116390](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116390). This is single cell RNA-seq data of CD8+ T cells from B16 melanoma tumours from mice pertaining to the experiments done by [S. Carmona et al.](https://doi.org/10.1080/2162402X.2020.1737369).

# Output

The workflow produces the result of KNN clustering on the scRNA-seq samples.

| Output File | Description |
| ----------- | ----------- |
| `cellClusters.csv` | A CSV file containing cluster labels of each cell identified through barcode and sample |
| `PCA_plot.png` | A plot of the first two principal components of the gene expression data, coloured by cluster |
| `TSNE_plot.png` | A t-SNE plot of the gene expression data, coloured by cluster |
| `UMAP_plot.png` | A UMAP plot of the gene expression data, coloured by cluster |

Using the sample data and the parameters specified in the Snakefile, here are the plot produced by the workflow.

### Principal Component Plot

![alt text](https://github.com/pattiey/CD8_clustering/blob/main/sample_output/PCA_plot.png)

### t-SNE Plot

![alt text](https://github.com/pattiey/CD8_clustering/blob/main/sample_output/TSNE_plot.png)

### UMAP plot

![alt text](https://github.com/pattiey/CD8_clustering/blob/main/sample_output/UMAP_plot.png)
