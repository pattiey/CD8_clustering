suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(Matrix))
suppressMessages(library(MAST))
suppressMessages(library(argparse))
suppressMessages(library(here))

parser <- ArgumentParser()

parser$add_argument("-d", "--data_dir", type = "character", default = here("data"),
                    help = "directory where data is stored")
parser$add_argument("-o", "--output_dir", type = "character", default = here("output"),
                    help = "destination directory for output")
parser$add_argument("-s", "--script_dir", type = "character", default = here("scripts"),
                    help = "directory where scripts are stored")
parser$add_argument("-p", "--project", type = "character", default = "test",
                    help = "name of project")
parser$add_argument("--reduction", type = "character", default = "PCA",
                    help = "method of dimensionality reduction for clustering, PCA, TSNE, or UMAP")
parser$add_argument("-k", "--k_param", type = "integer", default = 30,
                    help = "Number of neighbours to use for K nearest neighbours clustering")
parser$add_argument("-n", "--num_features", type = "integer", default = 1000,
                    help = "Number of features to use for SCTransform")
parser$add_argument("--resolution", type = "double", default = 0.5,
                    help = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.")
parser$add_argument("--dims", type = "integer", default = 10,
                    help = "Number of dimensions to use for clustering.")


args <- parser$parse_args()
print(args)

data_dir <- args$data_dir
output_dir <- args$output_dir
script_dir <- args$script_dir
project_name <- args$project
reduction <- args$reduction
k <- args$k_param
n <- args$num_features
resolution <- args$resolution
dims <- args$dims

data <- Read10X(data_dir)
data.seurat <- CreateSeuratObject(counts = data, project = project_name)

# Normalize count data with SCTransform
data.seurat <- SCTransform(data.seurat, variable.features.n = n, verbose = FALSE)

# Remove cell cycle, ribosomal, and mitochondrial genes from genes used for dimensionality reduction
# list of cell cycle genes from B16 Melanoma CD8 paper
cellCycle.genes <- read.csv(file.path(script_dir, "cellCycleGenes.csv"), as.is = TRUE)$x
mito.genes <- grep(pattern = "^mt-", x = rownames(data.seurat@assays$RNA@data), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[ls]", x = rownames(data.seurat@assays$RNA@data), value = TRUE)

SCT.features <- data.seurat@assays$SCT@var.features

print(paste("Number of features including cell cycle, ribosomal, and mitochondrial genes",
            length(SCT.features)))

data.seurat@assays$SCT@var.features <- SCT.features[!SCT.features %in% c(mito.genes, ribo.genes, cellCycle.genes)]

print(paste("Number of features excluding cell cycle, ribosomal, and mitochondrial genes",
            length(data.seurat@assays$SCT@var.features)))

# Run dimensionality reductions
set.seed(123)
DefaultAssay(object = data.seurat) <- "SCT"
data.seurat <- RunPCA(data.seurat, ndims.print = 1:5, nfeatures.print = 10)
data.seurat <- RunUMAP(data.seurat, reduction = "pca", dims = 1:10, seed.use = 123)
data.seurat <- RunTSNE(data.seurat, reduction = "pca", dims = 1:10, seed.use = 123, perplexity = 30)


# USE KNN to find neighbours/clusters based on selected dimensionality reduction method
data.seurat <- FindNeighbors(data.seurat, reduction = reduction, dims = 1:dims, k.param = k)
data.seurat <- FindClusters(data.seurat, resolution = resolution)

# Plot cell clusters
pca.plot <- DimPlot(data.seurat, reduction = "pca", group.by = "seurat_clusters") +
            ggtitle("TSNE plot of data coloured by cluster") +
            theme(aspect.ratio = 1, legend.position = "right")
ggsave(file.path(output_dir, 'PCA_plot.png'), pca.plot)
tsne.plot <- DimPlot(data.seurat, reduction = "tsne", group.by = "seurat_clusters") +
            ggtitle("TSNE plot of data coloured by cluster") +
            theme(aspect.ratio = 1, legend.position = "right")
ggsave(file.path(output_dir, 'TSNE_plot.png'), tsne.plot)
umap.plot <- DimPlot(data.seurat, reduction = "umap", group.by = "seurat_clusters") +
            ggtitle("UMAP plot of data coloured by cluster") +
            theme(aspect.ratio = 1, legend.position = "right")
ggsave(file.path(output_dir, 'UMAP_plot.png'), umap.plot)

# Save clusters to a csv identified by barcode
cluster.data <- data.frame(cell_id = colnames(data.seurat), cluster = data.seurat@meta.data$seurat_clusters)
write.csv(cluster.data, file = file.path(output_dir, "cellClusters.csv"), quote = FALSE, row.names = FALSE)
