suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(Matrix))
suppressMessages(library(MAST))
suppressMessages(library(argparse))
suppressMessages(library(here))

parser <- ArgumentParser()

parser$add_argument("-o", "--output_dir", type = "character", default = here("filtered"),
                    help = "destination directory for output")
parser$add_argument("-d", "--data_dir", type = "character", default = here("data"),
                    help = "directory where data is stored")
parser$add_argument("-p", "--project", type = "character", default = "test",
                    help = "name of project")
parser$add_argument("--mito", type = "double", default = 1.0,
                    help = "maximum percentage threshold of mitochondrial expression to filter")
parser$add_argument("--ribo", type = "double", default = 1.0,
                    help = "maximum percentage threshold of ribosomal expression to filter")
parser$add_argument("--nFeature_lo", type = "integer", default = 0,
                    help = "minimum number of features present in cells to keep")
parser$add_argument("--nFeature_hi", type = "integer", default = 10^6,
                    help = "maximum number of features present in cells to keep")
parser$add_argument("--nCount_lo", type = "integer", default = 0,
                    help = "minimum number of counts required to keep a cell")
parser$add_argument("--nCount_hi", type = "integer", default = 10^6,
                    help = "maximum number of counts required to keep a cell")

args <- parser$parse_args()

data_dir <- args$data_dir
output_dir <- args$output_dir
project_name <- args$project
mito <- args$mito
ribo <- args$ribo
nFeature_lo <- args$nFeature_lo
nFeature_hi <- args$nFeature_hi
nCount_lo <- args$nCount_lo
nCount_hi <- args$nCount_hi

data <- Read10X(data_dir)
data.seurat <- CreateSeuratObject(counts = data,
                                  project = project_name)

# remove non-expressed genes
print(paste("Number of features", nrow(data.seurat@assays$RNA)))
print(paste("Number of cells", ncol(data.seurat@assays$RNA)))
counts <- GetAssayData(data.seurat, assay = "RNA")
counts <- counts[-(which(rowSums(counts) == 0))]
data.seurat <- subset(data.seurat, features = rownames(counts))
print(paste("Number of feaures after filtering for non-expressed", nrow(data.seurat@assays$RNA)))

# Filter cells for mitochondrial percentage, and ribosomal percentage
data.seurat[["percent.mito"]] <- PercentageFeatureSet(data.seurat, pattern = "^mt-")
data.seurat[["percent.ribo"]] <- PercentageFeatureSet(data.seurat, pattern = "^Rp[ls]")

# VlnPlot(data.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), pt.size = 0.001)

# Filter cells for feature count, and UMI count
data.seurat <- subset(data.seurat, subset = nFeature_RNA > nFeature_lo & nFeature_RNA < nFeature_hi)
data.seurat <- subset(data.seurat, subset = nCount_RNA > nCount_lo & nCount_RNA < nCount_hi)

print(paste("Number of cells after filtering for feature and UMI count", ncol(data.seurat@assays$RNA)))

# Filter cells for mitochondrial and ribosomal percentage
data.seurat <- subset(data.seurat, subset = percent.mito < mito & percent.ribo < ribo)

print(paste("Number of cells after filtering for",
            "mitochondrial percentage and ribosomal percentage",
            ncol(data.seurat@assays$RNA)))

# Filter cells for expression of CD8 marker genes
data.seurat <- subset(data.seurat, subset = Cd2 > 0 & Cd8a > 0 & Cd8b1 > 0 & Cd4 < exp(-10))
print(paste("Number of cells after filtering for CD8 marker gene expression", ncol(data.seurat@assays$RNA)))

write10xCounts(output_dir, data.seurat@assays$RNA@counts, overwrite = TRUE)
