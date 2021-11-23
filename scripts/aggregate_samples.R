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
parser$add_argument("-o", "--output_dir", type = "character", default = here("aggregated"),
                    help = "destination directory for output")
parser$add_argument("-s", "--samples", nargs = "+", help = "sample names")
args <- parser$parse_args()
print(args)

data_dir <- args$data_dir
output_dir <- args$output_dir
samples <- args$samples
n_samples <- length(samples)

# process the first sample separately
sample1 <- samples[1]
fname <- file.path(data_dir, sample1)
data <- Read10X(fname)
data.seurat <- CreateSeuratObject(counts = data,
                                  project = sample1)
c <- 1

for (sample in samples[2:n_samples]) {
  c <- c + 1
  fname <- file.path(data_dir, sample)
  data <- Read10X(fname)
  sample_data <- CreateSeuratObject(counts = data,
                                    project = sample)
  sample_data <- RenameCells(sample_data, new.names = gsub("-1", paste0("-", c), colnames(sample_data)))
  data.seurat <- merge(data.seurat, y = sample_data)
}

write10xCounts(output_dir, data.seurat@assays$RNA@counts, overwrite = TRUE)
