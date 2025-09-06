library(Seurat)
library(Matrix)
library(tools)
library(openxlsx)
library(rstatix)
library(dplyr)
library(readr)
library(tidyverse)
library(tidyr)
library(stringr)
library(future)
library(tibble)
library(data.table)
library(GSEABase)
library(clusterProfiler)
library(entropy)
plan(sequential)



############################################################################
# Data load
############################################################################

# Collapse gene names into one row using median count values.
collapse_by_median <- function(mat) {
  base_names <- sub("\\..*", "", rownames(mat))
  u_base_names <- unique(base_names)
  
  if (length(u_base_names) == nrow(mat)) return(mat)
  
  dt <- as.data.table(as.matrix(mat))
  dt[, gene := base_names]
  
  dt_median <- dt[, lapply(.SD, median), by = gene]
  
  result <- as.matrix(dt_median[, -1])
  rownames(result) <- dt_median$gene
  colnames(result) <- colnames(mat)
  
  result <- Matrix::Matrix(result, sparse = TRUE)
  
  return(result)
}


# Load GSE161529
#---------------------------------------------------------------------------
# Files for samples of interest were downloaded only and then unarchived
# Make symbol gene names unique
features161529 <- read.delim('GSE161529/GSE161529_features.tsv', header = FALSE)
features161529$V2 <- make.unique(features161529$V2)
write.table(features161529, file = "GSE161529/GSE161529_features_unique.tsv",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Reformat sample names for futher match with clinical metadata
clean_sample_name_161529 <- function(name) {
  parts <- strsplit(name, "-", fixed = TRUE)[[1]]

  if (length(parts) == 3 && parts[1] == "TN" && parts[2] == "B1") {
    # TN-B1-XYZnnnn to TN-B1-nnnn
    cleaned_last <- sub("^[A-Za-z]+", "", parts[3])
    return(paste(parts[1], parts[2], cleaned_last, sep = "-"))

  } else if (length(parts) >= 2 && parts[1] == "TN") {
    # TN-XYZnnnn[-suffix] to TN-nnnn[-suffix]
    cleaned_second <- sub("^[A-Za-z]+", "", parts[2])
    suffix <- if (length(parts) > 2) paste(parts[3:length(parts)], collapse = "-") else NULL
    result <- paste("TN", cleaned_second, sep = "-")
    if (!is.null(suffix)) result <- paste(result, suffix, sep = "-")
    return(result)

  } else {
    # Change other types
    cleaned_second <- sub("^[A-Za-z]+", "", parts[2])
    suffix <- if (length(parts) > 2) paste(parts[3:length(parts)], collapse = "-") else NULL
    result <- paste(parts[1], cleaned_second, sep = "-")
    if (!is.null(suffix)) result <- paste(result, suffix, sep = "-")
    return(result)
  }
}

features161529 <- 'GSE161529/GSE161529_features_unique.tsv'
matrix_files161529 <- list.files("GSE161529/suppl", pattern = "*.mtx$", full.names = TRUE)
barcode_files161529 <- list.files("GSE161529/suppl", pattern = "barcodes\\.tsv$", full.names = TRUE)

matrix_files161529 <- sort(matrix_files161529)
barcode_files161529 <- sort(barcode_files161529)

sample_names161529 <- sub("^[^_]+_([^-]+-.*)-matrix\\.mtx$", "\\1", basename(matrix_files161529))
sample_names161529 <- sapply(sample_names161529, clean_sample_name_161529, USE.NAMES = FALSE)
sample_names161529[sample_names161529 == "ER-0040"] <- "ER-0040-T"

# Create matrix for each sample
create_matrix_161529 <- function(mtx_file, barcode_file, sample_name, features) {
  mat <- ReadMtx(
    mtx = mtx_file,
    features = features,
    cells = barcode_file,
    feature.column = 2
  )
  colnames(mat) <- paste(sample_name, colnames(mat), sep = "_")
  mat <- collapse_by_median(mat)
  print(sample_name)
  return(mat)
}

count_list_161529 <- mapply(
  create_matrix_161529,
  mtx_file = matrix_files161529,
  barcode_file = barcode_files161529,
  sample_name = sample_names161529,
  features = features161529,
  SIMPLIFY = FALSE
)

# Check that gene names are the same and in the same order
ref_rownames161529 <- rownames(count_list_161529[[1]])
all(sapply(count_list_161529, function(mat) {
  identical(rownames(mat), ref_rownames161529)
}))
# TRUE

# Make a single matrix
counts_161529 <- do.call(cbind, count_list_161529)
GSE161529 <- CreateSeuratObject(counts = counts_161529, project = "GSE161529")
GSE161529$sample_id <- sapply(strsplit(colnames(GSE161529), "_"), `[`, 1)
GSE161529$dataset_id <- "GSE161529"
rm(count_list_161529)
rm(counts_161529)
gc()


# Load GSE176078
#---------------------------------------------------------------------------

archive_paths176078 <- list.files("GSE176078/GSE176078_RAW/", pattern = "*.tar.gz$", full.names = TRUE)

excluded_samples176078 <- c("CID3963", "CID4066", "CID4398", "CID4513", "CID4523")
exclude_pattern176078 <- paste(excluded_samples176078, collapse = "|")

archive_paths176078 <- archive_paths176078[!grepl(exclude_pattern176078, archive_paths176078)]

count_list_176078 <- list()

for (archive in archive_paths176078) {
  sample_name <- file_path_sans_ext(basename(archive))
  sample_name <- file_path_sans_ext(sample_name)
  sample_name <- sub("^.*?_", "", sample_name)
  
  extract_dir <- file.path(tempdir(), sample_name)
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  untar(archive, exdir = extract_dir)
  
  inner_dir <- list.dirs(extract_dir, full.names = TRUE, recursive = FALSE)[1]
  matrix_file   <- file.path(inner_dir, "count_matrix_sparse.mtx")
  genes_file    <- file.path(inner_dir, "count_matrix_genes.tsv")
  barcodes_file <- file.path(inner_dir, "count_matrix_barcodes.tsv")
  
  mat <- ReadMtx(
    mtx = matrix_file,
    features = genes_file,
    cells = barcodes_file,
    feature.column = 1
  )
  
  mat <- collapse_by_median(mat)
  print(sample_name)
  
  count_list_176078[[sample_name]] <- mat
  
  unlink(extract_dir, recursive = TRUE)
}

ref_rownames176078 <- rownames(count_list_176078[[1]])
all(sapply(count_list_176078, function(mat) {
  identical(rownames(mat), ref_rownames176078)
}))
# TRUE

# Make a single matrix
counts_176078 <- do.call(cbind, count_list_176078)
GSE176078 <- CreateSeuratObject(counts = counts_176078, project = "GSE176078")
GSE176078$sample_id <- sapply(strsplit(colnames(GSE176078), "_"), `[`, 1)
GSE176078$dataset_id <- "GSE176078"
rm(count_list_176078)
rm(counts_176078)
gc()



# Load GSE167036
#---------------------------------------------------------------------------
  
counts167036 <- readMM("GSE167036/GSE167036_cell_counts_matrix.mtx")
features167036 <- read.csv("GSE167036/GSE167036_features.csv")
barcodes167036 <- read.csv("GSE167036/GSE167036_barcodes.csv")

# Reformat cell_id
barcodes167036$cell_id <- gsub("_", "-", barcodes167036$cell_id)

rownames(counts167036) <- features167036$x
colnames(counts167036) <- barcodes167036$cell_id

metadata167036 <- read.csv("GSE167036/GSE167036_meta_all.csv")
metadata167036$X <- NULL

metadata167036$cell_id <- gsub("_", "-", metadata167036$cell_id)
rownames(metadata167036) <- metadata167036$cell_id

# Exclude TCR cells
tcr167036 <- read.csv("GSE167036/GSE167036_meta_tcr.csv")
tcr167036$cell_id <- gsub("_", "-", tcr167036$cell_id)
cells_to_exclude <- tcr167036$cell_id

metadata167036 <- metadata167036[!(rownames(metadata167036) %in% cells_to_exclude), ]
counts167036 <- counts167036[, !(colnames(counts167036) %in% cells_to_exclude)]

all(rownames(metadata167036) == colnames(counts167036))
# TRUE

metadata167036$sample_id <- paste(metadata167036$patient_id, metadata167036$sample_type, sep = "-")
metadata167036$sample_id <- gsub("[ _]", "-", metadata167036$sample_id)

# Change cell_id to <sample_id>_<cell_id>
new_cell_ids <- paste(metadata167036$sample_id, metadata167036$cell_id, sep = "_")

colnames(counts167036) <- new_cell_ids
rownames(metadata167036) <- new_cell_ids

# Filter metadata
metadata167036 <- metadata167036[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample_id")]

all(rownames(metadata167036) == colnames(counts167036))
# TRUE

counts167036 <- collapse_by_median(counts167036)

GSE167036 <- CreateSeuratObject(
  counts = counts167036,
  project = "GSE167036"
)
GSE167036$sample_id <- sapply(strsplit(colnames(GSE167036), "_"), `[`, 1)
GSE167036$dataset_id <- "GSE167036"
rm(counts167036)
gc()


# Load Lambrecht_PD1
#---------------------------------------------------------------------------

lampd1 <- readRDS('Lambrecht_PD1/1863-counts_cells_cohort1.rds')
# Remove on-treatment samples
lampd1 <- lampd1[, grepl("Pre", colnames(lampd1))]

# Reformat cell names to distinguish sample_id
rename_id_lampd1 <- function(id) {
  
  last_underscore_pos <- max(gregexpr("_", id)[[1]])
  
  sample_id_part <- substr(id, 1, last_underscore_pos - 1)
  cell_id_part <- substr(id, last_underscore_pos + 1, nchar(id))
  
  sample_id_fixed <- gsub("_", "-", sample_id_part)
  
  paste0(sample_id_fixed, "_", cell_id_part)
}

colnames(lampd1) <- sapply(colnames(lampd1), rename_id_lampd1, USE.NAMES = FALSE)

lampd1 <- collapse_by_median(lampd1)

LambrechtPD1 <- CreateSeuratObject(
  counts = lampd1,
  project = "LambrechtPD1"
)

LambrechtPD1@meta.data$sample_id <- LambrechtPD1@meta.data$orig.ident
LambrechtPD1$dataset_id <- "LambrechtPD1"
rm(lampd1)
gc()


# Load GSE114727 samples BC09, BC10, BC11
#---------------------------------------------------------------------------
data_dir114727 <- "GSE114727/BC091011/"

matrix_files114727 <- list.files(data_dir114727, pattern = "_matrix.mtx.gz$", full.names = TRUE)
sample_names114727 <- str_replace(basename(matrix_files114727), "_matrix.mtx.gz", "")

read_sample_matrix_114727 <- function(sample) {
  matrix_path <- file.path(data_dir114727, paste0(sample, "_matrix.mtx.gz"))
  genes_path <- file.path(data_dir114727, paste0(sample, "_genes.tsv.gz"))
  barcodes_path <- file.path(data_dir114727, paste0(sample, "_barcodes.tsv.gz"))
  
  expression_matrix <- readMM(matrix_path)
  genes <- read.delim(genes_path, header = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE)
  
  rownames(expression_matrix) <- make.unique(genes$V2)
  colnames(expression_matrix) <- barcodes$V1
  colnames(expression_matrix) <- paste(sample, colnames(expression_matrix), sep = "_")
  
  expression_matrix <- collapse_by_median(expression_matrix)
  
  return(expression_matrix)
}

count_list_114727 <- lapply(sample_names114727, read_sample_matrix_114727)

ref_rownames114727 <- rownames(count_list_114727[[1]])
all(sapply(count_list_114727, function(mat) {
  identical(rownames(mat), ref_rownames114727)
}))
# TRUE

counts_114727 <- do.call(cbind, count_list_114727)

# Change cell names for consistency
clean_cell_name_114727 <- function(name) {
  name_no_prefix <- sub("^[^_]+_", "", name)
  sub("([^_]+)_", "\\1-", name_no_prefix)
}

colnames(counts_114727) <- unname(sapply(colnames(counts_114727), clean_cell_name_114727))

GSE114727 <- CreateSeuratObject(
  counts = counts_114727,
  project = "GSE114727"
)

# Add sample_id metadata 
GSE114727$sample_id <- sapply(strsplit(colnames(GSE114727), "_"), `[`, 1)
GSE114727$dataset_id <- "GSE114727"
rm(count_list_114727)
rm(counts_114727)
gc()


# Load GSE110686
#---------------------------------------------------------------------------

data_dir110686 <- "GSE110686"

matrix_files110686 <- list.files(data_dir110686, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
sample_names110686 <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files110686))

matrix_files110686 <- sort(matrix_files110686)
sample_names110686 <- sort(sample_names110686)

genes_paths <- file.path(data_dir110686, paste0(sample_names110686, "_genes.tsv.gz"))
barcodes_paths <- file.path(data_dir110686, paste0(sample_names110686, "_barcodes.tsv.gz"))

# Read count matrix
read_sample_matrix_110686 <- function(mtx_file, genes_file, barcodes_file, sample_name) {
  mat <- readMM(mtx_file)
  genes <- read.table(genes_file, header = FALSE, sep = "\t")
  barcodes <- read.table(barcodes_file, header = FALSE, sep = "\t")$V1
  
  rownames(mat) <- make.unique(genes$V2)
  colnames(mat) <- paste0(sample_name, "_", barcodes)
  
  mat <- collapse_by_median(mat)
  
  return(mat)
}

count_list_110686 <- mapply(
  read_sample_matrix_110686,
  mtx_file = matrix_files110686,
  genes_file = genes_paths,
  barcodes_file = barcodes_paths,
  sample_name = sample_names110686,
  SIMPLIFY = FALSE
)

all(rownames(count_list_110686[[1]]) == rownames(count_list_110686[[2]]))
# [1] TRUE

counts_110686 <- do.call(cbind, count_list_110686)

colnames(counts_110686) <- sub("^GSM[0-9]+_", "", colnames(counts_110686))


GSE110686 <- CreateSeuratObject(
  counts = counts_110686,
  project = "GSE110686"
)

GSE110686$sample_id <- sapply(strsplit(colnames(GSE110686), "_"), `[`, 1)
GSE110686$dataset_id <- "GSE110686"
rm(counts_110686)
rm(count_list_110686)
gc()


# Load GSE148673
#---------------------------------------------------------------------------
  
files148673 <- list.files('GSE148673/counts/', full.names = TRUE)

expr_list148673 <- list()

for (file in files148673) {
  sample_name <- sub(".*_([^_]+)\\.txt\\.gz$", "\\1", basename(file))
  
  expr <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  expr <- expr[3:nrow(expr), , drop = FALSE]
  expr[] <- lapply(expr, as.numeric)
  expr <- Matrix(as.matrix(expr), sparse = TRUE)
  
  colnames(expr) <- paste(sample_name, colnames(expr), sep = "_")
  
  expr <- collapse_by_median(expr)
  
  expr_list148673[[sample_name]] <- expr
}

all_genes148673 <- Reduce(union, lapply(expr_list148673, rownames))

# Fill in missing genes with zeros and reorder for matching
fill_missing_genes_148673 <- function(mat) {
  missing_genes <- setdiff(all_genes148673, rownames(mat))
  if (length(missing_genes) > 0) {
    filler <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                     dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, filler)
  }
  mat[all_genes148673, , drop = FALSE]
}

expr_list148673_filled <- lapply(expr_list148673, fill_missing_genes_148673)

ref_rownames148673 <- rownames(expr_list148673_filled[[1]])
all(sapply(expr_list148673_filled, function(mat) {
  identical(rownames(mat), ref_rownames148673)
}))
# TRUE

counts_148673 <- do.call(cbind, expr_list148673_filled)

GSE148673 <- CreateSeuratObject(counts = counts_148673, project = "GSE148673")
GSE148673$sample_id <- sapply(strsplit(colnames(GSE148673), "_"), `[`, 1)
GSE148673$dataset_id <- "GSE148673"
rm(counts_148673)
rm(expr_list148673)
rm(expr_list148673_filled)



# Make a single list of Seurat objects
#---------------------------------------------------------------------------

seurat_list_all <- list(GSE161529, GSE176078, GSE167036, LambrechtPD1, 
                        GSE114727, GSE110686, GSE148673)
names(seurat_list_all) <- c('GSE161529', 'GSE176078', 'GSE167036', 'LambrechtPD1', 
                            'GSE114727', 'GSE110686', 'GSE148673')

rm(GSE161529)
rm(GSE176078)
rm(GSE167036)
rm(LambrechtPD1)
rm(GSE114727)
rm(GSE110686)
rm(GSE148673)
gc()

# Add patient metadata
#---------------------------------------------------------------------------

# Metadata is taken from table S2 with columns names renamed in snake case
meta_full <- read.csv('meta_ts2.csv', sep=';', row.names = 1)
rownames(meta_full) <- gsub("_", "-", rownames(meta_full))
meta_full$Dataset_ID <- NULL

add_patient_metadata <- function(seu, metadata_df) {
  meta <- seu@meta.data
  meta$.__cell_barcode__ <- rownames(meta)
  
  metadata_df$sample_id <- rownames(metadata_df)
  
  merged_meta <- merge(meta, metadata_df, by = "sample_id", all.x = TRUE, sort = FALSE)
  
  rownames(merged_meta) <- merged_meta$.__cell_barcode__
  merged_meta$.__cell_barcode__ <- NULL
  
  seu@meta.data <- merged_meta
  return(seu)
}

seurat_list_all <- lapply(seurat_list_all, add_patient_metadata, meta_full)

############################################################################
# Quality control
############################################################################

# Filtering
#---------------------------------------------------------------------------

# Filtering cells
filter_min <- function(seu) {
  seu <- subset(seu, subset = nFeature_RNA >= 100)
  return(seu)
}

filter_mito <- function(seu) {
  if (!"percent.mt" %in% colnames(seu[[]])) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  }
  seu <- subset(seu, subset = percent.mt <= 10)
  return(seu)
}

filter_max <- function(seu) {
  dataset_id <- seu@project.name
  if (dataset_id == "GSE148673") {
    seu <- subset(seu, subset = nFeature_RNA <= 7000)
  } else {
    seu <- subset(seu, subset = nFeature_RNA <= 2500)
  }
  return(seu)
}

filter_2sd <- function(seu) {
  if (!"log10_nCount_RNA" %in% colnames(seu[[]])) {
    seu[["log10_nCount_RNA"]] <- log10(seu$nCount_RNA + 1)
  }
  log10_vals <- seu$log10_nCount_RNA
  upper_bound <- mean(log10_vals) + 2 * sd(log10_vals)
  lower_bound <- mean(log10_vals) - 2 * sd(log10_vals)
  
  seu <- subset(seu, subset = log10_nCount_RNA >= lower_bound & log10_nCount_RNA <= upper_bound)
  
  return(seu)
}

seurat_list_min <- lapply(seurat_list_all, filter_min)
seurat_list_mito <- lapply(seurat_list_min, filter_mito)
seurat_list_max <- lapply(seurat_list_mito, filter_max)
seurat_list_2sd <- lapply(seurat_list_max, filter_2sd)

# QC statistics
#---------------------------------------------------------------------------

# Calculate number of cells per dataset
get_cell_counts <- function(seu, colname) {
  df <- data.frame(dataset_name = seu@project.name)
  df[[colname]] <- ncol(seu)
  return(df)
}

cell_count_init <- do.call(rbind, lapply(seurat_list_all, get_cell_counts, colname = 'initial_cells'))
cell_count_min <- do.call(rbind, lapply(seurat_list_min, get_cell_counts, colname = 'cells_min_qc'))
cell_count_mito <- do.call(rbind, lapply(seurat_list_mito, get_cell_counts, colname = 'cells_mito_qc'))
cell_count_max <- do.call(rbind, lapply(seurat_list_max, get_cell_counts, colname = 'cells_max_qc'))
cell_count_2sd <- do.call(rbind, lapply(seurat_list_2sd, get_cell_counts, colname = 'cells_2sd_qc'))
cell_count_after_qc <- do.call(rbind, lapply(seurat_list_2sd, get_cell_counts, colname = 'cells_after_qc'))


cell_count_list = list(cell_count_init, cell_count_after_qc, cell_count_min, cell_count_mito,
                       cell_count_max, cell_count_2sd)
cell_count_merge <- Reduce(function(x, y) merge(x, y, by = "dataset_name", all = TRUE), cell_count_list)
cell_count_steps <- cell_count_merge[, -3]

# Add number of samples
get_sample_counts <- function(seu) {
  data.frame(
    dataset_name = seu@project.name,
    sample_number = length(unique(seu$sample_id))
  )
}
sample_num <- do.call(rbind, lapply(seurat_list_2sd, get_sample_counts))
stat_list <- list(sample_num, cell_count_merge)
cell_sample_count <- Reduce(function(x, y) merge(x, y, by = "dataset_name", all = TRUE), stat_list)

# Write table of interest
write.xlsx(cell_sample_count, file = "QC.xlsx", sheetName = "dataset_untreated", rowNames = FALSE)

rm(seurat_list_min)
rm(seurat_list_mito)
rm(seurat_list_max)
gc()


############################################################################
# Merge
############################################################################

merged_all <- merge(
  x = seurat_list_2sd[[1]],
  y = seurat_list_2sd[-1],
  add.cell.ids = names(seurat_list_2sd),
  project = "GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673"
)

rm(seurat_list_2sd)
gc()

merged_all <- NormalizeData(merged_all, verbose = T)
merged_all <- ScaleData(merged_all, verbose = T)
merged_all <- FindVariableFeatures(merged_all)
merged_all <- RunPCA(merged_all, reduction.name = "pca", verbose = T)


############################################################################
# Integration
############################################################################

# RPCA integration
#---------------------------------------------------------------------------

integrated_all_rpca <- IntegrateLayers(
  object = merged_all,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "rpca_integrated",
  verbose = T
)


# Harmony integration
#---------------------------------------------------------------------------

integrated_all_harmony <- IntegrateLayers(
  object = merged_all,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = T
)


############################################################################
# Clusterization and UMAP
############################################################################

# Merged object
#---------------------------------------------------------------------------

co3_merged <- 10
merged_all <- FindNeighbors(object = merged_all, dims = 1:co3_merged, verbose = F)
merged_all <- FindClusters(object = merged_all, verbose = F, resolution=0.2)
merged_all <- RunUMAP(object = merged_all, dims = 1:co3_merged, verbose = F)


# RPCA object
#---------------------------------------------------------------------------

co3_rpca <- 10
integrated_all_rpca <- FindNeighbors(integrated_all_rpca, dims = 1:co3_rpca, 
                                     verbose = F, reduction = "rpca_integrated")
integrated_all_rpca <- FindClusters(integrated_all_rpca, verbose = F, resolution = 0.2)
integrated_all_rpca <- RunUMAP(integrated_all_rpca, dims = 1:co3_rpca, 
                               verbose = F, reduction = "rpca_integrated")


# Harmony object
#---------------------------------------------------------------------------

co3_harmony <- 10
integrated_all_harmony <- FindNeighbors(integrated_all_harmony, dims = 1:co3_harmony, 
                                        verbose = F, reduction = "harmony")
integrated_all_harmony <- FindClusters(integrated_all_harmony, verbose = F, resolution = 0.2)
integrated_all_harmony <- RunUMAP(integrated_all_harmony, dims = 1:co3_harmony, 
                                  verbose = F, reduction = "harmony")


############################################################################
# Batch effect
############################################################################

# Cluster composition by dataset
#---------------------------------------------------------------------------
cluster_dataset_merged <- prop.table(table(merged_all@meta.data$seurat_clusters, 
                                           merged_all@meta.data$dataset_id), margin=1) * 100
cluster_dataset_merged_df <- as.data.frame(cluster_dataset_merged)
colnames(cluster_dataset_merged_df) <- c("Cluster", "Dataset", "Percentage")

cluster_dataset_rpca <- prop.table(table(integrated_all_rpca@meta.data$seurat_clusters, 
                                         integrated_all_rpca@meta.data$dataset_id), margin=1) * 100
cluster_dataset_rpca_df <- as.data.frame(cluster_dataset_rpca)
colnames(cluster_dataset_rpca_df) <- c("Cluster", "Dataset", "Percentage")


cluster_dataset_harmony <- prop.table(table(integrated_all_harmony@meta.data$seurat_clusters, 
                                            integrated_all_harmony@meta.data$dataset_id), margin=1) * 100
cluster_dataset_harmony_df <- as.data.frame(cluster_dataset_harmony)
colnames(cluster_dataset_harmony_df) <- c("Cluster", "Dataset", "Percentage")


# Entropy calculation
#---------------------------------------------------------------------------
  
entropy_meta_rpca <- integrated_all_rpca@meta.data %>%
  dplyr::select(seurat_clusters, dataset_id)

entropy_rpca <- entropy_meta_rpca %>%
  group_by(seurat_clusters) %>%
  summarise(entropy = {
    freq <- table(dataset_id)
    prob <- freq / sum(freq)
    entropy::entropy(prob, unit = "log2")
  })


entropy_meta_harmony <- integrated_all_harmony@meta.data %>%
  dplyr::select(seurat_clusters, dataset_id)

entropy_harmony <- entropy_meta_harmony %>%
  group_by(seurat_clusters) %>%
  summarise(entropy = {
    freq <- table(dataset_id)
    prob <- freq / sum(freq)
    entropy::entropy(prob, unit = "log2")
  })


# Mann-Whitney U test
entropy_harmony$method <- "Harmony"
entropy_rpca$method <- "RPCA"
entropy_all <- rbind(entropy_harmony, entropy_rpca)
wilcox.test(entropy ~ method, data = entropy_all)
# p-value = 0.6114

median_rpca <- median(entropy_rpca$entropy) # 1.781787
median_harmony <- median(entropy_harmony$entropy) # 1.754713


# Correspondence of cells between RPCA and Harmony clusters
#---------------------------------------------------------------------------
sankey_rpca <- integrated_all_rpca@meta.data["seurat_clusters"]
colnames(sankey_rpca) <- "seurat_clusters_rpca"

sankey_harmony <- integrated_all_harmony@meta.data["seurat_clusters"]
colnames(sankey_harmony) <- "seurat_clusters_harmony"

sankey_merged <- merge(
  sankey_rpca,
  sankey_harmony,
  by = "row.names",
  all = TRUE
)

rownames(sankey_merged) <- sankey_merged$Row.names
sankey_merged$Row.names <- NULL

sankey_all <- sankey_merged %>%
  count(seurat_clusters_rpca, seurat_clusters_harmony, name = "count")

# Distribution of datasets by clusters
#---------------------------------------------------------------------------

merged_selected_meta <- merged_all@meta.data %>%
  select(seurat_clusters, dataset_id)

merged_counts_dc <- merged_selected_meta %>%
  group_by(dataset_id, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop")

merged_counts_dc <- merged_counts_dc %>%
  group_by(dataset_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()



rpca_selected_meta <- integrated_all_rpca@meta.data %>%
  select(seurat_clusters, dataset_id)

rpca_counts_dc <- rpca_selected_meta %>%
  group_by(dataset_id, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop")

rpca_counts_dc <- rpca_counts_dc %>%
  group_by(dataset_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()



harmony_selected_meta <- integrated_all_harmony@meta.data %>%
  select(seurat_clusters, dataset_id)

harmony_counts_dc <- harmony_selected_meta %>%
  group_by(dataset_id, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop")

harmony_counts_dc <- harmony_counts_dc %>%
  group_by(dataset_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

rm(merged_all)
rm(integrated_all_harmony)
gc()


############################################################################
# Differentially expressed genes
############################################################################

joined_layers_rpca <- JoinLayers(integrated_all_rpca)
rm(integrated_all_rpca)
gc()


markers_rpca <- FindAllMarkers(joined_layers_rpca, only.pos = TRUE, logfc.threshold = 0.25)

# Write table of interest
write.xlsx(markers_rpca, file = "markers_rpca.xlsx", rowNames = T)

top_markers <- markers_rpca %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.xlsx(top_markers, file = "markers_rpca_top20.xlsx", rowNames = T)


############################################################################
# Cell cycle scores
############################################################################

# Calculate and assign cell scores
#---------------------------------------------------------------------------

# Get Seurat lists of S and G2/M phase genes
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# Assign cell cycle scores to normalized gene expression matrix
assign_cell_cycle_phase <- function(norm_data, s.genes, g2m.genes) {

  s.genes.use <- intersect(rownames(norm_data), s.genes)
  g2m.genes.use <- intersect(rownames(norm_data), g2m.genes)
  
  if (length(s.genes.use) == 0 || length(g2m.genes.use) == 0) {
    stop("No overlap between gene lists and expression matrix rows.")
  }
  
  s_score <- Matrix::colMeans(norm_data[s.genes.use, , drop = FALSE])
  g2m_score <- Matrix::colMeans(norm_data[g2m.genes.use, , drop = FALSE])
  
  return(data.frame(S.Score = s_score, G2M.Score = g2m_score))
}


# Assign cell cycle scores
norm_data_rpca <- joined_layers_rpca[["RNA"]]$data
cell_cycle_rpca <- assign_cell_cycle_phase(norm_data_rpca, s.genes, g2m.genes)

joined_layers_rpca$S.Score <- cell_cycle_rpca$S.Score
joined_layers_rpca$G2M.Score <- cell_cycle_rpca$G2M.Score

# Compare cell cycle scores for cluster 8, 10, 12, each vs pooled 1:7, 9, 11
#---------------------------------------------------------------------------
cell_cycle_long_rpca <- joined_layers_rpca@meta.data %>%
  select(seurat_clusters, S.Score, G2M.Score) %>%
  pivot_longer(cols = c(S.Score, G2M.Score),
               names_to = "CellCycleScore",
               values_to = "Score") %>%
  mutate(
    CellCycleScore = factor(CellCycleScore,
                            levels = c("S.Score", "G2M.Score"),
                            labels = c("S Score", "G2M Score")),
    seurat_clusters = factor(seurat_clusters)
  ) %>%
  group_by(seurat_clusters, CellCycleScore) %>%
  filter(
    Score >= mean(Score, na.rm = TRUE) - 2 * sd(Score, na.rm = TRUE),
    Score <= mean(Score, na.rm = TRUE) + 2 * sd(Score, na.rm = TRUE)
  ) %>%
  ungroup()

# Mann-Whitney tests for cluster 8, 10, 12, each vs pooled 1:7, 9, 11
comparison_clusters_mw <- as.character(c(1:7, 9, 11))

test_results_mw <- cell_cycle_long_rpca %>%
  filter(seurat_clusters %in% c(comparison_clusters_mw, "8", "10", "12")) %>%
  group_by(CellCycleScore) %>%
  summarise(
    p_cluster8 = {
      wt <- wilcox.test(
        Score[seurat_clusters == "8"],
        Score[seurat_clusters %in% comparison_clusters_mw],
        exact = TRUE
      )
      ifelse(wt$p.value < 2.2e-16, 2.2e-16, wt$p.value)
    },
    p_cluster10 = {
      wt <- wilcox.test(
        Score[seurat_clusters == "10"],
        Score[seurat_clusters %in% comparison_clusters_mw],
        exact = TRUE
      )
      ifelse(wt$p.value < 2.2e-16, 2.2e-16, wt$p.value)
    },
    p_cluster12 = {
      wt <- wilcox.test(
        Score[seurat_clusters == "12"],
        Score[seurat_clusters %in% comparison_clusters_mw],
        exact = TRUE
      )
      ifelse(wt$p.value < 2.2e-16, 2.2e-16, wt$p.value)
    },
    .groups = "drop"
  )

############################################################################
# GSEA
############################################################################


# Load immune markers
genesets <- getGmt('LITERATURE_immune_genesets_all.gmt')

genesets_df <- do.call(rbind, lapply(genesets, function(gs) {
  data.frame(
    gs_name = gs@setName,
    gene = gs@geneIds,
    stringsAsFactors = FALSE
  )
}))

genesets_df <- genesets_df[genesets_df$gene != '',]

clusters <- unique(markers_rpca$cluster)

gsea_results <- list()

for (cl in clusters) {
  cat("Running GSEA for cluster:", cl, "\n")
  
  markers_cl <- markers_rpca %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC)) %>%
    distinct(gene, .keep_all = TRUE)
  
  gene_list <- markers_cl$avg_log2FC
  names(gene_list) <- markers_cl$gene
  
  gene_list <- sort(na.omit(gene_list), decreasing = TRUE)
  
  gsea <- GSEA(
    geneList = gene_list,
    TERM2GENE = genesets_df,
    pvalueCutoff = 1,
    scoreType = "pos",
    eps = 0,
    verbose = FALSE
  )
  gsea_results[[as.character(cl)]] <- gsea
}

gsea_df <- do.call(rbind, lapply(names(gsea_results), function(cl) {
  res <- gsea_results[[cl]]@result
  if (nrow(res) > 0) {
    res$cluster <- cl
    return(res)
  } else {
    return(NULL)
  }
}))


# Filter significant GSEA results
gsea_filtered <- gsea_df[gsea_df$p.adjust < 0.25,]


############################################################################
# Gene marker expression
############################################################################

# Mean gene expression
#---------------------------------------------------------------------------
# Calculate mean gene expression per cell in cluster for each gene
mean_gene_expression <- function(expr_matrix, cluster_cells) {
  gene_expr_list <- list()
  
  for (cl in names(cluster_cells)) {
    cl_num <- as.numeric(cl)
    cl_cells <- cluster_cells[[cl]]
    
    gene_means <- round(rowMeans(expr_matrix[, cl_cells, drop = FALSE]), 4)
    
    df <- data.frame(
      gene = names(gene_means),
      cluster = cl_num,
      mean_expr = gene_means
    )
    
    gene_expr_list[[length(gene_expr_list) + 1]] <- df
  }
  
  gene_expr_df <- do.call(rbind, gene_expr_list)
  return(gene_expr_df)
}


data_matrix_rpca <- GetAssayData(joined_layers_rpca, layer = "data")

cell_clusters <- joined_layers_rpca$seurat_clusters
cluster_cells <- split(Cells(joined_layers_rpca), cell_clusters)

single_gene_expr_norm <- mean_gene_expression(data_matrix_rpca, cluster_cells)
single_gene_expr_norm <- single_gene_expr_norm %>%
  rename(mean_expr_norm = mean_expr)


# Percentage of cells in cluster expressing genes
#---------------------------------------------------------------------------
pct_expr_list <- list()

for (cl in names(cluster_cells)) {
  cl_num <- as.numeric(cl)
  cl_cells <- cluster_cells[[cl]]
  
  mat <- data_matrix_rpca[, cl_cells, drop = FALSE]
  
  pct_expressed <- round(rowSums(mat > 0) / length(cl_cells) * 100, 4)
  
  df_perc <- data.frame(
    gene = names(pct_expressed),
    cluster = cl_num,
    pct_cells_expressing = pct_expressed
  )
  
  pct_expr_list[[length(pct_expr_list) + 1]] <- df_perc
}

pct_expr_df <- do.call(rbind, pct_expr_list)

single_gene_all <- left_join(single_gene_expr_norm, pct_expr_df, by = c("gene", "cluster"))

# Expression of the 3 expert-curated gene sets
#---------------------------------------------------------------------------
markers1 <- read_csv('markers/markers1.txt', show_col_types = FALSE)
single_gene_markers1 <- merge(markers1[, c('cell_type', 'gene')], 
                              single_gene_all, by = "gene") %>%
  filter(!is.na(cell_type))


markers2 <- read_csv('markers/markers2.txt', show_col_types = FALSE)
single_gene_markers2 <- merge(markers2[, c('cell_type', 'gene')], 
                              single_gene_all, by = "gene") %>%
  filter(!is.na(cell_type))


markers3 <- read_csv('markers/markers3.txt', show_col_types = FALSE)
single_gene_markers3 <- merge(markers3[, c('cell_type', 'gene')], 
                              single_gene_all, by = "gene") %>%
  filter(!is.na(cell_type))


############################################################################
# Add cluster annotations
############################################################################

cluster_annotation <- read.xlsx("cluster_annot.xlsx")

cluster_annotation$cluster_number <- as.character(cluster_annotation$cluster_number)
joined_layers_rpca$seurat_clusters <- as.character(joined_layers_rpca$seurat_clusters)

matched_annotations <- cluster_annotation[match(joined_layers_rpca$seurat_clusters, 
                                                cluster_annotation$cluster_number), ]

joined_layers_rpca$cluster_annotation_main <- matched_annotations$cluster_annotation_main
joined_layers_rpca$Leukocyte_cluster <- matched_annotations$Leukocyte_cluster
joined_layers_rpca$Lymphocyte_cluster <- matched_annotations$Lymphocyte_cluster


############################################################################
# Leukocyte general investigation
############################################################################

# Leukocyte proportions
#---------------------------------------------------------------------------
# Leukocyte cell percentage of all cells
leuko_prop <- sum(joined_layers_rpca@meta.data$Leukocyte_cluster == TRUE) / nrow(joined_layers_rpca@meta.data) * 100
# 56.24

# Leukocyte percentage by dataset
leukocyte_summary <- joined_layers_rpca@meta.data %>%
  group_by(dataset_id) %>%
  summarise(
    total_cells = n(),
    leukocyte_cells = sum(Leukocyte_cluster == TRUE),
    proportion_leukocytes = leukocyte_cells / total_cells
  )

# Percentage of leukocyte types by dataset
cell_type_order <- c("Tcell cytotoxic&NK", "Tcell naive/memory", "Bcell naive" , "Plasma Bcells", 
                     "Cycling lymphocytes", "Cycling Neutrophils", "Myeloids")

leukocyte_data <- joined_layers_rpca@meta.data %>%
  filter(Leukocyte_cluster == TRUE) %>%
  mutate(
    cluster_annotation_main = factor(cluster_annotation_main, 
                                     levels = cell_type_order)
  ) %>%
  count(dataset_id, cluster_annotation_main) %>%
  group_by(dataset_id) %>%
  mutate(
    total_cells = sum(n),
    prop = n / total_cells * 100
  ) %>%
  ungroup()

############################################################################
# Leukocytes, leukocyte cell types and cell cycle by clinical subtypes
############################################################################

# Cell type by clinical subtype
#---------------------------------------------------------------------------
figure4_meta <- joined_layers_rpca@meta.data
figure4_meta <- figure4_meta %>%
  mutate(
    leukocyte_status = ifelse(Leukocyte_cluster == TRUE, "Leukocyte", "Non-Leukocyte")
  )

clinical_subtype_order <- c("TNBC", "HER2BC", "LBC")

leuko_prop_by_subtype <- figure4_meta %>%
  group_by(clinical_subtype, leukocyte_status) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(clinical_subtype) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  mutate(clinical_subtype = factor(clinical_subtype, 
                                   levels = clinical_subtype_order))

# Cell type proportions by clinical subtype
#---------------------------------------------------------------------------

leuko_meta <- joined_layers_rpca@meta.data %>%
  filter(Leukocyte_cluster == TRUE)

leuko_counts_s <- leuko_meta %>%
  group_by(clinical_subtype, cluster_annotation_main) %>%
  summarise(n = n(), .groups = "drop")

leuko_prop_s <- leuko_counts_s %>%
  group_by(clinical_subtype) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

stack_order <- c(
  "Myeloids", "Cycling Neutrophils", "Cycling lymphocytes", "Plasma Bcells",
  "Bcell naive", "Tcell naive/memory", "Tcell cytotoxic&NK", 
  "All B cells", "All T cells"
)

# Grouped data for all B cells and all T cells
leuko_grouped <- leuko_prop_s %>%
  mutate(
    cell_group = case_when(
      grepl("Tcell", cluster_annotation_main) ~ "All T cells",
      grepl("Bcell", cluster_annotation_main) ~ "All B cells",
      TRUE ~ as.character(cluster_annotation_main)
    ),
    cell_group = factor(cell_group, levels = stack_order)
  ) %>%
  group_by(clinical_subtype, cell_group) %>%
  summarise(proportion = sum(proportion), .groups = "drop") %>%
  mutate(bar_type = "Grouped")

# Detailed data for each cell type separately
leuko_detailed <- leuko_prop_s %>%
  mutate(
    cell_group = factor(cluster_annotation_main, levels = stack_order),
    bar_type = "Detailed"
  )


leuko_tot <- bind_rows(leuko_detailed, leuko_grouped) %>%
  mutate(
    clinical_subtype = factor(clinical_subtype, levels = clinical_subtype_order),
    x_pos = case_when(
      bar_type == "Detailed" ~ as.numeric(clinical_subtype) - 0.2,
      bar_type == "Grouped"  ~ as.numeric(clinical_subtype) + 0.2
    )
  )

# Examination of cycling neutrophils
cycling_neutro <- joined_layers_rpca@meta.data %>%
  dplyr::select(sample_id, dataset_id, cluster_annotation_main, clinical_subtype)

total_cycling_neutrophils <- nrow(cycling_neutro[
  cycling_neutro$cluster_annotation_main == 'Cycling Neutrophils', ])

cycling_neutrophil_summary <- cycling_neutro %>%
  group_by(sample_id, dataset_id, clinical_subtype) %>%
  summarise(
    total_cells = n(),
    cycling_neutrophils = sum(cluster_annotation_main == "Cycling Neutrophils"),
    proportion_sample = cycling_neutrophils / total_cells,
    .groups = "drop"
  ) %>%
  group_by(dataset_id) %>%
  mutate(
    total_cycling_in_dataset = sum(cycling_neutrophils),
    proportion_dataset = cycling_neutrophils / total_cycling_in_dataset,
    proportion_total = cycling_neutrophils / total_cycling_neutrophils
  ) %>%
  ungroup() %>%
  arrange(desc(proportion_sample))

cycling_neutrophil_median_dataset <- cycling_neutrophil_summary %>%
  group_by(dataset_id) %>%
  summarise(
    median_proportion_sample = median(proportion_sample, na.rm = TRUE),
    .groups = "drop"
  )


# Leukocyte cell type proportion by clinical subtype
#---------------------------------------------------------------------------
# Group data by patient

patient_fractions <- leuko_meta %>%
  group_by(patient_id, clinical_subtype, cluster_annotation_main) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

patient_fractions <- patient_fractions %>%
  mutate(
    cluster_annotation_main = factor(cluster_annotation_main, levels = cell_type_order),
    clinical_subtype = factor(clinical_subtype, levels = clinical_subtype_order)
  )

patient_fractions <- patient_fractions %>%
  mutate(cell_group = case_when(
    grepl("Tcell", cluster_annotation_main) ~ "All T cells",
    grepl("Bcell", cluster_annotation_main) ~ "All B cells",
    TRUE ~ cluster_annotation_main
  )) %>%
  mutate(cell_group = factor(cell_group, levels = c(
    "Myeloids", "Cycling Neutrophils", "Cycling lymphocytes", "Plasma Bcells",
    "Bcell naive", "Tcell naive/memory", "Tcell cytotoxic&NK", 
    "All B cells", "All T cells"
  )))

# Mann-Whitney U tests for each cell type separately
stat_tests_detailed <- patient_fractions %>%
  group_by(cluster_annotation_main) %>%
  wilcox_test(proportion ~ clinical_subtype, p.adjust.method = "BH") %>%
  add_xy_position(x = "cluster_annotation_main")

# Mann-Whitney U tests for all B cells and T cells
stat_tests_grouped <- patient_fractions %>%
  group_by(cell_group) %>%
  wilcox_test(proportion ~ clinical_subtype, p.adjust.method = "BH") %>%
  add_xy_position(x = "cell_group")


# Cell cycle scores by clinical subtype
#---------------------------------------------------------------------------

# Cell cycle scores grouped by patients
cell_cycle_clinical_patient <- leuko_meta %>%
  group_by(patient_id, clinical_subtype) %>%
  summarise(
    S.Score = mean(S.Score, na.rm = TRUE),
    G2M.Score = mean(G2M.Score, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(S.Score, G2M.Score), 
               names_to = "Score_Type", 
               values_to = "Score") %>%
  mutate(clinical_subtype = factor(clinical_subtype, 
                                   levels = c("TNBC", "HER2BC", "LBC")))

cell_cycle_clinical_patient <- cell_cycle_clinical_patient %>%
  mutate(Score_Group = paste0(Score_Type, "_", clinical_subtype))

# Mann-Whitney U tests for cell cycle scores
stat_tests_cc <- cell_cycle_clinical_patient %>%
  group_by(Score_Type) %>%
  wilcox_test(Score ~ clinical_subtype, p.adjust.method = "BH") %>%
  add_xy_position(x = "Score_Type")

# Leukocyte cell type proportion by sample
#---------------------------------------------------------------------------
heatmap_data <- joined_layers_rpca@meta.data %>%
  filter(Leukocyte_cluster == TRUE) %>%
  count(sample_id, cluster_annotation_main) %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(sample_id, cluster_annotation_main, prop) %>%
  tidyr::pivot_wider(
    names_from = sample_id,
    values_from = prop,
    values_fill = 0
  ) %>%
  as.data.frame()

heatmap_data <- heatmap_data %>%
  mutate(cluster_annotation_main = factor(cluster_annotation_main, levels = cell_type_order)) %>%
  arrange(cluster_annotation_main)

rownames(heatmap_data) <- heatmap_data$cluster_annotation_main

# Fisher's exact test was preferred as some classes contained 0-3 counts
heatmap_clusters <- read.csv("heatmap_clusters.txt")
heatmap_clusters <- heatmap_clusters %>%
  pivot_wider(names_from = subtype, values_from = count, values_fill = 0) %>%
  column_to_rownames("cluster")

heatmap_clusters_matrix <- as.matrix(heatmap_clusters)

fisher.test(heatmap_clusters_matrix, workspace = 2e8)
# p-value = 0.02006

# Validate Fisher's test result
cluster_fisher_results <- data.frame()

for(i in 1:nrow(heatmap_clusters_matrix)){
  
  cluster_counts <- heatmap_clusters_matrix[i, ]
  other_counts <- colSums(heatmap_clusters_matrix[-i, ])
  
  contingency <- rbind(cluster_counts, other_counts)
  rownames(contingency) <- c(paste0("cluster", i), "other_clusters")
  
  ft <- fisher.test(contingency)
  
  cluster_fisher_results <- rbind(cluster_fisher_results,
                                  data.frame(
                                    cluster = rownames(heatmap_clusters_matrix)[i],
                                    p_value = ft$p.value
                                  ))
}

cluster_fisher_results$adj_p <- p.adjust(cluster_fisher_results$p_value, method = "BH")
