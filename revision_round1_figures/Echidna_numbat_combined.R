# Tumor 2

library(dplyr)
count_matrix <- read.delim("numbat_files/RNA_gene_data", row.names=1)
barcodes <- read.table("numbat_files/ann_2c_no_dup.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- c("Cell", "Type")

count_matrix <- as.data.frame(count_matrix)
colnames(count_matrix) <- sub("\\d+$", "", colnames(count_matrix))

modified_rownames <- sub("^[^.]*\\.", "", rownames(count_matrix))
duplicated_rows <- modified_rownames[duplicated(modified_rownames) | duplicated(modified_rownames, fromLast = TRUE)]
count_matrix <- count_matrix[!modified_rownames %in% duplicated_rows, ]
rownames(count_matrix) <- sub("^[^.]*\\.", "", rownames(count_matrix))

valid_cells <- barcodes %>%
  filter(Type != "Tumor") %>%
  pull(Cell)

count_matrix_filtered <- count_matrix[, colnames(count_matrix) %in% valid_cells]

cell_annot <- read.delim("numbat_files/ann_2c_no_dup.txt")
cell_annot <- as.data.frame(cell_annot)
colnames(cell_annot) <- c("cell", "group")

library(numbat)
count_matrix_filtered <- as.matrix(count_matrix_filtered)
ref_internal = aggregate_counts(count_matrix_filtered, cell_annot)

barcodes <- read.delim("numbat_files/ann_all_no_dup.txt")
barcodes <- as.data.frame(barcodes)
colnames(barcodes) <- c("cell", "group")

tsv_data <- read.delim("numbat_files/numbat_t2_r1/tumor2_allele_counts.tsv.gz", header = TRUE)
tumor_cells <- barcodes$cell[barcodes$group == "Tumor"]
tsv_data_filtered <- tsv_data[tsv_data$cell %in% tumor_cells, ]

count_matrix_tumor <- count_matrix[, colnames(count_matrix) %in% tumor_cells]
count_matrix_tumor <- as.matrix(count_matrix_tumor)

out = run_numbat(
  count_mat = count_matrix_tumor, 
  lambdas_ref = ref_internal,
  df_allele = tsv_data_filtered,
  genome = "hg19",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = "numbat_files/numbat_run1"
)

#################

# Tumor 5

library(dplyr)
count_matrix <- read.delim("numbat5_files/RNA_gene_data_5", row.names=1)
barcodes <- read.table("numbat5_files/ann_ref_no_dup.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- c("Cell", "Type")

count_matrix <- as.data.frame(count_matrix)
colnames(count_matrix) <- sub("\\d+$", "", colnames(count_matrix))

modified_rownames <- sub("^[^.]*\\.", "", rownames(count_matrix))
duplicated_rows <- modified_rownames[duplicated(modified_rownames) | duplicated(modified_rownames, fromLast = TRUE)]
count_matrix <- count_matrix[!modified_rownames %in% duplicated_rows, ]
rownames(count_matrix) <- sub("^[^.]*\\.", "", rownames(count_matrix))

valid_cells <- barcodes %>%
  filter(Type != "Tumor") %>%
  pull(Cell)

count_matrix_filtered <- count_matrix[, colnames(count_matrix) %in% valid_cells]

cell_annot <- read.delim("numbat5_files/ann_ref_no_dup.txt", header=0)
cell_annot <- as.data.frame(cell_annot)
colnames(cell_annot) <- c("cell", "group")

library(numbat)
count_matrix_filtered <- as.matrix(count_matrix_filtered)
ref_internal = aggregate_counts(count_matrix_filtered, cell_annot)

barcodes <- read.delim("numbat5_files/ann_v5.txt")
barcodes <- as.data.frame(barcodes)
colnames(barcodes) <- c("cell", "group")

tsv_data <- read.delim("numbat5_files/numbat_t5_r1/tumor5_allele_counts.tsv.gz", header = TRUE)
tumor_cells <- barcodes$cell[barcodes$group == "Tumor"]
tsv_data_filtered <- tsv_data[tsv_data$cell %in% tumor_cells, ]

count_matrix_tumor <- count_matrix[, colnames(count_matrix) %in% tumor_cells]
count_matrix_tumor <- as.matrix(count_matrix_tumor)

out = run_numbat(
  count_mat = count_matrix_tumor, 
  lambdas_ref = ref_internal,
  df_allele = tsv_data_filtered,
  genome = "hg19",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = "numbat5_run2",
  max_entropy = 0.75
)

#################

# Defendseq data

count_matrix <- read.delim("~/numbat_defendseq/RNA_mat", row.names=1)
ref_barcodes <- read.table("~/numbat_defendseq/ref_barcodes.txt", header = FALSE, stringsAsFactors = FALSE)

colnames(ref_barcodes) <- c("cell", "group")

library(dplyr)
valid_cells <- ref_barcodes %>%
  pull(cell)

count_matrix_filtered <- count_matrix[, colnames(count_matrix) %in% valid_cells]

library(numbat)
count_matrix_filtered <- as.matrix(count_matrix_filtered)

ref_internal = aggregate_counts(count_matrix_filtered, ref_barcodes)
barcodes <- read.delim("~/numbat_defendseq/barcodes.txt")

barcodes <- as.data.frame(barcodes)

colnames(barcodes) <- c("cell", "group")
tsv_data <- read.delim("~/numbat_defendseq/defendseq_r1/defendseq_allele_counts.tsv.gz", header = TRUE)

tumor_cells <- barcodes$cell[grepl("^Tumor_", barcodes$group)]
tsv_data_filtered <- tsv_data[tsv_data$cell %in% tumor_cells, ]

count_matrix_tumor <- count_matrix[, colnames(count_matrix) %in% tumor_cells]
count_matrix_tumor <- as.matrix(count_matrix_tumor)

out = run_numbat(
  count_mat = count_matrix_tumor, 
  lambdas_ref = ref_internal,
  df_allele = tsv_data_filtered, 
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = "~/numbat_defendseq/numbat_defendseq_r1",
  max_entropy = 0.75
)
