# Tumor 2

RNA_gene_matrix <- read.delim("RNA_gene_matrix.txt", header=TRUE)

# make comb_barcodes file
RNA_Tcells <- read.table("RNA_Tcells.txt", quote="\"", comment.char="")
RNA_RNA2a <- read.table("RNA_RNA2a.txt", quote="\"", comment.char="")
RNA_RNA2b <- read.table("RNA_RNA2b.txt", quote="\"", comment.char="")

comb_RNA2 <- rbind(RNA_RNA2a, RNA_RNA2b)
comb_RNA2$cell_type <- "Tumor"

RNA_Tcells$cell_type <- "Immune"

comb_barcodes <- rbind(comb_RNA2, RNA_Tcells)

comb_barcodes$V1 <- substr(comb_barcodes$V1, 1, nchar(comb_barcodes$V1) - 2)
colnames(comb_barcodes) <- c("cell", "cell_type")
colnames(comb_barcodes) <- NULL 

cols_to_keep <- comb_barcodes$cell
RNA_gene_matrix <- RNA_gene_matrix[, colnames(RNA_gene_matrix) %in% cols_to_keep]

# create gene order file 
coords_ucsc <- read.delim("coords_ucsc.txt", header=FALSE)
kgxref <- read.delim("kgxref.txt", header=FALSE, comment.char="#")

colnames(kgxref) <- c("Old", "New")
colnames(coords_ucsc) <- c("Gene", "Prefix", "No1", "No2")

# perform the mapping
library(dplyr)
mapped_df <- coords_ucsc %>%
  left_join(kgxref, by = c("Gene" = "Old")) %>%
  mutate(Gene = coalesce(New, Gene)) %>%
  dplyr::select(-New)

# put together the prefix and gene name
View(mapped_df)
mapped_df$Gene_Name <- paste(mapped_df$Prefix, mapped_df$Gene, sep = ".")
View(mapped_df)

mapped_df <- mapped_df %>%
  dplyr::select(Gene_Name, everything())

mapped_df$Gene <- NULL
prefixes <- mapped_df$Prefix

mapped_df <- mapped_df %>%
  group_by(Gene_Name, Prefix) %>%
  summarise(
    Min_Value1 = min(No1),
    Max_Value2 = max(No2),
    .groups = "drop"  
  )

colnames(mapped_df) <- NULL

write.table(mapped_df, file = file.path("Inputs 2", "gene_order_file"), sep = "\t")
write.table(RNA_gene_matrix, file = file.path("Inputs 2", "RNA_gene_data"), sep = "\t")
write.table(comb_barcodes, file = file.path("Inputs 2", "comb_barcodes"), sep = "\t")

raw_counts <- read.delim("Inputs 2/RNA_gene_data", header=TRUE)
annotations <- read.delim("Inputs 2/comb_barcodes", header=FALSE)
gene_order <- read.delim("Inputs 2/gene_order_file", header=FALSE)

library(dplyr)

annotations <- annotations %>% slice(-1)
colnames(annotations) <- NULL
rownames(annotations) <- NULL

write.table(annotations, file = file.path("Inputs 2", "ann_v2.txt"), row.names = TRUE, col.names = FALSE, sep = "\t")

gene_order <- gene_order[,-1]
rownames(gene_order) <- gene_order$V2
gene_order$V2 <- NULL
colnames(gene_order) <- gene_order[1, ]

write.table(gene_order, file = file.path("Inputs 2", "go_v3.txt"), row.names = TRUE, col.names = FALSE, sep = "\t")

new_infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="Inputs 2/RNA_gene_data",
                                         annotations_file="Inputs 2/ann_v2.txt",
                                         delim="\t",
                                         gene_order_file="Inputs 2/go_v3.txt",
                                         ref_group_names=c("Immune")) 

new_infercnv_obj <- infercnv::run(new_infercnv_obj,
                                  cutoff=0.01, 
                                  out_dir="Output 1", 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  HMM=TRUE)

##################

# Tumor 5

RNA_gene_matrix_5 <- read.delim("RNA_gene_matrix_t5.txt", header=TRUE)

# merge barcodes
RNA_EC <- read.table("RNA_EC_5.txt", quote="\"", comment.char="")
RNA_Macro <- read.table("RNA_Macrophage_5.txt", quote="\"", comment.char="")
RNA_MKC <- read.table("RNA_Megakaryocyte_5.txt", quote="\"", comment.char="")
RNA_RNA5a <- read.table("RNA_RNA5a.txt", quote="\"", comment.char="")
RNA_RNA5b <- read.table("RNA_RNA5b.txt", quote="\"", comment.char="")
RNA_RNA5c <- read.table("RNA_RNA5c.txt", quote="\"", comment.char="")

comb_RNA5 <- rbind(RNA_RNA5a, RNA_RNA5b, RNA_RNA5c)
comb_RNA5$cell_type <- "Tumor"

RNA_EC$cell_type <- "EC"
RNA_Macro$cell_type <- "Macrophage"
RNA_MKC$cell_type <- "Megakaryocyte"

comb_barcodes5 <- rbind(comb_RNA5, RNA_EC, RNA_Macro, RNA_MKC)

comb_barcodes5$V1 <- substr(comb_barcodes5$V1, 1, nchar(comb_barcodes5$V1) - 2)
colnames(comb_barcodes5) <- c("cell", "cell_type")
colnames(comb_barcodes5) <- NULL 


# create gene order file
coords_ucsc <- read.delim("coords_ucsc.txt", header=FALSE)
kgxref <- read.delim("kgxref.txt", header=FALSE, comment.char="#")
colnames(kgxref) <- c("Old", "New")
colnames(coords_ucsc) <- c("Gene", "Prefix", "No1", "No2")

# mapping
library(dplyr)
mapped_df <- coords_ucsc %>%
  left_join(kgxref, by = c("Gene" = "Old")) %>%
  mutate(Gene = coalesce(New, Gene)) %>%
  dplyr::select(-New)

# put together the prefix and gene name
View(mapped_df)
mapped_df$Gene_Name <- paste(mapped_df$Prefix, mapped_df$Gene, sep = ".")

mapped_df <- mapped_df %>%
  dplyr::select(Gene_Name, everything())

mapped_df$Gene <- NULL
prefixes <- mapped_df$Prefix

mapped_df <- mapped_df %>%
  group_by(Gene_Name, Prefix) %>%
  summarise(
    Min_Value1 = min(No1),
    Max_Value2 = max(No2),
    .groups = "drop"  
  )

colnames(mapped_df) <- NULL

write.table(mapped_df, file = file.path("Final Inputs", "gene_order_file_5"), sep = "\t")
write.table(RNA_gene_matrix_5, file = file.path("Final Inputs", "RNA_gene_data_5"), sep = "\t")
write.table(comb_barcodes5, file = file.path("Final Inputs", "comb_barcodes_5"), sep = "\t")

library(infercnv)

raw_counts5 <- read.delim("Final Inputs/RNA_gene_data_5", header=TRUE)
annotations5 <- read.delim("Final Inputs/comb_barcodes_5", header=FALSE)
gene_order5 <- read.delim("Final Inputs/gene_order_file_5", header=FALSE)

library(dplyr)

annotations5 <- annotations5 %>% slice(-1)
colnames(annotations5) <- NULL
rownames(annotations5) <- NULL

annotations5 <- annotations5[-1]

write.table(annotations5, file = file.path("Final Inputs", "ann5_v2.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")

gene_order5 <- gene_order5[,-1]
rownames(gene_order5) <- gene_order5$V2
gene_order5$V2 <- NULL
colnames(gene_order5) <- gene_order5[1, ]

write.table(gene_order5, file = file.path("Final Inputs", "go5_v2.txt"), row.names = TRUE, col.names = FALSE, sep = "\t")

new_infercnv_obj_5 <- CreateInfercnvObject(raw_counts_matrix="Final Inputs/RNA_gene_data_5",
                                           annotations_file="Final Inputs/ann5_v2.txt",
                                           delim="\t",
                                           gene_order_file="Final Inputs/go5_v2.txt",
                                           ref_group_names=c("EC", "Macrophage", "Megakaryocyte")) 

new_infercnv_obj_5 <- infercnv::run(new_infercnv_obj_5,
                                    cutoff=0.001, 
                                    out_dir="Output 1", 
                                    cluster_by_groups=TRUE, 
                                    denoise=TRUE,
                                    HMM=TRUE)

##################

# Defendseq data

rna <- read.delim("defendseq_RNA_matrix.csv", sep=',')
rownames(rna) <- rna[,1]
rna <- rna[ , -1]

rna_mat <- as.matrix(rna)
rna_mat_t <- t(rna_mat)

barcodes <- read.delim("defendseq_barcodes.txt")

go <- read.delim("gene_order_merged2.txt")
go_mat <- as.matrix(go)

go_sorted <- go_mat[order(go_mat[,1]), ]

# clean go_sorted
library(stringr)
gene_suffixes <- str_extract(go_sorted[,1], "[^.]+$")

go_names <- go_sorted[, 1]
gene_suffixes <- sub("^[^.]*\\.", "", go_names)

rna_rownames <- rownames(rna_mat_t)

# map suffix -> full name with prefix
replacement_map <- setNames(go_sorted[,1], gene_suffixes)

# Replace matching rownames with prefixed names; keep others as is
new_rownames <- ifelse(rna_rownames %in% gene_suffixes, replacement_map[rna_rownames], rna_rownames)

# Update rownames of rna_mat_t
rna_mat_renamed <- rna_mat_t
rownames(rna_mat_renamed) <- new_rownames

rna_mat_renamed <- rna_mat_renamed[order(rownames(rna_mat_renamed)), ]

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
go_sorted_df <- as.data.frame(go_sorted)
colnames(go_sorted_df)[1:4] <- c("gene", "chrom", "start", "end")
go_filtered <- go_sorted_df[go_sorted_df$chrom %in% standard_chroms, ]

go_filtered <- go_filtered[!duplicated(go_filtered$gene), ]
colnames(go_filtered) <- NULL

write.table(rna_mat_t, file = file.path("InferCNV_defendseq", "RNA_mat"), sep = "\t")
write.table(barcodes, file=file.path("InferCNV_defendseq", "barcodes"), sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(go_filtered, file=file.path("InferCNV_defendseq", "gene_order"), sep = "\t", row.names=FALSE, col.names=FALSE)

library(infercnv)

defendseq_infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="InferCNV_defendseq/RNA_mat",
                                         annotations_file="InferCNV_defendseq/barcodes",
                                         delim="\t",
                                         gene_order_file="InferCNV_defendseq/gene_order",
                                         ref_group_names=c("Myeloid", "Oligodendrocyte"))

defendseq_infercnv_obj <- infercnv::run(defendseq_infercnv_obj,
                                  cutoff=0.01, 
                                  out_dir="InferCNV_defendseq/output1",
                                  cluster_by_groups=TRUE,
                                  denoise=TRUE,
                                  HMM=TRUE)
