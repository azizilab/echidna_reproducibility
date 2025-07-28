# Tumor 2 

library(devtools)
install_github("navinlabcode/copykat")
library(copykat)

rna_mat = read.table("Inputs 2/RNA_gene_data", quote="\"", comment.char="")
new_rownames <- sub("^[^.]*\\.", "", rownames(rna_mat))

unique_rows <- !duplicated(new_rownames)  
rna_mat <- rna_mat[unique_rows, ]     
rownames(rna_mat) <- new_rownames[unique_rows]  

copykat_obj <- copykat(rawmat=rna_mat,id.type="S", ngene.chr=1, sam.name="t2_r1")

####################

# Tumor 5

library(devtools)
install_github("navinlabcode/copykat")
library(copykat)
rna_mat_5 = read.table("Final Inputs/RNA_gene_data_5", quote="\"", comment.char="")

new_rownames_5 <- sub("^[^.]*\\.", "", rownames(rna_mat_5))

unique_rows_5 <- !duplicated(new_rownames_5)  
rna_mat_5 <- rna_mat_5[unique_rows_5, ]     
rownames(rna_mat_5) <- new_rownames_5[unique_rows_5]  

setwd("copykat_tumor5")
copykat_obj5 <- copykat(rawmat=rna_mat_5,id.type="S", ngene.chr=1, sam.name="t5_r1")

####################

# Defendseq data 

library(devtools)
install_github("navinlabcode/copykat")
library(copykat)
rna_mat_defendseq = read.table("InferCNV_defendseq/RNA_mat", quote="\"", comment.char="")

setwd("copykat_defendseq")
copykat_obj_defendseq <- copykat(rawmat=rna_mat_defendseq,id.type="S", ngene.chr=1, sam.name="defendseq_r1")

