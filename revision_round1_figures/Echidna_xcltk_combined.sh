#!/bin/bash

# Tumor 2
xcltk baf \
    --label t2_baf \
    --sam ~/xclone_files/short.bam \
    --barcode ~/xclone_files/ann_no_dup_unix.tsv \
    --snpvcf ~/xclone_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz \
    --region ~/xclone_files/region_file_reordered_t2.txt \
    --outdir ~/xclone_files/tumor2_baf \
    --gmap ~/xclone_files/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
    --eagle ~/xclone_files/Eagle_v2.4.1/eagle \
    --paneldir ~/bcf_files \
    --ncores 10 

xcltk basefc \
    --sam ~/xclone_files/short.bam \
    --barcode ~/xclone_files/tumor5_inputs/ann_all_no_dup.tsv \
    --region ~/xclone_files/region_file_reordered_t2.txt \
    --outdir ~/xclone_files/tumor2_rdr \
    --ncores 10

##########################

# Tumor 5
# Run cellsnp-lite manually to capture chr20:
cellsnp-lite -s ~/xclone_files/tumor5_inputs/short5.bam \
             -b ~/xclone_files/tumor5_inputs/ann_all_no_dup.tsv \
             -O ~/xclone_files/tumor5_outputs/tumor5_baf/1_pileup \
             -R ~/xclone_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz \
             -p 10 \
             --minCOUNT 5 \
             --minMAF 0.01 \
             --cellTAG CB \
             --UMItag UB \
             --gzip 

xcltk baf \
    --label t5_baf \
    --sam ~/xclone_files/tumor5_inputs/short5.bam \
    --barcode ~/xclone_files/tumor5_inputs/ann_all_no_dup.tsv \
    --snpvcf ~/xclone_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz \
    --region ~/xclone_files/tumor5_inputs/region_file_reordered_t5.txt \
    --outdir ~/xclone_files/tumor5_outputs/tumor5_baf \
    --gmap ~/xclone_files/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
    --eagle ~/xclone_files/Eagle_v2.4.1/eagle \
    --paneldir ~/bcf_files \
    --ncores 10 


xcltk basefc \
    --sam ~/xclone_files/tumor5_inputs/short5.bam \
    --barcode ~/xclone_files/tumor5_inputs/ann_all_no_dup.tsv \
    --region ~/xclone_files/tumor5_inputs/region_file_reordered_t5.txt \
    --outdir ~/xclone_files/tumor5_outputs/tumor5_rdr \
    --ncores 10 

##########################

# Defendseq data

xcltk baf \
    --label defendseq_baf \
    --sam ~/numbat_defendseq/tagged.bam \
    --barcode ~/numbat_defendseq/barcodes_cleaned.txt \
    --snpvcf ~/numbat_defendseq/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --region ~/xclone_defendseq/region_file_reordered_defendseq.txt \
    --outdir ~/xclone_defendseq/defendseq_baf \
    --gmap ~/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --eagle ~/Eagle_v2.4.1/eagle \
    --paneldir ~/numbat_defendseq/1000G_hg38/1000G_hg38 \
    --ncores 10 

xcltk basefc \
    --sam ~/numbat_defendseq/tagged.bam \
    --barcode ~/numbat_defendseq/barcodes_cleaned.txt \
    --region ~/xclone_defendseq/region_file_reordered_defendseq.txt \
    --outdir ~/xclone_defendseq/defendseq_rdr \
    --ncores 10 
