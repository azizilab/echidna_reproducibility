#!/bin/bash

# Tumor 2
Rscript numbat_files/numbat_package/numbat/bin/pileup_and_phase.R \
    --label tumor2_run1 \
    --samples tumor2 \
    --bams numbat_files/short.bam \
    --barcodes numbat_files/ann_no_dup.tsv \
    --outdir numbat_files/numbat_t2_r1 \
    --gmap numbat_files/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
    --eagle numbat_files/Eagle_v2.4.1/eagle \
    --snpvcf numbat_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf \
    --paneldir numbat_files/bcf_files \
    --ncores 20

##########################

# Tumor 5
Rscript numbat_files/numbat_package/numbat/bin/pileup_and_phase.R \
    --label tumor5_run1 \
    --samples tumor5 \
    --bams numbat5_files/short.bam \
    --barcodes numbat5_files/ann_all_no_dup.tsv \
    --outdir numbat5_files/numbat_t5_r1 \
    --gmap numbat_files/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
    --eagle numbat_files/Eagle_v2.4.1/eagle \
    --snpvcf numbat_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf \
    --paneldir numbat_files/bcf_files \
    --ncores 20

##########################

# Defendseq data
Rscript R/x86_64-pc-linux-gnu-library/4.4/numbat/bin/pileup_and_phase.R \
    --label defendseq_run1 \
    --samples defendseq\
    --bams ~/numbat_defendseq/tagged.bam \
    --barcodes ~/numbat_defendseq/barcodes_cleaned.txt \
    --outdir ~/numbat_defendseq/sims_r1 \
    --gmap ~/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --eagle ~/Eagle_v2.4.1/eagle \
    --snpvcf ~/numbat_defendseq/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir ~/numbat_defendseq/1000_hg38/ \
    --ncores 20
