#!/bin/bash

Rscript ~/ichorCNA_run1/runIchorCNA.R \
  --id PTO084_TB6913 \
  --WIG ~/ichorCNA_run1/PTO084_TB6913_500kb.wig \
  --ploidy 2 \
  --normal "c(0.1,0.2)" \
  --maxCN 50 \
  --gcWig ~/ichorCNA_run1/extdata/gc_hg38_500kb.wig \
  --mapWig ~/ichorCNA_run1/extdata/map_hg38_500kb.wig \
  --centromere ~/ichorCNA_run1/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --includeHOMD False \
  --chrs "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')" \
  --chrTrain "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')" \
  --estimateNormal True \
  --estimatePloidy False \
  --estimateScPrevalence True \
  --scStates "c(1,3)" \
  --txnE 0.9999 \
  --txnStrength 10000 \
  --outDir ~/ichorCNA_run1/output13 
