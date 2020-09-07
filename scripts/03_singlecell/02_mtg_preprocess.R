library(tidyverse)
library(SingleCellExperiment)
library(Seurat)


# set up directories
base_dir    = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
dropseq_dir = paste0(base_dir, '/data/singlecell')

homeFolder = paste0(base_dir,"/external/MTG_celltypes")
knitr::opts_knit$set(root.dir = homeFolder)

source("/gpfs/milgram/project/holmes/kma52/mdd_sst/external/MTG_celltypes/R/Support_code0_functions.R")
suppressPackageStartupMessages({
  library(feather)
  library(matrixStats)
  library(dplyr)
  library(edgeR)
  library(data.table)
})
options(stringsAsFactors=FALSE)


# read MTG data
tmp      = fread(paste0(dropseq_dir, '/ahba_mtg/human_MTG_2018-06-14_exon-matrix.csv'),header = T, sep = ',',verbose=FALSE)
exons    = as.matrix(tmp[,2:dim(tmp)[2]])
rownames(exons) = as.character(as.matrix(tmp[,1]))
tmp      = fread(paste0(dropseq_dir, '/ahba_mtg/human_MTG_2018-06-14_intron-matrix.csv'),header = T, sep = ',',verbose=FALSE)
introns  = as.matrix(tmp[,2:dim(tmp)[2]])
rownames(introns) = as.character(as.matrix(tmp[,1]))
geneInfo = read.csv(paste0(dropseq_dir, '/ahba_mtg/human_MTG_2018-06-14_genes-rows.csv'),row.names=1)
sampInfo = read.csv(paste0(dropseq_dir, '/ahba_mtg/human_MTG_2018-06-14_samples-columns.csv'),row.names=1)


# Omit cells with no class
kp = sampInfo$cluster!="no class"


# Format the cluster info
anno = auto_annotate(sampInfo[kp,])
anno$sample_id = anno$sample_name
# Update the correct cluster colors and ids
load(paste0(dropseq_dir, "/ahba_mtg/clusterInfoMTG.RData"), verbose=T)        #### May need to adjust location
anno$cluster_color = clusterInfoMTG$cluster_color[match(anno$cluster_label,clusterInfoMTG$cluster_label)]
anno$cluster_id    = clusterInfoMTG$cluster_id[match(anno$cluster_label,clusterInfoMTG$cluster_label)]



## Calculate CPM
CPM = cpm(introns[, kp] + exons[, kp])
rownames(CPM) = rownames(geneInfo)
colnames(CPM) = anno$sample_id

## Format appropriately
expr.data = transpose(as.data.frame(CPM))
colnames(expr.data) = rownames(CPM)
expr.data$sample_id = anno$sample_id



# Write annotation file
write_feather(anno,paste0(dropseq_dir, '/ahba_mtg/anno.feather'))
write_csv(anno, paste0(dropseq_dir, '/ahba_mtg/anno.csv'))

# Write data file
write_feather(expr.data,paste0(dropseq_dir, '/ahba_mtg/data.feather'))
write_csv(expr.data, paste0(dropseq_dir, '/ahba_mtg/data.csv'))










