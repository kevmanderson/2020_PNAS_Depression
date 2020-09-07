library(EWCE)
library(MAGMA.Celltyping)
library(seurat)
data(all_hgnc_wtEntrez)


# tmp
# tmp
mdd_gwas = read_delim('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/MDD2018_ex23andMe', delim='\t')
mdd_gwas$N = mdd_gwas$Nca + mdd_gwas$Nco
as.data.frame(head(mdd_gwas))
mdd_write = mdd_gwas[c('SNP','CHR','BP','P','N')]
write_delim(mdd_write, '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.txt', delim='\t')
# tmp
# tmp


#############
# set up dirs
#############
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
umi_dir  = file.path(base_dir, '/data/singlecell/')
ldsc_dir = paste0(base_dir, '/data/ldsc')



# MAGMA output
genome_ref_path    = paste0(base_dir, '/data/magma/g1000_eur')
gwas_sumstats_path = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.txt'


###############
# Lake 2018 DFC
###############
# Cell specific gene expression
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_lake_dfc_ctd_batch.rda'), verbose=T)
dfc_ctd = ctd
dfc_ctd = prepare.quantile.groups(dfc_ctd, specificity_species="human",numberOfBins=40)


# gene-set property analysis
cells = colnames(dfc_ctd[[2]]$mean_exp)
dfc_mean_exp         = as.data.frame(dfc_ctd[[2]]$mean_exp)
dfc_mean_exp$Average = rowMeans(dfc_mean_exp)
dfc_mean_exp$hgnc    = rownames(dfc_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
dfc_mean_exp_name = merge(x=dfc_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
dfc_gene_out   = dfc_mean_exp_name[c('entrezgene', cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/lake_dfc_expr_property.txt')
write_delim(x=dfc_gene_out, gene_cond_file, delim='\t')


# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
magma --gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_dfc_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_dfc_expr_property
# run in terminal

#
#
# more fine-grained categories
fine_cells = colnames(dfc_ctd[[1]]$mean_exp)
dfc_fine_mean_exp         = as.data.frame(dfc_ctd[[1]]$mean_exp)
dfc_fine_mean_exp$Average = rowMeans(dfc_fine_mean_exp)
dfc_fine_mean_exp$hgnc    = rownames(dfc_fine_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
dfc_fine_mean_exp_name = merge(x=dfc_fine_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
dfc_fine_gene_out   = dfc_fine_mean_exp_name[c('entrezgene', fine_cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/lake_dfc_fine_expr_property.txt')
write_delim(x=dfc_fine_gene_out, gene_cond_file, delim='\t')

# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/external/magma \
--gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_dfc_fine_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_dfc_fine_expr_property
# run in terminal






##################
# Lake 2018 - VIS
##################

# Cell specific gene expression
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_lake_vis_ctd_batch.rda'), verbose=T)
vis_ctd = ctd


# gene-set property analysis
cells = colnames(vis_ctd[[2]]$mean_exp)
vis_mean_exp         = as.data.frame(vis_ctd[[2]]$mean_exp)
vis_mean_exp$Average = rowMeans(vis_mean_exp)
vis_mean_exp$hgnc    = rownames(vis_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
vis_mean_exp_name = merge(x=vis_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
vis_gene_out   = vis_mean_exp_name[c('entrezgene', cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/lake_vis_expr_property.txt')
write_delim(x=vis_gene_out, gene_cond_file, delim='\t')


# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
magma --gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_vis_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_vis_expr_property
# run in terminal



#
#
# more fine-grained categories
fine_cells = colnames(vis_ctd[[1]]$mean_exp)
vis_fine_mean_exp         = as.data.frame(vis_ctd[[1]]$mean_exp)
vis_fine_mean_exp$Average = rowMeans(vis_fine_mean_exp)
vis_fine_mean_exp$hgnc    = rownames(vis_fine_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
vis_fine_mean_exp_name = merge(x=vis_fine_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
vis_fine_gene_out   = vis_fine_mean_exp_name[c('entrezgene', fine_cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/lake_vis_fine_expr_property.txt')
write_delim(x=vis_fine_gene_out, gene_cond_file, delim='\t')

# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/external/magma \
--gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_vis_fine_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/lake_vis_fine_expr_property
# run in terminal







##########
# ABA MTG
##########

# specific gene expression
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_aba_mtg_ctd.rda'), verbose=T)
mtg_ctd = ctd


# gene-set property analysis
cells = colnames(mtg_ctd[[2]]$mean_exp)
mtg_mean_exp         = as.data.frame(mtg_ctd[[2]]$mean_exp)
mtg_mean_exp$Average = rowMeans(mtg_mean_exp)
mtg_mean_exp$hgnc    = rownames(mtg_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
mtg_mean_exp_name = merge(x=mtg_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
mtg_gene_out   = mtg_mean_exp_name[c('entrezgene', cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/aba_mtg_expr_property.txt')
write_delim(x=mtg_gene_out, gene_cond_file, delim='\t')


# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
magma --gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/aba_mtg_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/aba_mtg_expr_property
# run in terminal




# gene-set property analysis
fine_cells = colnames(mtg_ctd[[1]]$mean_exp)
mtg_fine_mean_exp         = as.data.frame(mtg_ctd[[1]]$mean_exp)
mtg_fine_mean_exp$Average = rowMeans(mtg_fine_mean_exp)
mtg_fine_mean_exp$hgnc    = rownames(mtg_fine_mean_exp)

# map hgnc to entrez id (same as what is done in EWCE)
mtg_fine_mean_exp_name = merge(x=mtg_fine_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')

# format gene-level data and write for MAGMA analysis
mtg_fine_gene_out   = mtg_fine_mean_exp_name[c('entrezgene', fine_cells, 'Average')]
gene_cond_file = paste0(base_dir,'/data/magma/aba_mtg_fine_expr_property.txt')
write_delim(x=mtg_fine_gene_out, gene_cond_file, delim='\t')


# run in terminal
cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma
/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/external/magma \
--gene-results /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/magma_MDD2018_ex23andMe.genes.raw \
--gene-covar /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/aba_mtg_fine_expr_property.txt \
--model direction=pos condition='Average' \
--out /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/aba_mtg_fine_expr_property
# run in terminal











