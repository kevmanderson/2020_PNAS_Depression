library(tidyverse)
library(EWCE)
library(MAGMA.Celltyping)
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)


# This script will write text files of

####################
# set up directories
####################
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
umi_dir  = file.path(base_dir, '/data/singlecell/')
ldsc_dir = paste0(base_dir, '/data/ldsc')


# 37 gene mapping
ensembl37 = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")


# easy to map genes
hgnc_ens37_df = getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position'), mart = ensembl37)
hgnc_ens37_df = hgnc_ens37_df[which(hgnc_ens37_df$chromosome_name %in% 1:22),]
hgnc_ens37_df = hgnc_ens37_df[!is.na(hgnc_ens37_df$entrezgene_id),]
hgnc_ens37_df = hgnc_ens37_df[!duplicated(hgnc_ens37_df$hgnc_symbol),]



##################
# Lake 2018 - DFC
##################

# EWCE cell-specific gene expression
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_lake_dfc_ctd_batch.rda'), verbose=T)
dfc_ctd = ctd
dfc_ctd = prepare.quantile.groups(dfc_ctd, specificity_species="human",numberOfBins=10)


# load preprocessed data
region = 'FrontalCortex'

# current out dir
prefix      = 'lake_dfc_bigcat_ldsc'
fine_prefix = 'lake_dfc_superordinate_cat_ldsc'

# make if doesnt exist
current_dir_out     = paste0(ldsc_dir, '/', prefix, '_files')
current_finedir_out = paste0(ldsc_dir, '/', fine_prefix, '_files')


# write table with all genes and chromosomal positions
ldsc_out = hgnc_ens37_df[c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')]
colnames(ldsc_out) = c('GENE', 'CHR', 'START', 'END')
ldsc_out = ldsc_out[ldsc_out$CHR %in% as.character(1:22),]
ldsc_out = ldsc_out[order(ldsc_out$CHR),]
ldsc_out = ldsc_out[order(as.numeric(ldsc_out$CHR)),]
head(ldsc_out)


# write gene coordinates
out_path = paste0(current_dir_out, '/', prefix, '_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')
#
out_path = paste0(current_finedir_out, '/', fine_prefix, '_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')


# genes matched between hgnc coordinate data from and specificity data frame
hgnc_ens37_df_dfc   = hgnc_ens37_df %>% filter(hgnc_symbol %in% rownames(dfc_ctd[[1]]$specificity_quantiles))
hgnc_ens37_df_nodfc = hgnc_ens37_df %>% filter(!hgnc_symbol %in% rownames(dfc_ctd[[1]]$specificity_quantiles))


# subset specificity data
dfc_ctd[[1]]$specificity = dfc_ctd[[1]]$specificity[rownames(dfc_ctd[[1]]$specificity) %in% hgnc_ens37_df_dfc$hgnc_symbol,]
dfc_ctd[[2]]$specificity = dfc_ctd[[2]]$specificity[rownames(dfc_ctd[[2]]$specificity) %in% hgnc_ens37_df_dfc$hgnc_symbol,]

dfc_ctd = prepare.quantile.groups(dfc_ctd, specificity_species="human", numberOfBins=10)

# sanity check using SST
spec_check = as.data.frame(dfc_ctd[[1]]$specificity)
spec_check = spec_check[rev(order(spec_check$SST)),]
which(rownames(spec_check) == 'SST')
spec_check = spec_check[rev(order(spec_check$In6)),]
which(rownames(spec_check) == 'PVALB')
spec_check = as.data.frame(dfc_ctd[[2]]$specificity)
spec_check = spec_check[rev(order(spec_check$Inh)),]
which(rownames(spec_check) == 'GAD1')
# sanity check using SST


# FINE cell cat
nbins       = 10
genebin_num = round(nrow(hgnc_ens37_df_dfc)/nbins)
dfc_cells   = colnames(dfc_ctd[[1]]$specificity_quantiles)
quantile_df = as.data.frame(dfc_ctd[[1]]$specificity_quantiles)

for (cell in dfc_cells){
    write(cell,'')
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc     = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens37_df_dfc$ensembl_gene_id[hgnc_ens37_df_dfc$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_finedir_out, '/', fine_prefix, '_', cell, '_', str_pad(ct,2,pad="0"),'_',as.character(nbins),'_batch_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens37_df_nodfc$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_nomatch_batch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens37_df$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_allcontrol_batch_genes.txt'), quote=F, row.names=F, col.names=F)



# BIG cell cat
nbins       = 10
genebin_num = round(nrow(hgnc_ens37_df_dfc)/nbins)
dfc_cells   = colnames(dfc_ctd[[2]]$specificity_quantiles)
quantile_df = as.data.frame(dfc_ctd[[2]]$specificity_quantiles)

for (cell in dfc_cells){
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens37_df_dfc$ensembl_gene_id[hgnc_ens37_df_dfc$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_dir_out, '/', prefix, '_', cell, '_bigcat_', str_pad(ct,2,pad="0"),'_', as.character(nbins),'_batch_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens37_df_nodfc$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_nomatch_batch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens37_df$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_allcontrol_batch_genes.txt'), quote=F, row.names=F, col.names=F)







##################
# Lake 2018 - VIS
##################

# Lake VIS - specific gene expression
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_lake_vis_ctd_batch.rda'), verbose=T)
vis_ctd = ctd
vis_ctd = prepare.quantile.groups(vis_ctd, specificity_species="human",numberOfBins=10)

# load preprocessed data
region = 'VisualCortex'


# current out dir
prefix      = 'lake_vis_bigcat_ldsc'
fine_prefix = 'lake_vis_superordinate_cat_ldsc'

# make if doesnt exist
current_dir_out     = paste0(ldsc_dir, '/', prefix, '_files')
current_finedir_out = paste0(ldsc_dir, '/', fine_prefix, '_files')


# write table with all genes and chromosomal positions
ldsc_out = hgnc_ens37_df[c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')]
colnames(ldsc_out) = c('GENE', 'CHR', 'START', 'END')
ldsc_out = ldsc_out[ldsc_out$CHR %in% as.character(1:22),]
ldsc_out = ldsc_out[order(ldsc_out$CHR),]
ldsc_out = ldsc_out[order(as.numeric(ldsc_out$CHR)),]
head(ldsc_out)


# write gene coordinate file for both BIG/FINE cats
out_path = paste0(current_dir_out, '/', prefix, '_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')
#
out_path = paste0(current_finedir_out, '/', fine_prefix, '_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')



# genes matched between hgnc coordinate data from and specificity data frame
hgnc_ens37_df_vis   = hgnc_ens37_df %>% filter(hgnc_symbol %in% rownames(vis_ctd[[1]]$specificity_quantiles))
hgnc_ens37_df_novis = hgnc_ens37_df %>% filter(!hgnc_symbol %in% rownames(vis_ctd[[1]]$specificity_quantiles))


# subset specificity data
vis_ctd[[1]]$specificity = vis_ctd[[1]]$specificity[rownames(vis_ctd[[1]]$specificity) %in% hgnc_ens37_df_vis$hgnc_symbol,]
vis_ctd[[2]]$specificity = vis_ctd[[2]]$specificity[rownames(vis_ctd[[2]]$specificity) %in% hgnc_ens37_df_vis$hgnc_symbol,]

vis_ctd = prepare.quantile.groups(vis_ctd, specificity_species="human", numberOfBins=10)

# sanity check using SST
spec_check = as.data.frame(vis_ctd[[1]]$specificity)
spec_check = spec_check[rev(order(spec_check$SST)),]
which(rownames(spec_check) == 'SST')
spec_check = spec_check[rev(order(spec_check$In6)),]
which(rownames(spec_check) == 'PVALB')
spec_check = as.data.frame(vis_ctd[[2]]$specificity)
spec_check = spec_check[rev(order(spec_check$Inh)),]
which(rownames(spec_check) == 'GAD1')


# FINE cell cat
nbins       = 10
genebin_num = round(nrow(hgnc_ens37_df_vis)/nbins)
vis_cells   = colnames(vis_ctd[[1]]$specificity_quantiles)
quantile_df = as.data.frame(vis_ctd[[1]]$specificity_quantiles)

# for each cell, write gene lists for each decile (based on cell specificity)
for (cell in vis_cells){
    write(cell,'')
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc     = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens37_df_vis$ensembl_gene_id[hgnc_ens37_df_vis$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_finedir_out, '/', fine_prefix, '_', cell, '_', str_pad(ct,2,pad="0"),'_',as.character(nbins),'_batch_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens37_df_novis$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_nomatch_batch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens37_df$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_allcontrol_batch_genes.txt'), quote=F, row.names=F, col.names=F)



# BIG cell cat
nbins = 10
genebin_num = round(nrow(hgnc_ens37_df_vis)/nbins)
vis_cells   = colnames(vis_ctd[[2]]$specificity_quantiles)
quantile_df = as.data.frame(vis_ctd[[2]]$specificity_quantiles)

for (cell in vis_cells){
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens37_df_vis$ensembl_gene_id[hgnc_ens37_df_vis$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_dir_out, '/', prefix, '_', cell, '_bigcat_', str_pad(ct,2,pad="0"),'_', as.character(nbins),'_batch_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens37_df_novis$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_nomatch_batch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens37_df$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_allcontrol_batch_genes.txt'), quote=F, row.names=F, col.names=F)


# write table with all genes and chromosomal positions
ldsc_out = hgnc_ens37_df[c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')]
colnames(ldsc_out) = c('GENE', 'CHR', 'START', 'END')
ldsc_out = ldsc_out[ldsc_out$CHR %in% as.character(1:22),]
ldsc_out = ldsc_out[order(ldsc_out$CHR),]
ldsc_out = ldsc_out[order(as.numeric(ldsc_out$CHR)),]
out_path = paste0(current_dir_out, '/', prefix, '_gene_coords.txt')


head(ldsc_out)
write_delim(ldsc_out, out_path, delim='\t')
#
out_path = paste0(current_finedir_out, '/', fine_prefix, '_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')










# ABA MTG
# ABA MTG
# ABA MTG
# ABA MTG
# ABA MTG
# ABA MTG


load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_aba_mtg_ctd.rda'), verbose=T)
mtg_ctd = ctd
mtg_ctd = prepare.quantile.groups(mtg_ctd, specificity_species="human", numberOfBins=10)

# 37 gene mapping
ensembl38 = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# easy to map genes
hgnc_ens38_df = getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position'), mart = ensembl38)
hgnc_ens38_df = hgnc_ens38_df[which(hgnc_ens38_df$chromosome_name %in% 1:22),]
hgnc_ens38_df = hgnc_ens38_df[!is.na(hgnc_ens38_df$entrezgene_id),]
hgnc_ens38_df = hgnc_ens38_df[!duplicated(hgnc_ens38_df$hgnc_symbol),]


# load preprocessed data
load(paste0(umi_dir, '/ahba_mtg/mtg_singlecell_ahba_match_seurat_processed.Rdata'), verbose=T)
#mtg_genes = colnames(mtg_dat)[1:colnames(mtg_dat)-1]


# current out dir
prefix = 'aba_mtg_bigcat_ldsc'
fine_prefix = 'aba_mtg_finecat_ldsc'

# make if doesnt exist
current_dir_out = paste0(ldsc_dir, '/', prefix, '_files')
current_finedir_out = paste0(ldsc_dir, '/', fine_prefix, '_files')


ldsc_out = hgnc_ens38_df[c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')]
colnames(ldsc_out) = c('GENE', 'CHR', 'START', 'END')
ldsc_out = ldsc_out[ldsc_out$CHR %in% as.character(1:22),]
ldsc_out = ldsc_out[order(ldsc_out$CHR),]
ldsc_out = ldsc_out[order(as.numeric(ldsc_out$CHR)),]
out_path = paste0(current_dir_out, '/',prefix,'_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')
print(out_path)
#
out_path = paste0(current_finedir_out, '/',fine_prefix,'_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')


# genes matched between hgnc coordinate data from and specificity data frame
hgnc_ens38_df_mtg = hgnc_ens38_df[hgnc_ens38_df$hgnc_symbol %in% rownames(mtg_ctd[[1]]$specificity_quantiles),]
hgnc_ens38_df_non_mtg = hgnc_ens38_df[!hgnc_ens38_df$hgnc_symbol %in% rownames(mtg_ctd[[1]]$specificity_quantiles),]


# subset specificity data
mtg_ctd[[1]]$specificity = mtg_ctd[[1]]$specificity[rownames(mtg_ctd[[1]]$specificity) %in% hgnc_ens38_df_mtg$hgnc_symbol,]
mtg_ctd[[2]]$specificity = mtg_ctd[[2]]$specificity[rownames(mtg_ctd[[2]]$specificity) %in% hgnc_ens38_df_mtg$hgnc_symbol,]
#mtg_ctd = prepare.quantile.groups(mtg_ctd, specificity_species="human", numberOfBins=20)
mtg_ctd = prepare.quantile.groups(mtg_ctd, specificity_species="human", numberOfBins=10)


nbins = 10

genebin_num = round(nrow(hgnc_ens38_df)/nbins)
mtg_cells   = colnames(mtg_ctd[[1]]$specificity_quantiles)
quantile_df = as.data.frame(mtg_ctd[[1]]$specificity_quantiles)

for (cell in mtg_cells){
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens38_df_mtg$ensembl_gene_id[hgnc_ens38_df_mtg$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_finedir_out, '/', fine_prefix, '_', cell, '_', str_pad(ct,2,pad="0"),'_',as.character(nbins),'_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens38_df_non_mtg$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_nomatch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens38_df_mtg$ensembl_gene_id, paste0(current_finedir_out, '/', fine_prefix, '_allcontrol_genes.txt'), quote=F, row.names=F, col.names=F)




nbins = 10
genebin_num = round(nrow(hgnc_ens38_df)/nbins)
mtg_cells   = colnames(mtg_ctd[[2]]$specificity_quantiles)
quantile_df = as.data.frame(mtg_ctd[[2]]$specificity_quantiles)

for (cell in mtg_cells){
    cur_quantile = quantile_df[[cell]]
    for (bin in 0:nbins){
        ct = bin
        bin_hgnc = rownames(quantile_df)[which(cur_quantile == bin)]
        bin_ensembls = hgnc_ens38_df_mtg$ensembl_gene_id[hgnc_ens38_df_mtg$hgnc_symbol %in% bin_hgnc]
        print(length(bin_ensembls))

        out_file = paste0(current_dir_out, '/', prefix, '_', cell, '_bigcat_', str_pad(ct,2,pad="0"),'_', as.character(nbins),'_genes.txt')
        write.table(bin_ensembls, out_file, quote=F, row.names=F, col.names=F)
    }
}
write.table(hgnc_ens38_df_non_mtg$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_nomatch_genes.txt'), quote=F, row.names=F, col.names=F)
write.table(hgnc_ens38_df_mtg$ensembl_gene_id, paste0(current_dir_out, '/', prefix, '_allcontrol_genes.txt'), quote=F, row.names=F, col.names=F)



