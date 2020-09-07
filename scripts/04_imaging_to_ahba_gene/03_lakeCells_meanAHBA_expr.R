library(tidyverse)
library(ggcorrplot)
library(Cairo)

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
    dscalar_template = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/tmp/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotVolOnSurface'
    cmd     = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util\');',
                    matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                    write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}


# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'


# read gene-wise spatial correlations to MDD effect maps
all_ahba_cors = read_csv(paste0(base_dir, "/output/mdd_regr/all_mdd_ahba_genecors.csv"))


# reverse direction of RSFA to line up with other modalities
all_ahba_cors$ukbb_rsfa_cor_new = all_ahba_cors$ukbb_rsfa_cor*-1
all_ahba_cors$ukbb_rsfa_cor     = NULL
all_ahba_cors$gsp_rsfa_cor_new  = all_ahba_cors$gsp_rsfa_cor*-1
all_ahba_cors$gsp_rsfa_cor      = NULL


# AHBA expression data
ahba_in       = paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv')
ahba_expr_200 = read_csv(ahba_in)
parc_labels   = colnames(ahba_expr_200)[grepl('17Network', colnames(ahba_expr_200))]
ahba_expr_200 = ahba_expr_200[,c(order(parc_labels),201)]


# fix up the roi names for alingment to regressino data
ahba_expr_200_t           = as.data.frame(t(ahba_expr_200[,1:200]))
colnames(ahba_expr_200_t) = ahba_expr_200$gene
ahba_expr_200_t$roi       = rownames(ahba_expr_200_t)
ahba_expr_200_t$roi       = gsub('17Networks', 'net17', ahba_expr_200_t$roi)


# UKBB Thick
thick_file = paste0(base_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
ukbb_thick = read_csv(thick_file)

# UKBB RSFA
rsfa_file = paste0(base_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
ukbb_rsfa = read_csv(rsfa_file)
ukbb_rsfa$roi = gsub('_scale','',ukbb_rsfa$roi)

# UKBB GBC
gbc_file = paste0(base_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
ukbb_gbc = read_csv(gbc_file)

# Harvard RSFA
gsp_rsfa_file = paste0(base_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv')
gsp_rsfa = read_csv(gsp_rsfa_file)
gsp_rsfa$roi = gsub('_scale','',gsp_rsfa$roi)

# Harvard GBC
gsp_gbc_file = paste0(base_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv')
gsp_gbc = read_csv(gsp_gbc_file)


# format and combine
# ----------
thick_tmp = ukbb_thick[c('roi','cohens_d')]
colnames(thick_tmp) = c('roi','ukbb_thick_cohens_d')

ukbb_rsfa_tmp = ukbb_rsfa[c('roi','cohens_d')]
colnames(ukbb_rsfa_tmp) = c('roi','ukbb_rsfa_cohens_d')

ukbb_gbc_tmp = ukbb_gbc[c('roi','cohens_d')]
colnames(ukbb_gbc_tmp) = c('roi','ukbb_gbc_cohens_d')

gsp_rsfa_tmp = gsp_rsfa[c('roi','scaleneg_affect_mean_cohens_d')]
colnames(gsp_rsfa_tmp) = c('roi','gsp_rsfa_cohens_d')

gsp_gbc_tmp = gsp_gbc[c('roi','scaleneg_affect_mean_cohens_d')]
colnames(gsp_gbc_tmp) = c('roi','gsp_gbc_cohens_d')

mdd_pheno_df = merge(x=thick_tmp, y=ukbb_rsfa_tmp, by='roi')
mdd_pheno_df = merge(x=mdd_pheno_df, y=ukbb_gbc_tmp, by='roi')
mdd_pheno_df = merge(x=mdd_pheno_df, y=gsp_rsfa_tmp, by='roi')
mdd_pheno_df = merge(x=mdd_pheno_df, y=gsp_gbc_tmp, by='roi')


# --------
# ENIGMA
# --------
# read ENIGMA data and format desikan names for matching
enigma_mdd    = read_csv(paste0(base_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv'))

# gene by desikan ROI matrix
desikan_mat   = read_csv(paste0(base_dir, '/data/ahba_parcel/desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))
desikan_mat_t = as.data.frame(t(desikan_mat[,1:70]))
colnames(desikan_mat_t) = desikan_mat$gene
desikan_mat_t$roi = as.character(rownames(desikan_mat_t))
desikan_mat_t[1:5,1:5]




# load differential expression data
load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_frontal_diff_expr_batchcorrect.Rdata'))
load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_visual_diff_expr_batchcorrect.Rdata'))

# prepare single cell data
both_lake_cells = intersect(names(dfc_diff_expr_list), names(vis_diff_expr_list)) #names(dfc_diff_expr_list)
both_lake_sets  = list()
for (cell in both_lake_cells){
    # dfc
    dfc_dexpr_dat  = dfc_diff_expr_list[[cell]]
    dfc_dexpr_dat  = dfc_dexpr_dat[which(dfc_dexpr_dat$avg_logFC > 0),]
    dfc_dexpr_dat  = dfc_dexpr_dat[which(dfc_dexpr_dat$p_val_adj < .05),]

    # vis
    vis_dexpr_dat  = vis_diff_expr_list[[cell]]
    vis_dexpr_dat  = vis_dexpr_dat[which(vis_dexpr_dat$avg_logFC > 0),]
    vis_dexpr_dat  = vis_dexpr_dat[which(vis_dexpr_dat$p_val_adj < .05),]


    cell_genes = intersect(rownames(dfc_dexpr_dat), rownames(vis_dexpr_dat))
    print(cell)
    print(length(cell_genes))
    both_lake_sets[[cell]] = cell_genes
}


# get the average spatial expression of cell specific genes
# then correlate to each MDD imaging map
cell_cor_df = NULL
lake_cells  = names(both_lake_sets)
for (cell in lake_cells){

    lake_cell = both_lake_sets[[cell]]
    lake_cell = lake_cell[lake_cell %in% colnames(ahba_expr_200_t)]
    avg_expr  = rowMeans(ahba_expr_200_t[lake_cell])

    avg_expr_df = data.frame(roi=gsub('17Networks', 'net17', names(avg_expr)), expr=avg_expr)

    combined_df = merge(x=mdd_pheno_df, y=avg_expr_df, by='roi')
    cell_cors   = cor(combined_df$expr, combined_df[grep('cohens_d', colnames(combined_df))], use='pairwise.complete', method='spearman')
    cell_cors   = as.data.frame(cell_cors)
    cell_cors$cell = cell

    #desikan_weighted_expr  = apply(desikan_mat_t[lake_cell$Gene], 1, function(x) weighted.mean(x, weights))
    desikan_avg_expr  = rowMeans(desikan_mat_t[lake_cell])
    desikan_avg_expr_df = data.frame(roi=names(desikan_avg_expr), expr=desikan_avg_expr)

    enigma_combined_df = merge(x=enigma_mdd, y=desikan_avg_expr_df, by='roi')
    enigma_cor         = cor(enigma_combined_df$cohens_d, enigma_combined_df$expr, use='pairwise.complete', method='spearman')
    cell_cors$enigma_cohens_d = enigma_cor
    cell_cors$num_gene = nrow(lake_cell)
    cell_cor_df        = rbind(cell_cor_df, cell_cors)
}
cell_cor_df$ukbb_rsfa_cohens_d = cell_cor_df$ukbb_rsfa_cohens_d*-1
cell_cor_df$gsp_rsfa_cohens_d  = cell_cor_df$gsp_rsfa_cohens_d*-1

cell_cor_df$avg = rowMeans(cell_cor_df[grep('cohens_d', colnames(cell_cor_df))])
cell_cor_df = cell_cor_df[order(cell_cor_df$avg),]
head(cell_cor_df,10)

out_path = paste0(base_dir, '/output/tables/lake_diffexpr_weighted_mean_SUPERORDINATE.csv')
write_csv(x=cell_cor_df, path=out_path)




