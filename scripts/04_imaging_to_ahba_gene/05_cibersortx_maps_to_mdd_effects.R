library(tidyverse)
library(cifti)
library(Cairo)

plot_cell_cors = function(dfc_df_in, vis_df_in, pdf_out, ylim){

    both_cell_thick_cors = rbind(dfc_df_in, vis_df_in)

    plot_order = both_cell_thick_cors %>% group_by(cells) %>% summarize(m=mean(cor, na.rm=T))
    cell_order = as.character(plot_order$cells[order(plot_order$m)])
    both_cell_thick_cors$cells = factor(as.character(both_cell_thick_cors$cells), levels=cell_order)


    thick_plot_out = paste0(project_dir, '/figures/', pdf_out)
    thick_plot_out
    thick_cell_plot = ggplot(data=both_cell_thick_cors, aes(x=cells, y=cor, fill=cell_ref)) +
        geom_col(colour="black",width=0.7,
               position=position_dodge(0.7)) +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0), limits=c(ylim[1],ylim[2]), breaks = seq(ylim[1],ylim[2], by = ylim[3])) +
        scale_fill_manual(values=c('#6E727F','#BEC1CC')) +
                theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              plot.title = element_text(hjust = 0.5, size=12),
              axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black", size=8, angle = 90, hjust = 1),
              axis.text.y = element_text(colour="black", size=8),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),legend.position = "none")
    ggsave(thick_plot_out, plot=thick_cell_plot, height=1.25, width=2.50)
    print(thick_cell_plot)
    #dev.off()
    thick_plot_out
}

corr_eff_to_cells = function(df, cell_df, cell_names, region, csv_save){

    # combine effects with cells
    mdd_cell_df = merge(x=cell_df, y=df, by='roi')

    # correlate
    mdd_cors = cor(mdd_cell_df[cell_names], mdd_cell_df$cohens_d, use='pairwise.complete', method='spearman')
    mdd_df   = data.frame(cor=mdd_cors, cells=gsub(paste0(region, '_'),'',rownames(mdd_cors)))
    mdd_df$cell_ref = region
    mdd_df[order(mdd_df$cor),]

    write_csv(mdd_df, paste0(project_dir, '/output/mdd_regr/', csv_save))

    mdd_df = mdd_df[order(mdd_df$cor),]
    return(mdd_df)
}


# set up directories
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read Lake data averaged by parcel
nparcel  = '200'
net_num  = '17'
dfc_cell = read_csv(paste0(project_dir, '/data/ahba_cibersortx/schaeffer_LAKE_DFC_200_17Net_expr_matdfc_abs_orig_donorscaled_superordinate.csv'))

# transpose gene by parcel matrix
dfc_cell_t = as.data.frame(t(dfc_cell[,1:as.character(nparcel)]))
colnames(dfc_cell_t) = dfc_cell$gene
dfc_cell_t$roi = rownames(dfc_cell_t)
dfc_cells  = dfc_cell$gene
dfc_cell_t[1:5,1:5]


# Lake visual cortex
vis_cell = read_csv(paste0(project_dir, '/data/ahba_cibersortx/schaeffer_LAKE_VIS_200_17Net_expr_matvis_abs_orig_donorscaled_superordinate.csv'))

# transpose
vis_cell_t = as.data.frame(t(vis_cell[,1:as.character(nparcel)]))
colnames(vis_cell_t) = vis_cell$gene
vis_cell_t$roi = rownames(vis_cell_t)
vis_cells  = vis_cell$gene
vis_cell_t[1:5,1:5]


# cells that are present in both VIS/DFC
dfc_tmp   = gsub('DFC_', '', dfc_cells)
vis_tmp   = gsub('VIS_', '', vis_cells)
cells_use = intersect(dfc_tmp, vis_tmp)


# put the common cells back into DFC/VIS specific naming
dfc_cells = dfc_cells[dfc_cells %in% paste0('DFC_', cells_use)]
vis_cells = vis_cells[vis_cells %in% paste0('VIS_', cells_use)]


# combine both VIS/DFC data
cell_both_df = merge(x=dfc_cell_t, y=vis_cell_t, by='roi')
cell_both_df$roi = gsub('17Networks', 'net17', cell_both_df$roi)
cell_both_df[1:5,1:5]


# read parcel-wise regression estimates
# ---
# UKBB cortical thickness
mdd_thick        = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table.csv'))
dfc_ukb_thick_df = corr_eff_to_cells(df=mdd_thick, cell_df=cell_both_df, cell_names=dfc_cells, region='DFC', csv_save='ukbb_thick_to_lakeDFC_origlabels_cell.csv')
vis_ukb_thick_df = corr_eff_to_cells(df=mdd_thick, cell_df=cell_both_df, cell_names=vis_cells, region='VIS', csv_save='ukbb_thick_to_lakeVIS_origlabels_cell.csv')

# UKBB RSFA
mdd_rsfa     = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table.csv'))
mdd_rsfa$roi = gsub('_scale','',mdd_rsfa$roi)
dfc_ukb_rsfa_df = corr_eff_to_cells(df=mdd_rsfa, cell_df=cell_both_df, cell_names=dfc_cells, region='DFC', csv_save='ukbb_rsfa_to_lakeDFC_origlabels_cell.csv')
vis_ukb_rsfa_df = corr_eff_to_cells(df=mdd_rsfa, cell_df=cell_both_df, cell_names=vis_cells, region='VIS', csv_save='ukbb_rsfa_to_lakeVIS_origlabels_cell.csv')

# UKBB GBC
mdd_gbc   = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table.csv'))
dfc_ukb_gbc_df = corr_eff_to_cells(df=mdd_gbc, cell_df=cell_both_df, cell_names=dfc_cells, region='DFC', csv_save='ukbb_gbc_to_lakeDFC_origlabels_cell.csv')
vis_ukb_gbc_df = corr_eff_to_cells(df=mdd_gbc, cell_df=cell_both_df, cell_names=vis_cells, region='VIS', csv_save='ukbb_gbc_to_lakeVIS_origlabels_cell.csv')

# Harvard GSP RSFA
gsp_rsfa     = read_csv(paste0(project_dir, '/output/mdd_regr/gsp_rsfa_regr_table.csv'))
gsp_rsfa$cohens_d = gsp_rsfa$scaleneg_affect_mean_cohens_d
gsp_rsfa$roi = gsub('_scale', '', gsp_rsfa$roi)
dfc_gsp_rsfa_df = corr_eff_to_cells(df=gsp_rsfa, cell_df=cell_both_df, cell_names=dfc_cells, region='DFC', csv_save='gsp_rsfa_to_lakeDFC_origcell.csv')
vis_gsp_rsfa_df = corr_eff_to_cells(df=gsp_rsfa, cell_df=cell_both_df, cell_names=vis_cells, region='VIS', csv_save='gsp_rsfa_to_lakeDFC_origcell.csv')

# Harvard GSP GBC
gsp_gbc   = read_csv(paste0(project_dir, '/output/mdd_regr/gsp_gbc_regr_table.csv'))
gsp_gbc$cohens_d = gsp_gbc$scaleneg_affect_mean_cohens_d
dfc_gsp_gbc_df = corr_eff_to_cells(df=gsp_gbc, cell_df=cell_both_df, cell_names=dfc_cells, region='DFC', csv_save='gsp_gbc_to_lakeDFC_origcell.csv')
vis_gsp_gbc_df = corr_eff_to_cells(df=gsp_gbc, cell_df=cell_both_df, cell_names=vis_cells, region='VIS', csv_save='gsp_gbc_to_lakeDFC_origcell.csv')



# read ENIGMA data and format desikan names for matching
enigma_out = paste0(project_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv')
enigma_mdd = read_csv(enigma_out)


# read imputed cell density by desikan atlas
desikan_mat = read_csv(paste0(project_dir, '/data/ahba_cibersortx/desikan_lake_origlabels_absolute_superordinate_cell_mat.csv'))

# transpose data
desikan_mat_t = as.data.frame(t(desikan_mat[,1:70]))
colnames(desikan_mat_t) = desikan_mat$gene
desikan_mat_t$roi = as.character(rownames(desikan_mat_t))
desikan_mat_t[1:5,1:5]


# combine cibersortx with ENIGMA cohen's d
enigma_lake_df = merge(x=desikan_mat_t, y=enigma_mdd, by='roi')

# DFC; correlate cell types to ENIGMA cohen's d
enigma_dfc_mdd_cors = cor(enigma_lake_df[dfc_cells], enigma_lake_df$cohens_d, use='pairwise.complete', method='spearman')
enigma_dfc_mdd_cors = data.frame(cor=enigma_dfc_mdd_cors, cells=gsub('DFC_','',rownames(enigma_dfc_mdd_cors)))
enigma_dfc_mdd_cors$cell_ref = 'DFC'
enigma_dfc_mdd_cors = enigma_dfc_mdd_cors[order(enigma_dfc_mdd_cors$cor),]

# VIS; correlate cell types to ENIGMA cohen's d
enigma_vis_mdd_cors = cor(enigma_lake_df[vis_cells], enigma_lake_df$cohens_d, use='pairwise.complete', method='spearman')
enigma_vis_mdd_cors = data.frame(cor=enigma_vis_mdd_cors, cells=gsub('VIS_','',rownames(enigma_vis_mdd_cors)))
enigma_vis_mdd_cors$cell_ref = 'VIS'
enigma_vis_mdd_cors = as.data.frame(enigma_vis_mdd_cors[order(enigma_vis_mdd_cors$cor),])



# make some sanity check plots not used for publication
plot_cell_cors(dfc_df_in=dfc_ukb_thick_df, vis_df_in=vis_ukb_thick_df, pdf_out='ukbb_mdd_thick_lakecell_cors_frac.pdf', ylim=c(-.4,.4, .2))
plot_cell_cors(dfc_df_in=dfc_ukb_rsfa_df, vis_df_in=vis_ukb_rsfa_df, pdf_out='ukbb_mdd_rsfa_lakecell_cors_frac.pdf', ylim=c(-.5,.5,.25))
plot_cell_cors(dfc_df_in=dfc_ukb_gbc_df, vis_df_in=vis_ukb_gbc_df, pdf_out='ukbb_mdd_gbc_lakecell_cors_frac.pdf', ylim=c(-.4,.4,.2))
plot_cell_cors(dfc_df_in=dfc_gsp_rsfa_df, vis_df_in=vis_gsp_rsfa_df, pdf_out='gsp_mdd_rsfa_lakecell_cors_fracs.pdf', ylim=c(-.4,.4,.2))
plot_cell_cors(dfc_df_in=dfc_gsp_gbc_df, vis_df_in=vis_gsp_gbc_df, pdf_out='gsp_mdd_gbc_lakecell_cors_frac.pdf', ylim=c(-.5,.5,.25))
plot_cell_cors(dfc_df_in=enigma_dfc_mdd_cors, vis_df_in=enigma_vis_mdd_cors, pdf_out='enigma_mdd_thick_lakecell_cors_frac.pdf', ylim=c(-.7,.7, .35))



# combine VIS/DFC cell correlations
ukbb_mdd_thick = merge(dfc_ukb_thick_df, vis_ukb_thick_df, by='cells', all=T)
ukbb_mdd_rsfa  = merge(dfc_ukb_rsfa_df, vis_ukb_rsfa_df, by='cells', all=T)
ukbb_mdd_gbc   = merge(dfc_ukb_gbc_df, vis_ukb_gbc_df, by='cells', all=T)
gsp_mdd_rsfa   = merge(dfc_gsp_rsfa_df, vis_gsp_rsfa_df, by='cells', all=T)
gsp_mdd_gbc    = merge(dfc_gsp_gbc_df, vis_gsp_gbc_df, by='cells', all=T)
enigma_thick   = merge(enigma_dfc_mdd_cors, enigma_vis_mdd_cors, by='cells', all=T)


# average VIS/DFC correlations
ukbb_mdd_thick$avg_cor = rowMeans(cbind(ukbb_mdd_thick$cor.x, ukbb_mdd_thick$cor.y), na.rm=T)
ukbb_mdd_thick$pheno   = 'ukbb_thick'
ukbb_mdd_rsfa$avg_cor  = rowMeans(cbind(ukbb_mdd_rsfa$cor.x, ukbb_mdd_rsfa$cor.y), na.rm=T)*-1
ukbb_mdd_rsfa$pheno    = 'ukbb_rsfa'
ukbb_mdd_gbc$avg_cor   = rowMeans(cbind(ukbb_mdd_gbc$cor.x, ukbb_mdd_gbc$cor.y), na.rm=T)
ukbb_mdd_gbc$pheno     = 'ukbb_gbc'
gsp_mdd_rsfa$avg_cor   = rowMeans(cbind(gsp_mdd_rsfa$cor.x, gsp_mdd_rsfa$cor.y), na.rm=T)*-1
gsp_mdd_rsfa$pheno     = 'gsp_rsfa'
gsp_mdd_gbc$avg_cor    = rowMeans(cbind(gsp_mdd_gbc$cor.x, gsp_mdd_gbc$cor.y), na.rm=T)
gsp_mdd_gbc$pheno      = 'gsp_gbc'
enigma_thick$avg_cor   = rowMeans(cbind(enigma_thick$cor.x, enigma_thick$cor.y), na.rm=T)
enigma_thick$pheno     = 'enigma_thick'


# Combine all the cibersortX by MDD phenotype datasets
all_cibersort_cors_df = rbind(ukbb_mdd_thick, ukbb_mdd_rsfa, ukbb_mdd_gbc, gsp_mdd_rsfa, gsp_mdd_gbc, enigma_thick)

cor(cbind(ukbb_mdd_thick['cor.x'], -1*ukbb_mdd_rsfa['cor.x'], ukbb_mdd_gbc['cor.x'], -1*gsp_mdd_rsfa['cor.x'], gsp_mdd_gbc['cor.x'], enigma_thick['cor.x']))


# compute mean cell-type correlation to MDD phenotypes
plot_order = all_cibersort_cors_df %>% group_by(cells) %>% summarize(fullcor=mean(avg_cor))

# template dataframe to plug in averaged cell type correlations
avg_cor_df = ukbb_mdd_thick
avg_cor_df = merge(avg_cor_df, plot_order, by='cells')
avg_cor_df$avg_cor = avg_cor_df$fullcor
avg_cor_df$pheno   = 'average'
plot_order = plot_order[rev(order(plot_order$fullcor)),]

# combine the information across all modalities
all_cibersort_cors_df = rbind(all_cibersort_cors_df, avg_cor_df[colnames(all_cibersort_cors_df)])

# plot order
all_cibersort_cors_df$cells = factor(as.character(all_cibersort_cors_df$cells), levels=as.character(plot_order$cells))
all_cibersort_cors_df$pheno = factor(all_cibersort_cors_df$pheno, levels=c('ukbb_thick','enigma_thick','ukbb_rsfa','gsp_rsfa','ukbb_gbc','gsp_gbc','average'))

summary(all_cibersort_cors_df$avg_cor)
p = ggplot(data = all_cibersort_cors_df, aes(x=cells, y=pheno, fill=avg_cor)) +
            geom_tile() +
            coord_flip() +
            scale_fill_gradient2(low = "#4B5DAA", high = "#BD242F", mid = "white", midpoint = 0, limit = c(-.6,.6))+
            theme_classic() +
            scale_y_discrete(expand=c(0,0)) +
            scale_x_discrete(expand=c(0,0))
p

figure_out = paste0(project_dir, '/figures/lake_cibersortX_to_neuroimage_effects.pdf')
figure_out
ggsave(figure_out, plot=p, height=2.5, width=3.1)



