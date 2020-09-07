library(tidyverse)
library(Cairo)
library(RColorBrewer)


# directory
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read MDD regression tables
thick_file = paste0(base_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
ukbb_thick = read_csv(thick_file)

rsfa_file = paste0(base_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
ukbb_rsfa = read_csv(rsfa_file)
ukbb_rsfa$roi = gsub('_scale','',ukbb_rsfa$roi)

gbc_file = paste0(base_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
ukbb_gbc = read_csv(gbc_file)

gsp_rsfa_file = paste0(base_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv')
gsp_rsfa = read_csv(gsp_rsfa_file)

gsp_gbc_file = paste0(base_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv')
gsp_gbc = read_csv(gsp_gbc_file)
ukbb_rsfa$roi = gsub('_scale','',ukbb_rsfa$roi)



# Schaeffer parcellation
#schaef_dir    = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
#parc_label_in = read_csv(paste0(schaef_dir, '/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
#parc_labels   = parc_label_in$X1[grep('17Net', parc_label_in$X1)]
#parc_labels   = gsub('17Networks','net17',parc_labels)

info_file   = paste0(project_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt')
new_mapping = read.csv(header=F, info_file)
new_map_df  = data.frame(new_net=as.character(new_mapping$V1[grep('17Networks', new_mapping$V1)]),
                        new_info=as.character(new_mapping$V1[!grepl('17Networks', new_mapping$V1)]), stringsAsFactors=F)

parc_labels = gsub('17Networks','net17', new_map_df$new_net)


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


# merge MDD regression DF with AHBA expression
thick_ahba_df = merge(x=ukbb_thick, by.x='roi', y=ahba_expr_200_t, by.y='roi')


sst_marker_genes=c('SST', 'CORT', 'NPY')

# Correlate MDD THICKNESS with all AHBA genes
gene_names      = colnames(ahba_expr_200_t)[1:ncol(ahba_expr_200_t)-1]
thick_gene_cors = t(cor(thick_ahba_df$cohens_d, thick_ahba_df[gene_names], use='pairwise.complete', method='spearman'))


# Correlate MDD RSFA with all AHBA genes
rsfa_ahba_df = merge(x=ukbb_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
rsfa_gene_cors = t(cor(rsfa_ahba_df$cohens_d, rsfa_ahba_df[gene_names], use='pairwise.complete', method='spearman'))


# Correlate MDD GBC with all AHBA genes
gbc_ahba_df = merge(x=ukbb_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')
gbc_gene_cors = t(cor(gbc_ahba_df$cohens_d, gbc_ahba_df[gene_names], use='pairwise.complete', method='spearman'))



# SST expression to MDD cohen's D
summary(thick_ahba_df$cohens_d)
summary(thick_ahba_df$SST)
thick_ahba_df$sst_markers = rowMeans(thick_ahba_df[which(names(thick_ahba_df) %in% sst_marker_genes)])
summary(thick_ahba_df$sst_markers)

out_path = paste0(base_dir, '/figures/ukbb_thick_mdd_cohensD_corr_to_uniform_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')

ggplot(data=thick_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
        geom_point(size=1.5,stroke=.5) +
        geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
        geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
        geom_smooth(method='lm', linetype = 5, color='black', se=F) +
        scale_shape_manual(values=c(16, 6, 0)) +
        theme_classic() +
        theme(legend.position = 'none',
            text=element_text(size=16,  family='Helvetica', colour='black')) +
        scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
        scale_x_continuous(limits=c(-.05,.05), breaks=seq(-.05,.05,.025), expand=c(0,0)) +
        scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
        theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
                axis.ticks = element_line(colour = 'black'))
dev.off()
out_path


summary(rsfa_ahba_df$cohens_d)
summary(rsfa_ahba_df$SST)

out_path = paste0(base_dir, '/figures/ukbb_rsfa_mdd_cohensD_corr_to_uniform_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')
ggplot(rsfa_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
    geom_point(size=1.5,stroke=.5) +
    geom_smooth(method='lm', linetype = 5, color='black', se=F) +
    geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
    geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
    scale_shape_manual(values=c(16, 6, 0)) +
    theme_classic() +
    theme(legend.position = 'none',
        text=element_text(size=16,  family='Helvetica', colour='black')) +
    scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
    scale_x_continuous(limits=c(-.07,.07), breaks=seq(-.07,.07,.035), expand=c(0,0)) +
    scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'))
dev.off()
out_path



summary(gbc_ahba_df$cohens_d)
summary(gbc_ahba_df$SST)

out_path = paste0(base_dir, '/figures/ukbb_gbc_mdd_cohensD_corr_to_uniform_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')
ggplot(gbc_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
    geom_point(size=1.5,stroke=.5) +
    geom_smooth(method='lm', linetype = 5, color='black', se=F) +
    geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
    geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
    scale_shape_manual(values=c(16, 6, 0)) +
    theme_classic() +
    theme(legend.position = 'none',
        text=element_text(size=16,  family='Helvetica', colour='black')) +
    scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
    scale_x_continuous(limits=c(-.07,.07), breaks=seq(-.07,.07,.035), expand=c(0,0)) +
    scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'))
dev.off()
out_path





# Correlate
gsp_gbc_ahba_df = merge(x=gsp_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')
gsp_gbc_ahba_df$cohens_d = gsp_gbc_ahba_df$scaleneg_affect_mean_cohens_d
gsp_gbc_gene_cors = t(cor(gsp_gbc_ahba_df$cohens_d, gsp_gbc_ahba_df[gene_names], use='pairwise.complete', method='spearman'))

out_path = paste0(base_dir, '/figures/gsp_gbc_negaffect_COHENS_D_corr_to_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')
ggplot(gsp_gbc_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
    geom_point(size=1.5,stroke=.5) +
    geom_smooth(method='lm', linetype = 5, color='black', se=F) +
    geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
    geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
    scale_shape_manual(values=c(16, 6, 0)) +
    theme_classic() +
    theme(legend.position = 'none',
        text=element_text(size=16,  family='Helvetica', colour='black')) +
    scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
    scale_x_continuous(limits=c(-.20,.20), breaks=seq(-.20,.20,.1), expand=c(0,0)) +
    scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'))
dev.off()
out_path






# Correlate
gsp_rsfa$roi = gsub('_scale','', gsp_rsfa$roi)
gsp_rsfa_ahba_df = merge(x=gsp_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
gsp_rsfa_ahba_df$cohens_d = gsp_rsfa_ahba_df$scaleneg_affect_mean_cohens_d
gsp_rsfa_gene_cors = t(cor(gsp_rsfa_ahba_df$cohens_d, gsp_rsfa_ahba_df[gene_names], use='pairwise.complete', method='spearman'))

out_path = paste0(base_dir, '/figures/gsp_rsfa_negaffect_COHENS_D_corr_to_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')
ggplot(gsp_rsfa_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
    geom_point(size=1.5,stroke=.5) +
    geom_smooth(method='lm', linetype = 5, color='black', se=F) +
    geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
    geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
    scale_shape_manual(values=c(16, 6, 0)) +
    theme_classic() +
    theme(legend.position = 'none',
        text=element_text(size=16,  family='Helvetica', colour='black')) +
    scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
    scale_x_continuous(limits=c(-.20,.20), breaks=seq(-.20,.20,.1), expand=c(0,0)) +
    scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'))
dev.off()
out_path






# read ENIGMA data and format desikan names for matching
enigma_out = paste0(base_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv')
enigma_mdd = read_csv(enigma_out)


desikan_mat   = read_csv(paste0(base_dir, '/data/ahba_parcel/desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))
desikan_mat_t = as.data.frame(t(desikan_mat[,1:70]))
colnames(desikan_mat_t) = desikan_mat$gene
desikan_mat_t$roi = as.character(rownames(desikan_mat_t))

enigma_ahba_df = merge(x=desikan_mat_t, y=enigma_mdd, by='roi')

enigma_ahba_mdd_cors = cor(enigma_ahba_df[desikan_mat$gene], enigma_ahba_df$cohens_d, use='pairwise.complete', method='spearman')
enigma_ahba_mdd_cors = data.frame(cor=enigma_ahba_mdd_cors, genes=rownames(enigma_ahba_mdd_cors))
enigma_ahba_mdd_cors = enigma_ahba_mdd_cors[order(enigma_ahba_mdd_cors$cor),]




enigma_mdd_ahba_df = merge(x=enigma_mdd, by.x='roi', y=desikan_mat_t, by.y='roi')
enigma_mdd_ahba_df$cohens_d = enigma_mdd_ahba_df$cohens_d
enigma_gene_cors = t(cor(enigma_mdd_ahba_df$cohens_d, enigma_mdd_ahba_df[gene_names], use='pairwise.complete', method='spearman'))

out_path = paste0(base_dir, '/figures/enigma_thick_COHENS_D_corr_to_SST_markers.pdf')
out_path
CairoPDF(out_path, height=3, width=3, family='Helvetica')
ggplot(enigma_mdd_ahba_df, aes(x=cohens_d, y=SST, color=SST)) +
    geom_point(size=1.5,stroke=.5) +
    geom_smooth(method='lm', linetype = 5, color='black', se=F) +
    geom_smooth(aes(x=cohens_d, y=NPY), color='blue', method='lm', linetype = 5, se=F) +
    geom_smooth(aes(x=cohens_d, y=CORT), color='green', method='lm', linetype = 5, se=F) +
    scale_shape_manual(values=c(16, 6, 0)) +
    theme_classic() +
    theme(legend.position = 'none',
        text=element_text(size=16,  family='Helvetica', colour='black')) +
    scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2.5,2.5,1.25), expand=c(0,0)) +
    scale_x_continuous(limits=c(-.18,.18), breaks=seq(-.18,.18,.09), expand=c(0,0)) +
    scale_color_gradient(low = "#FEEFB7", high = "#BF423C", limits=c(-3,2)) +
    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'))
dev.off()
out_path














