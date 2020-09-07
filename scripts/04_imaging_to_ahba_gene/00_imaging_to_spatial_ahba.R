library(tidyverse)
library(Cairo)
library(RColorBrewer)


# distribution of gene correlations to MDD thickness
plot_genewise_density = function(cor_array, out_path, xlims, ylims, plot_gene){
    sst_val  = cor_array[ names(cor_array) == 'SST']
    mean_val = mean(cor_array)
    plot_df  = data.frame(plot_vals=as.numeric(cor_array), genes=names(cor_array))

    p = ggplot(plot_df, aes(x=plot_vals)) +
                geom_density(color="black", fill="lightgray", show.legend = FALSE) +
                geom_vline(aes(xintercept=sst_val)) +
                geom_vline(aes(xintercept=mean_val)) +
                theme_classic() +
                    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
                        axis.ticks = element_line(colour = 'black'),
                    axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = 'none', text=element_text(size=8,  family='Helvetica', colour='black')) +
                scale_x_continuous(limits=c(xlims[[1]],xlims[[2]]),breaks=seq(xlims[[1]],xlims[[2]],xlims[[3]]),expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0), limits=c(ylims[[1]],ylims[[2]]))
    ggsave(out_path, plot=p, height=.65, width=1.5)
    write(out_path,'')
}


# distribution of gene correlations to MDD thickness
plot_genewise_perm_density = function(plot_df, target_cor, out_path, xlims, ylims){
    mean_val = mean(plot_df$marker_perm)

    p = ggplot(plot_df, aes(x=marker_perm)) +
                geom_density(color="black", fill="lightgray", show.legend = FALSE) +
                geom_vline(aes(xintercept=target_cor)) +
                geom_vline(aes(xintercept=mean_val)) +
                theme_classic() +
                    theme(axis.text.x = element_text(colour = 'black'), axis.text.y = element_text(colour = 'black'),
                        axis.ticks = element_line(colour = 'black'),
                    axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = 'none', text=element_text(size=8,  family='Helvetica', colour='black')) +
                scale_x_continuous(limits=c(xlims[[1]],xlims[[2]]),breaks=seq(xlims[[1]],xlims[[2]],xlims[[3]]),expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0), limits=c(ylims[[1]],ylims[[2]]))
    ggsave(out_path, plot=p, height=.65, width=1.5)
    write(out_path,'')
    return(p)
}


marker_gene_permutation = function(gene_cors, genes, num_perm){

    target_val = mean(gene_cors$cor[gene_cors$genes %in% genes])

    perm_gene_arr = NULL
    perm_arr = NULL
    for (perm in 1:10000){
        if ((perm %% 1000 == 0)){
            print(perm)
        }

        random_genes = sample(1:length(gene_cors$genes),3)
        gene_3 = gene_cors$genes[random_genes]
        perm_gene_arr[[perm]] = gene_3
        random3  = mean(gene_cors$cor[random_genes])
        perm_arr = c(perm_arr, random3)
    }
    pval_perm_arr   = c(perm_arr, target_val)
    more_sig_combos = which(abs(pval_perm_arr) > abs(target_val))
    less_sig_combos = which(abs(pval_perm_arr) <= abs(target_val))

    sig_triplets  = perm_gene_arr[more_sig_combos]

    num_greater   = length(less_sig_combos)
    print(num_perm-num_greater)

    print(1-((num_greater+1)/length(pval_perm_arr)))
    return(list(perm_arr,sig_triplets))
}


cell_marker_gene_permutation = function(lake_dexpr, gene_cors, genes, num_perm){

    lake_dexpr_nosst = lake_dexpr[lake_dexpr$Cluster != 'In7|In8',]
    cell_types = unique(lake_dexpr_nosst$Cluster)
    target_val = mean(gene_cors$cor[gene_cors$genes %in% genes])

    all_perms = NULL
    for (cell in cell_types){
        write(cell,'')
        cell_dexpr = lake_dexpr_nosst[lake_dexpr_nosst$Cluster == cell,]
        combn      = t(combn(nrow(cell_dexpr), 3))
        print(dim(combn))
        cell_perms  = matrix(NA, 10000)
        random_idxs = sample(nrow(combn), 10000)
        row = 1
        for (rand_i in 1:length(random_idxs)){
            triplets = cell_dexpr$Gene[combn[rand_i,]]
            gene_3   = mean(gene_cors$cor[gene_cors$genes %in% triplets])
            cell_perms[row] = gene_3
            row = row + 1
        }
        all_perms = c(all_perms, cell_perms)
    }

    pval_perm_arr   = c(all_perms, target_val)
    more_sig_combos = which(abs(all_perms) > abs(target_val))
    less_sig_combos = which(abs(all_perms) <= abs(target_val))
    num_greater   = length(less_sig_combos)

    pval = 1-((num_greater+1)/length(pval_perm_arr))
    print(pval)
    return(list(pval, all_perms))
}


# Read Schaefer / AHBA / Regression tables
# --------------
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# Read Schaeffer parcellation
schaef_dir    = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
parc_label_in = read_csv(paste0(schaef_dir, '/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
parc_labels   = parc_label_in$X1[grep('17Net', parc_label_in$X1)]
parc_labels   = gsub('17Networks','net17',parc_labels)

# Read AHBA expression data
ahba_in       = paste0(project_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv')
ahba_expr_200 = read_csv(ahba_in)
parc_labels   = colnames(ahba_expr_200)[grepl('17Network', colnames(ahba_expr_200))]
ahba_expr_200 = ahba_expr_200[,c(order(parc_labels),201)]

# fix up the roi names for alignment to regression data
ahba_expr_200_t           = as.data.frame(t(ahba_expr_200[,1:200]))
colnames(ahba_expr_200_t) = ahba_expr_200$gene
ahba_expr_200_t$roi       = rownames(ahba_expr_200_t)
ahba_expr_200_t$roi       = gsub('17Networks', 'net17', ahba_expr_200_t$roi)


# "newmap" reflects an updated parcel labeling scheme from the Yeo Lab

# read MDD regression tables
# UKBB Thick
thick_file = paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
ukbb_thick = read_csv(thick_file)

# UKBB RSFA
rsfa_file     = paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
ukbb_rsfa     = read_csv(rsfa_file)
ukbb_rsfa$roi = gsub('_scale','',ukbb_rsfa$roi)

# UKBB GBC
gbc_file = paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
ukbb_gbc = read_csv(gbc_file)

# Harvard RSFA
gsp_rsfa_file = paste0(project_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv')
gsp_rsfa      = read_csv(gsp_rsfa_file)
gsp_rsfa$roi  = gsub('_scale','',gsp_rsfa$roi)

# Harvard GBC
gsp_gbc_file = paste0(project_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv')
gsp_gbc      = read_csv(gsp_gbc_file)

# merge MDD regression DF with AHBA expression
thick_expr_df = merge(x=ukbb_thick, by.x='roi', y=ahba_expr_200_t, by.y='roi')
rsfa_expr_df  = merge(x=ukbb_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
gbc_expr_df   = merge(x=ukbb_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')
#
gsp_rsfa_expr_df = merge(x=gsp_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
gsp_gbc_expr_df  = merge(x=gsp_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')




# Step 2: Read cell specific genes from Lake
# --------------
# read Lake differential expression
lake_dexpr = read_csv('/gpfs/milgram/project/holmes/kma52/mdd_sst/reference_files/lake_all_dexpr.csv')

# Correlate THICKNESS with AHBA marker genes
sst_marker_genes = c('SST', 'NPY', 'CORT')

# remove SST markers from lake data, so they aren't in the permuted gene triplets
lake_dexpr       = lake_dexpr[which(!lake_dexpr$Gene %in% sst_marker_genes),]
lake_dexpr       = lake_dexpr[!grepl('Cer|Purk|Gran', lake_dexpr$Cluster),]
lake_dexpr_nosst = lake_dexpr[!grepl('In7|In8', lake_dexpr$Cluster),]


set.seed(287934)




# Step 3: Spatial correlation of all genes to MDD Regression effects
# --------------
# UKB THICKNESS - genewise correlation
# ------
# ------
# ------
gene_names       = colnames(ahba_expr_200_t)[1:ncol(ahba_expr_200_t)-1]
thick_gene_cors  = t(cor(thick_expr_df$cohens_d, thick_expr_df[gene_names], use='pairwise.complete', method='spearman'))
thick_gene_cors  = thick_gene_cors[order(thick_gene_cors),]
head(thick_gene_cors)

# avg correlation of markers
sst_marker_cor     = thick_gene_cors[names(thick_gene_cors) %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)
sst_marker_cor_avg


# Step 3a: gene-based permutations
# randomly select three genes and see if any triplets are more associated
thick_gene_df    = data.frame(cor=thick_gene_cors, genes=names(thick_gene_cors))
out              = marker_gene_permutation(gene_cors=thick_gene_df, genes=sst_marker_genes, num_perm=10000)
thick_perm_arr   = out[[1]]
thick_perm_trips = out[[2]]

# # Step 3b: same, but only using cell makers
out_path = paste0(project_dir, '/figures/ukbb_thick_mdd_permutation_density_sst_markers.pdf')
out_path
thick_plot_df = data.frame(marker_perm = thick_perm_arr)
p = plot_genewise_perm_density(plot_df=thick_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.3,.3,.15), ylims=c(0,5))


# # Step 3c: same as above, but only select from genes that are significant markers
marker_thick_gene_df = thick_gene_df[thick_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm     = marker_gene_permutation(gene_cors=marker_thick_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/ukbb_thick_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
thick_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=thick_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.3,.3,.15), ylims=c(0,5))


# # Step 3c: now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr_nosst, gene_cors=marker_thick_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/ukbb_thick_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]] = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
thick_cell_marker_plot_df  = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=thick_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.45,.45,.15), ylims=c(0,5))
# ------
# ------
# ------




# UKB RSFA - genewise correlation
# ------
# ------
# ------
rsfa_gene_cors = t(cor(rsfa_expr_df$cohens_d, rsfa_expr_df[gene_names], use='pairwise.complete', method='spearman'))
rsfa_gene_cors = rsfa_gene_cors[order(rsfa_gene_cors),]
#rsfa_gene_cors[names(rsfa_gene_cors) %in% gene]

# avg correlation of markers
sst_marker_cor     = rsfa_gene_cors[names(rsfa_gene_cors) %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)

# randomly select three genes and see if any triplets are more associated
rsfa_gene_df    = data.frame(cor=rsfa_gene_cors, genes=names(rsfa_gene_cors))
out             = marker_gene_permutation(gene_cors=rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)
rsfa_perm_arr   = out[[1]]
rsfa_perm_trips = out[[2]]


# same, but only using cell makers
out_path = paste0(project_dir, '/figures/ukbb_rsfa_mdd_permutation_density_sst_markers.pdf')
out_path
rsfa_plot_df = data.frame(marker_perm = rsfa_perm_arr)
p = plot_genewise_perm_density(plot_df=rsfa_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.5,.5,.25), ylims=c(0,4))


# same as above, but only select from genes that are significant markers
marker_rsfa_gene_df  = rsfa_gene_df[rsfa_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm     = marker_gene_permutation(gene_cors=marker_rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/ukbb_rsfa_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
rsfa_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=rsfa_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.5,.5,.25), ylims=c(0,3))


# now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr_nosst, gene_cors=marker_rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)
cell_marker_gene_perm[[1]]

out_path = paste0(project_dir, '/figures/ukbb_rsfa_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]] = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
rsfa_cell_marker_plot_df = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=rsfa_cell_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.6,.6,.3), ylims=c(0,3))
# ------
# ------
# ------






# UKB GBC - genewise permutation
# ------
# ------
# ------
gbc_gene_cors = t(cor(gbc_expr_df$cohens_d, gbc_expr_df[gene_names], use='pairwise.complete', method='spearman'))
gbc_gene_cors = gbc_gene_cors[order(gbc_gene_cors),]
gbc_gene_cors[which(names(gbc_gene_cors) %in% sst_marker_genes)]
gbc_gene_cors[names(gbc_gene_cors) %in% sst_marker_genes]


# avg correlation of markers
sst_marker_cor     = gbc_gene_cors[names(gbc_gene_cors) %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)


# randomly select three genes and see if any triplets are more associated
gbc_gene_df     = data.frame(cor=gbc_gene_cors, genes=names(gbc_gene_cors))
out             = marker_gene_permutation(gene_cors=gbc_gene_df, genes=sst_marker_genes, num_perm=10000)
gbc_perm_arr    = out[[1]]
gbc_perm_trips  = out[[2]]


# same, but only using cell makers
out_path = paste0(project_dir, '/figures/ukbb_gbc_mdd_permutation_density_sst_markers.pdf')
out_path
gbc_plot_df = data.frame(marker_perm = gbc_perm_arr)
p = plot_genewise_perm_density(plot_df=gbc_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.4,.4,.2), ylims=c(0,4))


# same as above, but only select from genes that are significant markers
marker_gbc_gene_df  = gbc_gene_df[gbc_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm    = marker_gene_permutation(gene_cors=marker_gbc_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/ukbb_gbc_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
gbc_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=gbc_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.45,.45,.15), ylims=c(0,3))


# now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr_nosst, gene_cors=marker_gbc_gene_df, genes=sst_marker_genes, num_perm=10000)
cell_marker_gene_perm[[1]]

out_path = paste0(project_dir, '/figures/ukbb_gbc_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]] = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
gbc_cell_marker_plot_df    = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=gbc_cell_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.6,.6,.3), ylims=c(0,3))
# ------
# ------
# ------



# Harvard GBC
# ------
# ------
# ------
GSP_gbc_gene_cors = t(cor(gsp_gbc_expr_df$scaleneg_affect_mean_cohens_d, gsp_gbc_expr_df[gene_names], use='pairwise.complete', method='spearman'))
GSP_gbc_gene_cors = GSP_gbc_gene_cors[order(GSP_gbc_gene_cors),]
GSP_gbc_gene_cors[names(GSP_gbc_gene_cors) %in% sst_marker_genes]

# avg correlation of markers
sst_marker_cor     = GSP_gbc_gene_cors[names(GSP_gbc_gene_cors) %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)

# randomly select three genes and see if any triplets are more associated
GSP_gbc_gene_df = data.frame(cor=GSP_gbc_gene_cors, genes=names(GSP_gbc_gene_cors))
out             = marker_gene_permutation(gene_cors=GSP_gbc_gene_df, genes=sst_marker_genes, num_perm=10000)
GSP_gbc_perm_arr    = out[[1]]
GSP_gbc_perm_trips  = out[[2]]

# same, but only using cell makers
out_path = paste0(project_dir, '/figures/gsp_gbc_mdd_permutation_density_sst_markers.pdf')
out_path
GSP_gbc_plot_df = data.frame(marker_perm = GSP_gbc_perm_arr)
p = plot_genewise_perm_density(plot_df=GSP_gbc_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.4,.4,.2), ylims=c(0,4))

# same as above, but only select from genes that are significant markers
marker_GSP_gbc_gene_df = GSP_gbc_gene_df[GSP_gbc_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm       = marker_gene_permutation(gene_cors=marker_GSP_gbc_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/GSP_gbc_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
GSP_gbc_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=GSP_gbc_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.45,.45,.15), ylims=c(0,3))

# now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr, gene_cors=marker_GSP_gbc_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/GSP_gbc_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]] = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
GSP_gbc_cell_marker_plot_df    = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=GSP_gbc_cell_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.6,.6,.3), ylims=c(0,3))
# ------
# ------
# ------





# Harvard RSFA
# ------
# ------
# ------
GSP_rsfa_gene_cors = t(cor(gsp_rsfa_expr_df$scaleneg_affect_mean_cohens_d, gsp_rsfa_expr_df[gene_names], use='pairwise.complete', method='spearman'))
GSP_rsfa_gene_cors = GSP_rsfa_gene_cors[order(GSP_rsfa_gene_cors),]
GSP_rsfa_gene_cors[names(GSP_rsfa_gene_cors) %in% sst_marker_genes]


# avg correlation of markers
sst_marker_cor     = GSP_gbc_gene_cors[names(GSP_rsfa_gene_cors) %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)


# randomly select three genes and see if any triplets are more associated
GSP_rsfa_gene_df = data.frame(cor=GSP_rsfa_gene_cors, genes=names(GSP_rsfa_gene_cors))
out              = marker_gene_permutation(gene_cors=GSP_rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)
GSP_rsfa_perm_arr    = out[[1]]
GSP_rsfa_perm_trips  = out[[2]]


# same, but only using cell makers
out_path = paste0(project_dir, '/figures/gsp_rsfa_mdd_permutation_density_sst_markers.pdf')
out_path
GSP_rsfa_plot_df = data.frame(marker_perm = GSP_rsfa_perm_arr)
p = plot_genewise_perm_density(plot_df=GSP_rsfa_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.4,.4,.2), ylims=c(0,5))


# same as above, but only select from genes that are significant markers
marker_GSP_rsfa_gene_df = GSP_rsfa_gene_df[GSP_rsfa_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm       = marker_gene_permutation(gene_cors=marker_GSP_rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/GSP_rsfa_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
GSP_rsfa_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=GSP_rsfa_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.40,.40,.2), ylims=c(0,5))


# now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr_nosst, gene_cors=marker_GSP_rsfa_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/GSP_rsfa_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]]    = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
GSP_rsfa_cell_marker_plot_df  = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=GSP_rsfa_cell_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.4,.4,.2), ylims=c(0,5))
# ------
# ------
# ------





# ENIGMA
# ------
# ------
# ------
enigma_mdd  = read_csv(paste0(project_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv'))

# gene by desikan ROI matrix
desikan_mat   = read_csv(paste0(project_dir, '/data/ahba_parcel/desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))
desikan_mat_t = as.data.frame(t(desikan_mat[,1:70]))
colnames(desikan_mat_t) = desikan_mat$gene
desikan_mat_t$roi = as.character(rownames(desikan_mat_t))
desikan_mat_t[1:5,1:5]

# combine ENIGMA cohens d with gene matrix
enigma_ahba_df = merge(x=desikan_mat_t, y=enigma_mdd, by='roi')


# gene wise spatial correlation to cohen's d
enigma_ahba_mdd_cors = cor(enigma_ahba_df[desikan_mat$gene], enigma_ahba_df$cohens_d, use='pairwise.complete', method='spearman')
enigma_ahba_mdd_cors = data.frame(cor=enigma_ahba_mdd_cors, genes=rownames(enigma_ahba_mdd_cors))
enigma_ahba_mdd_cors = enigma_ahba_mdd_cors[order(enigma_ahba_mdd_cors$cor),]

# spatial cor of marker genes
enigma_ahba_mdd_cors[which(enigma_ahba_mdd_cors$genes %in% sst_marker_genes),]
sst_marker_cor     = enigma_ahba_mdd_cors$cor[enigma_ahba_mdd_cors$gene %in% sst_marker_genes]
sst_marker_cor_avg = mean(sst_marker_cor)

# triplet permutations
enigma_gene_df  = data.frame(cor=enigma_ahba_mdd_cors$cor, genes=enigma_ahba_mdd_cors$genes)
enigma_out      = marker_gene_permutation(gene_cors=enigma_gene_df, genes=sst_marker_genes, num_perm=10000)
enigma_gene_arr = enigma_out[[1]]

out_path       = paste0(project_dir, '/figures/enigma_mdd_permutation_density_sst_markers.pdf')
enigma_plot_df = data.frame(marker_perm = enigma_gene_arr)
p = plot_genewise_perm_density(plot_df=enigma_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.6,.6,.3), ylims=c(0,2.5))

# same as above, but only select from genes that are significant markers
marker_enigma_gene_df = enigma_gene_df[enigma_gene_df$genes %in% c(sst_marker_genes, lake_dexpr_nosst$Gene),]
marker_gene_perm      = marker_gene_permutation(gene_cors=marker_enigma_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/enigma_mdd_permutation_density_sst_markers_CELLMARKS.pdf')
out_path
enigma_marker_plot_df = data.frame(marker_perm = marker_gene_perm[[1]])
p = plot_genewise_perm_density(plot_df=enigma_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.6,.6,.3), ylims=c(0,2.5))


# now, make the permutation a bit more rigorous
# gene triplets will be selected from pools of cell specific genes
# e.g. three signficiant Ex1 marker genes will be grabbed (10k triplets per cell type)
cell_marker_gene_perm = cell_marker_gene_permutation(lake_dexpr=lake_dexpr_nosst, gene_cors=enigma_gene_df, genes=sst_marker_genes, num_perm=10000)

out_path = paste0(project_dir, '/figures/enigma_mdd_permutation_density_sst_markers_CELLMARKS_within.pdf')
out_path
cell_marker_gene_perm[[2]]    = cell_marker_gene_perm[[2]][!is.na(cell_marker_gene_perm[[2]])]
enigma_cell_marker_plot_df    = data.frame(marker_perm = cell_marker_gene_perm[[2]])
p = plot_genewise_perm_density(plot_df=enigma_cell_marker_plot_df, target_cor=sst_marker_cor_avg, out_path=out_path, xlims=c(-0.7,.7,.35), ylims=c(0,2))
# ------
# ------
# ------







# Step 4: Combine and Save
# ------------
# combine the gene-wise spatial correlations into a single dataframe
rsfa_gene_cors_df         = data.frame(genes=names(rsfa_gene_cors), ukbb_rsfa_cor=rsfa_gene_cors)
thick_gene_cors_df        = data.frame(genes=names(thick_gene_cors), ukbb_thick_cor=thick_gene_cors)
gbc_gene_cors_df          = data.frame(genes=names(gbc_gene_cors), ukbb_gbc_cor=gbc_gene_cors)
GSP_rsfa_gene_cors_df     = data.frame(genes=names(GSP_rsfa_gene_cors), gsp_rsfa_cor=GSP_rsfa_gene_cors)
GSP_gbc_gene_cors_df      = data.frame(genes=names(GSP_gbc_gene_cors), gsp_gbc_cor=GSP_gbc_gene_cors)
enigma_thick_gene_cors_df = data.frame(genes=enigma_ahba_mdd_cors$genes, enigma_thick_cor=enigma_ahba_mdd_cors$cor)


all_gene_df = merge(x=thick_gene_cors_df, y=rsfa_gene_cors_df, by='genes')
all_gene_df = merge(x=all_gene_df, y=gbc_gene_cors_df, by='genes')
all_gene_df = merge(x=all_gene_df, y=GSP_rsfa_gene_cors_df, by='genes')
all_gene_df = merge(x=all_gene_df, y=GSP_gbc_gene_cors_df, by='genes')
all_gene_df = merge(x=all_gene_df, y=enigma_thick_gene_cors_df, by='genes')

write_csv(x=all_gene_df, paste0(project_dir, "/output/mdd_regr/all_mdd_ahba_genecors.csv"))



