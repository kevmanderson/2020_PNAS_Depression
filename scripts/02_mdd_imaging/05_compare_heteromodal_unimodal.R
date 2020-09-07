library(tidyverse)
library(Cairo)
library(RColorBrewer)
library(corrplot)


# 1. plot uni-modal cohen's d for heteromodal and unimodal parcels
# 2. permute parcel labels to generate statistics


net_wise_effects = function(mdd_df, hex_df, plot_col, out_name, ylims){

    # get average RSFA for each networks
    plot_df = mdd_df
    plot_df$cohens_d = plot_df[[plot_col]]
    plot_df$net_names = unlist(lapply(plot_df$roi, function(x) strsplit(x,'_')[[1]][3]))
    plot_df$net_names[grep('Limbic_OFC', plot_df$roi)]='Limbic_OFC'
    plot_df$net_names[grep('Limbic_TempPole', plot_df$roi)]='Limbic_TempPole'
    plot_summary = plot_df %>%
                        group_by(net_names) %>%
                        summarise(mean=mean(cohens_d), sd=sd(cohens_d), se=sd(cohens_d)/sqrt(length(cohens_d)))

    plot_summary$net_names = factor(plot_summary$net_names, levels=colnames(hex_df))

    out_path = paste0(base_dir, '/figures/ukbb_rsfa_effect_by_network.pdf')
    out_path
    CairoPDF(out_path, width=8, height=3)
    ggplot(plot_summary, aes(x=net_names, y=mean, fill=net_names)) +
      geom_bar(stat="identity", color="black",
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4,
                     position=position_dodge(.9)) +
        scale_fill_manual(values=as.character(t(hex_df))) +
        theme_classic() + scale_y_continuous(limits = c(-.05, .05), breaks=seq(-.05,.05,.03), expand = c(0,0)) +
        theme(axis.text.x=element_text(angle=90), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
    dev.off()
}


het_uni_effects = function(plot_df, out_name, ylims){

    # get average RSFA for each networks
    plot_df$net_names = unlist(lapply(plot_df$roi, function(x) strsplit(x,'_')[[1]][3]))
    plot_df$uni_hetero = ifelse(grepl('Vis|SomMot', plot_df$net_names),'Unimodal', 'Heteromodal')
    plot_df %>% group_by(uni_hetero) %>%
                        summarise(mean=mean(cohens_d), sd=sd(cohens_d), se=sd(cohens_d)/sqrt(length(cohens_d)))

    out_path = paste0(base_dir, '/figures/',out_name,'_newmap.pdf')
    print(out_path)
    CairoPDF(out_path, width=1.25, height=1.9)
    p = ggplot(plot_df, aes(x=uni_hetero, y=cohens_d, fill=uni_hetero)) +
        geom_violin() +
        geom_boxplot(width=0.3) +
        theme_classic() +
        scale_fill_manual(values=c("#CE8196", "#9CC8C0")) +
        scale_y_continuous(limits = c(ylims[1], ylims[2]), breaks=round(seq(ylims[1], ylims[2], ylims[3]), 3), expand = c(0,0)) +
        theme(legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
    print(p)
    dev.off()
}

het_uni_perm_effects = function(plot_df){

    perm_df = plot_df
    # get average RSFA for each networks
    permutations_df = NULL
    for (perm in 1:1000){
        perm_df$uni_hetero_perm = perm_df$uni_hetero[sample(1:length(perm_df$uni_hetero))]
        uni_mean = mean(perm_df$cohens_d[perm_df$uni_hetero_perm == 'Unimodal'])
        het_mean = mean(perm_df$cohens_d[perm_df$uni_hetero_perm == 'Heteromodal'])
        out_row = data.frame(uni=uni_mean, het=het_mean, perm=perm)
        permutations_df = rbind(permutations_df, out_row)
    }
    obs_uni_mean = mean(plot_df$cohens_d[perm_df$uni_hetero == 'Unimodal'])
    obs_het_mean = mean(plot_df$cohens_d[perm_df$uni_hetero == 'Heteromodal'])
    obs_diff = obs_uni_mean - obs_het_mean
    perm_diffs = c(obs_diff, permutations_df$uni - permutations_df$het)
    print(length(which(abs(perm_diffs) >= abs(obs_diff)))/length(perm_diffs))
}


hex_df = data.frame(VisCent='#781286',
                    VisPeri='#ff0101',
                    SomMotA='#4682b2',
                    SomMotB='#2bcca2',
                    DorsAttnA='#4a9b3d',
                    DorsAttnB='#007610',
                    SalVentAttnA='#c43afb',
                    SalVentAttnB='#ff98d6',
                    LimbicB_OFC='#7a8733',
                    LimbicA_TempPole='#dcf8a5',
                    ContA='#e69423',
                    ContB='#87324c',
                    ContC='#778cb1',
                    DefaultA='#ffff02',
                    DefaultB='#cd3e50',
                    DefaultC='#000083',
                    TempPar='#0929fa')



# Step 1: Read parcel-wise estimates
# -----------
base_dir = "/gpfs/milgram/project/holmes/kma52/mdd_gene_expr"


# UKB thickness
ukbb_thick = read_csv(paste0(base_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv'))
colnames(ukbb_thick) = paste0('ukb_thick_', colnames(ukbb_thick))

# UKB RSFA
ukbb_rsfa  = read_csv(paste0(base_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv'))
ukbb_rsfa$roi = gsub('_scale','', ukbb_rsfa$roi)
colnames(ukbb_rsfa) = paste0('ukb_rsfa_', colnames(ukbb_rsfa))

# UKB GBC
ukbb_gbc   = read_csv(paste0(base_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv'))
colnames(ukbb_gbc) = paste0('ukb_gbc_', colnames(ukbb_gbc))

# GSP RSFA
gsp_rsfa   = read_csv(paste0(base_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv'))
gsp_rsfa$roi = gsub('_scale','', gsp_rsfa$roi)
colnames(gsp_rsfa) = paste0('gsp_rsfa_', colnames(gsp_rsfa))

# GSP GBC
gsp_gbc    = read_csv(paste0(base_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv'))
colnames(gsp_gbc) = paste0('gsp_gbc_', colnames(gsp_gbc))

# merge all regression tables
cohen_df = merge(x=ukbb_thick, y=ukbb_rsfa, by.x='ukb_thick_roi', by.y='ukb_rsfa_roi')
cohen_df = merge(x=cohen_df, y=ukbb_gbc, by.x='ukb_thick_roi', by.y='ukb_gbc_roi')
cohen_df = merge(x=cohen_df, y=gsp_rsfa, by.x='ukb_thick_roi', by.y='gsp_rsfa_roi')
cohen_df = merge(x=cohen_df, y=gsp_gbc, by.x='ukb_thick_roi', by.y='gsp_gbc_roi')


# Step 2: Cohen's D correlation matrix
# -----------
# correlation of MDD effects across ROIs
cor_cohens_df = cohen_df[grep('cohens_d', colnames(cohen_df))]
colnames(cor_cohens_df) = gsub('_cohens_d|_scaleneg_affect_mean_cohens_d','',colnames(cor_cohens_df))
cor_cohens_df = cor_cohens_df[c('ukb_thick','ukb_rsfa','gsp_rsfa','ukb_gbc','gsp_gbc')]
cor_mat  = cor(method='spearman', cor_cohens_df)

# set color scale
col      = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corr_out = paste0(base_dir, '/figures/cort_mdd_effects_multimodal_corrplot_newmap.pdf')
corr_out

# make correlation plot
corr_out
pdf(corr_out, height=3, width=3)
corrplot(cor_mat, method="color", col=rev(col(200)),
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=T)
dev.off()


# Step 3: GSP GBC between hetero/uni-modal cortex
# -----------
plot_df            = gsp_gbc
plot_df$cohens_d   = plot_df$gsp_gbc_scaleneg_affect_mean_cohens_d
plot_df$roi        = plot_df$gsp_gbc_roi
plot_df$net_names  = unlist(lapply(plot_df$gsp_gbc_roi, function(x) strsplit(x,'_')[[1]][3]))
plot_df$uni_hetero = ifelse(grepl('Vis|SomMot', plot_df$net_names),'Unimodal', 'Heteromodal')

het_uni_effects(plot_df=plot_df, out_name='gsp_gbc_effect_by_UnimodalHeteromodal', ylims=c(-.2,.2,.1))
# permutation based significance
het_uni_perm_effects(plot_df)
plot_df %>% group_by(uni_hetero) %>%
        summarize(m_cohens=mean(cohens_d), sd_cohens=sd(cohens_d))



# Step 4: GSP RSFA between hetero/uni-modal cortex
# -----------
plot_df            = gsp_rsfa
plot_df$cohens_d   = plot_df$gsp_rsfa_scaleneg_affect_mean_cohens_d
plot_df$roi        = plot_df$gsp_rsfa_roi
plot_df$net_names  = unlist(lapply(plot_df$roi, function(x) strsplit(x,'_')[[1]][3]))
plot_df$uni_hetero = ifelse(grepl('Vis|SomMot', plot_df$net_names),'Unimodal', 'Heteromodal')

het_uni_effects(plot_df=plot_df, out_name='gsp_rsfa_effect_by_UnimodalHeteromodal', ylims=c(-.2,.2,.1))
# permutation based significance
het_uni_perm_effects(plot_df)
plot_df %>% group_by(uni_hetero) %>%
        summarize(m_cohens=mean(cohens_d), sd_cohens=sd(cohens_d))



# Step 5: UKB RSFA between hetero/uni-modal cortex
# -----------
plot_df            = ukbb_rsfa
plot_df$cohens_d   = plot_df$ukb_rsfa_cohens_d
plot_df$roi        = plot_df$ukb_rsfa_roi
plot_df$net_names  = unlist(lapply(plot_df$roi, function(x) strsplit(x,'_')[[1]][3]))
plot_df$uni_hetero = ifelse(grepl('Vis|SomMot', plot_df$net_names),'Unimodal', 'Heteromodal')
het_uni_effects(plot_df=plot_df, out_name='ukbb_rsfa_effect_by_UnimodalHeteromodal', ylims=c(-.08,.08,.04))
# permutation based significance
het_uni_perm_effects(plot_df)
plot_df %>% group_by(uni_hetero) %>%
        summarize(m_cohens=mean(cohens_d), sd_cohens=sd(cohens_d))


# Step 6: UKB GBC between hetero/uni-modal cortex
# -----------
plot_df            = ukbb_gbc
plot_df$cohens_d   = plot_df$ukb_gbc_cohens_d
plot_df$roi        = plot_df$ukb_gbc_roi
plot_df$net_names  = unlist(lapply(plot_df$roi, function(x) strsplit(x,'_')[[1]][3]))
plot_df$uni_hetero = ifelse(grepl('Vis|SomMot', plot_df$net_names),'Unimodal', 'Heteromodal')
het_uni_effects(plot_df=plot_df, out_name='ukbb_gbc_effect_by_UnimodalHeteromodal', ylims=c(-.08,.08,.04))
# permutation based significance
het_uni_perm_effects(plot_df)
plot_df %>% group_by(uni_hetero) %>%
        summarize(m_cohens=mean(cohens_d), sd_cohens=sd(cohens_d))

