library(tidyverse)
library(ggcorrplot)
library(Cairo)
library(fgsea)
library(gridExtra)
library(grid)
library(RRHO)
library(RColorBrewer)


# 1. Rank-Rank Hypergeometric tests for consistent transcriptional signatures of MDD imaging phenotypes
# 2. Fast PreRanked Gene Set Enrichment Analysis (FGSEA)


####################
# Read AHBA/regr dat
####################
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read gene-wise spatial correlations to MDD effect maps
all_ahba_cors = read_csv(paste0(base_dir, "/output/mdd_regr/all_mdd_ahba_genecors.csv"))

# reverse direction of RSFA to line up with other modalities
all_ahba_cors$ukbb_rsfa_cor_new = all_ahba_cors$ukbb_rsfa_cor*-1
all_ahba_cors$ukbb_rsfa_cor     = NULL
all_ahba_cors$gsp_rsfa_cor_new  = all_ahba_cors$gsp_rsfa_cor*-1
all_ahba_cors$gsp_rsfa_cor      = NULL

# average together and sort
all_ahba_cors$avg = rowMeans(all_ahba_cors[grep('cor', colnames(all_ahba_cors))])
all_ahba_cors     = all_ahba_cors[order(all_ahba_cors$avg),]

# save
write_csv(all_ahba_cors, paste0(base_dir, "/output/mdd_regr/all_mdd_ahba_genecors_wavg.csv"))
write_csv(all_ahba_cors, paste0(base_dir, "/supp_data/SuppData_9_all_mdd_ahba_genecors_wavg.csv"))


# Define cell specific genes from Lake differential expression data



#####################
# Read Lake diffexpr
#####################
lake_dexpr = read_csv(paste0(base_dir, '/reference_files/lake_all_dexpr.csv'))

# get rid of cerebellar cell types
lake_dexpr = lake_dexpr[!grepl('Cer|Purk|Gran', lake_dexpr$Cluster),]
colnames(lake_dexpr) = gsub(' ','_',colnames(lake_dexpr))

# list of lake cells
lake_cells = unique(lake_dexpr$Cluster)


# load previously computed differential expression
load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_frontal_diff_expr_batchcorrect.Rdata'))
# dfc_diff_expr_list
load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_visual_diff_expr_batchcorrect.Rdata'))
# vis_diff_expr_list


# combine vis/dfc differentially expressed genes
both_lake_cells = intersect(names(dfc_diff_expr_list), names(vis_diff_expr_list))
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


# organize and write supplemental data
max_gene     = max(unlist(lapply(both_lake_sets, length)))
cell_gene_df = as.data.frame(matrix(NA, max_gene, length(both_lake_sets)))
colnames(cell_gene_df) = names(both_lake_sets)
for (cell in names(both_lake_sets)){
    write(cell,'')
    cell_genes  = both_lake_sets[[cell]]
    write_genes = c(cell_genes, rep(NA, max_gene - length(cell_genes)))
    cell_gene_df[cell] = write_genes
}
write_csv(cell_gene_df, paste0(base_dir, "/supp_data/SuppData_10_lake_cell_genes.csv"))



# Run FGSEA for the gene correlates of each MDD imaging phenotype
both_fgsea_all  = NULL
both_fgsea_list = NULL
fgsea_all_ahba_cors = all_ahba_cors
mri_phenos      = colnames(fgsea_all_ahba_cors)[grep('_cor', colnames(fgsea_all_ahba_cors))]
for (pheno in mri_phenos){
    write(pheno,'')
    ahba_gene_rank        = -1*fgsea_all_ahba_cors[[pheno]]
    names(ahba_gene_rank) = fgsea_all_ahba_cors$genes

    # run FGSEA
    fgseaRes = fgseaMultilevel(pathways = both_lake_sets,
                      stats = ahba_gene_rank,
                      minSize = 15,
                      maxSize = 1000)

    fgseaRes$pheno      = pheno
    fgseaRes$signlog10p = sign(fgseaRes$ES)*(1*log10(fgseaRes$padj))
    both_fgsea_all      = rbind(both_fgsea_all, fgseaRes)
    both_fgsea_list[[pheno]] = fgseaRes
}

# calc average FGSEA enrichment
fgsea_nes_df = cbind(both_fgsea_list[[1]]$NES,
                        both_fgsea_list[[2]]$NES,
                        both_fgsea_list[[3]]$NES,
                        both_fgsea_list[[4]]$NES,
                        both_fgsea_list[[5]]$NES,
                        both_fgsea_list[[6]]$NES)
fgsea_es_df = cbind(both_fgsea_list[[1]]$ES,
                        both_fgsea_list[[2]]$ES,
                        both_fgsea_list[[3]]$ES,
                        both_fgsea_list[[4]]$ES,
                        both_fgsea_list[[5]]$ES,
                        both_fgsea_list[[6]]$ES)

# construct the avg df
avg_NES            = rowMeans(fgsea_nes_df)
avg_fgseaRes       = fgseaRes
avg_fgseaRes$NES   = avg_NES
avg_fgseaRes$ES    = rowMeans(fgsea_es_df)
avg_fgseaRes$signlog10p = 0
avg_fgseaRes$pval  = 1
avg_fgseaRes$padj  = 1
avg_fgseaRes$log2err  = 1
avg_fgseaRes$pheno = 'avg'
both_fgsea_all     = rbind(both_fgsea_all, avg_fgseaRes)

# order to plot FGSEA values
pheno_order = c("ukbb_thick_cor",
                "enigma_thick_cor",
                "ukbb_rsfa_cor_new",
                "gsp_rsfa_cor_new",
                "ukbb_gbc_cor",
                "gsp_gbc_cor",
                'avg')

# start constructint the df to plot
both_fgsea_all$pheno   = factor(both_fgsea_all$pheno, levels=pheno_order)
plot_order             = both_fgsea_list[[1]]$pathway[rev(order(avg_NES))]
both_fgsea_all$pathway = factor(both_fgsea_all$pathway, levels=rev(plot_order))
both_fgsea_all$signlog10p_text = round(-1*log10(both_fgsea_all$padj), 2)

write_me = both_fgsea_all
write_me$leadingEdge = NULL
write.csv(as.data.frame(write_me), paste0(base_dir, "/supp_data/SuppData_11_AHBA_MDD_FGSEA.csv"))


both_fgsea_all$padj_text = formatC(both_fgsea_all$pval, format = "e", digits = 1)
both_fgsea_all$padj_text[both_fgsea_all$pval > .01] = round(both_fgsea_all$pval[both_fgsea_all$pval > .01], 2)

# FGSEA heatmap
p = ggplot(data = both_fgsea_all, aes(x=pathway, y=pheno, fill=NES)) +
            geom_tile() +
            geom_text(aes(x=pathway, y=pheno, label = padj_text), size=4, color = "black", size = 4) +
            coord_flip() +
            scale_fill_gradient2(low = "#5063AE", high = "#C32B2F", mid = "white", midpoint = 0, limit = c(-3,3))+
            theme_classic()

fgsea_out = paste0(base_dir, '/figures/big_celltype_ahba_fgsea.pdf')
fgsea_out
CairoPDF(fgsea_out, width=5)
print(p)
dev.off()



# FGSEA enrichment score plots for single cells
avg_ahba_gene_rank        = -1*fgsea_all_ahba_cors[['avg']]
names(avg_ahba_gene_rank) = fgsea_all_ahba_cors$genes

p = plotEnrichment(both_lake_sets[["Ast"]], avg_ahba_gene_rank) + labs(title="Astrocytes")
fgsea_out = paste0(base_dir, '/figures/astrocyte_ahba_avg_fgsea.pdf')
CairoPDF(fgsea_out, width=5, height=3)
print(p)
dev.off()
fgsea_out


p = plotEnrichment(both_lake_sets[["In6"]], avg_ahba_gene_rank) + labs(title="In6")
fgsea_out = paste0(base_dir, '/figures/PVALB_ahba_avg_fgsea.pdf')
CairoPDF(fgsea_out, width=5, height=3)
print(p)
dev.off()
fgsea_out




#######
# MTG
#######
load(verbose=T, paste0(base_dir, '/output/diff_expr/mtg_big_cat_ahba_diff_expr.Rdata'))
load(verbose=T, paste0(base_dir, '/output/diff_expr/mtg_cell_ALL_cell_fine_diff_expr_list.Rdata'))

fgsea_mtg = NULL
for (cell in names(mtg_fine_diff_expr_list)){

    cell_dexpr     = mtg_fine_diff_expr_list[[cell]]
    cell_dexpr_df  = cell_dexpr[which(cell_dexpr$avg_logFC > 0),]
    cell_dexpr_df  = cell_dexpr_df[which(cell_dexpr_df$p_val_adj < .05),]

    fgsea_mtg[[cell]] = rownames(cell_dexpr_df)
}


# Run FGSEA for the gene correlates of each MDD imaging phenotype
both_fgsea_all  = NULL
both_fgsea_list = NULL
fgsea_all_ahba_cors = all_ahba_cors
mri_phenos = colnames(fgsea_all_ahba_cors)[grep('_cor', colnames(fgsea_all_ahba_cors))]
for (pheno in mri_phenos){
    write(pheno,'')
    ahba_gene_rank        = -1*fgsea_all_ahba_cors[[pheno]]
    names(ahba_gene_rank) = fgsea_all_ahba_cors$genes

    # run FGSEA
    fgseaRes = fgseaMultilevel(pathways = fgsea_mtg,
                      stats = ahba_gene_rank,
                      minSize = 15,
                      maxSize = 1000)

    fgseaRes$pheno      = pheno
    fgseaRes$signlog10p = sign(fgseaRes$ES)*(1*log10(fgseaRes$padj))
    both_fgsea_all      = rbind(both_fgsea_all, fgseaRes)
    both_fgsea_list[[pheno]] = fgseaRes
}













