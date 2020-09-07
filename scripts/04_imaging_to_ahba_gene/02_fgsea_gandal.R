library(tidyverse)
library(ggcorrplot)
library(Cairo)
library(fgsea)
library(gridExtra)
library(grid)

# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'


# load Gandal disorder diff expr
table_dir  = '/gpfs/milgram/project/holmes/kma52/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master/results/results/tables'
gandal_mdd = read_csv(paste0(table_dir, '/Microarray_MDD_metaanalysis_092017_scale.csv'))
gandal_mdd$beta_scale = gandal_mdd$beta
gandal_scz = read_csv(paste0(table_dir, '/Microarray_SCZ_metaanalysis_092017_scale.csv'))
gandal_bd  = read_csv(paste0(table_dir, '/Microarray_BD_metaanalysis_092017_scale.csv'))
gandal_asd = read_csv(paste0(table_dir, '/Microarray_ASD_metaanalysis_092017_scale.csv'))
gandal_aad = read_csv(paste0(table_dir, '/Microarray_AAD_metaanalysis_092017_scale.csv'))


load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_frontal_diff_expr_batchcorrect.Rdata'))
# dfc_diff_expr_list
load(verbose=T, paste0(base_dir, '/output/diff_expr/lake_primary_cats_visual_diff_expr_batchcorrect.Rdata'))
# vis_diff_expr_list


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

both_fgsea_all  = NULL
both_fgsea_list = NULL
fgsea_all_ahba_cors = all_ahba_cors
mri_phenos = colnames(fgsea_all_ahba_cors)[grep('_cor', colnames(fgsea_all_ahba_cors))]



exvivo_fgsea = function(dexpr_df, cell_sets, pheno){

    dexpr_df     = dexpr_df[!duplicated(dexpr_df$symbol),]
    dexpr_df     = dexpr_df[!is.na(dexpr_df$beta_scale),]
    dexpr_df     = dexpr_df[order(dexpr_df$beta_scale),]
    dexpr_rank   = dexpr_df$beta_scale*-1

    names(dexpr_rank) = dexpr_df$symbol
    dexpr_fgsea = fgseaMultilevel(pathways = both_lake_sets,
                                  stats = dexpr_rank,
                                  minSize=15,
                                  maxSize=1500)
    dexpr_fgsea$signlog10p = sign(dexpr_fgsea$ES)*(1*log10(dexpr_fgsea$pval))
    dexpr_fgsea$log10p     = -1*log10(dexpr_fgsea$pval)
    dexpr_fgsea$pheno      = pheno
    return(dexpr_fgsea)
}


# run cell enrichment for all psych disorders
mdd_fgsea = exvivo_fgsea(dexpr_df=gandal_mdd, cell_sets=both_lake_sets, pheno='MDD')
scz_fgsea = exvivo_fgsea(dexpr_df=gandal_scz, cell_sets=both_lake_sets, pheno='SCZ')
bd_fgsea  = exvivo_fgsea(dexpr_df=gandal_bd, cell_sets=both_lake_sets, pheno='BD')
asd_fgsea = exvivo_fgsea(dexpr_df=gandal_asd, cell_sets=both_lake_sets, pheno='ASD')
aad_fgsea = exvivo_fgsea(dexpr_df=gandal_aad, cell_sets=both_lake_sets, pheno='AAD')



# prepare MDD FGSEA results for plotting
mdd_fgsea         = mdd_fgsea[order(mdd_fgsea$pval),]
mdd_fgsea$log10p  = -1*log10(mdd_fgsea$pval)
mdd_fgsea         = mdd_fgsea[order(mdd_fgsea$log10p),]
mdd_fgsea$pathway = factor(mdd_fgsea$pathway, levels=mdd_fgsea$pathway)

p = ggplot(data=mdd_fgsea, aes(x=pathway, y=log10p)) +
            geom_bar(position='dodge', stat='identity') +
            coord_flip() +
            theme_classic() +
            scale_y_continuous(expand=c(0,0), breaks=seq(0,12,4), limits=c(0,12)) +
            geom_hline(yintercept=-1*log10(.05))


fgsea_out = paste0(base_dir, '/figures/gandal_mdd_cellwise_fgsea.pdf')
fgsea_out
CairoPDF(fgsea_out, width=3.5)
print(p)
dev.off()


mdd_fgsea$leadingEdge = NULL
write.csv(as.data.frame(mdd_fgsea), paste0(base_dir, "/supp_data/SuppData_12_gandal_mdd_fgsea.csv"))





