library(tidyverse)
library(Seurat)
library(future)
library(EWCE)
library(MAGMA.Celltyping)
library(reshape2)


# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
umi_dir  = file.path(base_dir, 'data/singlecell/')

# load preprocessed data
region    = 'FrontalCortex'
load(file = paste0(umi_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_batchcorrect.Rdata'), verbose=T)


# read annotation file
lake_annot_df = read_csv(paste0(umi_dir, '/lake2018/lake_annot_df.csv'))

# parallelize diff expr calc
plan("multiprocess", workers = 18)


# collapse cell identities to overarching categories
cell_ident = as.character(seuset@ident)
cell_ident[grep('Ex3', cell_ident)] = 'Ex3'
cell_ident[grep('Ex5', cell_ident)] = 'Ex5'
cell_ident[grep('Ex6', cell_ident)] = 'Ex5'
cell_ident[grep('In1', cell_ident)] = 'In1'
cell_ident[grep('In4', cell_ident)] = 'In4'
cell_ident[grep('In6', cell_ident)] = 'In6'
cell_ident[grep('In7|In8', cell_ident)] = 'SST'

seuset    = SetIdent(seuset, ident.use = cell_ident)
dfc_cells = as.character(unique(seuset@ident))

# calc differential expression for each cell type, relative to all other
dfc_diff_expr_list = NULL
for (cell in dfc_cells){
    print(cell)
    cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
    dfc_diff_expr_list[[cell]] = cell_dexpr
}
diff_expr_save = paste0(base_dir, '/output/diff_expr/lake_primary_cats_frontal_diff_expr_batchcorrect.Rdata')
save(dfc_diff_expr_list, file=diff_expr_save)
load(verbose=T, file=diff_expr_save)


# Expression Weighted Cell Typing (EWCE) cell type specificity
level2 = as.character(seuset@ident)
level2[grep('In|SST', level2)] = 'Inh'
level2[grep('Ex', level2)] = 'Exc'


# drop uninformative genes (major cell cats)
exp_dfc_dropped = drop.uninformative.genes(exp=seuset@data,level2annot=level2)
annotLevels     = list(level1class=as.character(seuset@ident), level2class=level2)

# make "ctd" file with specificity information
setwd(paste0(base_dir, '/output/diff_expr/'))
fNames_CortexOnly = generate.celltype.data(exp=exp_dfc_dropped, annotLevels=annotLevels, groupName="EWCE_lake_dfc_ctd_batch")
print(fNames_CortexOnly)
load(paste0(base_dir, '/output/diff_expr/CellTypeData_EWCE_lake_dfc_ctd.rda'), verbose=T)
dfc_ctd = ctd
dfc_meanexpr = as.data.frame(dfc_ctd[[1]]$mean_exp)
dfc_meanexpr$genes = rownames(dfc_meanexpr)




# sanity check with marker genes
# -----------
spec = as.data.frame(ctd[[1]]$specificity)
spec = spec[rev(order(spec$Ex8)),]

genes = c("SST","PVALB","CNR1",'MOBP','CBLN2','NR4A2')
genes = c('CBLN2','NR4A2','POSTN','CALB2','CALB1')
exp   = melt(ctd[[1]]$mean_exp[genes,],id.vars="genes")
colnames(exp) = c("Gene","Cell","AvgExp")
p = ggplot(data=exp) +
            geom_bar(aes(x=Cell,y=AvgExp),stat="identity")+
            facet_grid(Gene~.)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(file='~/ctd_tmp.pdf', plot=p)






# load preprocessed data
region = 'VisualCortex'
#load(file = paste0(umi_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat.Rdata'), verbose=T)
load(file = paste0(umi_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_batchcorrect.Rdata'), verbose=T)

cell_ident = as.character(seuset@ident)
cell_ident[grep('Ex3', cell_ident)] = 'Ex3'
cell_ident[grep('Ex5', cell_ident)] = 'Ex5'
cell_ident[grep('Ex6', cell_ident)] = 'Ex5'
cell_ident[grep('In1', cell_ident)] = 'In1'
cell_ident[grep('In4', cell_ident)] = 'In4'
cell_ident[grep('In6', cell_ident)] = 'In6'
cell_ident[grep('In7|In8', cell_ident)] = 'SST'

seuset = SetIdent(seuset, ident.use = cell_ident)

vis_cells = as.character(unique(seuset@ident))
vis_diff_expr_list = NULL
for (cell in vis_cells){
    print(cell)
    cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
    vis_diff_expr_list[[cell]] = cell_dexpr
}
diff_expr_save = paste0(base_dir, '/output/diff_expr/lake_primary_cats_visual_diff_expr_batchcorrect.Rdata')
save(vis_diff_expr_list, file=diff_expr_save)


# EWCE cell type specificity
level2 = as.character(seuset@ident)
level2[grep('In|SST', level2)] = 'Inh'
level2[grep('Ex', level2)] = 'Exc'


exp_vis_dropped = drop.uninformative.genes(exp=seuset@data,level2annot=level2)
annotLevels     = list(level1class=as.character(seuset@ident), level2class=level2)


# make "ctd" file with specificity information
setwd('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/output/diff_expr/')
fNames_CortexOnly = generate.celltype.data(exp=exp_vis_dropped, annotLevels=annotLevels, groupName="EWCE_lake_vis_ctd_batch")
load('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/output/diff_expr/CellTypeData_EWCE_lake_vis_ctd_batch.rda', verbose=T)


vis_ctd = ctd

# sanity check with marker genes
# -----------
spec = as.data.frame(vis_ctd[[1]]$specificity)
spec = spec[rev(order(spec$Ex8)),]

genes = c("SST","PVALB","CNR1",'MOBP','CBLN2','NR4A2')
genes = c('CBLN2','NR4A2','POSTN','FAP','CALB1')
exp = melt(ctd[[1]]$mean_exp[genes,],id.vars="genes")
colnames(exp) = c("Gene","Cell","AvgExp")
p = ggplot(exp) +
    geom_bar(aes(x=Cell,y=AvgExp),stat="identity")+
    facet_grid(Gene~.)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
















