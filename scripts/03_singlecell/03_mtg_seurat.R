library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(future)
library(feather)
library(matrixStats)
library(dplyr)
library(edgeR)
library(data.table)
library(EWCE)
library(MAGMA.Celltyping)

options(stringsAsFactors=FALSE)



####################
# set up directories
####################
base_dir    = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
dropseq_dir = paste0(base_dir, '/data/singlecell')
homeFolder  = paste0(base_dir,"/external/MTG_celltypes")


# load some functions
source(paste0("/gpfs/milgram/project/holmes/kma52/mdd_sst/external/MTG_celltypes/R/Support_code0_functions.R"))



##################
# Read MTG dropseq
##################
anno      = read_csv(paste0(dropseq_dir, '/ahba_mtg/anno.csv'))
anno$sentinal = unlist(lapply(anno$cluster_label, function(x) strsplit(x, ' ')[[1]][[3]]))
Expr_data = read_csv(paste0(dropseq_dir, '/ahba_mtg/data.csv'))

# reformat MTG data
datIn      = as.matrix(Expr_data[,names(Expr_data) != "sample_id"])
rownames(datIn) = Expr_data$sample_id
datIn      = t(datIn)
rm(Expr_data)


# make an seuset expression object
seuset = CreateSeuratObject(
    counts = datIn,
    min.cells = 3,
    min.features = 200,
    assay = "RNA"
)
dim(seuset@data)

# cells with very clear outlier counts
seuset = FilterCells(
    object = seuset,
    subset.names = c("nGene"),
    high.thresholds = c(2e7)
)

# global-scaling normalization method, x10000
seuset = NormalizeData(
    object = seuset,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)
dim(seuset@data)


# regress number of detected genes per cell
seuset = ScaleData(
    object = seuset
)


# information about cells
anno = as.data.frame(read_csv(paste0(dropseq_dir, '/ahba_mtg/anno.csv')))

# process stringes to get rid of layer information/secondary marker info
str1 = unlist(lapply(anno$cluster_label, function(x) strsplit(x,' ')[[1]][[1]]))
str2 = unlist(lapply(anno$cluster_label, function(x) strsplit(x,' ')[[1]][[2]]))
str3 = unlist(lapply(anno$cluster_label, function(x) strsplit(x,' ')[[1]][[3]]))
str4 = unlist(lapply(anno$cluster_label, function(x) strsplit(x,' ')[[1]][[4]]))

anno$top_cluster_label = str1 #anno$cluster_label
anno$sentinal_marker = str3 #anno$cluster_label


# process a bit/gsub
replace_me  = names(table(paste0(' ', str2)))
replace_me1 = replace_me[grep('-', replace_me)]
for (var in replace_me1){ anno$top_cluster_label = gsub(var,'',anno$top_cluster_label) }

replace_me2 = replace_me[!grepl('-', replace_me)]
for (var in replace_me2){ anno$top_cluster_label = gsub(var,'',anno$top_cluster_label) }
anno$top_cluster_label = gsub(' ', '_', anno$top_cluster_label)


# place new cell categories into cell metadata
anno$cell_type_super
anno$cell_type_super = as.character(paste0(str1,'_',str3))

Idents(seuset) = anno$cell_type_super

# loop over each cell type
cell_types = as.character(unique(Idents(seuset)))
mtg_diff_expr_list = NULL
for (cell in cell_types){
    print(cell)
    cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
    mtg_diff_expr_list[[cell]] = cell_dexpr
}



#seuset = SetIdent(object=seuset, ident.use=anno$cell_type_super)
save(seuset, anno, file = paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_bigcat_singlecell_seurat.Rdata'))
# load previously processed data
load(paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_singlecell_seurat.Rdata'), verbose=T)

# parallelize diff expr calc
plan("multiprocess")


# loop over each cell type
cell_types = as.character(unique(seuset@ident))
mtg_diff_expr_list = NULL
for (cell in cell_types){
    print(cell)
    cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
    mtg_diff_expr_list[[cell]] = cell_dexpr
}
save(mtg_diff_expr_list, file=paste0(base_dir, '/output/diff_expr/mtg_big_cat_ahba_diff_expr.Rdata'))





seuset = SetIdent(object=seuset, ident.use=anno$top_cluster_label)
save(seuset, anno, file = paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_bigcat_singlecell_seurat.Rdata'))
# load previously processed data
load(paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_singlecell_seurat.Rdata'), verbose=T)

# parallelize diff expr calc
plan("multiprocess")


# loop over each cell type
cell_types = as.character(unique(seuset@ident))
mtg_diff_expr_list = NULL
for (cell in cell_types){
    print(cell)
    cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
    mtg_diff_expr_list[[cell]] = cell_dexpr
}
save(mtg_diff_expr_list, file=paste0(base_dir, '/output/diff_expr/mtg_big_cat_ahba_diff_expr.Rdata'))


# the same as above, except do it for the fine grained 75 cell categorizations
anno$cluster_fine_label = gsub(' |-', '_', anno$cluster_label)
table(anno$cluster_fine_label)
seuset = SetIdent(object=seuset, ident.use=anno$cluster_fine_label)
save(seuset, anno, file = paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_singlecell_finelevel_seurat.Rdata'))
load(paste0(dropseq_dir, '/ahba_mtg/ahba_mtg_singlecell_finelevel_seurat.Rdata'), verbose=T)

# parallelize diff expr calc
plan("multiprocess")

cell_types = as.character(unique(seuset@ident))
mtg_fine_diff_expr_list = NULL
for (cell in cell_types){
    print(cell)
    out_path   = paste0(base_dir, '/output/diff_expr/mtg_cell_', cell, '_cell_fine_diff_expr_list.Rdata')
    if ( file.exists(out_path) == F ) {
        cell_dexpr = FindMarkers(seuset, ident.1 = cell, ident.2 = NULL)
        save(x=cell_dexpr, file=out_path)
    }
    #mtg_fine_diff_expr_list[[cell]] = cell_dexpr
}
save(mtg_fine_diff_expr_list, file=paste0(base_dir, '/output/diff_expr/mtg_fine_diff_expr_list.Rdata'))


mtg_fine_diff_expr_list = NULL
for (cell in cell_types){
    print(cell)
    in_path  = paste0(base_dir, '/output/diff_expr/mtg_cell_', cell, '_cell_fine_diff_expr_list.Rdata')
    load(in_path, verbose=T)
    mtg_fine_diff_expr_list[[cell]] = cell_dexpr
}
out_path = paste0(base_dir, '/output/diff_expr/mtg_cell_ALL_cell_fine_diff_expr_list.Rdata')
save(x=mtg_fine_diff_expr_list, file=out_path)



# EWCE cell type specificity
level2 = as.character(seuset@ident)
level2[grep('In', level2)] = 'Inh'
level2[grep('Ex', level2)] = 'Exc'


# drop uninformative genes (major cell cats)
exp_mtg_dropped = drop.uninformative.genes(exp=seuset@data, level2annot=level2)
annotLevels     = list(level1class=as.character(anno$cell_type_super), level2class=level2)


# make "ctd" file with specificity information
setwd('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/output/diff_expr/')
fNames_CortexOnly = generate.celltype.data(exp=exp_mtg_dropped, annotLevels=annotLevels, groupName="EWCE_aba_mtg_ctd")
print(fNames_CortexOnly)

load('/gpfs/milgram/project/holmes/kma52/mdd_sst/output/diff_expr/CellTypeData_EWCE_aba_mtg_ctd.rda', verbose=T)


# calculate cell type specificity / quantiles
ctd = prepare.quantile.groups(ctd, specificity_species="human",numberOfBins=40)



# reformat gwas summary stats
#gwas_sumstats_path = '/gpfs/milgram/project/holmes/kma52/mdd_sst/data/ldsc_sc/Wray2018_PGC_UKB_depression_genome_wide.txt'
#gwas_sumstats_path = '/gpfs/milgram/project/holmes/kma52/mdd_sst/data/gwas/magma_Howard_PGC_UKB_depression_genome.txt'
#col_headers = format_sumstats_for_magma(gwas_sumstats_path)

genome_ref_path = '/gpfs/milgram/project/holmes/kma52/mdd_sst/external/magma/g1000_eur'
#genesOutPath = map.snps.to.genes(gwas_sumstats_path, genome_ref_path=genome_ref_path)

ctd[[1]]$quantiles = ctd[[1]]$specificity_quantiles
ctd[[2]]$quantiles = ctd[[2]]$specificity_quantiles
gwas_sumstats_path = '/gpfs/milgram/project/holmes/kma52/mdd_sst/data/gwas/magma_Howard_PGC_UKB_depression_genome.txt'
ctAssocsLinear = calculate_celltype_associations(ctd, analysis_name = "aba_mtg_MainRun", gwas_sumstats_path, genome_ref_path=genome_ref_path, specificity_species = "human")
#FigsLinear = plot_celltype_associations(ctAssocsLinear,ctd=ctd)


ctAssocsTop = calculate_celltype_associations(ctd,analysis_name = "aba_mtg_top10", gwas_sumstats_path,genome_ref_path=genome_ref_path,EnrichmentMode = "Top 10%", specificity_species = "human")
#FigsTopDecile = plot_celltype_associations(ctAssocsTop,ctd=ctd)





# read human genome reference
base_dir  = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ncbi_hgnc = read_delim(paste0(base_dir, '/ref_files/Homo_sapiens.gene_info'), delim='\t')


# create ensembl id array (useful for cross-build gene matching)
ens_array = unlist(lapply(ncbi_hgnc$dbXrefs, function(x) strsplit(x,'Ensembl:')[[1]][2]))
ncbi_hgnc$ensemblIDs = ens_array


# AHBA data
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)



# get rid of a few genes with duplicate ens/entrez maps
ncbi_hgnc = ncbi_hgnc[!duplicated(ncbi_hgnc$Symbol),]


ahba_mtg_seuset = seuset

# Gene symbols for Lake dataset
mtg_gene_names  = rownames(ahba_mtg_seuset@data)
mtg_gene_df     = data.frame(mtg_genes=mtg_gene_names)
mtg_gene_df$idx = 1:nrow(mtg_gene_df)
ncbi_hgnc$Symbol_merge = ncbi_hgnc$Symbol


# merge NCBI gene info with lake gene names
mtg_gene_df2 = merge(x=mtg_gene_df, y=ncbi_hgnc, by.x='mtg_genes', by.y='Symbol_merge', all.x=T, mult = "first", nomatch=0L)
mtg_gene_df2 = mtg_gene_df2[order(mtg_gene_df2$idx),]
mtg_gene_df2 = mtg_gene_df2[which(!duplicated(mtg_gene_df2$mtg_genes)),]


# genes that dont match the primary gene symbol in the database
nonmatch_names = mtg_gene_names[which(!mtg_gene_names %in% ncbi_hgnc$Symbol)]

# for the non-matches, check the gene aliases
ncbi_hgnc_v2 = ncbi_hgnc
ncbi_hgnc_v2$Symbol_v2 = NA

# prepare the synonym references array to match the genes to their respective entrez/ensembl/etc
synonym_ref = ncbi_hgnc$Synonyms
synonym_ref = gsub('[|]','_',synonym_ref)
synonym_ref = lapply(synonym_ref, function(x) paste0('_',x,'_'))
ct = 1
# for each gene that did not match intially, see if it matches a gene alias
for (s in unique(nonmatch_names)){
    # gene that didnt match
    s_to_match = paste0('_',s,'_') # the underscores help to prevent against substring matches
    # string match to the gene aliases
    idx = grep(s_to_match, synonym_ref)
    # if there is any match
    if (length(idx) == 1){
        ncbi_hgnc_v2$Symbol_v2[idx] = s
        print(s_to_match)
    }
    ct = ct + 1
}


# merge human genome build info with the lake gene info dataframe
keep_cols = c('mtg_genes','idx','GeneID','dbXrefs','Synonyms','ensemblIDs','Symbol_from_nomenclature_authority','type_of_gene')
mtg_gene_df_tomerge = mtg_gene_df2[keep_cols]
mtg_gene_df2 = merge(x=mtg_gene_df_tomerge, y=ncbi_hgnc_v2, by.x='mtg_genes', by.y='Symbol_v2', all.x=T)
mtg_gene_df2 = mtg_gene_df2[order(mtg_gene_df2$idx),]

# merge entrez id
mtg_gene_df2$GeneID = mtg_gene_df2$GeneID.x
mtg_gene_df2$GeneID[is.na(mtg_gene_df2$GeneID.x)] = mtg_gene_df2$GeneID.y[is.na(mtg_gene_df2$GeneID.x)]
# merge ensembl ids
mtg_gene_df2$ensemblIDs = mtg_gene_df2$ensemblIDs.x
mtg_gene_df2$ensemblIDs[is.na(mtg_gene_df2$ensemblIDs.x)] = mtg_gene_df2$ensemblIDs.y[is.na(mtg_gene_df2$ensemblIDs.x)]

# merge with probes
mtg_gene_dict = data.frame(mtg_orig_Symbol = mtg_gene_df2$mtg_genes,
                            ensembl_id = mtg_gene_df2$ensemblIDs,
                            entrez_id = mtg_gene_df2$GeneID,
                            idx=mtg_gene_df2$idx)
mtg_gene_dict = mtg_gene_dict[order(mtg_gene_dict$idx),]
length(which(mtg_gene_dict$mtg_orig_Symbol == rownames(seuset@data)))
mtg_seuset = seuset
save(mtg_seuset, mtg_gene_dict, file = paste0(dropseq_dir, '/ahba_mtg/mtg_singlecell_ahba_match_seurat_processed.Rdata'))
























