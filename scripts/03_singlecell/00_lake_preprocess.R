library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
library(GEOquery)
library(EWCE)
#library(MAGMA.Celltyping)


# set up directories
base_dir    = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
dropseq_dir = paste0(base_dir, '/data/singlecell')


# read annotation
GSE97930_soft = getGEO(filename=file.path(dropseq_dir, 'lake2018/GSE97930_family.soft.gz'))
# cell metadata
GSE97930_miseq = getGEO(filename=file.path(dropseq_dir, 'lake2018/GSE97930-GPL15520_series_matrix.txt'),GSEMatrix=TRUE)
GSE97930_miseq@phenoData@data
GSE97930_hiseq = getGEO(filename=file.path(dropseq_dir, 'lake2018/GSE97930-GPL16791_series_matrix.txt'),GSEMatrix=TRUE)
GSE97930_hiseq@phenoData@data


# re-format annotation dataframe
annot_df = data.frame()
for (id in names(GSE97930_soft@gsms)){

  sample_row = as.data.frame(GSE97930_soft@gsms[[id]]@header[c('geo_accession','title','platform_id','instrument_model')])
  keys = unlist(lapply(GSE97930_soft@gsms[[id]]@header$characteristics_ch1, function(x) strsplit(x,': ')[[1]][1]))
  vals = unlist(lapply(GSE97930_soft@gsms[[id]]@header$characteristics_ch1, function(x) strsplit(x,': ')[[1]][2]))

  sample_row[['age']]          = as.numeric(gsub('years','', vals[keys=='age']))
  sample_row[['race']]         = vals[keys=='race']
  sample_row[['brain_region']] = vals[grep('region',keys)]
  sample_row[['brodmann']]     = vals[keys=='brodmann area']
  sample_row[['tissue']]       = vals[keys=='tissue']
  sample_row[['experiment']]   = vals[keys=='experiment']
  sample_row[['patient']]      = vals[grep('patient',keys)]
  annot_df = rbind(annot_df, sample_row)
}
annot_df$title = gsub('Fcx', 'fcx', annot_df$title)
annot_df$title = gsub('Occ', 'occ', annot_df$title)
annot_df$title = gsub('Cbm', 'cbm', annot_df$title)


# save sample annotation data in dataframe format
write_csv(x=annot_df, paste0(dropseq_dir, '/lake2018/lake_annot_df.csv'))
annot_df = read_csv(paste0(dropseq_dir, '/lake2018/lake_annot_df.csv'))


# Perform Seurat style preprocessing on Lake data
for (region in c('FrontalCortex', 'VisualCortex'){

    # read UMI counts
    umi_in  = read.table(file.path(dropseq_dir, paste0('lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt')))
    umi_mat = as.matrix(umi_in)

    # extract cell types and create cell metadata df
    cell_type     = unlist(lapply(colnames(umi_in), function(x) strsplit(x,'_')[[1]][[1]]))
    sample_id     = unlist(lapply(colnames(umi_in), function(x) strsplit(x,'_')[[1]][[2]]))
    sample_df     = data.frame(sample_id)
    sample_df$idx = 1:nrow(sample_df)
    col_df        = data.frame(cell_type=cell_type, sample_id=sample_id)

    rownames(sample_df)    = colnames(umi_in)
    metadata_df            = merge(x=sample_df, y=annot_df, all.x=T, by.x='sample_id', by.y='title')
    metadata_df            = metadata_df[order(metadata_df$idx),]
    metadata_df$experiment = as.factor(metadata_df$experiment)
    metadata_df$patient    = as.factor(metadata_df$patient)
    metadata_df$instrument_model = as.factor(metadata_df$instrument_model)
    rownames(metadata_df)  = rownames(sample_df)
    length(which(metadata_df$sample_id != sample_id))

    # make an seuset expression object
    seuset = CreateSeuratObject(
        raw.data = umi_in,
        min.cells = 3,
        min.genes = 200,
        meta.data = metadata_df
    )
    dim(seuset@data)

    # cells with very clear outlier counts
    seuset = FilterCells(
        object = seuset,
        subset.names = c("nUMI"),
        high.thresholds = c(2e7)
    )
    dim(seuset@data)

    # global-scaling normalization method, x10000
    seuset = NormalizeData(
        object = seuset,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )

    # regress number of detected genes per cell
    if (region == 'FrontalCortex'){
        # regress batch cov
        seuset = ScaleData(
            object = seuset,
            vars.to.regress = c('nUMI','experiment','instrument_model')
        )
    else {
        seuset = ScaleData(
            object = seuset,
            vars.to.regress = c('nUMI','experiment')
        )
    }

    # collapse subtypes from particular cluster
    col_df$cell_type_super = as.character(col_df$cell_type)

    seuset = SetIdent(object=seuset, ident.use=col_df$cell_type_super)
    save(seuset, col_df, file = paste0(dropseq_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_batchcorrect.Rdata'))
}

# Next we want to align the gene names as much as possible between the AHBA and Lake datasets

# read human genome reference
tmp_dir  = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ncbi_hgnc = read_delim(paste0(tmp_dir, '/ref_files/Homo_sapiens.gene_info'), delim='\t')


# create ensembl id array (useful for cross-build gene matching)
ens_array = unlist(lapply(ncbi_hgnc$dbXrefs, function(x) strsplit(x,'Ensembl:')[[1]][2]))
ncbi_hgnc$ensemblIDs = ens_array


# AHBA data
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)


# get rid of a few genes with duplicate ens/entrez maps
ncbi_hgnc = ncbi_hgnc[!duplicated(ncbi_hgnc$Symbol),]


# read Seurat data
region = 'FrontalCortex'
load(verbose=T, file = paste0(dropseq_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_batchcorrect.Rdata'))
lake_dfc_seuset = seuset


# Gene symbols for Lake dataset
lake_gene_names  = rownames(lake_dfc_seuset@data)
lake_gene_df     = data.frame(lake_genes=lake_gene_names)
lake_gene_df$idx = 1:nrow(lake_gene_df)
ncbi_hgnc$Symbol_merge = ncbi_hgnc$Symbol


# merge NCBI gene info with lake gene names
lake_gene_df2 = merge(x=lake_gene_df, y=ncbi_hgnc, by.x='lake_genes', by.y='Symbol_merge', all.x=T, mult = "first", nomatch=0L)
lake_gene_df2 = lake_gene_df2[order(lake_gene_df2$idx),]
lake_gene_df2 = lake_gene_df2[which(!duplicated(lake_gene_df2$lake_genes)),]


# genes that dont match the primary gene symbol in the database
nonmatch_names = lake_gene_names[which(!lake_gene_names %in% ncbi_hgnc$Symbol)]


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
keep_cols     = c('lake_genes','idx','GeneID','dbXrefs','Synonyms','ensemblIDs','Symbol_from_nomenclature_authority','type_of_gene')
lake_gene_df_tomerge = lake_gene_df2[keep_cols]
lake_gene_df2 = merge(x=lake_gene_df_tomerge, y=ncbi_hgnc_v2, by.x='lake_genes', by.y='Symbol_v2', all.x=T)
lake_gene_df2 = lake_gene_df2[order(lake_gene_df2$idx),]

# merge entrez id
lake_gene_df2$GeneID = lake_gene_df2$GeneID.x
lake_gene_df2$GeneID[is.na(lake_gene_df2$GeneID.x)] = lake_gene_df2$GeneID.y[is.na(lake_gene_df2$GeneID.x)]
# merge ensembl ids
lake_gene_df2$ensemblIDs = lake_gene_df2$ensemblIDs.x
lake_gene_df2$ensemblIDs[is.na(lake_gene_df2$ensemblIDs.x)] = lake_gene_df2$ensemblIDs.y[is.na(lake_gene_df2$ensemblIDs.x)]

# merge with probes
lake_gene_dict = data.frame(lake_orig_Symbol = lake_gene_df2$lake_genes,
                            ensembl_id = lake_gene_df2$ensemblIDs,
                            entrez_id = lake_gene_df2$GeneID,
                            idx=lake_gene_df2$idx)
lake_dfc_gene_dict = lake_gene_dict[order(lake_gene_dict$idx),]
length(which(lake_dfc_gene_dict$lake_orig_Symbol == rownames(lake_dfc_seuset@data)))
save(lake_dfc_seuset, lake_dfc_gene_dict, file = paste0(dropseq_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed_bathcorrect.Rdata'))






region = 'VisualCortex'
load(verbose=T, file = paste0(dropseq_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_batchcorrect.Rdata'))
lake_vis_seuset = seuset

# Gene symbols for Lake dataset
lake_gene_names  = rownames(lake_vis_seuset@data)
lake_gene_df     = data.frame(lake_genes = lake_gene_names)
lake_gene_df$idx = 1:nrow(lake_gene_df)
ncbi_hgnc$Symbol_merge = ncbi_hgnc$Symbol

# merge NCBI gene info with lake gene names
lake_gene_df2 = merge(x=lake_gene_df, y=ncbi_hgnc, by.x='lake_genes', by.y='Symbol_merge', all.x=T, mult = "first", nomatch=0L)
# put the lake genes in their original order
lake_gene_df2 = lake_gene_df2[order(lake_gene_df2$idx),]
# remove genes with duplicated names
lake_gene_df2 = lake_gene_df2[which(!duplicated(lake_gene_df2$lake_genes)),]

# genes that dont match the primary gene symbol in the database
nonmatch_names = lake_gene_names[which(!lake_gene_names %in% ncbi_hgnc$Symbol)]

# for the non-matches, check the gene name aliases
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
keep_cols = c('lake_genes','idx','GeneID','dbXrefs','Synonyms','ensemblIDs','Symbol_from_nomenclature_authority','type_of_gene')
lake_gene_df_tomerge = lake_gene_df2[keep_cols]
lake_gene_df2 = merge(x=lake_gene_df_tomerge, y=ncbi_hgnc_v2, by.x='lake_genes', by.y='Symbol_v2', all.x=T)
lake_gene_df2 = lake_gene_df2[order(lake_gene_df2$idx),]

# merge entrez id
lake_gene_df2$GeneID = lake_gene_df2$GeneID.x
lake_gene_df2$GeneID[is.na(lake_gene_df2$GeneID.x)] = lake_gene_df2$GeneID.y[is.na(lake_gene_df2$GeneID.x)]
# merge ensembl ids
lake_gene_df2$ensemblIDs = lake_gene_df2$ensemblIDs.x
lake_gene_df2$ensemblIDs[is.na(lake_gene_df2$ensemblIDs.x)] = lake_gene_df2$ensemblIDs.y[is.na(lake_gene_df2$ensemblIDs.x)]

# merge with probes
lake_vis_gene_dict = data.frame(lake_orig_Symbol = lake_gene_df2$lake_genes,
                            ensembl_id = lake_gene_df2$ensemblIDs,
                            entrez_id = lake_gene_df2$GeneID,
                            idx=lake_gene_df2$idx)
lake_vis_gene_dict = lake_vis_gene_dict[order(lake_vis_gene_dict$idx),]
out_file =  paste0(dropseq_dir, '/lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed_batchcorrect.Rdata')
out_file
save(lake_vis_seuset, lake_vis_gene_dict, file = out_file)


