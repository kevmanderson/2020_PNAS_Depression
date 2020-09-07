library(tidyverse)
library(seuret)
library(XML)
library(dbparser)
library(gridExtra)
library(RColorBrewer)


source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/06_gene_interpretation/00_support_functions.R')

# set up directories
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



# read sample information
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)
colnames(ctx_data_scale) = paste0('wellid_', ctx_samp_all$well_id)
probes = read_csv(paste0(base_dir, '/data/ahba/donor10021/Probes.csv'))


# create ordered gene partitions for enrichment analysis
nbins   = 20
n_genes = length(all_ahba_cors$genes)
bin_len = round(n_genes/nbins)
partition_array = sort(rep(1:nbins, bin_len))
if (length(partition_array) < n_genes){
    partition_array = c(partition_array, rep(nbins, n_genes - length(partition_array)))
} else {
    partition_array = partition_array[1:n_genes]
}

toppgene_dir = paste0(base_dir, '/output/topp_gene')
enrich_list  = NULL
for (bin in 1:nbins){
    write(bin,'')

    bin_genes  = all_ahba_cors$genes[which(partition_array == bin)]
    print(length(bin_genes))

    bin_entrez = unique(donorDat$probes$entrez_id[donorDat$probes$gene_symbol %in% bin_genes])
    print(length(bin_entrez))

    enrich_list[[bin]] = as.character(bin_entrez)
}

enrich_df  = NULL
max_length = max(unlist(lapply(enrich_list, length)))
for (bin in 1:nbins){
    add_col   = c(enrich_list[[bin]], rep('NA', max_length - length(enrich_list[[bin]])))
    enrich_df = cbind(enrich_df, add_col)
}
enrich_df = as.data.frame(enrich_df)
out_file  = paste0(toppgene_dir, 'ahba_ukbb_order_bin_df_ENTREZ.csv')
write_csv(enrich_df, out_file)
out_file


#

# Submit deciles to ToppGene online, outside of this script

#


