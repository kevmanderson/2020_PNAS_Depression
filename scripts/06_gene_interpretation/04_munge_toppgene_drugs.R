library(tidyverse)
library(XML)
library(dbparser)

# set up directories
base_dir     = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
toppgene_dir = paste0(base_dir, '/output/topp_gene')

# read the enrichments
bin_enrich_list = NULL
for (bin in 1:10){
    toppgene_file = paste0(toppgene_dir, '/toppgene_ahba_mdd_bin', as.character(bin), '_10.txt')
    bin_enrich    = read_delim(toppgene_file, delim='\t')
    bin_enrich_list[[bin]] = bin_enrich
}

# get drug enrichments for top gene decile
toppgene_df   = bin_enrich_list[[1]]
colnames(toppgene_df) = gsub(' |&|-', '_', colnames(toppgene_df))
toppgene_drug = toppgene_df[toppgene_df$Category == 'Drug',]





