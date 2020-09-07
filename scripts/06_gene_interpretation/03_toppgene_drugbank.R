library(tidyverse)
library(seuret)
library(XML)
library(dbparser)
library(gridExtra)
library(RColorBrewer)


source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/06_gene_interpretation/00_support_functions.R')

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


for (bin in 1:10){
print(bin)
toppgene_df = bin_enrich_list[[bin]]
colnames(toppgene_df) = gsub(' |&|-', '_', colnames(toppgene_df))
toppgene_drug = toppgene_df[toppgene_df$Category == 'Drug',]

print(table(toppgene_drug$Source))
}

# make a similarity matrix of each drug in the enrichment table
# based on number of shared transcripts in the Hit list
drug_dist_mat =  make_go_dist(toppgene_drug)


# DrugBank Information
drugbank_path = paste0(base_dir, '/reference_files/full_database.xml')

# parse data from XML and save it to memory
get_xml_db_rows(drugbank_path)

# load drugs data
drugs = parse_drug()

# Give the input file name to the function.
result = xmlParse(file = drugbank_path)

# Exract the root node form the xml file.
rootnode = xmlRoot(result)


for (i in 1:10){
    toppgene_df = bin_enrich_list[[i]]
    colnames(toppgene_df) = gsub(' |&|-', '_', colnames(toppgene_df))
    toppgene_drug = toppgene_df[toppgene_df$Category == 'Drug',]

    toppgene_drug = toppgene_drug[grepl('Drug Bank', toppgene_drug$Source),]

    toppgene_drug$toppgene_clean = unlist(lapply(as.character(toppgene_drug$Name), function(x) strsplit(x,'[[]|;')[[1]][[1]]))
    drugs$name_match = gsub(' ','',drugs$name)
    toppgene_drug$toppgene_match = gsub(' ','', toppgene_drug$toppgene_clean)


    # match toppgene drugs to those in the DrugBank database
    matches       = toppgene_drug[which(tolower(toppgene_drug$toppgene_match) %in% tolower(drugs$name_match)),]
    non_matches   = toppgene_drug[which(!tolower(toppgene_drug$toppgene_match) %in% tolower(drugs$name_match)),]

    drugs_matches = drugs[which(tolower(drugs$name_match) %in% tolower(toppgene_drug$toppgene_match)),]




    # string search to identify pharmaceuticals tagged for MDD from DrugBank
    str_search = unlist(lapply(drugs$description, function(descr) grepl('anti-depressant|antidepressant|Major depressive disorder', descr, ignore.case=T)))
    mdd_total_drugs = drugs$name[str_search]
    mdd_total_drugs_df = drugs[str_search,]

    mdd_matches = NULL
    for (drug in mdd_total_drugs){
        #write(drug,'')
        idx = which(grepl(drug, toppgene_drug$toppgene_clean, ignore.case=T))
        if (length(idx) != 0){
            mdd_matches = c(mdd_matches, drug)
        }
    }
    print(length(mdd_matches))
    mdd_matches
}

# MDD drugs in the top MDD AHBA decile
str_search = unlist(lapply(drugs_matches$description, function(descr) grepl('MDD|Depressive Disorder, Major|anti-depressant|antidepressant|mood|Major depressive disorder', descr, ignore.case=T)))
drugs_matches$name[str_search]
mdd_drugs  = drugs_matches$name[str_search]
mdd_drugs

sort(drug_dist_mat[rownames(drug_dist_mat) == 'nortriptyline',])







