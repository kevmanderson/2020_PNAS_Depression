library(tidyverse)
library(WGCNA)
library(sva)


# This script will:
# 1. Read/Process AHBA microarray data
#    a. remove probes without entrez id
#    b. remove probes that are not expressed in more than 20% of cortical samples, on avarage across donors
#    c. select probe with highest mean expression
#    d. remove probes that are more than 4mm from closest surface vertex
#    e. z-transform expression matrix


# this function will map AHBA samples to their overarching regional category (e.g. mPFC >> CORTEX)
find_top_level = function(ref_df, struct_id, ontology, top_level) {
    ontology_row  = which(ontology$id == struct_id )
    if (length(ontology_row) > 1){ # print feedback if more than one match, this should never be the case
        write("ERROR")
    }
    ontology_info = ontology[ontology_row,]
    splits    = strsplit(as.character(ontology_info$structure_id_path), '/')[[1]]
    reg.match = ref_df$id[which(ref_df$id %in% splits)]
    if (length(reg.match) == 1){
        region_cat   = as.character(ref_df$top_level[which(ref_df$id %in% splits)])
        region_name  = ref_df[ref_df$id == reg.match,]$name
        return(c(reg.match, region_cat, region_name))
    } else {
        return(c(NA,NA,NA))
    }
}



# Step 1: Set up necessary paths
# -----------------
# Modify these filepaths for your local directory structure
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# AHBA subject list/directory
data_path  = paste(base_dir, '/data/ahba', sep = '') # path to AHBA microarray data
donor_nums = c('9861', '10021', '12876', '14380', '15496', '15697')

# Map ontology structure names to more general sample categories
# e.g. left motor gyrus --> CTX
top_level = read_csv(paste0(base_dir, '/reference_files/top_level_categories.csv'))



# Step 2: Read AHBA data from each donor
# -----------------
# Read data from each subject
ctx_data_micro   = matrix()
ctx_data_pa      = matrix()
ctx_data_sample  = NULL
ahba_data  = NULL
for ( donor in donor_nums ) {
    file = paste('donor', donor, sep='')
    write(paste('Reading and collapsing data for donor: ', donor, sep = ''),'')

    # create donor field
    ahba_data[[donor]] = NULL


    # read sample Information
    saname    = paste(data_path, file, 'SampleAnnot.csv', sep='/')
    samp_info = read_csv(saname)
    ahba_data[[donor]]$raw_samp = samp_info
    samp_info$brain = donor


    # Gene Ontology Info
    oname    = paste(data_path, file, 'Ontology.csv', sep='/')
    ont_data = read_csv(oname)
    ahba_data[[donor]]$raw_ont = ont_data
    ontology = ont_data

    # identify samples in the regions we want to analyze
    out         = do.call(rbind, lapply(samp_info$structure_id, find_top_level, ref_df=top_level, ontology=ontology))
    reg_info    = as.data.frame(out)
    reg_info$V1 = as.numeric(as.character(reg_info$V1))
    colnames(reg_info) = c('reg_num', 'top_level', 'region_clean')

    # Read Microarray Data
    fname                = paste(data_path, file, 'MicroarrayExpression.csv', sep='/')
    microdata            = read_csv(fname, col_names = F)
    micro_arr            = microdata

    # remove first column contains probe IDs
    micro_temp           = as.data.frame(microdata[,-which(colnames(microdata) == 'X1')]) # First column contains probe IDs
    rownames(micro_temp) = micro_arr$X1 # collapseRows requires rownames
    micro_df             = as.data.frame(micro_temp)
    #ahba_data[[donor]]$raw_micros = micro_df

    # Information about each Gene Probe (~50,000 probes for ~20,000 genes)
    pname     = paste(data_path, file, 'Probes.csv', sep='/')
    probes    = read_csv(pname)
    ahba_data[[donor]]$raw_probes = probes
    probeInfo = probes

    # discard probes without an entrez-id
    num_samples   = dim(ahba_data[[donor]]$raw_micros)[2]
    trash_me      = is.na(ahba_data[[donor]]$raw_probes$entrez_id) # & pasums > cutoff
    good_probes   = ahba_data[[donor]]$raw_probes[trash_me == FALSE,]

    # PA
    write('PA','')
    pname    = paste(data_path, file, 'PACall.csv', sep='/')
    pa_data  = read.csv(pname, header=F)
    rownames(pa_data) = pa_data$V1
    pa_data$V1 = NULL
    #ahba_data[[donor]]$raw_pas = pa_data


    # Select the probes with valid Entrez IDs
    ahba_data[[donor]]$probes_filter = good_probes
    micro_filter = micro_df[rownames(micro_df) %in% good_probes$probe_id,]
    ahba_data[[donor]]$micro_filter  = micro_filter
    pa_filter = pa_data[rownames(pa_data) %in% good_probes$probe_id,]
    ahba_data[[donor]]$pas_filter    = pa_filter


    ctx_samples = which(reg_info$top_level == 'CTX')
    n_samples   = length(ctx_samples)

    ctx_data_micro = cbind(ctx_data_micro, micro_filter[,ctx_samples])
    ctx_data_pa    = cbind(ctx_data_pa, pa_filter[,ctx_samples])
    ctx_data_sample = rbind(ctx_data_sample, samp_info[ctx_samples,])
}
ctx_data_micro$ctx_data_micro = NULL
ctx_data_pa$ctx_data_pa = NULL



# Step 3: Remove probes that are not strongly expressed above baseline in cortex
# -----------------
# probe-wise noise per donor
probe_noise = NULL
for ( donor in donor_nums ) {
    donor_idxs        = which(ctx_data_sample$brain == donor)
    donor_probe_noise = rowSums(ctx_data_pa[,donor_idxs])/length(donor_idxs)
    probe_noise       = cbind(probe_noise, donor_probe_noise)
}
avg_probe_noise = rowMeans(probe_noise)
head(avg_probe_noise)

# filter based on probes that, on average, are not expressed much in cortex
ctx_data_micro_EXPRESSED = ctx_data_micro[which(avg_probe_noise > .2),]
probes_EXPRESSED = good_probes[which(avg_probe_noise > .2),]



# Step 4: Collapse probes based on highest mean expression
# -----------------
rownames(ctx_data_micro_EXPRESSED) = as.character(probes_EXPRESSED$probe_name)
datExpr           = collapseRows(datET=ctx_data_micro_EXPRESSED, rowGroup=as.character(probes_EXPRESSED$gene_symbol), rowID=as.character(probes_EXPRESSED$probe_name), method="MaxMean", connectivityBasedCollapsing=FALSE)
donorDat = NULL
donorDat$micro    = datExpr$datETcollapsed
donorDat$samples  = ctx_data_sample
donorDat$ontology = ontology
probes2keep       = as.character(datExpr$group2row[,2])
probesKept        = probes_EXPRESSED[probes_EXPRESSED$probe_name %in% probes2keep,]
donorDat$probes   = probesKept[order(probesKept$gene_symbol),]
length(which(donorDat$probes$gene_symbol == rownames(datExpr$datETcollapsed)))


# For quicker loading after the above preprocessing has been run, save Rdata structure
save(x=donorDat, file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'))
load(verbose=T, paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'))




# Step 5: Within/Between samples normalization of expression
# -----------------
# sample information file with vertex mapping
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))

# remove sample greater than 4mm away from nearest vertex
usable_ctx_samples = sample_info[which(abs(sample_info$mm_to_surf) < 4),]
donor_nums     = c('9861', '10021', '12876', '14380', '15496', '15697')

# BOTH WITHIN-SAMPLE and ACROSS-SAMPLE z-transform cortical data separately for each donor
micro_sample_scale = scale(donorDat$micro)
ctx_data_scale = matrix()
ctx_samp_all   = NULL
for (donor in donor_nums){
    print(paste0('Mean and variance normalizing cortical data for : ', donor))

    ctx_samples = usable_ctx_samples %>% filter(brain == donor)
    ctx_idxs    = which(donorDat$samples$well_id %in% ctx_samples$well_id)
    ctx_micro   = micro_sample_scale[,ctx_idxs]
    ctx_samp    = donorDat$samples[ctx_idxs,]

    ctx_micro_scale = as.data.frame(t(apply(ctx_micro, 1, scale, center=TRUE, scale=TRUE)))
    rownames(ctx_micro_scale) = donorDat$probes$gene_symbol

    ctx_data_scale = cbind(ctx_data_scale, ctx_micro_scale)
    ctx_samp_all = rbind(ctx_samp_all, ctx_samp)
}
ctx_data_scale$ctx_data_scale = NULL
ctx_samp_all$ctx_samp_all     = NULL
out_file = paste0(base_dir, '/data/ahba/ahba_ctx_zWithinSubject_zWithinSample.Rdata')
save(ctx_data_scale, ctx_samp_all,file=out_file)








