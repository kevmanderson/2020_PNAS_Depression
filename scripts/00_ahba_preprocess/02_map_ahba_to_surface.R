library(tidyverse)
library(cifti)
library(gifti)
library(XML)
library(DescTools)
library(Cairo)


# This script will:
#   1. Map AHBA samples to Schaefer cortical parcels (average)
#   2. Similarly, map to Desikan parcels
#   3. Mid-way through the project, the Schaefer atlas was fixed to correct minor
#        errors in the naming convention of individual parcels. We map the old parcel
#        labels to the new ones to help with forward compatibility.


# function to plot parcel-wise metric files for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
    dscalar_template = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/tmp/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotVolOnSurface'
    cmd = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util\');',
                matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}



# Step 1: Read data
# -------------
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read sample information
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)

# read ctx expression data
expr_dat_in = paste0(base_dir, '/data/ahba/ahba_ctx_zWithinSubject_zWithinSample.Rdata')
load(verbose=T, expr_dat_in)

# read sample-to-vertex projection info
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))
sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]


# Step 2: Identify the closest surface vertex for each AHBA sample
# -------------
# cortical sample by normalized gene expression data frame
reg_micro_scale = as.data.frame(t(ctx_data_scale))
match_idxs      = match(ctx_samp_all$well_id, sample_dat$well_id)
reg_samples     = sample_dat[match_idxs[!is.na(match_idxs)],]
which(reg_samples$well_id != ctx_samp_all$well_id)


# stitch the left/right verticies together to match Schaeffer parcel cifti format
reg_samples$bihemi_vertex = reg_samples$vertex + 1 # cifti indices index at 0, R indexes at 1
right_ctx_idxs = intersect(grep('right', reg_samples$structure_name), which(reg_samples$top_level == 'CTX'))
reg_samples$bihemi_vertex[right_ctx_idxs] = reg_samples$bihemi_vertex[right_ctx_idxs] + 32492



# Step 3: Summarise expressino data within each Schaefer parcel
# -------------
# Read Schaeffer parcel info (for multiple parcel #'s, network assignments)
donor_arr =  c('9861','10021','12876','14380','15496','15697')
donor_specific_expression = NULL
for (parcels in c('200')){
    for (net in c('17')){

        # schaeffer parcellation by vertex
        schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
        schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)

        # corresponding parcel labels
        schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
        schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]


        # calculate the average gene expression within each parcel
        schaeffer_mat = matrix(NA, ncol=as.numeric(parcels), nrow=ncol(reg_micro_scale))
        for (donor in donor_arr){
            donor_specific_expression[[donor]] = as.data.frame(schaeffer_mat)
        }

        # loop over each parcel, find matching samples, calculate average expression
        for (idx in 1:length(schaef_labels)){
            write(idx,'')
            parcel_idxs = which(schaef_cii$data == idx) # schaeffer indices
            match_idxs  = which(reg_samples$bihemi_vertex %in% parcel_idxs) # foci within this parcel

            if (length(match_idxs) < 2){
                next
            }
            # expr data for this parcel
            match_samples  = reg_samples[match_idxs,]
            schaeffer_expr = apply(reg_micro_scale[match_idxs,], 2, mean)

            # plug in values to the pre-allocated matrix
            schaeffer_mat[,idx] = schaeffer_expr

            #
            for (donor in donor_arr){
                donor_idxs = which(reg_samples$brain == donor)
                donor_matches = intersect(donor_idxs, match_idxs)
                if (length(donor_matches) < 2){
                    next
                }

                # expr data for this parcel
                donor_match_samples  = reg_samples[donor_matches,]
                donor_schaeffer_expr = apply(reg_micro_scale[donor_matches,], 2, mean)
                donor_specific_expression[[donor]][,idx] = donor_schaeffer_expr
            }
        }
        # add column/row names
        schaef_out = as.data.frame(schaeffer_mat)
        rownames(schaef_out) = colnames(reg_micro_scale)
        colnames(schaef_out) = schaef_labels

        for (donor in donor_arr){
            rownames(donor_specific_expression[[donor]]) = colnames(reg_micro_scale)
            colnames(donor_specific_expression[[donor]]) = schaef_labels
            donor_specific_expression[[donor]]$gene = colnames(reg_micro_scale)

            out_path = paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_donor',as.character(donor),'_',parcels,'_',net,'Net_expr_mat.csv')
            write_csv(donor_specific_expression[[donor]], path=out_path)
        }

        schaef_out$gene = colnames(reg_micro_scale)
        out_path = paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_',parcels,'_',net,'Net_expr_mat.csv')
        write_csv(schaef_out, path=out_path)
    }
}




# Step 4: Make sure the Scheafer parcel labels reflect the latest/"fixed" naming conventions
# -------------
# read the ROI summarized expression data from above
# plot interneuron marker distributions, summarized across scheaffer parcels
schaeffer_mat = read_csv(paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))

# Fix schaefer mapping
info_file   = paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt')
new_mapping = read.csv(header=F, info_file)
new_map_df = data.frame(new_net=as.character(new_mapping$V1[grep('17Networks', new_mapping$V1)]),
                        new_info=as.character(new_mapping$V1[!grepl('17Networks', new_mapping$V1)]), stringsAsFactors=F)


old_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_200Parcels_17Networks_order_info.txt')
old_mapping  = read.csv(header=F, old_info_file)
old_map_df = data.frame(old_net=as.character(old_mapping$V1[grep('17Networks', old_mapping$V1)]),
                        old_info=as.character(old_mapping$V1[!grepl('17Networks', old_mapping$V1)]), stringsAsFactors=F)

both_map_df = cbind(new_map_df, old_map_df)
remap_df = both_map_df[both_map_df$old_net != both_map_df$new_net,]

old_df_cols = colnames(schaeffer_mat)
new_df_cols = old_df_cols
for (replace_row in 1:nrow(remap_df)){
    cur_fix_row  = remap_df[replace_row,]
    cur_fix_idxs = grep(cur_fix_row$old_net, old_df_cols)

    cur_old_names = old_df_cols[cur_fix_idxs]
    cur_new_names = gsub(cur_fix_row$old_net, cur_fix_row$new_net, cur_old_names)

    new_df_cols[cur_fix_idxs] = cur_new_names

    print(old_df_cols[cur_fix_idxs])
    print(new_df_cols[cur_fix_idxs])
    print('')
}
colnames(schaeffer_mat) = new_df_cols

write_csv(schaeffer_mat, paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv'))
write_csv(schaeffer_mat, paste0(base_dir, '/supp_data/SuppData_1_schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv'))


for (donor in donor_arr){
    donor_mat = read_csv(paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_donor',as.character(donor),'_200_17Net_expr_mat.csv'))
    colnames(donor_mat) = new_df_cols
    write_csv(donor_mat, paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_donor',as.character(donor),'_200_17Net_expr_mat_NEWMAP.csv'))
}


# Step 5: Summarise AHBA expression within each Desikan ROI
# -------------
# Desikan atlas
desikan     = paste0(base_dir, '/reference_files/desikan_atlas_32k.dlabel.nii')
desikan_cii = read_cifti(desikan, drop_data = TRUE, trans_data = TRUE)

rh_idxs = ((length(desikan_cii$data)/2)+1):length(desikan_cii$data)
desikan_cii$data[rh_idxs] = desikan_cii$data[rh_idxs]+35

# corresponding parcel labels
desikan_path   = paste0(base_dir, '/reference_files/desikan_atlas_32k.txt')
desikan_labels_in = read_csv(col_names=F, file=desikan_path)
desikan_labels = desikan_labels_in$X1[seq(1, length(desikan_labels_in$X1), by=2)]
desikan_labels = c(paste0('lh_',desikan_labels), paste0('rh_',desikan_labels))

# calculate the average gene expression within each parcel
desikan_mat = matrix(NA, ncol=length(desikan_labels), nrow=ncol(reg_micro_scale))

# calculate the average gene expression within each parcel
donor_desikan_expression = NULL
for (donor in donor_arr){
    donor_desikan_expression[[donor]] = as.data.frame(desikan_mat)
}


for (idx in 1:length(desikan_labels)){
    write(idx,'')
    parcel_idxs   = which(desikan_cii$data == idx) # schaeffer indices
    match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs) # foci within this parcel
    match_samples = reg_samples[match_idxs,] # data for this parcel

    desikan_expr = colMeans(reg_micro_scale[match_idxs,]) # average expressino of every gene, across samples in this parcel
    desikan_mat[,idx] = desikan_expr # plug in values to the pre-allocated matrix

    for (donor in donor_arr){
        donor_idxs = which(reg_samples$brain == donor)
        donor_matches = intersect(donor_idxs, match_idxs)

        # expr data for this parcel
        donor_match_samples  = reg_samples[donor_matches,]
        # average expression of every gene, across samples in this parcel
        desikan_expr = colMeans(reg_micro_scale[donor_matches,])
        donor_desikan_expression[[donor]][,idx] = desikan_expr
    }
}
desikan_mat = as.data.frame(desikan_mat)
colnames(desikan_mat) = desikan_labels
desikan_mat$gene = colnames(reg_micro_scale)


write_csv(desikan_mat, path=paste0(base_dir, '/data/ahba/desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))
write_csv(desikan_mat, paste0(base_dir, '/supp_data/SuppData_2_desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))


for (donor in donor_arr){
    colnames(donor_desikan_expression[[donor]]) = desikan_labels
    donor_desikan_expression[[donor]]$gene = colnames(reg_micro_scale)
    write_csv(donor_desikan_expression[[donor]], path=paste0(base_dir, '/data/ahba/desikan_ahba_ctx_zWithinSubject_zWithinSample_donor',as.character(donor),'_200_17Net_expr_mat.csv'))
}


