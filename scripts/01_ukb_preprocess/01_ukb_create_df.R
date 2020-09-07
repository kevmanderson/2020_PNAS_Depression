library(feather)

# This script will combined raw UKB data into one large dataframe


# Step 1: Read individual UKB data buckets
# ----------------
# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read raw data (separated in different UKB bucket releases)
source(paste0(base_dir, '/data/ukb/enc/ukb26541.r'))
bd_26541 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb27118.r'))
bd_27118 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb27119.r'))
bd_27119 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb27121.r'))
bd_27121 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb27123.r'))
bd_27123 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb27124.r'))
bd_27124 = bd

source(paste0(base_dir, '/data/ukb/enc/ukb29196.r'))
bd_29196 = bd


# merge data
ukb_df = merge(x=bd_26541, y=bd_27118, by='f.eid', all.x=T, all.y=T)
rm(bd_26541)
ukb_df = merge(x=ukb_df, y=bd_27119, by='f.eid', all.x=T, all.y=T)
ukb_df = merge(x=ukb_df, y=bd_27121, by='f.eid', all.x=T, all.y=T)
ukb_df = merge(x=ukb_df, y=bd_27123, by='f.eid', all.x=T, all.y=T)
ukb_df = merge(x=ukb_df, y=bd_27124, by='f.eid', all.x=T, all.y=T)
ukb_df = merge(x=ukb_df, y=bd_29196, by='f.eid', all.x=T, all.y=T)


# Step 2: Remove duplicate columns present in more than one bucket
# ----------------
duplicate_cols = colnames(ukb_df)[grep('.x', colnames(ukb_df))]
for (col in duplicate_cols){
    new_col = gsub('.x','',col)
    ukb_df[[new_col]] = ukb_df[[col]]
}
rm_cols = colnames(ukb_df)[grep('.x|.y', colnames(ukb_df))]
ukb_df = ukb_df[!colnames(ukb_df) %in% rm_cols]



# Step 3: Save
# -----------
# save multiple forms (csv/feather) of the data
write_feather(ukb_df, paste0(args[1], '/data/ukb/enc/ukb_pheno.feather'))
write_csv(ukb_df, paste0(base_dir, '/data/ukb/enc/ukb_pheno.csv'))


# subset to freesurfer completed subjects
fs_subs    = read_csv('/gpfs/milgram/project/holmes/kma52/mdd_sst/data/freesurfer/fs_list.txt', col_names=FALSE)
ukb_mri_df = ukb_df[ukb_df$f.eid %in% fs_subs$X1,]
write_feather(ukb_mri_df, paste0(base_dir, '/data/ukb/enc/ukb_pheno_mri.feather'))
write_csv(ukb_mri_df, paste0(base_dir, '/data/ukb/enc/ukb_pheno_mri.csv'))



