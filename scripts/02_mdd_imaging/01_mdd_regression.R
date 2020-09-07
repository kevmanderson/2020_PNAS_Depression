library(tidyverse)
library(Cairo)
library(RColorBrewer)


# 1. QC / quantitative checks on UKB MRI data and covariates
# 2. parcel-wise regression of UKB data
# 3. save and plot maps


cohens_d = function(t, df){
    (2*t)/sqrt(df)
}

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
    dscalar_template = paste0(base_dir, '/reference_files/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/reference_files/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0(base_dir, '/reference_files/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotVolOnSurface'
    cmd     = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util\');',
                    matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                    write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}


# run parcel-wise regression in UKB
mdd_regr_function = function(regr_df, dv_vars, iv_var, descriptor){

    covariate = c('white_mixed_other',
                    'genetic_ethnicity',
                    'scanner_X_pos_2_0_scale',
                    'scanner_Y_pos_2_0_scale',
                    'scanner_Z_pos_2_0_scale',
                    'sex.0.0',
                    'age_at_scan_scale',
                    'age_at_scan_scale*sex.0.0',
                    'age_square_scale*sex.0.0',
                    'MRI_T1_invSNR_2_0_scale',
                    'lesion_vol_scale',
                    'brain_size_scale',
                    'MRI_REST_motion_2_0_scale',
                    'MRI_REST_invSNR_fullPreproc_2_0_scale',
                    'diastolic_2_0_scale',
                    'systolic_2_0_scale',
                    'UK_Biobank_assessment_centre.2.0')

    covars = paste(covariate, collapse=' + ')
    out_df = NULL
    for (c in dv_vars){

        formula  = paste0('scale(', c,') ~ ',iv_var,' + ', covars)
        lm_model = lm(data=regr_df, formula=formula)
        out      = summary(lm_model)
        out_row  = out$coefficients[rownames(out$coefficients) == iv_var,]

        out_row = as.data.frame(t(out_row))
        colnames(out_row) = c('est','se','tval','pval')
        out_row$roi = c
        out_row$df  = out$df[2]
        out_df = rbind(out_df, out_row)
    }
    head(out_df,2)

    out_df$roi = gsub(paste0('_', descriptor), '', out_df$roi)

    # convert to Cohen's D
    out_df$cohens_d = as.numeric(lapply(out_df$tval, function(x) cohens_d(t=x, df=out_df$df[1])))
    summary(out_df$cohens_d)
    head(out_df[order(out_df$cohens_d),])

    return(out_df)
}


# Step 1: Read data
# --------------
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
pheno_path  = paste0(project_dir, '/data/ukb/enc/ukb_pheno_processed_mri_FS_RSFA_GBC.csv')
image_df    = read_csv(pheno_path)
dim(image_df)

# read self-reported ethnicity
source('/gpfs/milgram/project/holmes/kma52/mdd_sst/data/ukb/enc/ethnicity.r')
ethnicity_df = bd
ethnicity_df$ethnicity = ethnicity_df$f.21000.0.0

# merge the two
image_df = merge(x=image_df, by.x='UKB_ID', y=ethnicity_df, by.y='f.eid')




# Step 2: Fix Schaefer naming
# --------------
# Fixed schaefer mapping. The Yeo Lab CBIG repository found an error in the naming conventions of their label file halfway through analyses).
# This doesn't affect any reported results, but we correct minor ROI naming issues to help with forward compatability
info_file   = paste0(project_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt')
new_mapping = read.csv(header=F, info_file)
new_map_df  = data.frame(new_net=as.character(new_mapping$V1[grep('17Networks', new_mapping$V1)]),
                        new_info=as.character(new_mapping$V1[!grepl('17Networks', new_mapping$V1)]), stringsAsFactors=F)


old_info_file = paste0(project_dir, '/data/Schaefer/Schaefer2018_200Parcels_17Networks_order_info.txt')
old_mapping  = read.csv(header=F, old_info_file)
old_map_df = data.frame(old_net=as.character(old_mapping$V1[grep('17Networks', old_mapping$V1)]),
                        old_info=as.character(old_mapping$V1[!grepl('17Networks', old_mapping$V1)]), stringsAsFactors=F)

both_map_df = cbind(new_map_df, old_map_df)
remap_df    = both_map_df[both_map_df$old_net != both_map_df$new_net,]

old_df_cols = colnames(image_df)
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
colnames(image_df) = new_df_cols



# Step 3: Create/format variables
# -------------
# make some phenotypes
image_df$brain_size = image_df$MRI_T1_volGM_WM.2.0 + image_df$MRI_T1_ventCSF.2.0
image_df$UK_Biobank_assessment_centre.2.0 = as.character(image_df$UK_Biobank_assessment_centre.2.0)
image_df$genetic_ethnicity = image_df$genetic_ethnicity.0.0
image_df$genetic_ethnicity[is.na(image_df$genetic_ethnicity)] = 'NotCaucasian'



# Step 4: Remove imaging outliers
# -------------
# imaging columns
thick_cols = colnames(image_df)[grep('_thickness', colnames(image_df))]
gbc_cols   = colnames(image_df)[grep('_gbc', colnames(image_df))]
rsfa_cols  = colnames(image_df)[grep('_rsfa', colnames(image_df))]

# get rid of subjects missing any imaging phenotype
dim(image_df)
image_df = image_df[!is.na(rowMeans(image_df[gbc_cols])),]
image_df = image_df[!is.na(rowMeans(image_df[rsfa_cols])),]
image_df = image_df[!is.na(rowMeans(image_df[thick_cols])),]
dim(image_df)

# global outliers greater than 3sd
dim(image_df)
image_df = image_df[which(abs(scale(rowMeans(image_df[gbc_cols]))) < 3),]
image_df = image_df[which(abs(scale(rowMeans(image_df[rsfa_cols]))) < 3),]
image_df = image_df[which(abs(scale(rowMeans(image_df[thick_cols]))) < 3),]
dim(image_df)

# ROI-based outliers
# thickness
thick_idxs = NULL
for (thick in thick_cols){
    roi_idxs   = which(abs(scale(image_df[[thick]])) > 3)
    thick_idxs = c(thick_idxs, roi_idxs)
}

# RSFA
rsfa_idxs = NULL
for (rsfa in rsfa_cols){
    roi_idxs   = which(abs(scale(image_df[[rsfa]])) > 3)
    rsfa_idxs = c(rsfa_idxs, roi_idxs)
}

# GBC
gbc_idxs = NULL
for (gbc in gbc_cols){
    roi_idxs   = which(abs(scale(image_df[[gbc]])) > 3)
    gbc_idxs = c(gbc_idxs, roi_idxs)
}

# remove subjects with at least 10 outlier ROIs
thick_rm_idxs = names(table(thick_idxs)[table(thick_idxs) > 10])
rsfa_rm_idxs  = names(table(rsfa_idxs)[table(rsfa_idxs) > 10])
gbc_rm_idxs   = names(table(gbc_idxs)[table(gbc_idxs) > 10])

rm_idxs  = as.numeric(unique(c(thick_rm_idxs, rsfa_rm_idxs, gbc_rm_idxs)))
dim(image_df)
image_df = image_df[!1:nrow(image_df) %in% rm_idxs,]
dim(image_df)

# rename imaging features so they are compatible with lm()
for (feat in c(thick_cols, gbc_cols, rsfa_cols)){
    new_feat = paste0(gsub('rh_|lh_','', gsub('17Networks', 'net17', feat)), '_use')
    image_df[[new_feat]] = image_df[[feat]]
}
thick_cols = colnames(image_df)[grep('_thickness_use', colnames(image_df))]
gbc_cols   = colnames(image_df)[grep('_gbc_use', colnames(image_df))]
rsfa_cols  = colnames(image_df)[grep('_rsfa_use', colnames(image_df))]


# filter on ethnicity
image_df = image_df[!is.na(image_df$ethnicity),]

# remove brain size/lesion outliers
dim(image_df)
image_df = image_df[which(abs(scale(image_df$brain_size)) < 3),]
image_df = image_df[which(abs(scale(image_df$lesion_vol)) < 3),]
dim(image_df)


# save
write_csv(image_df, paste0(project_dir, '/data/ukb/ukb_all_mdd_levels_ctl_FIX.csv'))



# binarize mdd -- only look at control vs severe/moderate
image_df$mdd_status.2.0_cc = image_df$mdd_status.2.0
severe_ctl  = image_df[image_df$mdd_status.2.0_cc %in% c(0,2,3),]
severe_ctl$mdd_status.2.0_cc[severe_ctl$mdd_status.2.0_cc > 0]=1
regr_data   = severe_ctl

# reformat self-reported ethnicity to white/mixed/other
regr_data$white_mixed_other = 'Other'
regr_data$white_mixed_other[grep('mixed|and|Mixed', regr_data$ethnicity, ignore.case=T)] = 'mixed'
regr_data$white_mixed_other[grep('British|White|Irish', regr_data$ethnicity, ignore.case=T)] = 'white'
regr_data$white_mixed_other[grep('and|or', regr_data$ethnicity, ignore.case=T)] = 'mixed'
table(regr_data$ethnicity,regr_data$white_mixed_other)
table(regr_data$white_mixed_other,regr_data$genetic_ethnicity)




# Step 4: Scale quantitative covariates
# -------------
quant_covars = c('scanner_X_pos.2.0',
                 'scanner_Y_pos.2.0',
                 'scanner_Z_pos.2.0',
                 'age_at_scan',
                 'MRI_T1_invSNR.2.0',
                 'lesion_vol',
                 'brain_size',
                 'MRI_REST_motion.2.0',
                 'MRI_REST_invSNR_fullPreproc.2.0',
                 'diastolic.2.0',
                 'systolic.2.0')
for (covar in quant_covars){
    print(covar)
    new_covar = paste0(gsub('[.]','_',covar), '_scale')
    regr_data[[new_covar]] = scale(regr_data[[covar]])
}
regr_data$age_square_scale = regr_data$age_at_scan^2


# save checkpoint
severe_ctl = regr_data
write_csv(regr_data, paste0(project_dir, '/data/ukb/ukb_severe_ctl_FIX.csv'))
severe_ctl = read_csv(paste0(project_dir, '/data/ukb/ukb_severe_ctl_FIX.csv'))



# Step 5: Regression, MDD predicting imaging
# -----------
# Thickness
# ----------
thick_df = mdd_regr_function(regr_df=severe_ctl, dv_vars=thick_cols, iv_var='mdd_status.2.0_cc', descriptor='thickness_use')
head(thick_df[order(thick_df$cohens_d),])
thick_df[p.adjust(thick_df$pval, method='BH') < .05,]

# save regression dat
thick_write_file = paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
write_csv(x=thick_df, path=thick_write_file)

thick_write_file = paste0(project_dir, '/supp_data/SuppData_3_ukbb_thick_regr_table_newmap.csv')
write_csv(x=thick_df, path=thick_write_file)


# RSFA
# ----------------
rsfa_scale = apply(regr_data[rsfa_cols], 1, scale)
rsfa_scale = t(rsfa_scale)
colnames(rsfa_scale) = paste0(rsfa_cols,'_scale')

# add back to data frame
tmp_regr_data = cbind(regr_data, rsfa_scale)

# run regression
rsfa_df = mdd_regr_function(regr_df=tmp_regr_data, dv_vars=paste0(rsfa_cols,'_scale'), iv_var='mdd_status.2.0_cc', descriptor='rsfa_use')
rsfa_df[p.adjust(rsfa_df$pval, method='BH') < .05,]

# save
rsfa_write_file = paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
write_csv(x=rsfa_df, path=rsfa_write_file)

thick_write_file = paste0(project_dir, '/supp_data/SuppData_4_ukbb_rsfa_regr_table_newmap.csv')
write_csv(x=thick_df, path=thick_write_file)



# GBC
# ----------------
gbc_df = mdd_regr_function(regr_df=regr_data, dv_vars=gbc_cols, iv_var='mdd_status.2.0_cc', descriptor='gbc_use')
gbc_df[p.adjust(gbc_df$pval, method='BH') < .05,]

# save regression dat
gbc_write_file = paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
write_csv(x=gbc_df, path=gbc_write_file)

gbc_write_file = paste0(project_dir, '/supp_data/SuppData_5_ukbb_gbc_regr_table_newmap.csv')
write_csv(x=gbc_df, path=gbc_write_file)



# Step 6: Merge Thickness/RSFA/GBC regressino data
# ----------------
rsfa_write_file  = paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
rsfa_df          = read_csv(rsfa_write_file)
thick_write_file = paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
thick_df         = read_csv(thick_write_file)
gbc_write_file   = paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
gbc_df           = read_csv(gbc_write_file)

thick_merge = thick_df
gbc_merge   = gbc_df
rsfa_merge  = rsfa_df
rsfa_merge$roi = gsub('_scale','', rsfa_merge$roi)

colnames(gbc_merge)   = paste0('gbc_', colnames(gbc_merge))
colnames(rsfa_merge)  = paste0('rsfa_', colnames(rsfa_merge))
rsfa_merge$rsfa_roi   = gsub('_scale','',rsfa_merge$rsfa_roi)
colnames(thick_merge) = paste0('thick_', colnames(thick_merge))

all_regr = merge(x=thick_merge, y=rsfa_merge, by.x='thick_roi', by.y='rsfa_roi')
all_regr = merge(x=all_regr, y=gbc_merge, by.x='thick_roi', by.y='gbc_roi')





# Step 7: Prep data for cortical surface plotting
# ----------
parcel_num  = '200'
net_num = '17'
#net_name_in = read.csv(header=F, paste0(project_dir, '/data/Schaefer/Schaefer2018_200Parcels_17Networks_order_info.txt'))

# Fixed schaefer mapping
info_file = paste0(project_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt')
net_name_in = read.csv(header=F, info_file)
net_names   = as.character(net_name_in$V1[grep('Net', net_name_in$V1)])
net17       = data.frame(net=gsub('17Networks','net17',net_names), idx=1:as.numeric(parcel_num))


# merge regression data with schaefer plot order
plot_df = merge(x=net17, y=all_regr, by.x='net', by.y='thick_roi')
plot_df = plot_df[order(plot_df$idx),]
all_regr$thick_roi[!all_regr$thick_roi %in% net17$net]



# Step 8: Plot regression Cohen's D values
# ----------
# create averaged plots
nparcel = '200'
net_num = '17'
head(plot_df)

# thickness
out_path = paste0(project_dir, '/figures/surface_plots/ukbb_thick200_net17_ModSevMDD_v_ctl_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(plot_df$thick_cohens_d), out_path=out_path , parcel_num=nparcel, net_num=net_num)

# rsfa
out_path = paste0(project_dir, '/figures/surface_plots/ukbb_rsfa200_net17_ModSevMDD_v_ctl_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(plot_df$rsfa_cohens_d),out_path=out_path , parcel_num=nparcel, net_num=net_num)

# gbc
out_path = paste0(project_dir, '/figures/surface_plots/ukbb_gbc200_net17_ModSevMDD_v_ctl_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(plot_df$gbc_cohens_d),out_path=out_path , parcel_num=nparcel, net_num=net_num)

# average
avg_cohens = rowMeans(cbind(plot_df$thick_cohens_d, plot_df$gbc_cohens_d, plot_df$rsfa_cohens_d*-1))
out_path = paste0(project_dir, '/figures/surface_plots/ukbb_avg200_net17_ModSevMDD_v_ctl_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(avg_cohens),out_path=out_path , parcel_num=nparcel, net_num=net_num)






# Control analysis
# ----------------
# See if medication status changes reported imaging effects

# UKBB medication code
rx_table        = read_delim(paste0(project_dir, '/reference_files/ukbb_rx_codes.tsv'), delim='\t')
wray_rx_table   = read.csv(paste0(project_dir, '/reference_files/ukbb_wray_med_class.csv'))
antidepressants = wray_rx_table[grep('N06A', wray_rx_table$Medication_ATC_code),]
rx_cols = colnames(severe_ctl)[grep('Rx_code.2', colnames(severe_ctl))]


# determine if subject taking antidepressant
rx_grep  = lapply(rx_cols, function(x) severe_ctl[[x]] %in% antidepressants$Coding)
rx_table = do.call(cbind, rx_grep)
antidepr_ct = rowSums(rx_table)

# binary - taking anti-depressant or not
severe_ctl$antidepr_ct = antidepr_ct
severe_ctl$antidepr_bin = severe_ctl$antidepr_ct
severe_ctl$antidepr_bin[severe_ctl$antidepr_bin > 0] = 1

# run parcel-wise regression in UKB
mdd_regr_function = function(regr_df, dv_vars, iv_var, descriptor){

    covariate = c('antidepr_bin',
                    'white_mixed_other',
                    'genetic_ethnicity',
                    'scanner_X_pos_2_0_scale',
                    'scanner_Y_pos_2_0_scale',
                    'scanner_Z_pos_2_0_scale',
                    'sex.0.0',
                    'age_at_scan_scale',
                    'age_at_scan_scale*sex.0.0',
                    'age_square_scale*sex.0.0',
                    'MRI_T1_invSNR_2_0_scale',
                    'lesion_vol_scale',
                    'brain_size_scale',
                    'MRI_REST_motion_2_0_scale',
                    'MRI_REST_invSNR_fullPreproc_2_0_scale',
                    'diastolic_2_0_scale',
                    'systolic_2_0_scale',
                    'UK_Biobank_assessment_centre.2.0')

    covars = paste(covariate, collapse=' + ')
    out_df = NULL
    for (c in dv_vars){

        formula  = paste0('scale(', c,') ~ ',iv_var,' + ', covars)
        lm_model = lm(data=regr_df, formula=formula)
        out      = summary(lm_model)
        out_row  = out$coefficients[rownames(out$coefficients) == iv_var,]

        out_row = as.data.frame(t(out_row))
        colnames(out_row) = c('est','se','tval','pval')
        out_row$roi = c
        out_row$df  = out$df[2]
        out_df = rbind(out_df, out_row)
    }
    head(out_df,2)

    out_df$roi = gsub(paste0('_', descriptor), '', out_df$roi)

    # convert to Cohen's D
    out_df$cohens_d = as.numeric(lapply(out_df$tval, function(x) cohens_d(t=x, df=out_df$df[1])))
    summary(out_df$cohens_d)
    head(out_df[order(out_df$cohens_d),])

    return(out_df)
}


# Thickness
rx_thick_df = mdd_regr_function(regr_df=severe_ctl, dv_vars=thick_cols, iv_var='mdd_status.2.0_cc', descriptor='thickness_use')
thick_df    = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv'))
thick_merge_df = merge(x=thick_df, y=rx_thick_df, by='roi')
cor.test(thick_merge_df$cohens_d.x, thick_merge_df$cohens_d.y, method='spearman')

# RSFA
rsfa_scale = apply(severe_ctl[rsfa_cols], 1, scale)
rsfa_scale = t(rsfa_scale)
colnames(rsfa_scale) = paste0(rsfa_cols,'_scale')

# add back to data frame
tmp_regr_data = cbind(severe_ctl, rsfa_scale)

rx_rsfa_df = rx_mdd_regr_function(regr_df=tmp_regr_data, dv_vars=paste0(rsfa_cols,'_scale'), iv_var='mdd_status.2.0_cc', descriptor='rsfa_use')
rsfa_df    = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv'))
rsfa_merge_df = merge(x=rsfa_df, y=rx_rsfa_df, by='roi')
cor.test(rsfa_merge_df$cohens_d.x, rsfa_merge_df$cohens_d.y, method='spearman')


# GBC
rx_gbc_df = rx_mdd_regr_function(regr_df=severe_ctl, dv_vars=gbc_cols, iv_var='mdd_status.2.0_cc', descriptor='gbc_use')
gbc_df    = read_csv(paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv'))
gbc_merge_df = merge(x=gbc_df, y=rx_gbc_df, by='roi')
cor.test(gbc_merge_df$cohens_d.x, gbc_merge_df$cohens_d.y, method='spearman')





