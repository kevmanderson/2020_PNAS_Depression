library(tidyverse)
library(Cairo)


# 1. prepare and QC GSP data
# 2. parcel-wise regression with trait negative affect against RSFA and GBC
# 3. compare to UKB estimates


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

# gsp_rsfa_df  = mdd_regr_function(regr_dat=tmp_gsp_filt, dv_vars=paste0(rsfa_cols,'_scale'), iv_var='neg_affect_mean', descriptor='rsfa_')

mdd_regr_function = function(regr_dat, dv_vars, iv_var, descriptor){

    covariate = c('ht',
                    'wt',
                    'EstIQ_Ship',
                    'MF',
                    'Age_A',
                    'MF*Age_A',
                    'age_square',
                    'MF*age_square',
                    'Educ_A',
                    'Scanner',
                    'Console',
                    'ICV')
    covars = paste(covariate, collapse=' + ')

    regr_df = NULL
    for (roi in dv_vars){

        formula  = paste0('scale(', roi,') ~ scale(',iv_var,') + ', covars)
        lm_out   = summary(lm(data=regr_dat, formula=formula))
        out_row  = lm_out$coefficients[rownames(lm_out$coefficients) == paste0('scale(',iv_var,')'),]
        out_row  = as.data.frame(t(out_row))
        colnames(out_row) = c('est','se','tval','pval')

        colnames(lm_out$coefficients) = c('est','se','tval','pval')
        rownames(lm_out$coefficients) = gsub('[(]|[)]', '', rownames(lm_out$coefficients))

        # reformat coefficient matrix
        out_mat = NULL
        for (r in 1:nrow(lm_out$coefficients) ){
            row = lm_out$coefficients[r,]
            names(row) = paste(rownames(lm_out$coefficients)[r], names(row), sep='_')
            out_mat = c(out_mat, row)
        }
        out_row  = as.data.frame(t(out_mat))
        out_row$roi = roi
        out_row$df  = lm_out$df[[2]]
        regr_df = rbind(regr_df, out_row)
    }
    head(regr_df,2)
    regr_df$roi = gsub(descriptor, '', regr_df$roi)

    # convert to Cohen's D
    regr_df[[paste0('scale', iv_var,'_cohens_d')]] = as.numeric(lapply(regr_df[[paste0('scale', iv_var,'_tval')]], function(x) cohens_d(t=x, df=regr_df$df[1])))
    summary(regr_df$cohens_d)

    return(regr_df)
}


# Step 1: Read data
# --------------
base_dir     = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
gsp_dir      = paste0(base_dir, '/data/gsp')
gsp_cbig_dir = paste0(base_dir, '/data/gsp/cbig/')

# read processed RSFA data and GSP phenotype information
gsp_rsfa  = read_csv(paste0(gsp_dir, '/gsp_rsfa_net17_200.csv'))
gsp_pheno = read_csv(paste0(gsp_dir, '/GSP_MegaOverview_140124.csv'))

# merge
colnames(gsp_rsfa) = paste0('rsfa_', colnames(gsp_rsfa))
gsp_rsfa_df = merge(x=gsp_pheno, y=gsp_rsfa, by.x='Label', by.y='rsfa_ID')

# parcellation info
parc_200_17   = read_csv(paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=F)
parcel_labels = parc_200_17$X1[grep('Net', parc_200_17$X1)]

# make a list of subject IDs with GBC data
pconn_files = list.files(gsp_cbig_dir, pattern='*txt')
sub_arr     = unlist(lapply(as.list(pconn_files), function(x) paste0(strsplit(x[1], '_')[[1]][1:2], collapse='_')))
unique_subs = unique(sub_arr)


# Step 2: Calculate average GBC for each subject
# --------------
gsp_gbc_df = NULL
for (sub in unique_subs){
    print(sub)

    # connectivity data for this subject
    sub_pconn_files = list.files(gsp_cbig_dir, pattern=paste0(sub))
    sub_pconn_files = sort(sub_pconn_files)

    # select one file
    use_pconn       = sub_pconn_files[1]
    use_pconn_file  = paste0(gsp_cbig_dir, '/', use_pconn)

    # connectivity matrix
    cur_conn = read_delim(use_pconn_file, delim='\t', col_names=F)

    # parcelwise GBC
    sub_gbc = NULL
    for (parc in parcel_labels){
        roi_conn = cur_conn[which(parcel_labels != parc), which(parcel_labels == parc)]
        roi_mean = mean(roi_conn[[1]])
        sub_gbc  = c(sub_gbc, roi_mean)
    }
    sub_gbc_row = data.frame(t(c(sub_gbc, sub)))
    colnames(sub_gbc_row) = c(parcel_labels, 'UKB_ID')
    gsp_gbc_df = rbind(gsp_gbc_df, sub_gbc_row)
}
write_csv(x=gsp_gbc_df, paste0(gsp_dir, '/gbc_Schaefer2018_200parc_17net_df.csv'))
gsp_gbc_df = read_csv(paste0(gsp_dir, '/gbc_Schaefer2018_200parc_17net_df.csv'))


# merge
gbc_names = paste0('gbc_', colnames(gsp_gbc_df))
colnames(gsp_gbc_df) = gbc_names
gbc_names = gbc_names[1:200]
gsp_df    = merge(x=gsp_rsfa_df, y=gsp_gbc_df, by.x='Label', by.y='gbc_UKB_ID')


old_df_cols = colnames(gsp_df)
new_df_cols = old_df_cols



# Step 3: Fix schaefer mapping
# --------------
info_file   = paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt')
new_mapping = read.csv(header=F, info_file)
new_map_df  = data.frame(new_net=as.character(new_mapping$V1[grep('17Networks', new_mapping$V1)]),
                        new_info=as.character(new_mapping$V1[!grepl('17Networks', new_mapping$V1)]), stringsAsFactors=F)


old_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_200Parcels_17Networks_order_info.txt')
old_mapping  = read.csv(header=F, old_info_file)
old_map_df = data.frame(old_net=as.character(old_mapping$V1[grep('17Networks', old_mapping$V1)]),
                        old_info=as.character(old_mapping$V1[!grepl('17Networks', old_mapping$V1)]), stringsAsFactors=F)

both_map_df = cbind(new_map_df, old_map_df)
remap_df = both_map_df[both_map_df$old_net != both_map_df$new_net,]

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
colnames(gsp_df) = new_df_cols
colnames(gsp_df) = gsub('17Networks', 'net17', colnames(gsp_df))


# schaeffer column names
roi_labels = colnames(gsp_df)[grep('net17', colnames(gsp_df))]

# convert GBC to numeric
for (label in roi_labels){
    gsp_df[[label]] = as.numeric(as.character(gsp_df[[label]]))
}

# age-squared
gsp_df$age_square = gsp_df$Age_A^2


# Step 4: Calculate trait negative affect
# --------------
neg_aff = cbind(scale(gsp_df$POMS_TotMdDisturb), scale(gsp_df$TCI_HarmAvoidance), scale(gsp_df$BISBAS_BIS), scale(gsp_df$NEO_N), scale(gsp_df$STAI_tAnxiety))



# Step 5: Subset data/QC/remove imaging outliers
# --------------
# identify RSFA outliers
dat_only  = gsp_df[roi_labels]
dat_names = roi_labels

# get rid of subjects with zero in any parcel
dat_only = dat_only[(rowSums(dat_only[dat_names]) == 0) == 0,]

# GBC ROI-based outliers
gbc_rois = roi_labels[grep('gbc',roi_labels)]
all_outliers = NULL
for (gbc in gbc_rois){
  scale_gbc    = scale(dat_only[[gbc]])
  outlier_idxs = which(abs(scale_gbc) > 3)
  all_outliers = c(all_outliers, outlier_idxs)
}

rsfa_rois = roi_labels[grep('rsfa',roi_labels)]
for (rsfa in rsfa_rois){
  scale_rsfa    = scale(dat_only[[rsfa]])
  outlier_idxs = which(abs(scale_rsfa) > 3)
  all_outliers = c(all_outliers, outlier_idxs)
}

# add negative affect as columns
gsp_df$neg_affect_mean = rowMeans(neg_aff)
remove_idxs            = as.numeric(names(table(all_outliers))[table(all_outliers) > 10])
gsp_filt               = gsp_df[!1:nrow(gsp_df) %in% remove_idxs,]

# subset to Tim12 scanner
gsp_filt = gsp_filt[gsp_filt$Coil == 'Tim_12',]

# subjects with behavioral data passing QC
gsp_filt  = gsp_filt[which(gsp_filt$SM_Comp >= -8),]
gsp_filt  = gsp_filt[which(gsp_filt$iCog_Comp >= -8),]

# global RSFA outliers
rsfa_cols = colnames(gsp_filt)[grep('rsfa_net', colnames(gsp_filt))]
gsp_filt  = gsp_filt[abs(scale(rowMeans(gsp_filt[rsfa_cols]))) < 4,]

# global GBC outliers
gbc_cols = colnames(gsp_filt)[grep('gbc_net', colnames(gsp_filt))]
gsp_filt = gsp_filt[abs(scale(rowMeans(gsp_filt[gbc_cols]))) < 4,]

# plot average values
rsfa_mean = colMeans(gsp_filt[rsfa_cols])
out_path  = paste0(base_dir, '/figures/surface_plots/gsp_rsfa_mean_newmap.dscalar.nii')
plot_matlab(values=as.numeric(rsfa_mean),out_path=out_path , parcel_num='200', net_num='17')
out_path



# Step 6: RSFA regression by negative affect
# ----------------
rsfa_cols  = colnames(gsp_filt)[grep('rsfa_net', colnames(gsp_filt))]
rsfa_scale = apply(gsp_filt[rsfa_cols], 1, scale)
rsfa_scale = t(rsfa_scale)
colnames(rsfa_scale) = paste0(rsfa_cols,'_scale')

# run ROI regression
tmp_gsp_filt = cbind(gsp_filt, rsfa_scale)
gsp_rsfa_df  = mdd_regr_function(regr_dat=tmp_gsp_filt, dv_vars=paste0(rsfa_cols,'_scale'), iv_var='neg_affect_mean', descriptor='rsfa_')
gsp_rsfa_df$scaleneg_affect_mean_qval = p.adjust(gsp_rsfa_df$scaleneg_affect_mean_pval, method='BH')
head(gsp_rsfa_df[rev(order(gsp_rsfa_df$scaleneg_affect_mean_cohens_d)),])

#gsp_rsfa_df = mdd_regr_function(regr_dat=gsp_filt, dv_vars=rsfa_cols, iv_var='neg_affect_mean', descriptor='rsfa_')

# save regression datgsp_filt = cbind(gsp_filt, x)
rsfa_write_file = paste0(base_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv')
write_csv(x=gsp_rsfa_df, path=rsfa_write_file)

rsfa_write_file = paste0(base_dir, '/supp_data/SuppData_7_gsp_rsfa_regr_table_newmap.csv')
write_csv(x=gsp_rsfa_df, path=rsfa_write_file)
# ----------------

# create averaged plots
nparcel = '200'
net_num = '17'

# rsfa
out_path = paste0(base_dir, '/figures/surface_plots/gsp_negaff_scale_rsfa200_net17_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(gsp_rsfa_df$scaleneg_affect_mean_cohens_d),out_path=out_path , parcel_num=nparcel, net_num=net_num)
out_path



# Step 7: GBC regression by negative affect
# ----------------
gbc_cols = colnames(gsp_filt)[grep('gbc_net', colnames(gsp_filt))]
gsp_gbc_df = mdd_regr_function(regr_dat=gsp_filt, dv_vars=gbc_cols, iv_var='neg_affect_mean', descriptor='gbc_')

gsp_gbc_df$scaleneg_affect_mean_qval = p.adjust(gsp_gbc_df$scaleneg_affect_mean_pval, method='BH')
head(gsp_gbc_df[order(gsp_gbc_df$scaleneg_affect_mean_cohens_d),])

# save regression dat
gbc_write_file = paste0(base_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv')
write_csv(x=gsp_gbc_df, path=gbc_write_file)

gbc_write_file = paste0(base_dir, '/supp_data/SuppData_8_gsp_gbc_regr_table_newmap.csv')
write_csv(x=gsp_gbc_df, path=gbc_write_file)


# gbc
out_path = paste0(base_dir, '/figures/surface_plots/gsp_negaff_gbc200_net17_COHENSD_newmap.dscalar.nii')
plot_matlab(values=as.numeric(gsp_gbc_df$scaleneg_affect_mean_cohens_d),out_path=out_path , parcel_num=nparcel, net_num=net_num)
out_path




# Step 8: GSP to UKB Global Brain Connectivity consistency
# ----------------
ukb_gbc = read_csv(paste0(base_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv'))
colnames(ukb_gbc) = paste0('ukb_', colnames(ukb_gbc))

# merge and corr
ukb_gsp_gbc = merge(x=ukb_gbc, y=gsp_gbc_df, by.x='ukb_roi', by.y='roi')
cor(ukb_gsp_gbc$scaleneg_affect_mean_cohens_d, ukb_gsp_gbc$ukb_cohens_d)
cor(ukb_gsp_gbc$scaleneg_affect_mean_cohens_d, ukb_gsp_gbc$ukb_cohens_d, method='spearman')


# plot
out_path = paste0(base_dir, '/figures/ukbb_gbc_to_gsp_gbc_corrplot_newmap.pdf')
out_path
p = ggplot(ukb_gsp_gbc, aes(x=ukb_cohens_d, y=scaleneg_affect_mean_cohens_d)) +
      geom_point(size=2.5, fill='#D6D6D6', color='black', stroke=.25, shape=21) +
      geom_smooth(method=lm, fullrange=TRUE, linetype='dashed', color='black', se=F) +
          theme_classic() +
            theme(
            text=element_text(size=8,  family='Helvetica', colour='black')) +
        scale_y_continuous(limits = c(-.20, .20), labels=c('-.20','-.10','0','.10','.20'), breaks=c(-.20,-.10,0,.10,.20), expand = c(0, 0)) +
        scale_x_continuous(limits = c(-.08, .08), breaks=seq(-.08,.08,.04), expand = c(0, 0)) +
        theme(text=element_text(colour="black"), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))

CairoPDF(out_path, width=2.5, height=2.5)
p
dev.off()



# Step 9: GSP to UKB RSFA consistency
# ----------------
ukb_rsfa           = read_csv(paste0(base_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv'))
colnames(ukb_rsfa) = paste0('ukb_', colnames(ukb_rsfa))
ukb_rsfa$ukb_roi   = gsub('_scale','',ukb_rsfa$ukb_roi)
gsp_rsfa_df$roi    = gsub('_scale','',gsp_rsfa_df$roi)

# combined
ukb_gsp_rsfa = merge(x=ukb_rsfa, y=gsp_rsfa_df, by.x='ukb_roi', by.y='roi')

cor(ukb_gsp_rsfa$scaleneg_affect_mean_cohens_d, ukb_gsp_rsfa$ukb_cohens_d)
cor(ukb_gsp_rsfa$scaleneg_affect_mean_cohens_d, ukb_gsp_rsfa$ukb_cohens_d, method='spearman')

# plot
out_path = paste0(base_dir, '/figures/ukbb_gbc_to_gsp_rsfa_corrplot_newmap.pdf')
out_path

p = ggplot(ukb_gsp_rsfa, aes(x=ukb_cohens_d, y=scaleneg_affect_mean_cohens_d)) +
      geom_point(size=2.5, fill='#D6D6D6', color='black', stroke=.25, shape=21) +
      geom_smooth(method=lm, fullrange=TRUE, linetype='dashed', color='black', se=F) +
          theme_classic() +
            theme(
            text=element_text(size=8,  family='Helvetica', colour='black')) +
        scale_y_continuous(limits = c(-.18, .18), breaks=seq(-.18,.18,.09), expand = c(0, 0)) +
        scale_x_continuous(limits = c(-.08, .08), breaks=seq(-.08,.08,.04), expand = c(0, 0)) +
        theme(text=element_text(colour="black"), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))

CairoPDF(out_path, width=2.5, height=2.5)
p
dev.off()




