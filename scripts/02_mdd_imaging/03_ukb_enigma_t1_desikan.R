library(tidyverse)
library(Cairo)
library(extrafont)
loadfonts()


# 1. parcel-wise UKB regression using desikan atlas data
# 2. compare to ENIGMA cohen's d estimates from Schmaall meta analysis

cohens_d = function(t, df){
    (2*t)/sqrt(df)
}

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab_desikan = function(values, out_path){
    base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
    dscalar_template = '/gpfs/milgram/project/holmes/kma52/ukb_pymood/ref_files/desikan_atlas_32k.dscalar.nii'
    write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/tmp/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotDesikanVolOnSurface'
    cmd     = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/ukb_pymood/analyses/multimodal_mood/scripts/01_prep_data\');',
                    matfunc, '(\'', dscalar_template, '\',\'', write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}

mdd_regr_function = function(regr_dat, dv_vars, iv_var, descriptor){

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

    regr_df = NULL
    for (c in dv_vars){
        formula  = paste0('scale(', c,') ~ ',iv_var,' + ', covars)

        lm_model = lm(data=regr_dat, formula=formula)
        out      = summary(lm_model)
        out_row  = out$coefficients[rownames(out$coefficients) == iv_var,]

        out_row = as.data.frame(t(out_row))
        colnames(out_row) = c('est','se','tval','pval')
        out_row$roi = c
        out_row$df  = out$df[2]
        regr_df = rbind(regr_df, out_row)
    }
    head(regr_df,2)
    regr_df$roi = gsub(paste0('_', descriptor), '', regr_df$roi)

    # convert to Cohen's D
    regr_df$cohens_d = as.numeric(lapply(regr_df$tval, function(x) cohens_d(t=x, df=regr_df$df[1])))
    summary(regr_df$cohens_d)

    return(regr_df)
}



# Step 1: Read Data
# ------------
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read UKB thickness, summarised to desikan atlas
ukb_desikan_path = paste0(project_dir, '/data/ukb/dataframes/desikan_stats.csv')
ukb_desikan = read_csv(ukb_desikan_path)
thick_cols  = colnames(ukb_desikan)[grep('thickness', colnames(ukb_desikan))]
thick_cols  = thick_cols[!grepl('Mean', thick_cols)]

# merge with UKB phenotype
file_in = paste0(project_dir, '/data/ukb/ukb_severe_ctl.csv')
ukb_df  = read_csv(file_in)
ukb_desikan_df = merge(x=ukb_df, y=ukb_desikan, by.x='UKB_ID', by.y='UKB_ID')


# Step 2: Regression Desikan Thickness by MDD
# ----------------
thick_df = mdd_regr_function(regr_dat=ukb_desikan_df, dv_vars=thick_cols, iv_var='mdd_status.2.0_cc', descriptor='thickness')

# save regression dat
thick_write_file = paste0(project_dir, '/output/mdd_regr/ukbb_thick_desikan_regr_table.csv')
write_csv(x=thick_df, path=thick_write_file)

enigma_out = paste0(project_dir, '/supp_data/SuppData_6_ukbb_thick_desikan_regr_table.csv')
write_csv(thick_df, enigma_out)




# Step 3: Read ENIGMA data
# -----------
# read ENIGMA data and format desikan names for matching
enigma_in  = paste0(project_dir, '/reference_files/ENIGMA_MDD_thickness_desikan_multiple.csv')
enigma_mdd = read_csv(enigma_in)
enigma_mdd$roi = gsub('Left ','lh_',enigma_mdd$roi)
enigma_mdd$roi = gsub('Right ','rh_',enigma_mdd$roi)
enigma_mdd$roi = gsub(' ','',enigma_mdd$roi)
enigma_mdd$roi = gsub('cortex','',enigma_mdd$roi)
enigma_mdd$roi = gsub('gyrus','',enigma_mdd$roi)
enigma_mdd$roi = gsub('lobule','',enigma_mdd$roi)
enigma_mdd$roi = gsub('sulcus','',enigma_mdd$roi)
enigma_mdd$roi = gsub('superiortemporalsulcus','sts',enigma_mdd$roi)
enigma_mdd$roi = gsub('bankssuperiortemporal','bankssts',enigma_mdd$roi)
enigma_mdd     = enigma_mdd[!grepl('average', enigma_mdd$roi),]

enigma_out = paste0(project_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv')
write_csv(enigma_mdd, enigma_out)


# plot the desikan data for viewing in workbench
# -------------
ref_in      = paste0(project_dir, '/reference_files/desikan_colors.txt')
desikan_ref = read.table(ref_in, header=T)
desikan_ref$LabelName = gsub('-','_',gsub('ctx-','',desikan_ref$LabelName))

write_array = NULL
for (name in desikan_ref$LabelName){
    match_idx = which(enigma_mdd$roi == name)
    if (length(match_idx) == 1){
        write_val = as.numeric(enigma_mdd[match_idx,'cohens_d'])
    } else {
        write_val = 0
    }
    write_array = c(write_array, write_val)
}
out_file = paste0(project_dir, '/figures/surface_plots/enigma_mdd_desikan_cohens_d.dscalar.nii')
plot_matlab_desikan(values=write_array, out_path=out_file)
conv_cmd = paste('/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command -cifti-change-mapping ',out_file,'ROW',out_file,'-scalar')
system(conv_cmd)


# plot the UKBB data
# -------------
thick_df$roi = gsub('_thickness', '', thick_df$roi)
write_array = NULL
for (name in desikan_ref$LabelName){
    match_idx = which(thick_df$roi == name)
    if (length(match_idx) == 1){
        write_val = as.numeric(thick_df[match_idx,'cohens_d'])
    } else {
        write_val = 0
    }
    write_array = c(write_array, write_val)
}
out_file = paste0(project_dir, '/figures/surface_plots/ukb_mdd_desikan_cohens_d.dscalar.nii')
plot_matlab_desikan(values=write_array, out_path=out_file)
conv_cmd = paste('/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command -cifti-change-mapping ',out_file,'ROW',out_file,'-scalar')
system(conv_cmd)
# -------------



# Step 4: Compare UKB / ENIGMA effects
# -------------
# merge ukbb and enigma estimates
colnames(thick_df) = paste0('ukb_', colnames(thick_df))
ukb_enigma_df = merge(x=enigma_mdd, by.x='roi', y=thick_df, by.y='ukb_roi')
cor.test(ukb_enigma_df$ukb_cohens_d, ukb_enigma_df$cohens_d)
cor.test(ukb_enigma_df$ukb_cohens_d, ukb_enigma_df$cohens_d, method='spearman')



# plot
out_path = paste0(project_dir, '/figures/ukbb_enigma_thick_corrplot.pdf')
out_path
p = ggplot(ukb_enigma_df, aes(x=ukb_cohens_d, y=cohens_d)) +
      geom_point(size=2.5, fill='#D6D6D6', color='black', stroke=.25, shape=21) +
      geom_smooth(method=lm, fullrange=TRUE, linetype='dashed', color='black', se=F) +
          theme_classic() +
            theme(
            text=element_text(size=8,  family='Helvetica', colour='black')) +
        scale_y_continuous(limits = c(-.18, .12), breaks=seq(-.18,.12,.06), expand = c(0, 0)) +
        scale_x_continuous(limits = c(-.04, .04), breaks=seq(-.04,.04,.02), expand = c(0, 0)) +
        theme(text=element_text(colour="black"), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))

CairoPDF(out_path, width=2.5, height=2.5, family="Arial")
p
dev.off()










