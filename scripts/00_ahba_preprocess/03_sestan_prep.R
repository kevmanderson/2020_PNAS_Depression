library(tidyverse)


# Read and prerocess PsychENCODE/BrainSpan developmental expression data
#   1. remove low expressed transcripts
#   2. normalized expression within each donor
#   3. identify two Schaefer parcels that align to Brainspan regions


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


# Step 1: Read PsychENCODE Brainspan data
# --------------
# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
pec_dir  = '/gpfs/milgram/project/holmes/kma52/devgen/data/PsychENCODE_dev'

# RPKM
rpkm         = read.table(header=T, paste0(pec_dir, '/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt'))
gene_symbols = unlist(lapply(as.character(rpkm$Geneid), function(x) strsplit(x,'[|]')[[1]][[2]]))
ensembl_arr  = unlist(lapply(as.character(rpkm$Geneid), function(x) strsplit(x,'[|]')[[1]][[1]]))
rpkm$symbol  = gene_symbols
rpkm$ensembl = ensembl_arr
subj_columns = colnames(rpkm)[grep('HSB', colnames(rpkm))]

# QC
mrna_qc    = read.csv(paste0(pec_dir, '/mRNA-seq_QC.csv'))
mrna_qc$id = paste0(mrna_qc$Braincode,'.', mrna_qc$Regioncode)
mrna_qc$long_id = paste0(mrna_qc$Braincode,'.', mrna_qc$Regioncode, '.', mrna_meta$Hemisphere)

# Sample metadata
mrna_meta  = read_csv(paste0(pec_dir, '/mRNA-seq_Sample_metadata.csv'))
mrna_meta  = mrna_meta[!is.na(mrna_meta$Braincode),]




# Step 2: Subset to adult samples/expressed transcripts
# -------------
adult_meta    = mrna_meta[mrna_meta$Age %in% c('18Y','19Y','21Y','23Y','30Y','36Y','37Y','40Y'),]
adult_mrna_qc = mrna_qc[mrna_qc$Braincode %in% adult_meta$Braincode,]
cort_regions  = c('V1C','M1C','S1C','DFC','VFC','A1C','IPC','MFC','OFC','ITC','STC')
adult_mrna_qc = adult_mrna_qc[adult_mrna_qc$Regioncode %in% cort_regions,]
adult_rpkm    = rpkm[,which(colnames(rpkm) %in% c('ensembl', 'symbol', adult_mrna_qc$id))]
subj_columns  = colnames(adult_rpkm)

mean((adult_meta$Days-270)/365)
sd((adult_meta$Days-270)/365)


# remove transcripts with median RPKM >= 1 across all samples
median_expr     = apply(adult_rpkm[adult_mrna_qc$id], 1, median)
brain_mrna_expr = adult_rpkm[which(median_expr >= 1),]


# Step 3: Scale within each subject
# -------------
rpkm_scale = NULL
sub_ids = unique(as.character(adult_mrna_qc$Braincode))
for (id in sub_ids){
    write(id,'')
    subj_expr = brain_mrna_expr[grep(id, colnames(brain_mrna_expr))]
    subj_expr_scale = t(apply(subj_expr, 1, scale))
    colnames(subj_expr_scale) = colnames(subj_expr)
    rpkm_scale = cbind(rpkm_scale, subj_expr_scale)
}
rpkm_scale_df = as.data.frame(rpkm_scale)
rpkm_scale_df$symbol  = brain_mrna_expr$symbol
rpkm_scale_df$ensembl = brain_mrna_expr$ensembl

write.csv(rpkm_scale_df, paste0(base_dir, '/data/brainspan/bspan_rpkm_scaled_ctx.csv'), row.names=F)
rpkm_scale_df = read.csv(paste0(base_dir,'/data/brainspan/bspan_rpkm_scaled_ctx.csv'))


# add hemisphere to sample names
new_col_names = colnames(rpkm_scale_df)
for (brain in adult_meta$Braincode){
    write(brain,'')
    hemi = adult_meta$Hemisphere[adult_meta$Braincode == brain]
    new_col_names = gsub(brain, paste0(brain, '.', hemi), new_col_names)
}
colnames(rpkm_scale_df) = new_col_names


# Step 4: average expression within each region/hemi
# -------------
col_names = NULL
avg_expr  = NULL
cort_regions = c('V1C','M1C','S1C','DFC','VFC','A1C','IPC','MFC','OFC','ITC','STC')
for (hemi in c('L','R')){
    for (cort in cort_regions){
        write(cort,'')
        col_names = c(col_names, paste0(cort, '_', hemi))
        reg_idxs  = grep(cort, colnames(rpkm_scale_df))
        hemi_idxs = grep(hemi, colnames(rpkm_scale_df))
        use_idxs  = intersect(reg_idxs, hemi_idxs)
        print(head(rpkm_scale_df[use_idxs]))
        reg_avg   = rowMeans(rpkm_scale_df[use_idxs])
        avg_expr  = cbind(avg_expr, reg_avg)
    }
}
avg_expr_df = as.data.frame(avg_expr)
colnames(avg_expr_df) = col_names
avg_expr_df$symbol  = rpkm_scale_df$symbol
avg_expr_df$ensembl = rpkm_scale_df$ensembl




# Step 5: Map Schaefer Parcels to Brainspan Regions
# -------------
# Read Schaeffer parc
schaef_dir    = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
#parc_label_in = read_csv(paste0(schaef_dir, '/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
parc_label_in = read_csv(paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
parc_labels   = parc_label_in$X1[grep('17Net', parc_label_in$X1)]
parc_labels   = gsub('17Networks','net17',parc_labels)

schaefer_df = data.frame(schaefer=parc_labels, parcel_num=1:length(parc_labels))

# identify approximate homologs to the Sestan regions
schaefer_df$sestan_roi = 'NA'

# OFC
ofc_schaef_names = c('net17_LH_LimbicB_OFC_1', 'net17_LH_LimbicB_OFC_2','net17_RH_LimbicB_OFC_1', 'net17_RH_LimbicB_OFC_2')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% ofc_schaef_names] = 'OFC'

# MFC
mfc_schaef_names = c('net17_LH_DefaultA_PFCm_1', 'net17_LH_DefaultA_PFCm_3','net17_RH_LimbicB_OFC_3', 'net17_RH_DefaultA_PFCm_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% mfc_schaef_names] = 'MFC'

# VFC
vfc_schaef_names = c('net17_LH_DefaultB_PFCv_4', 'net17_LH_SalVentAttnA_FrOper_2','net17_RH_DefaultB_PFCv_1', 'net17_RH_ContA_PFCl_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% vfc_schaef_names] = 'VFC'

# DFC
dfc_schaef_names = c('net17_LH_ContB_PFCl_1', 'net17_LH_SalVentAttnB_PFCl_1','net17_RH_ContB_PFCld_2', 'net17_RH_ContB_PFCld_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% dfc_schaef_names] = 'DFC'

# M1C
m1c_schaef_names = c('net17_LH_SomMotA_4', 'net17_LH_SomMotB_Cent_2','net17_RH_SomMotA_6', 'net17_RH_SomMotA_2')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% m1c_schaef_names] = 'M1C'

# S1C
s1c_schaef_names = c('net17_LH_SomMotA_3', 'net17_LH_SomMotA_2','net17_RH_SomMotA_4', 'net17_RH_SomMotA_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% s1c_schaef_names] = 'S1C'

# ITC
itc_schaef_names = c('net17_LH_LimbicA_TempPole_1', 'net17_LH_LimbicA_TempPole_2','net17_RH_LimbicA_TempPole_2', 'net17_RH_LimbicA_TempPole_3')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% itc_schaef_names] = 'ITC'

# A1C
a1c_schaef_names = c('net17_LH_SomMotB_Aud_3', 'net17_LH_SomMotB_Aud_2', 'net17_RH_SomMotB_S2_3', 'net17_RH_SomMotB_S2_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% a1c_schaef_names] = 'A1C'

# IPC
ipc_schaef_names = c('net17_LH_SalVentAttnB_IPL_1', 'net17_LH_SalVentAttnA_ParOper_1', 'net17_RH_SalVentAttnB_IPL_1', 'net17_RH_SalVentAttnA_ParOper_1')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% ipc_schaef_names] = 'IPC'

# STC
stc_schaef_names = c('net17_LH_TempPar_1', 'net17_LH_TempPar_2', 'net17_RH_TempPar_3', 'net17_RH_TempPar_2')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% stc_schaef_names] = 'STC'

# V1C
v1c_schaef_names = c('net17_LH_VisCent_ExStr_3', 'net17_LH_VisPeri_StriCal_1', 'net17_RH_VisCent_ExStr_3', 'net17_RH_VisPeri_ExStrSup_3')
schaefer_df$sestan_roi[schaefer_df$schaefer %in% v1c_schaef_names] = 'V1C'

write.csv(schaefer_df, paste0(base_dir, '/reference_files/bspan_schaefer_17net_200.csv'), row.names=F)
write.csv(avg_expr_df, paste0(base_dir, '/reference_files/bspan_avg_expr_df.csv'), row.names=F)


# colors for each Brainspan region
hex_df = data.frame(V1C=c('120', '18', '134'),
            OFC=c('122', '135', '51'),
            MFC=c('205','62','80'),
            ITC=c('220','248','165'),
            VFC=c('230','148','35'),
            DFC=c('135','50','76'),
            M1C=c('74','131','177'),
            S1C=c('9','41','250'),
            A1C=c('43','204','162'),
            IPC=c('0','118','16'),
            STC=c('255','255','2'))


# make a dlabel file to show which regions are being analyzed
labels = read.csv(header=F, paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), stringsAsFactors=F)
nets   = labels$V1[grep('Network', labels$V1)]

label_df = schaefer_df[schaefer_df$sestan_roi != 'NA',]
not_used = nets[!nets %in% gsub('net17','17Networks', label_df$schaefer)]

labels_new = labels
for (reg in unique(label_df$sestan_roi)){
    write(reg,'')
    reg_info     = label_df[label_df$sestan_roi == reg,]
    match_labels = gsub('net17', '17Networks', reg_info$schaefer)
    label_rows   = which(labels_new$V1 %in% match_labels)
    color_rows   = label_rows + 1
    for (row in color_rows){
        val_first = strsplit(labels_new$V1[row], ' ' )[[1]][[1]]
        new_row   = paste(val_first, paste0(as.character(hex_df[[reg]]), collapse=' '), '255')
        labels_new$V1[row] = new_row
    }
    col_dat = labels[color_rows,]
    unique(label_df$reg)
}
label_rows = which(labels_new$V1 %in% not_used)
color_rows = label_rows + 1
for (row in color_rows){
    val_first = strsplit(labels_new$V1[row], ' ' )[[1]][[1]]
    new_row = paste(val_first, '234 234 234 255')
    labels_new$V1[row] = new_row
}
write_delim(labels_new, paste0(base_dir, '/reference_files/SESTAN_MATCH_Schaefer2018_200Parcels_17Networks_order_info.txt'), delim='\t', col_names=F)


# run this command in terminal
# wb_command -cifti-label-import ./Schaefer2018_200Parcels_17Networks_order.dlabel.nii ./SESTAN_MATCH_Schaefer2018_200Parcels_17Networks_order_info.txt ./SESTAN_MATCH_Schaefer2018_200Parcels_17Networks_order.dlabel.nii


label_df$match_roi = paste0(label_df$sestan_roi, '_', ifelse(grepl('LH', label_df$schaefer), 'L', 'R'))
expr_array = rep(0, 200)

######
# SST
######
sst_expr    = t(avg_expr_df[which(avg_expr_df$ensembl == 'ENSG00000157005'),1:22])
sst_df      = data.frame(match_roi=rownames(sst_expr), expr=as.numeric(sst_expr))
sst_plot_df = merge(x=label_df, y=sst_df, by='match_roi')
expr_array[sst_plot_df$parcel_num] = sst_plot_df$expr

out_path = paste0(base_dir, '/figures/surface_plots/sst_sestan_expr.dscalar.nii')
plot_matlab(values=expr_array, out_path=out_path, parcel_num='200', net_num='17')


######
# CORT
######
expr_array = rep(0, 200)
cort_expr    = t(avg_expr_df[which(avg_expr_df$ensembl == 'ENSG00000241563'),1:22])
cort_df      = data.frame(match_roi=rownames(cort_expr), expr=as.numeric(cort_expr))
cort_plot_df = merge(x=label_df, y=cort_df, by='match_roi')
expr_array[cort_plot_df$parcel_num] = cort_plot_df$expr

out_path = paste0(base_dir, '/figures/surface_plots/cort_sestan_expr.dscalar.nii')
plot_matlab(values=expr_array, out_path=out_path, parcel_num='200', net_num='17')


######
# NPY
######
expr_array = rep(0, 200)
npy_expr    = t(avg_expr_df[which(avg_expr_df$ensembl == 'ENSG00000122585'),1:22])
npy_df      = data.frame(match_roi=rownames(npy_expr), expr=as.numeric(npy_expr))
npy_plot_df = merge(x=label_df, y=npy_df, by='match_roi')
expr_array[npy_plot_df$parcel_num] = npy_plot_df$expr

out_path = paste0(base_dir, '/figures/surface_plots/npy_sestan_expr.dscalar.nii')
plot_matlab(values=expr_array, out_path=out_path, parcel_num='200', net_num='17')







# Step 6: Map Desikan ROIs to Sestan Regions
# -----------------
ref_in      = paste0(base_dir, '/reference_files/desikan_colors.txt')
desikan_ref = read.table(ref_in, header=T)
desikan_ref$LabelName = gsub('-','_',gsub('ctx-','',desikan_ref$LabelName))

desikan_df = data.frame(desikan=desikan_ref$LabelName , parcel_num=1:length(desikan_ref$LabelName ))
desikan_df$desikan = as.character(desikan_df$desikan)
desikan_df$sestan_roi = 'NA'

# OFC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_lateralorbitofrontal', 'rh_lateralorbitofrontal')] = 'OFC'

# MFC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_rostralanteriorcingulate', 'rh_rostralanteriorcingulate')] = 'MFC'

# VFC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_parsopercularis', 'rh_parsopercularis')] = 'VFC'

# DFC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_rostralmiddlefrontal', 'rh_rostralmiddlefrontal')] = 'DFC'

# M1C
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_precentral', 'rh_precentral')] = 'M1C'

# S1C
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_postcentral', 'rh_postcentral')] = 'S1C'

# ITC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_fusiform', 'rh_fusiform')] = 'ITC'

# A1C
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_transversetemporal', 'rh_transversetemporal')] = 'A1C'

# IPC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_supramarginal', 'rh_supramarginal')] = 'IPC'

# STC
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_bankssts', 'rh_bankssts')] = 'STC'

# V1C
desikan_df$sestan_roi[desikan_df$desikan %in% c('lh_pericalcarine', 'rh_pericalcarine')] = 'V1C'


write.csv(desikan_df, paste0(base_dir, '/reference_files/bspan_desikan_200.csv'), row.names=F)




















