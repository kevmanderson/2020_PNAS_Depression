library(tidyverse)
library(gifti)
library(XML)


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util/brain_plot_functions.R')
'
'
# read sample information
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_zWithinSubject_zWithinSample.Rdata'), verbose=T)
colnames(ctx_data_scale) = paste0('wellid_', ctx_samp_all$well_id)

sub_list = unique(donorDat$samples$brain)
donor    = sub_list[1]


# read sample-to-vertex projection info
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))
sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]
sample_dat$well_id_match = paste0('wellid_', sample_dat$well_id)


# cortical sample by normalized gene expression data frame
reg_micro_scale = as.data.frame(t(ctx_data_scale))
reg_samples     = sample_dat


# stitch the left/right verticies together to match Schaeffer parcel cifti format
reg_samples$bihemi_vertex = reg_samples$vertex + 1 # cifti indices index at 0, R indexes at 1
right_ctx_idxs = intersect(grep('right', reg_samples$structure_name), which(reg_samples$top_level == 'CTX'))
reg_samples$bihemi_vertex[right_ctx_idxs] = reg_samples$bihemi_vertex[right_ctx_idxs] + 32492


# vis data
vis_path = paste0(base_dir, '/data/ahba_cibersortx/schaeffer_LAKE_VIS_200_17Net_expr_matvis_abs_orig_donorscaled_superordinate.csv')
vis_df   = read_csv(vis_path)


# dfc data
dfc_path = paste0(base_dir, '/data/ahba_cibersortx/schaeffer_LAKE_DFC_200_17Net_expr_matdfc_abs_orig_donorscaled_superordinate.csv')
dfc_df = read_csv(dfc_path)


# prepare data for plotting on cortical surface
parcel_num  = '200'
net_num     = '17'
net_name_in = read.csv(header=F, '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/Schaefer/Schaefer2018_200Parcels_17Networks_order_info.txt')
net_names   = as.character(net_name_in$V1[grep('Net', net_name_in$V1)])
net17       = data.frame(net=gsub('17Networks','net17',net_names), idx=1:as.numeric(parcel_num))


vis_df_t           = as.data.frame(t(vis_df[,1:200]))
colnames(vis_df_t) = vis_df$gene
vis_df_t$roi       = gsub('17Networks_','net17_',rownames(vis_df_t))


# merge regression data with schaefer plot order
plot_df = merge(x=net17, y=vis_df_t, by.x='net', by.y='roi')
plot_df = plot_df[order(plot_df$idx),]
vis_plot_df = plot_df

for (cell in colnames(plot_df)[grep('VIS', colnames(plot_df))]){
    write(cell,'')
    out_path = paste0(base_dir, '/figures/surface_plots/lake_vis_',cell,'_abs_orig_labels_donorscaled.dscalar.nii')
    plot_matlab(values=as.numeric(scale(plot_df[[cell]])), out_path=out_path , parcel_num=parcel_num, net_num=net_num)
}

dfc_df_t           = as.data.frame(t(dfc_df[,1:200]))
colnames(dfc_df_t) = dfc_df$gene
dfc_df_t$roi       = gsub('17Networks_','net17_',rownames(dfc_df_t))


# merge regression data with schaefer plot order
plot_df     = merge(x=net17, y=dfc_df_t, by.x='net', by.y='roi')
plot_df     = plot_df[order(plot_df$idx),]
dfc_plot_df = plot_df

for (cell in colnames(plot_df)[grep('DFC', colnames(plot_df))]){
    write(cell,'')
    out_path = paste0(base_dir, '/figures/surface_plots/lake_dfc_',cell,'_abs_orig_labels_donorscaled.dscalar.nii')
    plot_matlab(values=as.numeric(scale(plot_df[[cell]])), out_path=out_path , parcel_num=parcel_num, net_num=net_num)
}

# in1 averaged
plot_arr = rowMeans(cbind(dfc_plot_df$DFC_In1, vis_plot_df$VIS_In1))
out_path = paste0(base_dir, '/figures/surface_plots/lake_avg_In1_abs_orig_labels_donorscaled.dscalar.nii')
plot_matlab(values=plot_arr, out_path=out_path , parcel_num=parcel_num, net_num=net_num)

# sst averaged
plot_arr = rowMeans(cbind(dfc_plot_df$DFC_SST, vis_plot_df$VIS_SST))
out_path = paste0(base_dir, '/figures/surface_plots/lake_avg_SST_abs_orig_labels_donorscaled.dscalar.nii')
plot_matlab(values=plot_arr, out_path=out_path , parcel_num=parcel_num, net_num=net_num)

# Ast averaged
plot_arr = rowMeans(cbind(dfc_plot_df$DFC_Ast, vis_plot_df$VIS_Ast))
out_path = paste0(base_dir, '/figures/surface_plots/lake_avg_Ast_abs_orig_labels_donorscaled.dscalar.nii')
plot_matlab(values=plot_arr, out_path=out_path , parcel_num=parcel_num, net_num=net_num)





