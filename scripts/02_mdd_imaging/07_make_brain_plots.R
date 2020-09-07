library(tidyverse)
library(gifti)
library(XML)


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util/brain_plot_functions.R')

# read sample information
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_zWithinSubject_zWithinSample.Rdata'), verbose=T)
colnames(ctx_data_scale) = paste0('wellid_', ctx_samp_all$well_id)

sub_list = unique(donorDat$samples$brain)
donor = sub_list[1]

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


# create individual expression plots for each gene of interest
for (gene in c('SST','CORT','NPY','CBLN2','CNR1')){
  for (hemi in c('lh','rh')){
    for ( donor in sub_list ){
        if (hemi=='lh'){
            name_string   = paste0(gene, '_LH')
            vert_ctx_hemi = sample_dat[intersect(which(sample_dat$brain == as.numeric(donor)), grep('left', sample_dat$structure_name)),]
            ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'
            projection_structure = 'CORTEX_LEFT'

        } else if (hemi=='rh') {
            name_string   =  paste0(gene, '_RH')
            vert_ctx_hemi = sample_dat[intersect(which(sample_dat$brain == as.numeric(donor)), grep('right', sample_dat$structure_name)),]
            ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'
            projection_structure = 'CORTEX_RIGHT'
        if (!donor %in% c('9861','10021')){
          next
        }
        }
        if (gene == 'SST'){
            mypal = colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D') )( 100 )
        } else if (gene == 'PVALB'){
            mypal = colorRampPalette( c("white", '#EBF0B5', '#46B6C2', '#283B8E')  )( 100 )
        } else if (gene == 'CORT'){
            mypal = colorRampPalette( c("white", '#BBDC9F', '#4DB56F', '#006C42') )( 100 )
        } else if (gene == 'NPY'){
            mypal = colorRampPalette( c("white", '#DFE2EE', '#9388BA', '#3C1F5D') )( 100 )
        }


        out_name = paste(base_dir, '/figures/surface_plots/PaperFoci_donor', donor, '_', name_string, '_', projection_structure, '.foci', sep='')
        print(out_name)
        # subset LH and donor sample information

        # vertices and MNI coords
        merge_df = data.frame(expr=reg_micro_scale[[gene]], well_id_match=rownames(reg_micro_scale))
        vert_ctx_hemi = merge(x=vert_ctx_hemi, y=merge_df, by='well_id_match')
        mni_coords = vert_ctx_hemi[c('mni_x','mni_y','mni_z','R','A','S')]
        vertices   = vert_ctx_hemi$vertex


        # convert expression values to rgb
        rgb_matrix = t(col2rgb(map2color(vert_ctx_hemi$expr, mypal, limits=c(-2.0, 2.0))))/256
        rgba_color = cbind(rgb_matrix, rep(1,dim(rgb_matrix)[1]))
        rgba_color = as.data.frame(rgba_color)
        colnames(rgba_color) = c('r', 'g', 'b', 'a')
        rgba_color = as.data.frame(rgba_color)

        # Do the plotting
        plot_foci(mni.coords=mni_coords,
                  vertices=vertices,
                  foci.names=vert_ctx_hemi$well_id,
                  ref.surface=ref_surface,
                  rgba.color=rgba_color,
                  wb.path=wb_path,
                  projection.structure=projection_structure,
                  out.name=out_name)
    }
  }
}


# color scale for figure 1
for (gene in c('SST','CORT','NPY')){
    if (gene == 'SST'){
        mypal = colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D') )( 100 )
    } else if (gene == 'CORT'){
        mypal = colorRampPalette( c("white", '#EBF0B5', '#46B6C2', '#283B8E')  )( 100 )
    } else if (gene == 'NPY'){
        mypal = colorRampPalette( c("white", '#BBDC9F', '#4DB56F', '#006C42') )( 100 )
    } else if (gene == 'CNR1'){
        mypal = colorRampPalette( c("white", '#DFE2EE', '#9388BA', '#3C1F5D') )( 100 )
    }
    out_name = paste0(base_dir, '/figures/surface_plots/PaperFoci_palette_', gene, '.pdf')
    print(out_name)
    pdf(out_name)
    z=matrix(1:100,nrow=1)
    x=1
    y=seq(-3,3,len=100)
    image(x,y,z,col=mypal,axes=FALSE,xlab="",ylab="")
    dev.off()
}


# Schaeffer parcellation
parc_label_in  = read_csv(paste0(base_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
parc_labels    = parc_label_in$X1[grep('17Net', parc_label_in$X1)]
parc_labels    = gsub('17Networks','net17',parc_labels)
parc_labels_df = data.frame(roi=parc_labels, idx=1:length(parc_labels))

# AHBA expression data
ahba_in       = paste0(base_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv')
ahba_expr_200 = read_csv(ahba_in)
parc_labels   = colnames(ahba_expr_200)[grepl('17Network', colnames(ahba_expr_200))]
#ahba_expr_200 = ahba_expr_200[,c(order(parc_labels),201)]

# fix up the roi names for alingment to regressino data
ahba_expr_200_t           = as.data.frame(t(ahba_expr_200[,1:200]))
colnames(ahba_expr_200_t) = ahba_expr_200$gene
ahba_expr_200_t$roi       = rownames(ahba_expr_200_t)
ahba_expr_200_t$roi       = gsub('17Networks', 'net17', ahba_expr_200_t$roi)


roi_intensity = colMeans(ahba_expr_200[,1:200])
out_path = paste0(base_dir, '/figures/surface_plots/ahba_roi_intensity_newmap.dscalar.nii')
plot_matlab(values=as.numeric(roi_intensity), out_path=out_path , parcel_num=nparcel, net_num=net_num)



gene='SST'
mypal       = colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D') )( 100 )
plot_dlabel = paste0(base_dir, '/figures/surface_plots/',gene,'_expr_parc200_net17_newmap.dlabel.nii')
plot_label(plot_vals=as.numeric(ahba_expr_200_t[[gene]]), limits=c(-2,2), mypal=mypal, parcel_num='200', net_num='17', out_path=plot_dlabel)
plot_dlabel

gene = 'CORT'
mypal = colorRampPalette( c("white", '#BBDC9F', '#4DB56F', '#006C42') )( 100 )
plot_dlabel = paste0(base_dir, '/figures/surface_plots/',gene,'_expr_parc200_net17.dlabel.nii')
plot_label(plot_vals=as.numeric(ahba_expr_200_t[[gene]]), limits=c(-2,2), mypal=mypal, parcel_num='200', net_num='17', out_path=plot_dlabel)
plot_dlabel

gene = 'NPY'
mypal = colorRampPalette( c("white", '#DFE2EE', '#9388BA', '#3C1F5D') )( 100 )
plot_dlabel = paste0(base_dir, '/figures/surface_plots/',gene,'_expr_parc200_net17.dlabel.nii')
plot_label(plot_vals=as.numeric(ahba_expr_200_t[[gene]]), limits=c(-2,2), mypal=mypal, parcel_num='200', net_num='17', out_path=plot_dlabel)
plot_dlabel



gene='DUSP6'
roi_intensity = colMeans(ahba_expr_200[,1:200])
out_path      = paste0(base_dir, '/figures/surface_plots/',gene,'_expr_parc200_net17_newmap.dscalar.nii')
gene_expr = as.numeric(ahba_expr_200_t[[gene]])
gene_expr[is.na(gene_expr)] = NaN
plot_matlab(values=gene_expr, out_path=out_path , parcel_num='200', net_num='17')



gene='EMX1'
roi_intensity = colMeans(ahba_expr_200[,1:200])
out_path      = paste0(base_dir, '/figures/surface_plots/',gene,'_expr_parc200_net17.dscalar.nii')
gene_expr = as.numeric(ahba_expr_200_t[[gene]])
gene_expr[is.na(gene_expr)] = NaN
plot_matlab(values=gene_expr, out_path=out_path , parcel_num='200', net_num='17')

