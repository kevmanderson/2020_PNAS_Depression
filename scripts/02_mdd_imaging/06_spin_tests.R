library(tidyverse)
library(cifti)
library(gifti)
library(matrixStats)

source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/external/rotate_parcellation/R/rotate.parcellation.R')

# 1. run spatial permutation spin-based tests

perm.sphere.p = function(x,y,perm.id,corr.type='spearman') {

    nroi  = dim(perm.id)[1]  # number of regions
    nperm = dim(perm.id)[2] # number of permutations

    rho.emp = cor(x,y,method=corr.type, use='pairwise.complete')  # empirical correlation

    # permutation of measures
    x.perm = y.perm = array(NA,dim=c(nroi,nperm))
    for (r in 1:nperm) {
        for (i in 1:nroi) {
            x.perm[i,r] = x[perm.id[i,r]]
            y.perm[i,r] = y[perm.id[i,r]]
        }
    }

    # correlation to unpermuted measures
    rho.null.xy = rho.null.yx = vector(length=nperm)
    for (r in 1:nperm) {
        rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type,use='pairwise.complete')
        rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type, use='pairwise.complete')
    }

    # p-value definition depends on the sign of the empirical correlation
    if (rho.emp>0) {
        p.perm.xy = sum(rho.null.xy>rho.emp)/nperm
        p.perm.yx = sum(rho.null.yx>rho.emp)/nperm
    } else {
        p.perm.xy = sum(rho.null.xy<rho.emp)/nperm
        p.perm.yx = sum(rho.null.yx<rho.emp)/nperm
    }

    # return average p-value
    return(list((p.perm.xy+p.perm.yx)/2, rho.null.xy, rho.null.yx))
}


# Step 1: Read data/surface giftis
# ------------
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read surface information
sphere_dir = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/Pipelines/global/templates/standard_mesh_atlases'
r_sphere = paste0(sphere_dir, '/R.sphere.32k_fs_LR.surf.gii')
l_sphere = paste0(sphere_dir, '/L.sphere.32k_fs_LR.surf.gii')

l_sphere_dat = readgii(l_sphere)
r_sphere_dat = readgii(r_sphere)


# schaefer parcellation
parcellation  = read_cifti(paste0(project_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order.dlabel.nii'))
parcel_labels = parcellation$data
n_verts       = length(parcellation$data)/2

# split labels by hemisphere
lh_parcel_labels = parcel_labels[1:n_verts]
rh_parcel_labels = parcel_labels[(n_verts+1):length(parcel_labels)]


# Step 2: Identify average vertex location for each parcel
# ------------
centroid_l = NULL
for ( parcel in 1:100 ){
    which(parcel_labels == parcel)
    avg_coord  = colMeans(l_sphere_dat$data$pointset[which(lh_parcel_labels == parcel),])
    centroid_l = rbind(centroid_l, avg_coord)
}

centroid_r = NULL
for ( parcel in 101:200 ){
    which(parcel_labels == parcel)
    avg_coord  = colMeans(r_sphere_dat$data$pointset[which(rh_parcel_labels == parcel),])
    centroid_r = rbind(centroid_r, avg_coord)
}


# Step 3: Do the rotation
perm.id = rotate.parcellation(centroid_l, centroid_r, 10000)
save(perm.id, file='/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ukb/perm_id.Rdata')
load('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ukb/perm_id.Rdata')


# Step 4: Permutation Spin Tests
# ------------
# RSFA
rsfa_gsp = read.csv(paste0(project_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv'))
rsfa_ukb = read.csv(paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv'))
cor.test(rsfa_gsp$scaleneg_affect_mean_cohens_d, rsfa_ukb$cohens_d, method='spearman')
perm_rsfa = perm.sphere.p(x=rsfa_ukb$cohens_d, y=rsfa_gsp$scaleneg_affect_mean_cohens_d, perm.id=perm.id,corr.type='spearman')

# GBC
# ------------
gbc_gsp = read.csv(paste0(project_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv'))
gbc_ukb = read.csv(paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv'))
cor.test(gbc_gsp$scaleneg_affect_mean_cohens_d, gbc_ukb$cohens_d, method='spearman')
perm_gbc = perm.sphere.p(x=gbc_ukb$cohens_d, y=gbc_gsp$scaleneg_affect_mean_cohens_d, perm.id=perm.id,corr.type='spearman')



# Step 5: Desikan Spin Tests
# ------------
sphere_dir = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/Pipelines/global/templates/standard_mesh_atlases'
r_sphere = paste0(sphere_dir, '/R.sphere.32k_fs_LR.surf.gii')
l_sphere = paste0(sphere_dir, '/L.sphere.32k_fs_LR.surf.gii')

l_sphere_dat = readgii(l_sphere)
r_sphere_dat = readgii(r_sphere)

parcel_names_in = read_csv(paste0(project_dir, '/reference_files/desikan_atlas_32k.txt'), col_names=F)
parcel_names = parcel_names_in$X1[!grepl(' ', parcel_names_in$X1)]
parcel_idxs = which(parcel_names != 'corpuscallosum')

parcellation  = read_cifti(paste0(project_dir, '/reference_files/desikan_atlas_32k.dlabel.nii'))
parcel_labels = parcellation$data
n_verts = length(parcellation$data)/2

lh_parcel_labels = parcel_labels[1:n_verts]
rh_parcel_labels = parcel_labels[(n_verts+1):length(parcel_labels)]

centroid_l = NULL
for ( parcel in parcel_idxs ){
    which(parcel_labels == parcel)
    avg_coord  = colMeans(l_sphere_dat$data$pointset[which(lh_parcel_labels == parcel),])
    centroid_l = rbind(centroid_l, avg_coord)
}

centroid_r = NULL
for ( parcel in parcel_idxs ){
    which(parcel_labels == parcel)
    avg_coord  = colMeans(r_sphere_dat$data$pointset[which(rh_parcel_labels == parcel),])
    centroid_r = rbind(centroid_r, avg_coord)
}

desikan_perm_id = rotate.parcellation(centroid_l, centroid_r, 10000)
save(desikan_perm_id, file='/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ukb/desikan_perm_id.Rdata')
load('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ukb/desikan_perm_id.Rdata', verbose=T)


# T1 SPIN
# read ENIGMA data and format desikan names for matching
enigma_out = paste0(project_dir, '/reference_files/formated_ENIGMA_MDD_thickness_desikan_multiple.csv')
enigma_mdd = read_csv(enigma_out)

t1_desikan_ukb  = read.csv(paste0(project_dir, '/output/mdd_regr/ukbb_thick_desikan_regr_table.csv'))
colnames(t1_desikan_ukb) = paste0('ukb_', colnames(t1_desikan_ukb))
t1_desikan_comb = merge(x=enigma_mdd, y=t1_desikan_ukb, by.x='roi', by.y='ukb_roi')
cor.test(t1_desikan_comb$cohens_d, t1_desikan_comb$ukb_cohens_d, method='spearman')

perm_gbc = perm.sphere.p(x=t1_desikan_comb$cohens_d, y=t1_desikan_comb$ukb_cohens_d, perm.id=desikan_perm_id, corr.type='spearman')







# Schaeffer parcellation
parc_label_in = read_csv(paste0(project_dir, '/reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), col_names=FALSE)
parc_labels   = parc_label_in$X1[grep('17Net', parc_label_in$X1)]
parc_labels   = gsub('17Networks','net17',parc_labels)

# AHBA expression data
ahba_in       = paste0(project_dir, '/data/ahba_parcel/schaeffer_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat_NEWMAP.csv')
ahba_expr_200 = read_csv(ahba_in)
parc_labels   = colnames(ahba_expr_200)[grepl('17Network', colnames(ahba_expr_200))]
ahba_expr_200 = ahba_expr_200[,c(order(parc_labels),201)]

# fix up the roi names for alingment to regressino data
ahba_expr_200_t           = as.data.frame(t(ahba_expr_200[,1:200]))
colnames(ahba_expr_200_t) = ahba_expr_200$gene
ahba_expr_200_t$roi       = rownames(ahba_expr_200_t)
ahba_expr_200_t$roi       = gsub('17Networks', 'net17', ahba_expr_200_t$roi)




# UKBB Thick
thick_file = paste0(project_dir, '/output/mdd_regr/ukbb_thick_regr_table_newmap.csv')
ukbb_thick = read_csv(thick_file)

# UKBB RSFA
rsfa_file = paste0(project_dir, '/output/mdd_regr/ukbb_rsfa_regr_table_newmap.csv')
ukbb_rsfa = read_csv(rsfa_file)
ukbb_rsfa$roi = gsub('_scale','',ukbb_rsfa$roi)

# UKBB GBC
gbc_file = paste0(project_dir, '/output/mdd_regr/ukbb_gbc_regr_table_newmap.csv')
ukbb_gbc = read_csv(gbc_file)

# Harvard RSFA
gsp_rsfa_file = paste0(project_dir, '/output/mdd_regr/gsp_rsfa_regr_table_newmap.csv')
gsp_rsfa = read_csv(gsp_rsfa_file)
gsp_rsfa$roi = gsub('_scale','',gsp_rsfa$roi)

# Harvard GBC
gsp_gbc_file = paste0(project_dir, '/output/mdd_regr/gsp_gbc_regr_table_newmap.csv')
gsp_gbc = read_csv(gsp_gbc_file)


sst_genes     = c('SST','CORT','NPY')
thick_expr_df = merge(x=ukbb_thick, by.x='roi', y=ahba_expr_200_t, by.y='roi')
for (gene in sst_genes){
    write(gene,'')
    ukbb_thick_gene = perm.sphere.p(x=thick_expr_df$cohens_d, y=thick_expr_df[[gene]], perm.id=perm.id, corr.type='spearman')
    corr = cor.test(x=thick_expr_df$cohens_d, y=thick_expr_df[[gene]], method='spearman')
    print(corr)
    print(ukbb_thick_gene[[1]])
}


rsfa_expr_df  = merge(x=ukbb_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
for (gene in sst_genes){
    write(gene,'')
    ukbb_rsfa_gene = perm.sphere.p(x=rsfa_expr_df$cohens_d, y=rsfa_expr_df[[gene]], perm.id=perm.id, corr.type='spearman')
    corr = cor.test(x=rsfa_expr_df$cohens_d, y=rsfa_expr_df[[gene]], method='spearman')
    print(corr)
    print(ukbb_rsfa_gene[[1]])
}


gbc_expr_df   = merge(x=ukbb_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')
for (gene in sst_genes){
    write(gene,'')
    ukbb_gbc_gene = perm.sphere.p(x=gbc_expr_df$cohens_d, y=gbc_expr_df[[gene]], perm.id=perm.id, corr.type='spearman')
    corr = cor.test(x=gbc_expr_df$cohens_d, y=gbc_expr_df[[gene]], method='spearman')
    print(corr)
    print(ukbb_gbc_gene[[1]])
}


#
gsp_rsfa_expr_df = merge(x=gsp_rsfa, by.x='roi', y=ahba_expr_200_t, by.y='roi')
for (gene in sst_genes){
    write(gene,'')
    gsp_rsfa_gene = perm.sphere.p(x=gsp_rsfa_expr_df$scaleneg_affect_mean_cohens_d, y=gsp_rsfa_expr_df[[gene]], perm.id=perm.id, corr.type='spearman')
    corr = cor.test(x=gsp_rsfa_expr_df$scaleneg_affect_mean_cohens_d, y=gsp_rsfa_expr_df[[gene]], method='spearman')
    print(corr)
    print(gsp_rsfa_gene[[1]])
}


gsp_gbc_expr_df  = merge(x=gsp_gbc, by.x='roi', y=ahba_expr_200_t, by.y='roi')
for (gene in sst_genes){
    write(gene,'')
    gsp_gbc_gene = perm.sphere.p(x=gsp_gbc_expr_df$scaleneg_affect_mean_cohens_d, y=gsp_gbc_expr_df[[gene]], perm.id=perm.id, corr.type='spearman')
    corr = cor.test(x=gsp_gbc_expr_df$scaleneg_affect_mean_cohens_d, y=gsp_gbc_expr_df[[gene]], method='spearman')
    print(corr)
    print(gsp_gbc_gene[[1]])
}





# gene by desikan ROI matrix
desikan_mat = read_csv(paste0(project_dir, '/data/ahba/desikan_ahba_ctx_zWithinSubject_zWithinSample_200_17Net_expr_mat.csv'))
desikan_mat_t = as.data.frame(t(desikan_mat[,1:70]))
colnames(desikan_mat_t) = desikan_mat$gene
desikan_mat_t$roi = as.character(rownames(desikan_mat_t))
desikan_mat_t[1:5,1:5]

# combine ENIGMA cohens d with gene matrix
enigma_ahba_df = merge(x=desikan_mat_t, y=enigma_mdd, by='roi')
for (gene in sst_genes){
    write(gene,'')
    enigma_gene = perm.sphere.p(x=enigma_ahba_df$cohens_d, y=enigma_ahba_df[[gene]], perm.id=desikan_perm_id, corr.type='spearman')
    corr = cor.test(x=enigma_ahba_df$cohens_d, y=enigma_ahba_df[[gene]], method='spearman')
    print(corr)
    print(enigma_gene[[1]])
}








