#!/bin/python

import os
import glob
import pandas as pd
import numpy as np
from sklearn import preprocessing
from scipy import stats
import statsmodels.regression as smr
import statsmodels.formula.api as smf
import cifti
import copy
from nibabel import cifti2



# Stitch all of the imaging data / UKB phenotypes together


# Step 1: Read UKB phenotype
# -----------
# project dir
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'

# read phenotype data
pheno_path = os.path.join(base_dir, 'data/ukb/enc/ukb_pheno_processed_mri.csv')
pheno_df   = pd.read_csv(pheno_path)
pheno_df['UKB_ID'] = pheno_df['f.eid'].astype(int).astype(str)



# Step 2: Read imaging dataframes
# -----------
# read lesion data
lesion_df = pd.read_csv(os.path.join(base_dir, 'data/ukb/dataframes/lesion_df.csv'))
lesion_df['UKB_ID'] = lesion_df['UKB_ID'].astype(int).astype(str)

# read RSFA data
rsfc_df_1  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n0-5000_gbc.csv'))
rsfc_df_2  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n10000-15000_gbc.csv'))
rsfc_df_3  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n15000-end_gbc.csv'))
rsfc_df_4  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n5000-10000_gbc.csv'))
rsfc_df = rbind(rsfc_df_1, rsfc_df_2, rsfc_df_3, rsfc_df_4)


# read anatomical data
rsfc_df_1  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n0-5000_gbc.csv'))
rsfc_df_2  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n10000-15000_gbc.csv'))
rsfc_df_3  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n15000-end_gbc.csv'))
rsfc_df_4  = read_csv(paste0(base_dir, '/data/ukb_imaging/rsfc_conn_Schaef2018_200_17Net_df_n5000-10000_gbc.csv'))
rsfc_df = rbind(rsfc_df_1, rsfc_df_2, rsfc_df_3, rsfc_df_4)


# read thickness
lh_thick  = read_delim(paste0(base_dir, '/data/ukb_imaging/lh_t1_Schaefer2018_200Parcels_17Networks_sum_table_thickness.txt'), delim='\t')
lh_thick['UKB_ID'] = lh_thick['lh.t1_Schaefer2018_200Parcels_17Networks.sum.thickness']
lh_thick = lh_thick[,grep('_LH_|UKB_ID', colnames(lh_thick))]

rh_thick  = read_delim(paste0(base_dir, '/data/ukb_imaging/rh_t1_Schaefer2018_200Parcels_17Networks_sum_table_thickness.txt'), delim='\t')
rh_thick['UKB_ID'] = rh_thick['rh.t1_Schaefer2018_200Parcels_17Networks.sum.thickness']
rh_thick  = rh_thick[,grep('_RH_|UKB_ID', colnames(rh_thick))]

# combine LH/RH
thick_df   = lh_thick_df.merge(rh_thick_df, on='UKB_ID')
thick_df['UKB_ID'] = thick_df['UKB_ID'].astype(int).astype(str)



# Step 3: Merge everything
# -----------
# merge everything
pheno_subset_df = pheno_df.loc[pheno_df['UKB_ID'].isin(thick_df['UKB_ID'])]
pheno_mri_df    = pheno_subset_df.merge(thick_df, on='UKB_ID')
pheno_mri_df    = pheno_mri_df.merge(gbc_df, on='UKB_ID')
pheno_mri_df    = pheno_mri_df.merge(rsfa_df, on='UKB_ID')
pheno_mri_df    = pheno_mri_df.merge(lesion_df, on='UKB_ID')


rsfa_cols  = pheno_mri_df.columns[pheno_mri_df.columns.str.contains('rsfa')]
thick_cols = pheno_mri_df.columns[pheno_mri_df.columns.str.contains('thickness')]
gbc_cols   = pheno_mri_df.columns[pheno_mri_df.columns.str.contains('gbc')]



# Step 4: Prepare covariates
# ------------
# prepare covariates
pheno_mri_df['bmi_imaging'] = pheno_mri_df['weight_preimaging.2.0']/((pheno_mri_df['height_12144.2.0']/100)**2)

# median value replacement for missing diastolic/systolic BP
pheno_mri_df.loc[np.isnan(pheno_mri_df['systolic.2.0']), 'systolic.2.0'] = pheno_mri_df.loc[np.isnan(pheno_mri_df['systolic.2.0']), 'systolic.0.0']
pheno_mri_df.loc[np.isnan(pheno_mri_df['systolic.2.0']),'systolic.2.0'] = np.nanmedian(pheno_mri_df['systolic.2.0'])

pheno_mri_df.loc[np.isnan(pheno_mri_df['diastolic.2.0']), 'diastolic.2.0'] = pheno_mri_df.loc[np.isnan(pheno_mri_df['diastolic.2.0']), 'diastolic.0.0']
pheno_mri_df.loc[np.isnan(pheno_mri_df['diastolic.2.0']), 'diastolic.2.0'] = np.nanmedian(pheno_mri_df['diastolic.2.0'])


# inspect covariates
thick_quant_covars = ['age_at_scan',
                      'bmi_imaging',
                      'height_12144.2.0',
                      'weight_preimaging.2.0',
                      'lesion_vol',
                      'MRI_T1_invSNR.2.0',
                      'diastolic.2.0',
                      'systolic.2.0',
                      'MRI_T1_vetntCSF_norm.2.0',
                      'MRI_T1_volGM_WM_headNorm.2.0',
                      'scanner_X_pos.2.0',
                      'scanner_Z_pos.2.0',
                      'scanner_Y_pos.2.0']

# merge
pheno_mri_regr_df = pheno_mri_df.loc[~pheno_mri_df['lesion_vol'].isna()]
np.sum(pheno_mri_regr_df[thick_quant_covars].isna())
#pheno_mri_regr_df = pheno_mri_regr_df.reindex(np.arange(pheno_mri_regr_df.shape[0]))


# recode some vars
pheno_mri_regr_df['sex_0_0'] = pheno_mri_regr_df['sex.0.0']
pheno_mri_regr_df['UK_Biobank_assessment_centre_2_0'] = pheno_mri_regr_df['UK_Biobank_assessment_centre.2.0'].astype(int).astype(str)

# get MDD status at visit 1
pheno_mri_regr_df['mdd_status.0.0'] = copy.deepcopy(pheno_mri_regr_df['bp_mdd_status.0.0'])
pheno_mri_regr_df.loc[pheno_mri_regr_df['mdd_status.0.0'] == 4, 'mdd_status.0.0'] = np.nan
pheno_mri_regr_df.loc[pheno_mri_regr_df['mdd_status.0.0'] == 5, 'mdd_status.0.0'] = np.nan


# binarize CTL vs moderate/severe MDD
pheno_mri_regr_df['mdd_status_2_0_bin'] = 0
pheno_mri_regr_df.loc[pheno_mri_regr_df['mdd_status.2.0'] > 1, 'mdd_status_2_0_bin'] = 1



# Step 5: Save
# ----------
pheno_path = os.path.join(base_dir, 'data/ukb/enc/ukb_pheno_processed_mri_FS_RSFA_GBC.csv')
pheno_mri_regr_df.to_csv(pheno_path, index=False)
pheno_mri_regr_df = pd.read_csv(os.path.join(base_dir, 'data/ukb/enc/ukb_pheno_processed_mri_FS_RSFA_GBC.csv'))



