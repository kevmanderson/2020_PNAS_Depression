#!/bin/python

import os
import glob
import numpy as np
import pandas as pd
import progressbar
import multiprocessing as mp
from functools import partial
from itertools import product
from itertools import combinations


# this script will read previously derived RSFA/GBC data into a single larger dataframe


# function to read RSFC data
def read_rsfc(rest_path, scan_string):
    sub = rest_path.split('/')[-1]
    print(sub)
    rsfc_path = os.path.join(rest_path, 'surf', sub + scan_string)
    if os.path.exists(rsfc_path):
        rsfa_dat  = pd.read_csv(rsfc_path, header=None, sep='\t')
        conn_flat = pd.DataFrame([rsfa_dat.values.flatten()])
        conn_flat['UKB_ID'] = sub
        return(conn_flat)



# Step 1: Set up paths/data lists/parcel info
# ------------
# project dir
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'


# dir with CBIG processed rest
proj_dir  = '/gpfs/milgram/data/UKB/REPOSITORY/cbig_rest'
rest_subs = glob.glob(os.path.join(proj_dir, '*'))


# parcellation
schaef_dir    = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
parc_label_in = pd.read_csv(os.path.join(schaef_dir, 'Schaefer2018_200Parcels_17Networks_order_info.txt'), header=None)
parc_labels   = [x for x in parc_label_in[0] if 'Net' in x]
print(parc_labels[1:5])



# Step 2: Read RSFA data from subject's directory
# --------------
check_subs = []
rsfa_subs  = []
rsfa_df    = pd.DataFrame()
pbar       = progressbar.ProgressBar(widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()], term_width=50)
for rest_path in pbar(rest_subs):
    sub       = rest_path.split('/')[-1]
    rsfa_path = os.path.join(rest_path, 'surf', sub + '_bld001_rest_mc_residc_interp_FDRMS0.3_DVARS75_fs6_sm2_fs5_censor_Schaef2018_200_17Net_fsavg5_rsfa.txt')
    if os.path.exists(rsfa_path):
        rsfa_dat = pd.read_csv(rsfa_path, header=None)[0]
        rsfa_row = pd.DataFrame([rsfa_dat.values.tolist()], columns=parc_labels)
        rsfa_row['UKB_ID'] = sub
        rsfa_df  = pd.concat([rsfa_df,rsfa_row])
    else:
        check_subs.append(rest_path)

# save RSFA dataframe, summarized by Schaefer parcellation
rsfa_df.to_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfa_Schaef2018_200_17Net_df.csv'), index=False)



# Step 3: Read RSFC data
# ---------------
# the data are large, so we read them in parallel and by chunks
rest_subs.sort()
scan_string = '_bld001_rest_mc_residc_interp_FDRMS0.3_DVARS75_bp_0.009_0.08_fs6_sm2_fs5_censor_Schaefer2018_200parc_17net_fsavg5_pconn.txt'
pool     = mp.Pool(24)
rsfc_out = pool.map(partial(read_rsfc, scan_string=scan_string), rest_subs[0:5000])
rsfc_mat = pd.concat(rsfc_out)
rsfc_mat.to_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n0-5000.csv'), index=False)


rest_subs.sort()
pool     = mp.Pool(24)
rsfc_out = pool.map(partial(read_rsfc, scan_string=scan_string), rest_subs[5000:10000])
rsfc_mat = pd.concat(rsfc_out)
rsfc_mat.to_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n5000-10000.csv'), index=False)


rest_subs.sort()
pool     = mp.Pool(24)
rsfc_out = pool.map(partial(read_rsfc, scan_string=scan_string), rest_subs[10000:15000])
rsfc_mat = pd.concat(rsfc_out)
rsfc_mat.to_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n10000-15000.csv'), index=False)


rest_subs.sort()
pool     = mp.Pool(24)
rsfc_out = pool.map(partial(read_rsfc, scan_string=scan_string), rest_subs[15000:])
rsfc_mat = pd.concat(rsfc_out)
rsfc_mat.to_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n15000-end.csv'), index=False)


# Step 4: Map matrix indices to flattened values
# ---------
# read schaefer parcellation order
Schaefer200_17    = pd.read_csv(os.path.join(base_dir, 'reference_files/Schaefer2018_200Parcels_17Networks_order_info.txt'), header=None)
net17_200_ordered = np.array(Schaefer200_17[0].loc[Schaefer200_17[0].str.contains('17Net')])

# create a dataframe to map parcels to flattened connectivity matrix
net_matrix    = np.chararray((len(net17_200_ordered), len(net17_200_ordered)))
net_matrix[:] = 'na'

#
x, y = net_matrix.shape
x_, y_ = zip(*product(range(x), range(y)))

# create a dataframe with mapping info
df        = pd.DataFrame(net_matrix.flatten()).assign(x=x_, y=y_)
df['idx'] = np.arange(0, df.shape[0])
index_arr = [z for z in zip(df['x'], df['y'])]
df['idx_zip'] = index_arr
df['new_col_name'] = ['_'.join(map(str,x)) for x in zip(df['x'],df['y'])]

uniq_combos = [comb for comb in combinations(np.arange(200),2)]
use_indices = df.loc[df['idx_zip'].isin(uniq_combos),'idx']
new_names   = np.array(df.loc[use_indices,'new_col_name'].astype(str))

parcel_coords = list()
for edge in np.arange(0, 200):
    print(edge)
    zips = [x for x in uniq_combos if edge in x]
    idxs = df.loc[df['idx_zip'].isin(zips), 'idx']
    parcel_coords.append( idxs )


# identify which columns to take
rsfc_paths = []
rsfc_paths.append(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n0-5000.csv'))
rsfc_paths.append(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n5000-10000.csv'))
rsfc_paths.append(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n10000-15000.csv'))
rsfc_paths.append(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n15000-end.csv'))


# Step 5: Calculate Global Brain Connectivity for each subject
# ------------
parcel_headers = ','.join(net17_200_ordered)+',UKB_ID'
for rsfc_path in rsfc_paths:
    print(rsfc_path)

    # open file for smaller csv
    rsfc_small_path = rsfc_path.replace('.csv', '_gbc.csv')
    f_out = open(rsfc_small_path, 'w')
    f_out.write(parcel_headers+'\n')

    with open(rsfc_path, 'r') as fp:
        cnt  = 1
        line = fp.readline().split(',')
        while line:
            print(cnt)
            line  = fp.readline()
            if line != '':
                edges = np.array(line.split(','))
                ukb_id = edges[len(edges) - 1]
                avg_conn_arr = []
                for edge in np.arange(0, 200):
                    idxs     = parcel_coords[edge]
                    edge_arr = edges[df.loc[idxs, 'idx']]
                    edge_arr[edge_arr == ''] = np.nan
                    parcel_conn = np.float32(edge_arr)
                    avg_conn_arr.append(np.nanmean(parcel_conn) )
                sub_gbc_est  = ','.join(np.array(avg_conn_arr).astype(str)) + ',' + ukb_id
                f_out.write(sub_gbc_est)
                cnt+=1
    f_out.close()



# make the RSFC files smaller
for rsfc_path in rsfc_paths:
    print(rsfc_path)

    # open file for smaller csv
    rsfc_small_path = rsfc_path.replace('.csv', '_small.csv')
    f_out = open(rsfc_small_path, 'w')
    with open(rsfc_path) as fp:

        # write new header
        line   = fp.readline().split(',')
        header_write = ','.join(new_names) + ',UKB_ID\n'
        f_out.write(header_write)
        cnt = 1
        while line:
            print(cnt)
            line     = fp.readline()
            if line != '':
                rsfc     = np.array(line.split(','))
                ukb_id   = rsfc[len(rsfc) - 1]
                use_rsfc = rsfc[use_indices]
                use_rsfc[np.where(use_rsfc == '')] = np.nan # if there are missing values
                use_rsfc = np.float32(use_rsfc).astype(str)
                rsfc_write = ','.join(use_rsfc) + ',' + ukb_id
                f_out.write(rsfc_write)
                cnt += 1

    f_out.close()


# open an interactive session with lots of RAM and read all the data
rsfc_1 = pd.read_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n0-5000.csv'))
rsfc_2 = pd.read_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n5000-10000.csv'))
rsfc_3 = pd.read_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n10000-15000.csv'))
rsfc_4 = pd.read_csv(os.path.join(base_dir, 'data/ukb/dataframes/rsfc_conn_Schaef2018_200_17Net_df_n15000-end.csv'))


