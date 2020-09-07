#!/bin/python

import os
import shutil
import feather
import glob
import zipfile
import pandas as pd
import numpy as np
import nibabel
import scipy.stats as stats
import multiprocessing as mp
from itertools import repeat
import subprocess
import progressbar
import sys
from collections import defaultdict

sys.path.append('/gpfs/milgram/project/holmes/kma52/ukb_pymood/scripts')
from utilities.slurm_util import writeSlurm, submitSlurm


def force_symlink(src, dest):
    try:
        os.symlink(src, dest)
    except:
        os.remove(dest)
        os.symlink(src, dest)


def mri_anatomical_stats(lh_annot, rh_annot, out_file):
    '''

    :param parcellation:
    :param measure:
    :return:
    '''

    # load packages
    base_cmd = 'module load FreeSurfer/6.0.0 \n'
    base_cmd = base_cmd + '. /gpfs/milgram/apps/hpc.rhel7/software/FreeSurfer/6.0.0/SetUpFreeSurfer.sh \n'
    base_cmd = base_cmd + 'module load Python/Anaconda2 \n\n'
    base_cmd = base_cmd + 'module load FSL/5.0.10 \n'
    base_cmd = base_cmd + '. /gpfs/milgram/apps/hpc.rhel7/software/FSL/5.0.10-centos7_64/etc/fslconf/fsl.sh \n\n'
    base_cmd = base_cmd + 'export SUBJECTS_DIR='+str(fs_dir)+' \n\n'

    subject  = lh_annot.split('/')[-3]
    seg_cmd  = 'mris_anatomical_stats -a ' + lh_annot + ' -f ' + sum_out.replace('t1','lh.t1') + ' ' + subject + ' lh '
    base_cmd = base_cmd + seg_cmd + '\n\n'
    seg_cmd  = 'mris_anatomical_stats -a ' + rh_annot + ' -f ' + sum_out.replace('t1','rh.t1') + ' ' + subject + ' rh '
    base_cmd = base_cmd + seg_cmd + '\n\n'

    return(base_cmd)


def combine_sum_files(sum_file_list, output_table):
    # volumetric aseg table
    print('volumetric stats')
    slurm_cmd = os.path.join(out_dir, 'slurm_aseg_table.txt')
    hemi_aseg_table = os.path.join(out_dir, 'aseg_table.txt')
    convert_cmd = cmd_base + '/gpfs/milgram/apps/hpc.rhel7/software/FreeSurfer/6.0.0/bin/asegstats2table --subjectsfile ' + fs_sub_list + ' --meas volume --tablefile ' + hemi_aseg_table
    print(convert_cmd)
    cmd_path = writeSlurm(slurm_file=slurm_cmd, partition='short', nthreads=1, cmd=convert_cmd, stime='01:00:00', jobName=hemi + '_aparc_convert')
    job1_id  = submitSlurm(cmd_path, dependencies=None)
    job_ids.append(job1_id)


def main(sub):

    # Step 1: Set up paths
    # ----------
    proj_dir   = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
    mri_dir    = '/gpfs/milgram/data/UKB/REPOSITORY/MRI'
    fs_dir     = '/gpfs/milgram/data/UKB/REPOSITORY/FS_v6_T1only'
    fs_out     = os.path.join(proj_dir, 'data/ukb/freesurfer/surf_files')
    schaef_dir = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6'
    out_dir    = fs_out

    # identify subjects with freesurfer data to process
    fs_sub_list = os.path.join(proj_dir, 'data/ukb/freesurfer/fs_list.txt')
    fs_subs     = pd.read_csv(fs_sub_list, header=None)
    fs_subs     = fs_subs[0]

    for parcnum in [200]:
        for sub in fs_subs:
            print(str(sub) + ' : ' + str(parcnum))

            # current parcellation
            parc_string = 'Schaefer2018_' + str(parcnum) + 'Parcels_17Networks'

            # sub MRI dir
            sub_mri = os.path.join(mri_dir, str(sub))
            if not os.path.exists(sub_mri):
                os.mkdir(sub_mri)

            slurm_mri = os.path.join(sub_mri, 'slurm')
            if not os.path.exists(slurm_mri):
                os.mkdir(slurm_mri)

            # where to place the summary stats
            sum_dir = os.path.join(out_dir, str(sub), 'anat_sum_stats')
            if not os.path.exists(sum_dir):
                os.mkdir(sum_dir)

            # subject specific parcellation
            lh_annot = os.path.join(fs_dir, str(sub), 'label', 'lh.' + parc_string + '.annot')
            rh_annot = os.path.join(fs_dir, str(sub), 'label', 'rh.' + parc_string + '.annot')
            cur_parc = os.path.join(fs_dir, str(sub), 'label', parc_string + '.nii.gz')
            if os.path.exists(cur_parc):
                # mri_segstats(parc=cur_parc, measure)

                # use freesurfer utilities to summarized stats according to Scheafer parcellation
                sum_out   = os.path.join(sum_dir, 't1_' + parc_string + '.sum')
                seg_cmd   = mri_anatomical_stats(lh_annot=lh_annot, rh_annot=rh_annot, out_file=sum_out)
                slurm_out = os.path.join(out_dir, str(sub), 'slurm', 't1_' + str(sub) + '_' + parc_string)
                cmd_path  = writeSlurm(slurm_file=slurm_out, partition='short', nthreads=1, cmd=seg_cmd,
                                      stime='06:00:00', jobName='t2_flair_' + str(sub))
                submitSlurm(cmd_path, dependencies=None)


list_dict = defaultdict(list)
ct = 0
for parcnum in [200]:
    for sub in fs_subs:
        print('ct='+str(ct)+' '+str(sub) + ' : ' + str(parcnum))
        ct+=1
        # current parcellation
        parc_string = 'Schaefer2018_' + str(parcnum) + 'Parcels_17Networks'

        # sub MRI dir
        sub_mri = os.path.join(mri_dir, str(sub))

        # where to place the summary stats
        sum_dir = os.path.join(out_dir, str(sub), 'anat_sum_stats')

        # collect flair sum files
        for hemi in ['lh','rh']:
            dict_key = hemi + '_t1_' + parc_string + '.sum'
            sum_out  = os.path.join(sum_dir, hemi + '.t1_' + parc_string + '.sum')
            if os.path.exists(sum_out):
                list_dict[dict_key].append(sum_out)


# prepare T1 data file structure
sub_list = pd.DataFrame([x.split('/')[-3] for x in list_dict['lh_t1_Schaefer2018_200Parcels_17Networks.sum']])
sub_arr  = np.array(sub_list[0])
sub_arr.sort()
for sub in sub_arr:
    print(sub)

    # the aparc summation requires the data to be in a {sub}/stats dir. do that
    sub_stats_dir = os.path.join(proj_dir, 'data/freesurfer/surf_files/', str(sub), 'stats')
    if not os.path.exists(sub_stats_dir):
        os.mkdir(sub_stats_dir)

    cp_files = glob.glob(os.path.join(proj_dir, 'data/freesurfer/surf_files/',str(sub),'anat_sum_stats/*t1*'))
    for f in cp_files:
        new_file = f.replace('anat_sum_stats', 'stats') + '.stats'
        shutil.copy(f, new_file)


table_out_dir = os.path.join(proj_dir, 'data/freesurfer/summary_stats/anat_features')
slurm_dir     = os.path.join(proj_dir, 'data/freesurfer/summary_stats/slurm')

for key in ['lh_t1_Schaefer2018_200Parcels_17Networks.sum', 'rh_t1_Schaefer2018_200Parcels_17Networks.sum']:
    sum_files = list_dict[key]
    sub_list  = pd.DataFrame([x.split('/')[-3] for x in sum_files])
    print(key + ' : ' + str(len(sum_files)))
    sum_file_path = os.path.join(table_out_dir, key.replace('.sum', '_sum_filepaths.txt'))
    sub_list.to_csv(sum_file_path, index=None, header=None)

    # where to write data
    base_cmd = 'export SUBJECTS_DIR=/gpfs/milgram/project/holmes/kma52/mdd_sst/data/freesurfer/surf_files/ \n\n'
    base_cmd = base_cmd + 'module load Python/Anaconda2 \n\n'
    slurm_file = os.path.join('/gpfs/milgram/project/holmes/kma52/mdd_sst/data/freesurfer/slurm', key.replace('.sum',''))

    # lh/rh signifies whether it is T1
    hemi = key.split('_')[0]
    for meas in ['area','volume','thickness','meancurv']:
        slurm_file = os.path.join('/gpfs/milgram/project/holmes/kma52/mdd_sst/data/freesurfer/slurm', meas + '_' +
                                  key.replace('.sum', ''))
        table_out   = os.path.join(table_out_dir, key.replace('.sum', '_sum_table_' + meas + '.txt'))
        convert_cmd = base_cmd + '/gpfs/milgram/apps/hpc.rhel7/software/FreeSurfer/6.0.0/bin/aparcstats2table --parc=' + key.replace('lh_','').replace('rh_','') + ' --hemi=' + hemi + ' --subjectsfile=' + sum_file_path + ' --tablefile ' + table_out + ' --meas ' + meas
        print(convert_cmd)
        cmd_path    = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=1, cmd=convert_cmd, stime='06:00:00', jobName=meas+'_'+key)
        #os.system(convert_cmd)
        submitSlurm(cmd_path, dependencies=None)

    print(convert_cmd)


if __name__ == "__main__":
    main()




