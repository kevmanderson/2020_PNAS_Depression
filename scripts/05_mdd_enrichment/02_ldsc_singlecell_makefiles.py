#!/bin/python

import os
import sys
import numpy as np
import scipy.io as io
import pandas as pd
import nibabel as nb
import subprocess
import time
import glob

def submitSlurm(cmd, dependencies=None):
    if dependencies != None:
        # execute
        p = subprocess.Popen(['sbatch', '--dependency=afterany:' + ':'.join(dependencies), cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    else:
        p = subprocess.Popen(['sbatch', cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    # get slurm job id to set up job submission dependency
    job_id = str(out).split(' ')[-1].replace("\\n'", '')
    return(job_id)



def writeSlurm(slurm_file, partition, cmd, jobName, stime='6:00:00', nthreads=None, mem=None):
    '''
    Submit batch job to Yale Milgram cluster (SLURM)

    required arguments:
        - slurm_file        base filepath string for writing slurm command and output
                                (e.g. /gpfs/milgram/project/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/MRI/100024/slurm/100024_hello_world_cmd)
        - partition         short/long/scavenge
                                (e.g. short)
        - nthreads          up to 28 cpus on a node
                                (e.g. 4)
        - cmd           full string of the command to be written/run
                            (e.g. module load Apps/FREESURFER/5.3.0\\n print('Hello World'))
    '''
    slurm_name = slurm_file + '_slurm.txt'
    slurm_out  = slurm_file + '_slurmOut.txt'
    slurm_file = open(slurm_name, "w")
    slurm_file.write('#!/bin/bash\n')
    slurm_file.write('#SBATCH --partition=' +  partition + '\n')
    slurm_file.write('#SBATCH --output=' + slurm_out + '\n')
    slurm_file.write('#SBATCH --nodes=1\n')
    if mem != None:
        slurm_file.write('#SBATCH --mem=' + str(mem) + '\n')
    if nthreads != None:
        slurm_file.write('#SBATCH --ntasks=1 --cpus-per-task=' + str(nthreads) + '\n')
    slurm_file.write('#SBATCH --job-name=' + jobName + '\n')
    slurm_file.write('#SBATCH --time=' + stime + '\n')
    slurm_file.write(str(cmd))
    slurm_file.close()
    subprocess.call(['chmod', '0770', slurm_name])
    return(slurm_name)


def annot_function(annot_txt_list, prefix, ldsc_dir, nbins):
    '''
    Make annotation files for LDSC single cell enrichment.
    Submits to SLURM automatically

    required arguments:
        - annot_txt_list    list of txt files that to annotate. each line is an ensembl ID

        - prefix            name of the folder for this enrichment run
                                (e.g. lake_dfc_bigcat_ldsc_files)

        - ldsc_dir          analysis directory
                                (e.g. /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc)

        - nbins             ...figure it out...
    '''


    # create directories if they dont exist
    folder_name = prefix + '_files'
    out_dir     = os.path.join(ldsc_dir, folder_name)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.join(out_dir, 'slurm'):
        os.mkdir(os.path.join(out_dir, 'slurm'))

    # iterate through each annotation file
    ct = 1
    for gene_set_file in annot_txt_list:
        print(gene_set_file)
        for chr in np.arange(1, 23):
            print(chr)
            chr = str(chr)

            # paths to relevant files
            gene_coord_file = os.path.join(ldsc_dir, folder_name, prefix + '_gene_coords.txt')
            eur_1000G_file  = os.path.join(ldsc_dir, '1000G_EUR_Phase3_plink/1000G.EUR.QC.' + chr)
            lake_annot_file = gene_set_file.replace('_genes.txt', '.') + chr
            hapmap_file     = os.path.join(ldsc_dir, 'hapmap3_snps/hm.' + chr + '.snp')

            # commands for LDSC annotation file generation
            ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/mdd_sst/external/ldsc/ \n'
            ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
            ldsc_cmd = ldsc_cmd + f'''python make_annot.py --gene-set-file {gene_set_file} --gene-coord-file {gene_coord_file} --windowsize 100000 --bimfile {eur_1000G_file}.bim --annot-file {lake_annot_file}.annot.gz \n\n'''
            ldsc_cmd = ldsc_cmd + f'''python ldsc.py --l2 --bfile {eur_1000G_file} --ld-wind-cm 1 --annot {lake_annot_file}.annot.gz --thin-annot --out {lake_annot_file} --print-snps {hapmap_file} \n\n'''

            slurm_file = os.path.join(ldsc_dir, folder_name, 'slurm', prefix + '_' + str(chr) + '_' + str(nbins))
            cmd_path = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=1, cmd=ldsc_cmd, stime='06:00:00',
                                  jobName=str(ct) + '_' + str(chr))
            submitSlurm(cmd_path, dependencies=None)


def make_ldcts_files(ldsc_dir, prefix, cell_list, control_files,nbins,mtg_bin=False):
    folder_name = prefix + '_files'
    ldcts_list = []
    for cell in cell_list:
        ldcts_path = os.path.join(ldsc_dir, folder_name, prefix + '_' + cell + '.ldcts')
        print(cell)
        with open(ldcts_path, 'w') as the_file:
            ldcts_list.append((cell, ldcts_path)  )
            if mtg_bin == True:
                path_to_cell_files = os.path.join(ldsc_dir, folder_name, prefix + '_nomatch.')
            else:
                path_to_cell_files = os.path.join(ldsc_dir, folder_name, prefix + '_nomatch_batch.')
            bin_name = 'nomatch_genes'
            the_file.write(bin_name + '\t' + path_to_cell_files + ',' + control_files + '\n')

            for ct in np.arange(1, nbins + 1):
                print(ct)
                ct = str(ct).zfill(2)
                if mtg_bin == True:
                    path_to_cell_files = os.path.join(ldsc_dir, folder_name,
                                                      prefix + '_' + cell + '_{}_{}.'.format(ct, nbins))
                else:
                    path_to_cell_files = os.path.join(ldsc_dir, folder_name,
                                                  prefix + '_' + cell + '_{}_{}_batch.'.format(ct, nbins))

                bin_name = '{}_bin_{}_{}'.format(cell, ct, nbins)
                print((bin_name + '\t' + path_to_cell_files + ',' + control_files + '\n'))
                the_file.write(bin_name + '\t' + path_to_cell_files + ',' + control_files + '\n')
    return(ldcts_list)



def run_ldsc_enrich(ldcts_list, ldsc_dir, prefix):

    folder_name = prefix + '_files'
    for row in ldcts_list:
        print(row)
        cell  = row[0]
        ldcts_path = row[1]

        #for cell in cell_list:

        # Wray
        name_out      = 'Wray_MDD'
        sumstats_file = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc/Wray2018_PGC_UKB_depr.sumstats.gz'

        # Howard
        #name_out      = 'Howard_PGC_UKB_depression_genome'
        #sumstats_file = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc/2019_PGC_UKB_depression_genome_wide_munge.sumstats.gz'

        baseline_path = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc/1000G_EUR_Phase3_baseline/baseline.'
        out_path      = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc/' + folder_name + '/' + prefix + '_' + name_out + '_' + cell + '_enrich'
        weights       = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc/weights_hm3_no_hla/weights.'

        ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/mdd_gene_expr/external/ldsc/ \n'
        ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
        ldsc_cmd = ldsc_cmd + f'''python ldsc.py --h2-cts {sumstats_file} --ref-ld-chr {baseline_path} --out {out_path} --ref-ld-chr-cts {ldcts_path} --w-ld-chr {weights}\n\n'''

        slurm_file = os.path.join(ldsc_dir, folder_name, 'slurm', prefix + '_' + cell + '_enrichment_run')
        cmd_path = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=1, cmd=ldsc_cmd, stime='06:00:00',
                              jobName=cell)
        submitSlurm(cmd_path, dependencies=None)


# set up directories
nbins          = 10
ldsc_dir       = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc'


# ----------
# DFC - big cat
# ----------
prefix         = 'lake_dfc_bigcat_ldsc'
folder_name    = prefix + '_files'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*bigcat*batch*txt')).sort()
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

# list of cells to run enrichment on
cell_list = list(np.unique([x.split('_ldsc_')[2].split('_')[0] for x in annot_txt_list]))
cell_list = [x for x in cell_list if 'nomatch' not in x and 'gene' not in x and 'bigcat' not in x]
cell_list = [x + '_bigcat' for x in cell_list if x != 'bigcat']
control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol.')

# create ldcts files needed for LDSC and run enrichments
ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------


# ----------
# DFC - fine cat
# ----------
prefix         = 'lake_dfc_superordinate_cat_ldsc'
folder_name    = prefix + '_files'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*superordinate*batch*txt'))
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

cell_list     = list(np.unique([x.split('_ldsc_')[2].split('_')[0] for x in annot_txt_list]))
cell_list     = [x for x in cell_list if 'nomatch' not in x and 'gene' not in x and 'superordinate' not in x and 'control' not in x]
control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol_batch.')

ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------




# ----------
# VIS - big cat
# ----------
prefix         = 'lake_vis_bigcat_ldsc'
folder_name    = prefix + '_files'
ldsc_dir       = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*bigcat*batch*txt'))
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

cell_list = list(np.unique([x.split('_ldsc_')[2].split('_')[0] for x in annot_txt_list]))
cell_list = [x for x in cell_list if 'nomatch' not in x and 'gene' not in x and 'bigcat' not in x and 'control' not in x]
cell_list = [x + '_bigcat' for x in cell_list if x != 'bigcat']
control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol_batch.')

ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------



# ----------
# VIS - fine categories
# ----------
prefix         = 'lake_vis_superordinate_cat_ldsc'
folder_name    = prefix + '_files'
ldsc_dir       = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/ldsc'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*superordinate*batch*txt'))
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

cell_list = list(np.unique([x.split('_ldsc_')[2].split('_')[0] for x in annot_txt_list]))
cell_list = [x for x in cell_list if 'nomatch' not in x and 'gene' not in x and 'superordinate' not in x and 'control' not in x]
control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol_batch.')

ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------




# ----------
# ABA MTG - big categories
# ----------
prefix         = 'aba_mtg_bigcat_ldsc'
folder_name    = prefix + '_files'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*bigcat*txt'))
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

annot_txt_list = [x for x in annot_txt_list if 'control' not in x and 'nomatch' not in x and 'gene_coords' not in x]
cell_list = list(np.unique([x.split('_bigcat_')[2].split('ldsc_')[1] for x in annot_txt_list]))
cell_list = [x + '_bigcat' for x in cell_list]

control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol.')

ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins, mtg_bin=True)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------





# ----------
# ABA MTG - fine categories
# ----------
prefix         = 'aba_mtg_finecat_ldsc'
folder_name    = prefix + '_files'
annot_txt_list = glob.glob(os.path.join(ldsc_dir, folder_name, '*finecat*txt'))
annot_txt_list.sort()

# annotate gene bins
annot_function(annot_txt_list=annot_txt_list, prefix=prefix, ldsc_dir=ldsc_dir, nbins=nbins)

annot_txt_tmp = [x for x in annot_txt_list if 'control' not in x and 'nomatch' not in x and 'gene_coords' not in x]
cell_list = list(np.unique([x.split('_finecat_')[2].replace('ldsc_','') for x in annot_txt_tmp]))
cell_list = list(np.unique(['_'.join(x.split('_')[0:2]) for x in cell_list]))

control_files = os.path.join(ldsc_dir, folder_name, prefix + '_allcontrol.')

ldcts_list = make_ldcts_files(ldsc_dir=ldsc_dir, prefix=prefix, cell_list=cell_list, control_files=control_files, nbins=nbins, mtg_bin=True)
run_ldsc_enrich(ldcts_list=ldcts_list, ldsc_dir=ldsc_dir, prefix=prefix)
# ----------





