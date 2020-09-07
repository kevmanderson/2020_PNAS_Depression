#!/bin/bash


# set up directories
base_dir=/gpfs/milgram/project/holmes/kma52/mdd_gene_expr
gene_dir=/gpfs/milgram/data/UKB/REPOSITORY/GWAS/NON_IMPUTED
plink_bin=${base_dir}/external/plink_1.09/plink

# subject file
sub_file=${base_dir}/data/prsice/ukb_subj_list.txt

# quality control variants
bfile_list_file=${base_dir}/data/prsice/plink_bfile_list.txt
rm ${bfile_list_file}
touch ${bfile_list_file}
for chr in {1..22}
do
    echo ${chr}
    out_file=${base_dir}/data/prsice/ukb_cal_chr${chr}_v2_maf05_hwe1e6_geno1_mind1
    echo ${out_file} >> ${bfile_list_file}
    ${plink_bin} --bfile ${gene_dir}/ukb_cal_chr${chr}_v2 \
                    --keep ${sub_file} \
                    --maf .05 \
                    --hwe 1e-6 \
                    --geno 0.1 \
                    --mind 0.1 \
                    --make-bed --out ${out_file}
    echo ${out_file} >> ${bed_list_file}
done

all_file=${base_dir}/data/prsice/ukb_cal_chrALL_v2_maf05_hwe1e6_geno1_mind1
${plink_bin} --merge-list ${bfile_list_file} --make-bed --out ${all_file}


#
grm_list_file=${base_dir}/data/prsice/plink_grm_list.txt
rm ${grm_list_file}
touch ${grm_list_file}

for chr in {1..22}
do
echo ${chr}

# set up paths
grm_slurm_file=${base_dir}/data/prsice/slurm/gcta_make_grm_chr${chr}.txt
grm_slurmOut_file=${base_dir}/data/prsice/slurm/gcta_make_grm_chr${chr}_out.txt
out_file=${base_dir}/data/prsice/ukb_cal_chr${chr}_v2_maf05_hwe1e6_geno1_mind1

echo ${out_file} >> ${grm_list_file}

cat >${grm_slurm_file} <<EOL
#!/bin/bash
#SBATCH --partition=short
#SBATCH --output=${grm_slurmOut_file}
#SBATCH --nodes=1
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --job-name=chr${chr}_gcta
#SBATCH --time=06:00:00

${base_dir}/external/gcta_1.91.1beta/gcta64 --bfile ${out_file} --autosome --make-grm --out ${out_file} --thread-num 10
EOL

sbatch ${grm_slurm_file}

done

grm_all=${base_dir}/data/prsice/ukb_cal_chrALL_v2_maf05_hwe1e6_geno1_mind1
${base_dir}/external/gcta_1.91.1beta/gcta64 --mgrm ${grm_list_file} --make-grm --out ${grm_all}


# calculate principle components
${base_dir}/external/gcta_1.91.1beta/gcta64 --grm ${grm_all} --pca 20 --out ${grm_all}_eigen --thread-num 10


# prsice calculation

gwas_base=${base_dir}/data/ldsc/2019_PGC_UKB_depression_genome_wide.txt
target=${base_dir}/data/prsice/ukb_cal_chrALL_v2_maf05_hwe1e6_geno1_mind1

${base_dir}/external/prsice/PRSice_linux \
    --base ${gwas_base} \
    --target ${target} \
    --thread 1 \
    --stat OR \
    --binary-target T \
    --fastscore --no-regress --missing CENTER


gwas_base=/gpfs/milgram/project/holmes/kma52/ukb_pymood/data/gwas_sumstats/PRSice_PGC_MDD_ex23andme2018.txt
target=${base_dir}/data/prsice/ukb_cal_chrALL_v2_maf05_hwe1e6_geno1_mind1

${base_dir}/external/prsice/PRSice_linux \
    --base ${gwas_base} \
    --target ${target} \
    --thread 1 \
    --stat OR \
    --binary-target T \
    --fastscore --no-regress --missing CENTER



