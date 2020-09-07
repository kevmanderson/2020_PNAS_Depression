#!/bin/bash



##############
# set up paths
##############
base_dir=/gpfs/milgram/project/holmes/kma52/mdd_gene_expr
magma=${base_dir}/external/magma

#scz_gwas_snps=/gpfs/milgram/project/holmes/kma52/mdd_sst/data/ldsc_sc/Wray2018_PGC_UKB_depression_genome_wide.txt
#scz_gwas_snps=/gpfs/milgram/project/holmes/kma52/mdd_sst/data/gwas/magma_Howard_PGC_UKB_depression_genome.txt



# set up paths
#snp_loc_file=/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/data/magma/MDD2018_ex23andMe
#snp_loc_file=${base_dir}/data/magma/Wray2018_PGC_UKB_depression_genome_wide.txt
#genome_build=${base_dir}/data/magma/NCBI37.3.gene.loc
#annot_out=${base_dir}/data/magma/Wray2018_PGC_UKB_depression_genome_wide
# Howard MDD GWAS
#gwas_snps=${base_dir}/data/magma/magma_Howard_PGC_UKB_depression_genome.txt
#annot_out=${base_dir}/data/magma/magma_Howard_PGC_UKB_depression_genome


# GWAS related files
snp_loc_file=${base_dir}/data/magma/MDD2018_ex23andMe
gwas_snps=${base_dir}/data/magma/magma_MDD2018_ex23andMe.txt
annot_out=${base_dir}/data/magma/magma_MDD2018_ex23andMe
ref_bfile=${base_dir}/data/magma/g1000_eur


##########################
# MAGMA #1: annotate genes
##########################
${magma} --annotate window=5,5 \
            --snp-loc ${gwas_snps} \
            --gene-loc ${genome_build} \
            --out ${annot_out}


# Run MAGMA
#scz_gwas=/gpfs/milgram/project/holmes/kma52/mdd_sst/data/gwas/magma_Howard_PGC_UKB_depression_genome.txt
#gwas_snps=${base_dir}/data/magma/Wray2018_PGC_UKB_depr_magma_USE.txt


##############################
# MAGMA #2: Association Tests
##############################
${magma} --bfile ${ref_bfile} \
            --pval ${gwas_snps}  ncol=N \
            --gene-annot ${annot_out}.genes.annot \
            --out ${annot_out}





