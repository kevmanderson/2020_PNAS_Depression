library(tidyverse)
library(RColorBrewer)
library(Cairo)
library(effsize)

# 1. UKB depression status compared to neuroticism, online MDD, antidepressant count, polygenic MDD PRS


# Step 1: Read/Prepare data
# -----------
project_dir = '/gpfs/milgram/project/holmes/kma52/mdd_sst'
pheno_path  = paste0(project_dir, '/data/ukb/enc/ukb_pheno_processed_mri_FS_RSFA_GBC.csv')
image_df    = read_csv(pheno_path)
dim(image_df)
#image_df = read_csv(paste0(project_dir, '/data/ukb/ukb_severe_ctl.csv'))


# remove subjects without MDD status
image_df = image_df[which(!is.na(image_df$mdd_status.2.0)),]

# remove subjects with missing genetic ethnicity
image_df = image_df %>% filter(!is.na(genetic_ethnicity.0.0))

# calculate some covariates
image_df$brain_size = image_df$MRI_T1_volGM_WM.2.0 + image_df$MRI_T1_ventCSF.2.0
image_df$age_square_scale = image_df$age_at_scan^2
image_df$UK_Biobank_assessment_centre.2.0 = as.character(image_df$UK_Biobank_assessment_centre.2.0)

# imaging columns
thick_cols = colnames(image_df)[grep('_thickness_use', colnames(image_df))]
gbc_cols   = colnames(image_df)[grep('_gbc_use', colnames(image_df))]
rsfa_cols  = colnames(image_df)[grep('_rsfa_use', colnames(image_df))]


# frequency of MDD diagnoses by severity
mdd_freq     = data.frame(table(image_df$mdd_status.2.0)/length(image_df$mdd_status.2.0), mdd=c('CTL','Single','Moderate','Severe'))
mdd_freq$freq_rnd = round(mdd_freq$Freq,2)
mdd_freq$mdd = factor(mdd_freq$mdd, levels=c('CTL','Single','Moderate','Severe'))




# Step 2: Plot Frequency of MDD status
# -------------
freq_plot = paste0(project_dir, '/figures/ukbb_mdd_freq_by_category_190824.pdf')
freq_plot
CairoPDF(freq_plot, width=3, height=2)
ggplot(data=mdd_freq[2:4,], aes(x=mdd, y=freq_rnd, fill=mdd)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=freq_rnd), vjust=-0.3, size=3.5)+
    scale_fill_manual(values=c('#ADD3E6', "#2786C2", "#1B447D")) +
    theme_classic() +
    scale_y_continuous(limits = c(0, .12), breaks=seq(0,.12,by=.03), expand=c(0,0)) +
    theme(axis.text.x = element_text(color='black'), axis.text.y = element_text(color='black'))
dev.off()



# Step 3: Plot Neuroticism by MDD status
# -------------
image_df$mdd_status_2_0_fac = as.factor(image_df$mdd_status.2.0)
neu_plot = paste0(project_dir, '/figures/ukbb_mdd_by_neuroticism_190824.pdf')
neu_plot
CairoPDF(neu_plot, width=2, height=2)
ggplot(image_df, aes(x=mdd_status_2_0_fac, y=neuroticism.2.0, fill=mdd_status_2_0_fac)) +
    geom_boxplot(fill='#A4A4A4') +
    theme_classic() +
    scale_fill_manual(values=c("#FEFBDE", "#999999", "#E69F00", "#56B4E9")) +
    scale_y_continuous(limits = c(0, 15), breaks=seq(0,15,by=5), expand=c(0,0)) +
    theme(axis.text.x = element_text(color='black'), axis.text.y = element_text(color='black'))
dev.off()

# stats for neu group differences
summary(lm(neuroticism.2.0 ~ mdd_status_2_0_fac + sex_0_0 + age_at_scan + age_square + age_at_scan*sex_0_0 + age_square*sex_0_0, data=image_df))

# means/sd/se for neuroticism
mdd_neu_df = image_df %>%
                group_by(mdd_status.2.0) %>%
                summarize(neu_mean=mean(neuroticism.2.0,
                          neu_df=sd(neuroticism.2.0),
                          neu_se=sd(neuroticism.2.0)/sqrt(length(neuroticism.2.0))))


# Step 3: Plot UKB online MDD symptom
# -------------
omdd_plot = paste0(project_dir, '/figures/ukbb_mdd_by_onlineMDD_190824.pdf')
omdd_plot
CairoPDF(omdd_plot, width=2, height=2)
ggplot(image_df, aes(x=mdd_status_2_0_fac, y=mdd_online_symptom_sum.0.0, fill=mdd_status_2_0_fac)) +
    geom_boxplot(fill='#A4A4A4') +
    theme_classic() +
    scale_fill_manual(values=c("#FEFBDE", "#999999", "#E69F00", "#56B4E9")) +
    scale_y_continuous(limits = c(0, 10), breaks=seq(0,10,by=2), expand=c(0,0)) +
    theme(axis.text.x = element_text(color='black'), axis.text.y = element_text(color='black'))
dev.off()

# stats for omdd group differences
summary(lm(mdd_online_symptom_sum.0.0 ~ mdd_status_2_0_fac + sex_0_0 + age_at_scan + age_square + age_at_scan*sex_0_0 + age_square*sex_0_0, data=image_df))




# Step 4: Plot MDD Polygenic risk score by MDD status
# -------------
hyde_gwas = read_delim('/gpfs/milgram/project/holmes/kma52/ukb_pymood/data/prsice/PRSice_PGC_MDD_ex23andme2018/PRSice.all.score', delim=' ')
colnames(hyde_gwas) = gsub('[.]','_', paste0('hyde_', colnames(hyde_gwas)))
image_df_hyde = merge(x=image_df, y=hyde_gwas, by.x='UKB_ID', by.y='hyde_IID')
image_df_hyde$hyde_0_100000 = as.numeric(image_df_hyde$hyde_0_100000)

freq_plot = paste0(project_dir, '/figures/mdd_by_HydePRS.pdf')
freq_plot
CairoPDF(freq_plot, width=2, height=2)
ggplot(image_df_hyde, aes(x=mdd_status_2_0_fac, y=hyde_0_100000, fill=mdd_status_2_0_fac)) +
    geom_boxplot(fill='#A4A4A4') +
    theme_classic() +
    scale_fill_manual(values=c("#FEFBDE", "#999999", "#E69F00", "#56B4E9")) +
    scale_y_continuous(limits = c(-2e-4, 2e-4), breaks=seq(-2e-4,2e-4,by=1e-4), expand=c(0,0)) +
    theme(axis.text.x = element_text(color='black'), axis.text.y = element_text(color='black'))
dev.off()
summary(lm(hyde_0_100000~ mdd_status_2_0_fac + sex_0_0 + age_at_scan + age_square + age_at_scan*sex_0_0 + age_square*sex_0_0, data=image_df_hyde))




# Step 5: Antidepressant Rx by MDD status
# -------------
# UKBB medication code
rx_table        = read_delim(paste0(project_dir, '/reference_files/ukbb_rx_codes.tsv'), delim='\t')
wray_rx_table   = read.csv(paste0(project_dir, '/reference_files/ukbb_wray_med_class.csv'))
antidepressants = wray_rx_table[grep('N06A', wray_rx_table$Medication_ATC_code),]
rx_cols = colnames(image_df)[grep('Rx_code.2', colnames(image_df))]

# determine if subject taking antidepressant
rx_grep  = lapply(rx_cols, function(x) image_df[[x]] %in% antidepressants$Coding)
rx_table = do.call(cbind, rx_grep)
antidepr_ct = rowSums(rx_table)

# binary - taking anti-depressant or not
image_df$antidepr_ct = antidepr_ct
image_df$antidepr_bin = image_df$antidepr_ct
image_df$antidepr_bin[image_df$antidepr_bin > 0] = 1

# calc antidepr by group
antidepr_by_group = image_df %>% group_by(mdd_status_2_0_fac) %>% summarise(rx_use=mean(antidepr_bin))
antidepr_by_group$freq_rnd = round(antidepr_by_group$rx_use,2)

rx_plot = paste0(project_dir, '/figures/mdd_by_antidepr_freq_190824.pdf')
rx_plot
CairoPDF(rx_plot, width=4, height=2)
ggplot(data=antidepr_by_group, aes(x=mdd_status_2_0_fac, y=rx_use, fill=mdd_status_2_0_fac)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=freq_rnd), vjust=-0.3, size=3.5)+
    scale_fill_manual(values=c("#FEFBDE",'#ADD3E6', "#2786C2", "#1B447D")) +
    theme_classic() +
    scale_y_continuous(limits = c(0, .3), breaks=seq(0,.3,by=.1), expand=c(0,0)) +
    theme(axis.text.x = element_text(color='black'), axis.text.y = element_text(color='black'))
dev.off()

# write dataframe
#pheno_out_path = paste0(project_dir, '/data/ukb/enc/ukb_pheno_processed_mri_FS_RSFA_GBC_step2.csv')
#write_csv(image_df, pheno_out_path)





