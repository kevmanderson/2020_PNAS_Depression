library(tidyverse)
library(ggcorrplot)
library(Cairo)
library(psych)


# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'


# read gene-wise spatial correlations to MDD effect maps
all_ahba_cors = read_csv(paste0(base_dir, "/output/mdd_regr/all_mdd_ahba_genecors.csv"))


# reverse direction of RSFA to line up with other modalities
all_ahba_cors$ukbb_rsfa_cor_new = all_ahba_cors$ukbb_rsfa_cor*-1
all_ahba_cors$ukbb_rsfa_cor     = NULL
all_ahba_cors$gsp_rsfa_cor_new  = all_ahba_cors$gsp_rsfa_cor*-1
all_ahba_cors$gsp_rsfa_cor      = NULL


# average together and sort
all_ahba_cors$avg = rowMeans(all_ahba_cors[grep('cor', colnames(all_ahba_cors))])
all_ahba_cors     = all_ahba_cors[order(all_ahba_cors$avg),]


# load diff expr
table_dir  = '/gpfs/milgram/project/holmes/kma52/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master/results/results/tables'
gandal_mdd = read_csv(paste0(table_dir, '/Microarray_MDD_metaanalysis_092017_scale.csv'))
gandal_scz = read_csv(paste0(table_dir, '/Microarray_SCZ_metaanalysis_092017_scale.csv'))
gandal_bd  = read_csv(paste0(table_dir, '/Microarray_BD_metaanalysis_092017_scale.csv'))
gandal_asd = read_csv(paste0(table_dir, '/Microarray_ASD_metaanalysis_092017_scale.csv'))
gandal_aad = read_csv(paste0(table_dir, '/Microarray_AAD_metaanalysis_092017_scale.csv'))


# gene-wise correlation to AHBA spatial and downreg

# MDD
ahba_mdd_gandal = merge(x=all_ahba_cors, y=gandal_mdd, by.x='genes', by.y='symbol')
mdd_cor         = cor.test(ahba_mdd_gandal$beta_scale, ahba_mdd_gandal$avg)


# SCZ
ahba_scz_gandal = merge(x=all_ahba_cors, y=gandal_scz, by.x='genes', by.y='symbol')
scz_cor         = cor.test(ahba_scz_gandal$beta_scale, ahba_scz_gandal$avg)


# BD
ahba_bd_gandal = merge(x=all_ahba_cors, y=gandal_bd, by.x='genes', by.y='symbol')
bd_cor         = cor.test(ahba_bd_gandal$beta_scale, ahba_bd_gandal$avg)


# ASD
ahba_asd_gandal = merge(x=all_ahba_cors, y=gandal_asd, by.x='genes', by.y='symbol')
asd_cor         = cor.test(ahba_asd_gandal$beta_scale, ahba_asd_gandal$avg)


# AAD
ahba_aad_gandal = merge(x=all_ahba_cors, y=gandal_aad, by.x='genes', by.y='symbol')
aad_cor         = cor.test(ahba_aad_gandal$beta_scale, ahba_aad_gandal$avg)


cor.test(abs(ahba_mdd_gandal$beta_scale), ahba_mdd_gandal$avg)
cor.test(abs(ahba_scz_gandal$beta_scale), ahba_scz_gandal$avg)
cor.test(abs(ahba_bd_gandal$beta_scale), ahba_bd_gandal$avg)
cor.test(abs(ahba_asd_gandal$beta_scale), ahba_asd_gandal$avg)
cor.test(abs(ahba_aad_gandal$beta_scale), ahba_aad_gandal$avg)



# organize for plotting
genewise_cor = rbind(data.frame(cor=mdd_cor$estimate, pheno='MDD', l95=mdd_cor$conf.int[1], u95=mdd_cor$conf.int[2], pval=mdd_cor$p.value),
                        data.frame(cor=bd_cor$estimate, pheno='BD', l95=bd_cor$conf.int[1], u95=bd_cor$conf.int[2], pval=bd_cor$p.value),
                        data.frame(cor=asd_cor$estimate, pheno='ASD', l95=asd_cor$conf.int[1], u95=asd_cor$conf.int[2], pval=asd_cor$p.value),
                        data.frame(cor=aad_cor$estimate, pheno='AAD', l95=aad_cor$conf.int[1], u95=aad_cor$conf.int[2], pval=aad_cor$p.value),
                        data.frame(cor=scz_cor$estimate, pheno='SCZ', l95=scz_cor$conf.int[1], u95=scz_cor$conf.int[2], pval=scz_cor$p.value))


# barplot
p = ggplot(genewise_cor, aes(x=pheno, y=cor, fill=pheno)) +
            geom_bar(stat="identity", color="black", position=position_dodge(), width=0.85) +
            geom_errorbar(aes(ymin=l95, ymax=u95), width=.4, position=position_dodge(.9)) +
            theme_classic() +
            scale_fill_manual(values=c('red','grey69','grey69','grey69','grey69')) +
            scale_y_continuous(limit=c(-0.07,0.07), breaks=seq(-0.07,0.07, 0.035), expand=c(0,0))

p

figure_out = paste0(base_dir, '/figures/invivo_gene_corr_to_gandal_downreg.pdf')
figure_out
CairoPDF(figure_out, height=2, width=3)
print(p)
dev.off()



# summarize data by bins
ahba_mdd_gandal = ahba_mdd_gandal[!is.na(ahba_mdd_gandal$beta_scale),]
ahba_mdd_gandal = ahba_mdd_gandal[order(ahba_mdd_gandal$beta_scale),]

ahba_mdd_gandal$group = as.numeric(cut_number(1:nrow(ahba_mdd_gandal), 40))
by_decile             = ahba_mdd_gandal %>% group_by(group) %>% summarise(ahba_cor=mean(avg), downreg=mean(beta_scale))
mdd_decile_cor        = cor.test(by_decile$ahba_cor, by_decile$downreg, method='spearman')

mdd_decile_cor = cor.ci(cbind(by_decile$downreg, by_decile$ahba_cor), method='spearman', plot=F)
mdd_decile_cor

plot(by_decile$downreg, by_decile$ahba_cor)

# ---------
summary(by_decile$ahba_cor)
summary(by_decile$downreg)
# plot MDD correlation
plot_point = ggplot(data=by_decile, aes(x=ahba_cor, y=downreg)) +
                      geom_point(shape=21, show.legend = FALSE, type=21) +
                      geom_smooth(size=.5, fill='gray69', color='black', show.legend = FALSE, linetype = 'dashed', fullrange = TRUE, method = lm, se = FALSE) +
                      theme_classic() +
                      scale_x_continuous(limits=c(-.04,.04), expand=c(0,0), breaks=seq(-.04,.04,.02)) +
                      scale_y_continuous(limits=c(-.5,.5), expand=c(0,0), breaks=seq(-.5,.5,.25)) +
                      scale_color_manual(values='red') +
                      theme(legend.position = "none", axis.text=element_text(color='black'), axis.ticks=element_line(color='black'))

figure_out = paste0(base_dir, '/figures/gene_mddcor_by_mdd_downgreg.pdf')
figure_out
CairoPDF(figure_out, height=2.5, width=2.5)
print(plot_point)
dev.off()



# SCZ
# ----
ahba_scz_gandal = ahba_scz_gandal[!is.na(ahba_scz_gandal$beta_scale),]
ahba_scz_gandal = ahba_scz_gandal[order(ahba_scz_gandal$beta_scale),]
ahba_scz_gandal$group = as.numeric(cut_number(1:nrow(ahba_scz_gandal), 40))
scz_by_decile   = ahba_scz_gandal %>% group_by(group) %>% summarise(ahba_cor=mean(avg), downreg=mean(beta_scale))
scz_decile_cor  = cor.test(scz_by_decile$ahba_cor, scz_by_decile$downreg, method='spearman')
scz_decile_cor  = cor.ci(cbind(scz_by_decile$downreg, scz_by_decile$ahba_cor), method='spearman', plot=F)


# BD
# ----
ahba_bd_gandal = ahba_bd_gandal[!is.na(ahba_bd_gandal$beta_scale),]
ahba_bd_gandal = ahba_bd_gandal[order(ahba_bd_gandal$beta_scale),]
ahba_bd_gandal$group = as.numeric(cut_number(1:nrow(ahba_bd_gandal), 40))
bd_by_decile   = ahba_bd_gandal %>% group_by(group) %>% summarise(ahba_cor=mean(avg), downreg=mean(beta_scale))
bd_decile_cor  = cor.test(bd_by_decile$ahba_cor, bd_by_decile$downreg, method='spearman')
bd_decile_cor  = cor.ci(cbind(bd_by_decile$downreg, bd_by_decile$ahba_cor), method='spearman', plot=F)


# ASD
# ----
ahba_asd_gandal = ahba_asd_gandal[!is.na(ahba_asd_gandal$beta_scale),]
ahba_asd_gandal = ahba_asd_gandal[order(ahba_asd_gandal$beta_scale),]
ahba_asd_gandal$group = as.numeric(cut_number(1:nrow(ahba_asd_gandal), 40))
asd_by_decile   = ahba_asd_gandal %>% group_by(group) %>% summarise(ahba_cor=mean(avg), downreg=mean(beta_scale))
asd_decile_cor  = cor.test(asd_by_decile$ahba_cor, asd_by_decile$downreg, method='spearman')
asd_decile_cor  = cor.ci(cbind(asd_by_decile$downreg, asd_by_decile$ahba_cor), method='spearman', plot=F)


# AAD
# ----
ahba_aad_gandal = ahba_aad_gandal[!is.na(ahba_aad_gandal$beta_scale),]
ahba_aad_gandal = ahba_aad_gandal[order(ahba_aad_gandal$beta_scale),]
ahba_aad_gandal$group = as.numeric(cut_number(1:nrow(ahba_aad_gandal), 40))
aad_by_decile   = ahba_aad_gandal %>% group_by(group) %>% summarise(ahba_cor=mean(avg), downreg=mean(beta_scale))
aad_decile_cor  = cor.test(aad_by_decile$ahba_cor, aad_by_decile$downreg, method='spearman')
aad_decile_cor  = cor.ci(cbind(aad_by_decile$downreg, aad_by_decile$ahba_cor), method='spearman', plot=F)



decile_genewise_cor = rbind(data.frame(cor=mdd_decile_cor$rho[1,2], pheno='MDD', l95=mdd_decile_cor$ci$lower, u95=mdd_decile_cor$ci$upper, pval= mdd_decile_cor$p),
                            data.frame(cor=bd_decile_cor$rho[1,2], pheno='BD', l95=bd_decile_cor$ci$lower, u95=bd_decile_cor$ci$upper, pval=bd_decile_cor$p),
                            data.frame(cor=asd_decile_cor$rho[1,2], pheno='ASD', l95=asd_decile_cor$ci$lower, u95=asd_decile_cor$ci$upper, pval=asd_decile_cor$p),
                            data.frame(cor=scz_decile_cor$rho[1,2], pheno='SCZ', l95=scz_decile_cor$ci$lower, u95=scz_decile_cor$ci$upper, pval=scz_decile_cor$p),
                            data.frame(cor=aad_decile_cor$rho[1,2], pheno='AAD', l95=aad_decile_cor$ci$lower, u95=aad_decile_cor$ci$upper, pval=aad_decile_cor$p))

decile_genewise_cor$pheno = factor(decile_genewise_cor$pheno, levels=rev(decile_genewise_cor$pheno))


p = ggplot(decile_genewise_cor, aes(x=pheno, y=cor, fill=pheno)) +
            geom_bar(stat="identity", color="black", position=position_dodge(), width=0.85) +
            geom_errorbar(aes(ymin=l95, ymax=u95), width=.3, position=position_dodge(.9)) +
            theme_classic() +
            coord_flip() +
            scale_fill_manual(values=rev(c('red','grey69','grey69','grey69','grey69'))) +
            scale_y_continuous(limit=c(-1,1), breaks=seq(-1,1,.25), expand=c(0,0))


figure_out = paste0(base_dir, '/figures/decile_cors_by_pheno.pdf')
figure_out
CairoPDF(figure_out, height=2.5, width=3)
print(p)
dev.off()


