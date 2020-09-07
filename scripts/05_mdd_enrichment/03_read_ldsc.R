library(tidyverse)


parse_cell_type_results = function(enrich_outs, nbins,  descriptor){
    nbins = 10
    # higher bin (e.g. bin_20_20) means more enriched in cell type
    out_df = NULL
    cor_out_df = NULL
    for (e in enrich_outs){
        write(e,'')
        cell_type = strsplit(strsplit(as.character(e),'_ldsc_')[[1]][[2]],'_enrich')[[1]][[1]]
        e_path    = paste0(cur_out, '/', e)
        e_df      = read.table(e_path, header=T)
        e_tmp_df  = e_df[grep('bin', e_df$Name),]
        e_tmp_df  = e_tmp_df[rev(order(e_tmp_df$Name)),]
        e_tmp_df$idx = 1:nbins
        e_tmp_df$cell = cell_type
        e_tmp_df2 = e_df[grep('nomatch', e_df$Name),]
        e_tmp_df2$idx = nbins+1

        cell_cor = cor.test(e_tmp_df$idx, e_tmp_df$Coefficient)
        cor_out_df = rbind(cor_out_df, data.frame(cell=cell_type, bin_cor = cell_cor$estimate))
        out_df   = rbind(out_df, e_tmp_df)
    }
    out_df$category = descriptor
    return(out_df)
}

#############
# set up dirs
#############
base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'


##################
# Lake 2018 - DFC
# Big cat
##################
cur_out     = paste0(base_dir, '/data/ldsc/lake_dfc_bigcat_ldsc_files')
enrich_outs = list.files(cur_out, pattern=paste0('cell_type_results.txt'))
enrich_outs = enrich_outs[grep('Wray', enrich_outs)]
enrich_outs = enrich_outs[!grepl('filter', enrich_outs)]

# filter to top most bin
dfc_bigcat_df = parse_cell_type_results(enrich_outs=enrich_outs, nbins=nbins,  descriptor='dfc_bigcat')
dfc_bigcat_df = dfc_bigcat_df %>% filter(idx==1)
dfc_bigcat_df$q = p.adjust(dfc_bigcat_df$Coefficient_P_value, method='BH')
dfc_bigcat_df$log10q = -1*log10(dfc_bigcat_df$q)



##################
# Lake 2018 - VIS
# Big cat
##################
cur_out     = paste0(base_dir, '/data/ldsc/lake_vis_bigcat_ldsc_files')
enrich_outs = list.files(cur_out, pattern=paste0('cell_type_results.txt'))
enrich_outs = enrich_outs[grep('Wray', enrich_outs)]
enrich_outs = enrich_outs[!grepl('filter', enrich_outs)]

# filter to top bin
vis_bigcat_df = parse_cell_type_results(enrich_outs=enrich_outs, nbins=nbins,  descriptor='vis_bigcat')
vis_bigcat_df = vis_bigcat_df %>% filter(idx==1)
vis_bigcat_df$q = p.adjust(vis_bigcat_df$Coefficient_P_value, method='BH')
vis_bigcat_df$log10q = -1*log10(vis_bigcat_df$q)




##########
# ABA MTG
# Big cat
##########
# LDSC ABA MTG - Big cat
cur_out = paste0(base_dir, '/data/ldsc/aba_mtg_bigcat_ldsc_files')
enrich_outs = list.files(cur_out, pattern=paste0('cell_type_results.txt'))
enrich_outs = enrich_outs[grep('Wray', enrich_outs)]
enrich_outs = enrich_outs[!grepl('filter', enrich_outs)]

mtg_bigcat_df = parse_cell_type_results(enrich_outs=enrich_outs, nbins=nbins,  descriptor='mtg_bigcat')
mtg_bigcat_df = mtg_bigcat_df %>% filter(idx==1)
mtg_bigcat_df$q = p.adjust(mtg_bigcat_df$Coefficient_P_value, method='BH')
mtg_bigcat_df$log10q = -1*log10(mtg_bigcat_df$q)


###################
# LDSC - BIG
# combine and plot
###################
all_enrich_df      = rbind(dfc_bigcat_df, vis_bigcat_df, mtg_bigcat_df)
all_enrich_df$cell = gsub('Astro','Ast',all_enrich_df$cell)
all_enrich_df$cell = gsub('Oligo','Oli',all_enrich_df$cell)
all_enrich_df$cell = gsub('Micro','Mic',all_enrich_df$cell)
all_enrich_df$cell = gsub('Endo','End',all_enrich_df$cell)


write_me = all_enrich_df
write_me$category = gsub('mtg_bigcat', 'ABA_MTG', write_me$category)
write_me$category = gsub('vis_bigcat', 'Lake_VIS', write_me$category)
write_me$category = gsub('dfc_bigcat', 'Lake_DFC', write_me$category)
cur_out = paste0(base_dir, '/supp_data/SuppData_13_ldsc_enrichment.csv')
write_csv(x=write_me, cur_out)




# prepare for plotting
top_bin_df        = all_enrich_df %>% filter(idx==1)
top_bin_df$log10p = -1*log10(top_bin_df$Coefficient_P_value)
# MTG data doesnt have perinchymal cells, add filler row to make plot look uniform
fill_row          = data.frame(Name='Per_bigcat_bin_10_10',
                                Coefficient=0,
                                Coefficient_std_error=0,
                                Coefficient_P_value=0,
                                idx=1,q=1, log10q=0,
                                cell='Wray_MDD_Per_bigcat',
                                category='mtg_bigcat',
                                log10p=0)
top_bin_df = rbind(top_bin_df, fill_row)


# order
plot_df      = top_bin_df[c('log10q','category','cell','Name')]
plot_df$cell = gsub('_bigcat', '', gsub('Wray_MDD_', '', plot_df$cell))
cell_cats    = rev(c('Inh','Exc','OPC','Ast','Mic','Oli','Per','End'))


#avg_plot_df = plot_df %>% group_by(cell) %>% summarize(meanlog10 = mean(log10p))
#cell_order   = gsub('_bigcat', '', gsub('Wray_MDD_', '', cell_order))

plot_df$cell     = factor(plot_df$cell, levels=cell_cats)
plot_df$category = gsub('_bigcat', '', plot_df$category)
plot_df$category = factor(plot_df$category, levels=c('mtg','vis','dfc'))

outpath = paste0(base_dir, '/figures/ldsc_bigcat_enrichments.pdf')
colors  = c('#DEA255', '#98C3CF', '#65A8BC')
p = ggplot(plot_df, aes(y=log10q, x=cell, fill=category))+
            geom_bar(stat='identity', position='dodge') +
            scale_y_continuous(limits=c(0,1.5), expand=c(0,0)) +
            coord_flip() +
            geom_hline(yintercept=-1*log10(.05)) +
            theme_classic() +
            scale_fill_manual(values=colors)
p
ggsave(outpath, plot=p, height=3, width=3)
outpath



# Read magma results
load(paste0(base_dir, '/output/magma/aba_mtg_magma_results.Rdata'), verbose=T)
mtg_ctAssocsLinear = ctAssocsLinear
mtg_ctAssocsLinear[[2]]$results$category = 'aba_bigcat'

load(paste0(base_dir, '/output/magma/lake_dfc_magma_results.Rdata'), verbose=T)
dfc_ctAssocsLinear = ctAssocsLinear
dfc_ctAssocsLinear[[2]]$results$category = 'dfc_bigcat'

load(paste0(base_dir, '/output/magma/lake_vis_magma_results.Rdata'), verbose=T)
vis_ctAssocsLinear = ctAssocsLinear
vis_ctAssocsLinear[[2]]$results$category = 'vis_bigcat'


all_magma_bigcat = rbind(mtg_ctAssocsLinear[[2]]$results, dfc_ctAssocsLinear[[2]]$results, vis_ctAssocsLinear[[2]]$results)




# MAGMA enrichment
lake_vis = read.table(paste0(base_dir, '/data/magma/lake_vis_expr_property.gsa.out'), header=T)
lake_vis = lake_vis[lake_vis$VARIABLE != 'Average',]
lake_vis$cell_source = 'vis'
lake_vis$q = p.adjust(lake_vis$P, method='BH')

lake_dfc = read.table(paste0(base_dir, '/data/magma/lake_dfc_expr_property.gsa.out'), header=T)
lake_dfc = lake_dfc[lake_dfc$VARIABLE != 'Average',]
lake_dfc$cell_source = 'dfc'
lake_dfc$q = p.adjust(lake_dfc$P, method='BH')

magma_mtg = read.table(paste0(base_dir, '/data/magma/aba_mtg_expr_property.gsa.out'), header=T)
magma_mtg = magma_mtg[magma_mtg$VARIABLE != 'Average',]
magma_mtg$cell_source = 'mtg'
magma_mtg$q = p.adjust(magma_mtg$P, method='BH')

# MTG data doesnt have perinchymal cells, add filler row to make plot look uniform
fill_row          = data.frame(VARIABLE='Per',
                                TYPE='COVAR',
                                MODEL=0,
                                NGENES=0,
                                BETA_STD=0,
                                BETA=0,
                                SE=0,
                                P=0,
                                q=0,
                                cell_source = 'mtg')
magma_mtg = rbind(magma_mtg, fill_row)

magma_mtg$VARIABLE = gsub('Astro', 'Ast', magma_mtg$VARIABLE)
magma_mtg$VARIABLE = gsub('Endo', 'End', magma_mtg$VARIABLE)
magma_mtg$VARIABLE = gsub('Micro', 'Mic', magma_mtg$VARIABLE)
magma_mtg$VARIABLE = gsub('Oligo', 'Oli', magma_mtg$VARIABLE)


magma_plot_df = rbind(lake_dfc, lake_vis, magma_mtg)
cell_cats = rev(c('Inh','Exc','OPC','Ast','Mic','Oli','Per','End'))
magma_plot_df$VARIABLE = factor(magma_plot_df$VARIABLE, levels=cell_cats)
magma_plot_df$log10q = -1*log10(magma_plot_df$q)
magma_plot_df$log10q[is.infinite(magma_plot_df$log10q)] = 0
magma_plot_df$category = factor(magma_plot_df$cell_source, levels=c('mtg','vis','dfc'))

outpath = paste0(base_dir, '/figures/magma_gene_property_enrichments.pdf')
colors = c('#DEA255', '#98C3CF', '#65A8BC')
p = ggplot(magma_plot_df, aes(y=log10q, x=VARIABLE, fill=category))+
        geom_bar(stat='identity', position='dodge') +
        scale_y_continuous(limits=c(0,6), expand=c(0,0)) +
        coord_flip() +
        geom_hline(yintercept=-1*log10(.05)) +
        theme_classic() +
        scale_fill_manual(values=colors)
ggsave(outpath, plot=p, height=3, width=3)
outpath

write_csv(x=magma_plot_df, paste0(base_dir, '/supp_data/SuppData_14_magma_enrichment.csv'))
