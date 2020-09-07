library(tidyverse)
library(seuret)
library(XML)
library(dbparser)
library(gridExtra)
library(RColorBrewer)


source('/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/06_gene_interpretation/00_support_functions.R')

plot_bin_enrichments = function(enrich_list, cat, nbin, ylim){
    color_arr = rev(colorRampPalette(brewer.pal(9,"Blues"))(10))

    bin_df = NULL
    for (bin in 1:nbin){
        bin_len = length(which(bin_enrich_list[[bin]]$Category == cat))
        bin_df  = rbind(bin_df, data.frame(bin=bin, cat=cat, n_enrich=bin_len))
    }
    bin_df$bin = as.factor(bin_df$bin)
    p = ggplot(bin_df, aes(x=bin, y=n_enrich, fill=bin)) +
                geom_bar(stat='identity', colour="black", width=.75) +
                scale_fill_manual(values=color_arr) +
                theme_classic() +
                scale_y_continuous(expand=c(0,0), limits=ylim[1:2], breaks=seq(ylim[1], ylim[2], ylim[3])) +
                theme(legend.position = "none", axis.title.y=element_blank())

    return(p)
}


# set up directories
base_dir     = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
toppgene_dir = paste0(base_dir, '/output/topp_gene')

# read the enrichments
nbins = 10
bin_enrich_list = NULL
for (bin in 1:nbins){
    toppgene_file = paste0(toppgene_dir, '/toppgene_ahba_mdd_bin', as.character(bin), '_10.txt')
    bin_enrich    = read_delim(toppgene_file, delim='\t')
    bin_enrich_list[[bin]] = bin_enrich
}


# Plot the number of enrichment categories for each decile
cat_arr = c('GO: Molecular Function',
            'GO: Biological Process',
            'GO: Cellular Component',
            'Pathway')
enrich_list = bin_enrich_list
ylim_arr = list(c(0,60,20), c(0,400,100), c(0,200,50), c(0,200,50))
ct = 1
plot_arr = NULL
for (cat in cat_arr){
    print(cat)
    for (bin in 1:10){
        print(length(which(bin_enrich_list[[bin]]$Category == cat)))
    }
    p = plot_bin_enrichments(enrich_list=bin_enrich_list, cat=cat, nbin=10, ylim=ylim_arr[[ct]])
    plot_arr[[ct]] = p
    ct = ct + 1
    out_plot = paste0(base_dir, '/figures/toppgene_topdecile_',cat,'_enrichment.pdf')
    ggsave(out_plot, plot=p, height=2, width=1.5)
}


# plot the disease enrichment for the top MDD gene decile
top_bin         = bin_enrich_list[[1]]
top_bin_disease = top_bin[which(top_bin$Category == 'Disease'),]
colnames(top_bin_disease) = gsub('&|-| ','_',colnames(top_bin_disease))
top_bin_disease$Name      = factor(top_bin_disease$Name, levels=rev(top_bin_disease$Name))
top_bin_disease$log10_p   = -1*log10(top_bin_disease$p_value)

threshold = -1*log10(.05)

p = ggplot(data=top_bin_disease, aes(x=Name, y=log10_p)) +
            geom_bar(stat='identity') +
            coord_flip() +
            geom_hline(yintercept=threshold) +
            theme_classic() +
            scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=seq(0,8,2))
out_plot = paste0(base_dir, '/figures/toppgene_disease_enrich_barplot.pdf')
out_plot
ggsave(out_plot, plot=p, height=2.5, width=3)



check_top_genes = function(top_genes, go_enrich){

  for (hit_list in go_enrich$Hit_in_Query_List){
      term_genes = strsplit(as.character(hit_list),',')[[1]]
      int_genes  = intersect(term_genes, top_genes)
      diff_genes = term_genes[!term_genes %in% top_genes]

      # select random
      if (length(int_genes) < 5){
          print('heyoo')
          print(length(int_genes))
          rand_gene = diff_genes[sample(length(diff_genes),4)]
          print(rand_gene)
          top_genes = c(top_genes, rand_gene)
      }
  }
  return(top_genes)
}



get_gene_hit_matrix = function(go_enrich, genes){

  n_row  = nrow(go_enrich)
  n_gene = length(genes)
  gene_enrich_mat = matrix(0, n_row, n_gene)
  rownames(gene_enrich_mat) = go_enrich$Name
  colnames(gene_enrich_mat) = genes
  col = 1
  for (gene in genes){
    row = 1
    for (hit_list in go_enrich$Hit_in_Query_List){
      gene_query = gene %in% strsplit(hit_list,',')[[1]]
      if (gene_query == T){
        gene_enrich_mat[row,col] = 1
      }
    row = row + 1
  }
  col = col + 1
  }
  return(gene_enrich_mat)
}


# Plot a select subset of the hundreds of enrichment terms from the top MDD gene decile
toppgene_plot = read.csv(paste0(toppgene_dir, '/mdd_toppgene_10_1_forPLOT.csv'))
colnames(toppgene_plot) = gsub('[.]','_',colnames(toppgene_plot))


go_all_genes = unlist(lapply(toppgene_plot$Hit_in_Query_List, function(x) strsplit(as.character(x),',')[[1]]))
go_gene_cts  = sort(table(go_all_genes))
top25        = rev(sort(go_gene_cts))[1:40]

# add gabra genes
top25_genes  = c(names(top25), c('GABRA5', 'GABRA3'))
top_genes    = check_top_genes(top_genes=top25_genes, go_enrich=toppgene_plot)
gene_hit_mat = get_gene_hit_matrix(toppgene_plot, genes=top_genes)


dist_mat    =  make_go_dist(toppgene_plot)
go_dist_mat = dist(dist_mat)
go_hc       = hclust(1-go_dist_mat, "ward.D2")
go_hcdata   = dendro_data_k(go_hc, 5)

p1 = plot_ggdendro(go_hcdata,
                    direction   = "lr",
                    expand.y    = 0.15) +
                    theme(axis.text.x      = element_text(color = "#ffffff"),
                          panel.background = element_rect(fill  = "#ffffff"),
                          axis.ticks       = element_blank()) +
                    scale_color_brewer(palette = "Paired") +
                    ylab('dendro') +
                    theme(axis.text.y = element_text(angle = 90)) +
                    xlab(NULL)
p1



gene_dist_mat = dist(t(gene_hit_mat))
gene_hc       = hclust(1-gene_dist_mat, "ward.D2")
gene_hcdata   = dendro_data_k(gene_hc, 6)
gene_order    = as.character(gene_hcdata$labels$label[order(gene_hcdata$labels$clust)])

gene_melt = reshape2::melt(gene_hit_mat)
gene_melt = data.table(gene_melt)
gene_melt[, Var1 := as.factor(Var1)]
idx       = match(rownames(gene_hit_mat), go_hcdata$labels$label)
gene_melt[, go_idx := rep(idx, ncol(gene_hit_mat))]

gene_melt$Var2 = factor(gene_melt$Var2, levels=gene_order)

# 'bubblemap'
gene_melt$value = as.character(gene_melt$value)
p2 = ggplot(gene_melt, aes(x=Var2, y=go_idx, fill=value)) +
        geom_tile(color='gray') +
        scale_color_viridis_c(direction = -1) +
        theme_minimal() +
        theme(axis.text.y = element_blank()) +
        xlab(NULL) +
        ylab(NULL) +
        theme(axis.text.x = element_text(angle = 90, size=4)) +
        scale_fill_manual(values=c('#f2f1ef','#5492d8'))

p2
library(gridExtra)
library(grid)
out_path =  paste0(base_dir, '/figures/toppgene_enrich_table.pdf')
g = grid.arrange(p1, p2, ncol = 2, widths = c(1,2))
ggsave(out_path, plot=g, height=4, width=10)







