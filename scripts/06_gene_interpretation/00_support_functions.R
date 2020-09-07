library(ggdendro)
library(ggplot2)
library(Cairo)
library(tidyverse)
library(data.table)
library(gridExtra)

dendro_data_k <- function(hc, k) {

  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust

  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {

  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)

  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)

  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }

  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle

  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)

  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }

  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)

  p
}



make_go_dist = function(df){
  n_go = nrow(df)
  dist_mat = matrix(NA, n_go, n_go)
  colnames(dist_mat) = df$Name
  rownames(dist_mat) = df$Name

  for (i in 1:n_go){
      print(paste0(as.character(i),'/',as.character(n_go)))
    for (j in i:n_go){
      i_genes = as.character(df$Hit_in_Query_List[i])
      i_gene_arr = strsplit(i_genes,',')[[1]]
      j_genes = as.character(df$Hit_in_Query_List[j])
      j_gene_arr = strsplit(j_genes,',')[[1]]

      n_intersect =  length(intersect(i_gene_arr, j_gene_arr))
      all_hits = length(union(i_gene_arr, j_gene_arr))

      sim_val = ( n_intersect / all_hits )
      dist_mat[i,j]=sim_val
      dist_mat[j,i]=sim_val
    }
  }
  return(dist_mat)
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

check_top_genes = function(top_genes, go_enrich){

  for (hit_list in go_enrich$Hit_in_Query_List){
    term_genes = strsplit(hit_list,',')[[1]]
    int_genes = intersect(term_genes, top_genes)
    # select random
    if (length(int_genes) == 0){
      print('heyoo')
      rand_gene = term_genes[sample(length(term_genes),1)]
      top_genes = c(top_genes, rand_gene)
    }
  }


  return(top_genes)
}
