library(tidyverse)
library(fgsea)
library(glue)
library(ComplexHeatmap)
library(circlize)


read_pseudobulk_res <- function(path){
  read_csv(path) %>%
    dplyr::select(featurekey, contains('t_stat')) %>%
    dplyr::select(!contains('vs'))
}


get_gsea_res <- function(pb_res, custom_pathway_genes, save_name=NA){
  # get clusters in resolution
  clusters <- str_split(colnames(pb_res[-1]), ':') %>% 
    map(1) %>% 
    unlist() %>% 
    as.numeric() %>% 
    sort() %>%
    as.character()
  
  # create dataframe to store results
  res_table <- tibble(cluster = rep(NA, length(clusters)), 
                      res = rep(NA, length(clusters)))
  
  # run gsea for each cluster
  for (i in seq_along(clusters)){
    # print current cluster
    print(glue('cluster: {clusters[i]}'))
    
    # sort genes for gsea
    ranks <- pb_res %>%
      rename('t_stat' = glue('{clusters[i]}:pseudobulk_t_stat'),
             'gene_name' = 'featurekey') %>%
      dplyr::select('gene_name','t_stat') %>%
      na.omit() %>%
      distinct() %>%
      arrange(desc(t_stat)) %>%
      deframe()
    
    fgseaRes <- fgsea(pathways=custom_pathway_genes, 
                      stats=ranks,
                      minSize=10,
                      nPermSimple=10000)
    
    res_table$cluster[i] <- clusters[i]
    res_table$res[i] <- list(fgseaRes)
  }  
  
  res_table <- res_table %>% unnest(cols = c(res))
  
  if (!is.na(save_name)){
    # save gsea res
    write_csv(res_table %>%
                  rowwise() %>% 
                  mutate_if(is.list, ~paste(unlist(.), collapse = '|')),
                glue('{save_name}.csv'))
  }
  
  return(res_table)
}

plot_gsea_heatmap <- function(gsea_df){
  es_mat <- gsea_df %>%
    dplyr::select(pathway, cluster, ES) %>%
    pivot_wider(names_from = cluster, values_from = ES)
  es_mat_rownames <- es_mat$pathway
  es_mat <- es_mat %>% 
    dplyr::select(!pathway) %>% 
    as.matrix()
  rownames(es_mat) <- es_mat_rownames
  
  # get pvalue matrix
  pval_mat <- gsea_df %>%
    dplyr::select(pathway, cluster, pval) %>%
    pivot_wider(names_from = cluster, values_from = pval)
  pval_mat_rownames <- pval_mat$pathway
  pval_mat <- pval_mat %>% 
    dplyr::select(!pathway) %>% 
    as.matrix()
  rownames(pval_mat) <- pval_mat_rownames
  
  # check matrices the same order
  stopifnot(colnames(pval_mat) == colnames(es_mat))
  stopifnot(rownames(pval_mat) == rownames(es_mat))
  
  # define cell color range
  col_fun <- colorRamp2(c(floor(min(es_mat)), 0, ceiling(max(es_mat))),
                        c("blue", "white", "red"))
  
  ht <- Heatmap(es_mat, 
                name = "GSEA Enrichment Score",
                col = col_fun,
                cluster_columns = F,
                cluster_rows = F,
                row_title = NULL,
                column_names_side = "top",
                column_names_rot = 45,
                column_title = 'Top 25 Tstat',
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if((pval_mat[i, j] <= 0.1) & !is.na(pval_mat[i,j]))
                    grid.text(sprintf("%.2f", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                },
                width = ncol(es_mat)*unit(15, "mm"), 
                height = nrow(es_mat)*unit(15, "mm"))
  
  return(ht)
}

