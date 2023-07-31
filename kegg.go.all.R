library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(gprofiler2)
library(ggplot2)
library(readxl)
library(enrichR)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(ggraph)
library(RColorBrewer)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was provided
if (length(args) == 0) {
  stop("No Excel file provided. Please provide the path to an Excel file as a command line argument.")
}

# Read the Excel file
kegg_go <- read_excel(args[1])

background_genes <- kegg_go$background

background_genes_entrez <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform enrichment analysis
perform_analysis <- function(gene_list, list_name) {
  # Convert gene symbols to ENTREZ IDs
  gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Check if any genes were mapped
  if (is.null(gene_list_entrez$ENTREZID)) {
    message("No genes could be mapped in list ", list_name, ". Skipping enrichment analysis.")
    return(NULL)
  }
  
  # Perform GO enrichment analysis
  go_result <- enrichGO(gene           = gene_list_entrez$ENTREZID,
                        OrgDb          = org.Mm.eg.db,
                        ont            = "BP",  # change to "CC" or "MF" for cellular component or molecular function
                        pAdjustMethod  = "BH",
                        pvalueCutoff   = 0.05,
                        qvalueCutoff   = 0.2,
                        universe       = background_genes_entrez$ENTREZID,
                        readable       = TRUE)
  
  # Check if any GO enrichment was found
  if (is.null(go_result) || nrow(go_result) == 0) {
    message("No GO enrichment found in list ", list_name, ". Skipping plot.")
    return(NULL)
  }
  
  # To visualize the results, you can use dotplot
  p <- dotplot(go_result, showCategory=20)
  ggsave(file=paste0("go_result_", list_name, ".png"), plot=p)
  
  write.csv(go_result, file=paste0("go_result_", list_name, ".csv"))
  
  # Perform KEGG enrichment analysis
  kegg_result <- enrichr(genes = gene_list , databases = "KEGG_2019_Mouse")[[1]]
  
  # Check if any KEGG enrichment was found
  if (is.null(kegg_result) || nrow(kegg_result) == 1) {
    message("No KEGG enrichment found in list ", list_name, ". Skipping plot.")
    return(NULL)
  }

  # Check the number of elements in kegg_results$P.value
  n <- length(kegg_result$P.value)

  
  # If there are less than 20 elements, use all of them; otherwise, use the last 20 elements
  if (n < 20) {
    values_to_plot <- -log10(kegg_result$P.value)[n:1]
    names_to_plot <- kegg_result$Term[n:1]
  } else {
    values_to_plot <- -log10(kegg_result$P.value)[20:1]
    names_to_plot <- kegg_result$Term[20:1]
  }
  
  png(paste0("kegg_result_", list_name, ".png"))
  par(mfrow = c(1, 1), mar = c(3, 25, 2, 1))
  barplot(height = values_to_plot, names.arg = names_to_plot, horiz = T, las = 1, border = F, cex.names = 1)
  abline(v = c(-log10(0.05)), lty = 2)
  abline(v = 0, lty = 1)
  dev.off()
  
  write.csv(kegg_result, file=paste0("kegg_result_", list_name, ".csv"))
}

string_db <- STRINGdb$new(version = "11.5", species = 10090, score_threshold = 400, input_directory = "")

# Perform network analysis
perform_network_analysis <- function(gene_list, list_name) {
  
  if (length(gene_list) > 2000) {
    message("Gene list ", list_name, " is too long. Skipping network analysis.")
    return(NULL)
  }
  
  # Convert gene symbols to ENTREZ IDs
  gene_list_entrez <- gene_list %>% bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db",drop = T)
  
  # Map the ENTREZ IDs to STRING IDs
  data_mapped <- gene_list_entrez %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                   removeUnmappedRows = TRUE)
  string_db$plot_network(data_mapped$STRING_id )

  data_links <- data_mapped$STRING_id %>% string_db$get_interactions()
  links <- data_links %>%
    mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
    mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
    dplyr::select(from, to , last_col()) %>% 
    dplyr::rename(weight = combined_score)
  nodes <- links %>% { data.frame(bulk.DEGgenes = c(.$from, .$to)) } %>% distinct()
  net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
  igraph::V(net)$deg <- igraph::degree(net) 
  igraph::V(net)$size <- igraph::degree(net)/5 #
  igraph::E(net)$width <- igraph::E(net)$weight/10
  
  p <- ggraph(net,layout = "linear", circular = TRUE)+
    geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
    geom_node_point(aes(size=size), color="orange", alpha=0.7)+
    geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
    scale_edge_width(range = c(0.2,1))+
    scale_size_continuous(range = c(1,10) )+
    guides(size=F)+
    theme_graph()
  
  
  # Save the plot
  ggsave(plot = p, file = paste0("network_", list_name, ".png"), width = 15, height = 7)
}


# A helper function to apply perform_analysis to a single column
perform_analysis_on_column <- function(i) {
  
  # Perform enrichment analysis
  perform_analysis(kegg_go[[i]], names(kegg_go)[i])
  
  # Perform protein network analysis
  perform_network_analysis(kegg_go[[i]], names(kegg_go)[i])
  
}

# Apply the function to each column of the data frame (excluding the last column)
lapply(seq_along(kegg_go)[-ncol(kegg_go)], perform_analysis_on_column)

