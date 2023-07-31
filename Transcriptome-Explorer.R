# Load necessary libraries
# Install libraries if not installed

packages <- c("DESeq2", "pheatmap", "tidyverse", "cowplot", "Matrix.utils", 
              "edgeR", "dplyr", "magrittr", "Matrix", "purrr", "reshape2", 
              "S4Vectors", "tibble", "apeglm", "png", "RColorBrewer", "optparse", 
              "limma", "edgeR", "openxlsx", "ggplot2", "ggrepel", "clusterProfiler", 
              "enrichplot", "gprofiler2", "ggplot2", 
              "enrichR", "tidyverse", "STRINGdb", "igraph", "ggraph", "RColorBrewer","readxl")
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))

# get options
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="counts matrix file", metavar="character"),
  make_option(c("-f", "--samplefile"), type="character", default=NULL, 
              help="sample file", metavar="character"),
  make_option(c("-b", "--background"), type="character", default=NULL, 
              help="background gene list file (e.g. all/no/a list of genes)", metavar="character"),
  make_option(c("-t", "--tool"), type="character", default="DESeq2", 
              help="tool selection (DESeq2/Limma), default=DESeq2 ", metavar="character"),
  make_option(c("-d", "--datatype"), type="character", default="raw", 
              help="data type (raw/norm), default=raw ", metavar="character"),
  make_option(c("-g", "--groupnames"), type="character", default=NULL, 
              help="group names (e.g., condition,wt,ko)", metavar="character"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05, 
              help="p-value cutoff, default=0.05", metavar="numeric"),
  make_option(c("-q", "--fdr"), type="numeric", default=0.05, 
              help="FDR cutoff, default=0.05", metavar="numeric"),
  make_option(c("-F", "--foldchange"), type="numeric", default=0.58, 
              help="fold change cutoff, default=0.58", metavar="numeric"),
  make_option(c("-T", "--top"), type="integer", default=20, 
              help="number of top DEGs, default=20", metavar="integer"),
  make_option(c("-D", "--formula"), type="character", default=NULL, 
              help="design formula, e.g. ~sex + age + treatment", metavar="character"),
  make_option(c("-P", "--plottype"), type="character", default="png", 
              help="plot type (png/pdf/svg)", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default="mouse",
  help="Reference genome, can be 'mouse' or 'human'"),
  make_option(c("-m", "--missing_option")) 
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

check_file_exists <- function(filename) {
  print(paste("Checking file:", filename))  # Add this line
  if (!file.exists(filename)) {
    stop(paste("File does not exist:", filename))
  }
}

check_file_exists(opt$counts)
check_file_exists(opt$samplefile)

# Check if necessary options are provided
if (is.null(opt$groupnames)) {
  stop("You must provide the --groupnames option.")
}

if (is.null(opt$formula)) {
  stop("You must provide the --formula option.")
}

check_tool = function(tool) {
  if (!(tool %in% c("DESeq2", "Limma"))) {
    stop(paste("Invalid tool:", tool))
  }
}

if (opt$reference == "mouse") {
  library(org.Mm.eg.db)
  orgdb <- org.Mm.eg.db
} else if (opt$reference == "human") {
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
} else {
  stop("Invalid reference genome. Must be 'mouse' or 'human'.")
}


check_tool(opt$tool)

counts = read_xlsx(opt$counts)
counts = as.data.frame(counts)
counts_select = counts[,-c(1,2)]

if (opt$tool == "DESeq2" && opt$datatype != "raw") {
  counts_select <- round(counts_select)
}

rownames(counts_select) <- counts$ID

sample = read_xlsx(opt$samplefile)
sample = as.data.frame(sample)
samples_data = sample



if (nchar(opt$background) > 0 && !(opt$background %in% c('all', 'no'))) {
  background <- readxl::read_xlsx(opt$background)
  background <- as.data.frame(background)
}
  
genename <- cbind(counts$ID,counts$geneid)
colnames(genename) <- c("ID","geneid")
background_all = as.data.frame(genename)

formula <- as.formula(opt$formula)
groupnames <- strsplit(opt$groupnames, ",")[[1]]
coldata = data.frame(row.names=colnames(counts_select), samples_data[,-1])

filetype <- opt$plottype
#-------------------------------Plots--------------------------------

save_ggplot <- function(plot, filetype, filename) {
  if (filetype == "png") {
    ggsave(paste0(filename, ".png"), plot)
  } else if (filetype == "pdf") {
    ggsave(paste0(filename, ".pdf"), plot)
  } else if (filetype == "svg") {
    ggsave(paste0(filename, ".svg"), plot)
  } else {
    print(paste("Unsupported file type:", filetype))
  }
}

save_heatmap <- function(rld_cor, filetype, filename) {
  # Define the color scheme
  color_scheme <- colorRampPalette(c("blue", "white", "red"))
  
  if (filetype == "png") {
    png(paste0(filename, ".png"))
    pheatmap::pheatmap(rld_cor, color = color_scheme(100))
    dev.off()
  } else if (filetype == "pdf") {
    pdf(paste0(filename, ".pdf"))
    pheatmap::pheatmap(rld_cor, color = color_scheme(100))
    dev.off()
  } else if (filetype == "svg") {
    svg(paste0(filename, ".svg"))
    pheatmap::pheatmap(rld_cor, color = color_scheme(100))
    dev.off()
  } else {
    print(paste("Unsupported file type:", filetype))
  }
}

# Run pheatmap with the annotation
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
# Run pheatmap with the annotation
save_heatmap1 <- function(sig_norm, heat_colors, filetype, filename, annotation, annotation_colors, cluster_cols = TRUE) {
  if (nrow(sig_norm) > 30) {
    show_rownames = FALSE
  } else {
    show_rownames = TRUE
  }
  
   if (filetype == "png") {
    png(paste0(filename, ".png"))
    pheatmap(sig_norm, 
             color = heat_colors, 
             cluster_rows = T, 
             show_rownames = show_rownames, 
             cluster_cols = cluster_cols,
             scale = "row",
             border_color = NA, 
             fontsize = 10, 
             fontsize_row = 10, 
             height = 20,
             annotation_col = annotation,
             annotation_colors = annotation_colors)
    dev.off()
  } else if (filetype == "pdf") {
    pdf(paste0(filename, ".pdf"))
    pheatmap(sig_norm, 
             color = heat_colors, 
             cluster_rows = T, 
             show_rownames = show_rownames,   
             cluster_cols = cluster_cols,
             scale = "row",
             border_color = NA, 
             fontsize = 10, 
             fontsize_row = 10, 
             height = 20,
             annotation_col = annotation,
             annotation_colors = annotation_colors)
    dev.off()
  } else if (filetype == "svg") {
    svg(paste0(filename, ".svg"))
    pheatmap(sig_norm, 
             color = heat_colors, 
             cluster_rows = T, 
             show_rownames = show_rownames,   
             cluster_cols = cluster_cols,
             scale = "row",
             border_color = NA, 
             fontsize = 10, 
             fontsize_row = 10, 
             height = 20,
             annotation_col = annotation,
             annotation_colors = annotation_colors)
    dev.off()
  } else {
    print(paste("Unsupported file type:", filetype))
  }
}





#-------------------------------DEseq2-------------------------------
  
  run_deseq2 <- function(counts_select, coldata, design_formula) {
    
    dds = DESeqDataSetFromMatrix(countData = counts_select, colData = coldata, design = design_formula)
    dds <- DESeq(dds) 
    dds
    rld <- rlog(dds, blind=TRUE)
    
    # PCA
    plot <- DESeq2::plotPCA(rld, intgroup=c("condition"))
    save_ggplot(plot, filetype, "PCA")
    
    
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    # sample correlation heatmap
    save_heatmap(rld_cor, filetype, "heatmap_sample")
    
    res = results(dds, contrast=c(groupnames))
    res = res[order(res$pvalue),]
    head(res)
    summary(res)
    res_tbl <- res %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>%
      as_tibble()
    res_tbl
    res_tbl <- merge(res_tbl,genename,by.x = "gene", by.y = "ID" )
    res_tbl$regulation <- ifelse(res_tbl$log2FoldChange > opt$foldchange , "up", 
                                 ifelse(res_tbl$log2FoldChange < -opt$foldchange, "down", "not"))
    write.csv(res_tbl, file="all_results.csv")
    
    
    # Set thresholds
    padj_cutoff <- opt$fdr
    
    # Subset the significant results
    sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
    
    # Check significant genes output
    sig_res
    
    # Write significant results to file
    write.csv(sig_res, file = "sig_genes.csv",
              quote = FALSE, 
              row.names = FALSE)
    
    # Subset the data frame to only include rows where regulation is "up"
    upregulated_genes = subset(sig_res, regulation == "up")
    
    # Write upregulated genes to file
    write.csv(upregulated_genes, file = "upregulated_genes.csv",
              quote = FALSE, 
              row.names = FALSE)
    
    # Subset the data frame to only include rows where regulation is "up"
    downregulated_genes = subset(sig_res, regulation == "down")
    
    # Write upregulated genes to file
    write.csv(downregulated_genes, file = "downregulated_genes.csv",
              quote = FALSE, 
              row.names = FALSE)
    
    ## dots of top genes
    normalized_counts <- counts(dds, 
                                normalized = TRUE)
    
    ## Order results by padj values
    top_sig_genes <- sig_res %>%
      dplyr::arrange(padj) %>%
      dplyr::pull(gene) %>%
      head(n=opt$top)
    
    top_sig_norm <- data.frame(normalized_counts) %>%
      rownames_to_column(var = "gene") %>%
      dplyr::filter(gene %in% top_sig_genes)
    
    gathered_top_sig <- top_sig_norm %>%
      gather(colnames(top_sig_norm)[2:length(colnames(top_sig_norm))], key = "samplename", value = "normalized_counts")
    gathered_top_sig <- merge(gathered_top_sig,genename,by.x = "gene", by.y = "ID" )
    
    
    ## plot using ggplot2
    top = ggplot(gathered_top_sig) +
      geom_point(aes(x = geneid, 
                     y = normalized_counts, 
                     color = samplename), 
                 position=position_jitter(w=0.1,h=0)) +
      scale_y_log10() +
      xlab("Genes") +
      ylab("log10 Normalized Counts") +
      ggtitle("Top Significant DE Genes") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    save_ggplot(top, filetype, "topDEGs")
    
    ##topDEGs_heatmap
    # Prepare the annotation data frame
    samples_data$condition <- as.factor(samples_data$condition)
    annotation <- data.frame(
      Group = factor(samples_data$condition)
    )
    rownames(annotation) <- samples_data$samplename
    groups <- unique(annotation$Group) 
    colors <- c("lightblue", "pink") 
    names(colors) <- groups 
    annotation_colors <- list(Group = colors)

    # First, you need to sort your genes based on log fold change
    sig_res <- sig_res %>% arrange(desc(log2FoldChange))
    
    # Then you can select the top 10 up-regulated and down-regulated genes
    top_up <- head(sig_res, 10)
    top_down <- tail(sig_res, 10)
    
    # Now you can extract the expression data for these genes from your normalized counts
    selected_genes <- rbind(top_up, top_down)
    selected_data <- normalized_counts[rownames(normalized_counts) %in% selected_genes$gene,]
    selected_data <- as.data.frame(selected_data)
    selected_data$gene <- rownames(selected_data)
    selected_data <- merge(selected_data,genename,by.x = "gene", by.y = "ID" )
    rownames(selected_data) <- selected_data$geneid
    selected_data <- subset(selected_data, select = -c(gene,geneid))

    # Now you can create the heatmap using the save_heatmap1 function
    save_heatmap1(selected_data, 
                  heat_colors, 
                  filetype, 
                  "heatmap_top10_up_down_genes", 
                  annotation,
                  annotation_colors, 
                  cluster_cols = T)
    
    ##DEGs_heatmap
    sig_norm <- data.frame(normalized_counts) %>%
      rownames_to_column(var = "gene") %>%
      dplyr::filter(gene %in% sig_res$gene)
    rownames(sig_norm) <- sig_norm$gene
    sig_norm <- dplyr::select(sig_norm, -gene)
    
    save_heatmap1(sig_norm, heat_colors, filetype, "heatmap_DEGs", annotation, annotation_colors, cluster_cols = T)
    
    ## volcano plot
    res_table_thres <- res_tbl %>% 
      mutate(significant = padj < opt$fdr & abs(log2FoldChange) >= opt$foldchange, 
             color = case_when(
               log2FoldChange < 0 & significant ~ 'down', 
               log2FoldChange > 0 & significant ~ 'up', 
               TRUE ~ 'not'
             ))
    
    # Get the top 5 up and down genes
    top_up <- res_table_thres %>% filter(color == 'up') %>% arrange(desc(log2FoldChange)) %>% head(5)
    top_down <- res_table_thres %>% filter(color == 'down') %>% arrange(log2FoldChange) %>% head(5)
    
    # Count the number of up, down, and not changed genes
    n_up <- sum(res_table_thres$color == 'up')
    n_down <- sum(res_table_thres$color == 'down')
    n_not <- sum(res_table_thres$color == 'not')
    
    vo = ggplot(res_table_thres) +
      geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = color, size = -log10(padj))) +
      geom_text(data = top_up, aes(x = log2FoldChange, y = -log10(padj), label = geneid), vjust = -1) +
      geom_text(data = top_down, aes(x = log2FoldChange, y = -log10(padj), label = geneid), vjust = 1) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 3, 
               label = paste("Up: ", n_up, " Down: ", n_down, " Not changed: ", n_not), size = 5) +
      scale_color_manual(values = c("down" = "blue", "up" = "red", "not" = "grey")) +
      scale_size(range = c(1, 5)) +
      labs(caption = paste("Differential Expression Analysis:", groupnames[2], "vs", groupnames[3])) +
      xlab("log2 fold change") + 
      ylab("-log10 adjusted p-value") +
      scale_y_continuous(limits = c(0,50)) +
      theme(legend.position = "none",
            plot.caption = element_text(hjust = 0.5, size = rel(1.5)),
            axis.title = element_text(size = rel(1.25)))  
    
    save_ggplot(vo, filetype, "volcano")
  }


#---------------------------------------Limma-----------------------------------
  
  run_limma <- function(counts_select, coldata, design_formula) {
    
    dge <- DGEList(counts = counts_select)
    
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    
    variables <- all.vars(design_formula)
    coldata[variables] <- lapply(coldata[variables], as.factor)
    design_formula <- "~ 0"
    for (var in variables) {
      design_formula <- paste(design_formula, "+", var)
    }
    design <- model.matrix(as.formula(design_formula), data = coldata)
    
    
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    summary(decideTests(fit))
    output <- topTable(fit, coef=ncol(design), n=Inf)
    sum(output$adj.P.Val<0.05)
    res_tbl = na.omit(output) 
    write.csv(res_tbl, "all_results.csv") 
    
    
    
    foldChange = opt$foldchange
    padj = opt$fdr
    sig_res <- res_tbl[(res_tbl$adj.P.Val < padj & (res_tbl$logFC>foldChange | res_tbl$logFC < (-foldChange))),]
    dim(sig_res)
    sig_res <- merge(sig_res,genename,by.x = "gene", by.y = "ID" )
    write.csv(sig_res, "sig_genes.csv")
    
    upregulated_genes <-  sig_res[(sig_res$P.Value < padj & (sig_res$logFC > foldChange)),]
    write.csv(diffup, "upregulated_genes.csv")
    #
    downregulated_genes <- sig_res[(sig_res$P.Value < padj & (sig_res < -foldChange)),]
    write.csv(diffdown, "downregulated_genes.csv")
    
    
    logFC <- res_tbl$logFC
    deg.padj <- res_tbl$P.Value
    data <- data.frame(logFC = logFC, padj = deg.padj)
    data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "not"
    data$group[(data$padj <= 0.05 & data$logFC > 1)] <-  "up"
    data$group[(data$padj <= 0.05 & data$logFC < -1)] <- "down"
    x_lim <- max(logFC,-logFC)
    
    # volcano Plot
    pdf('volcano.pdf',width = 7,height = 6.5)  
    label = subset(res_tbl,P.Value <0.05 & abs(logFC) > 0.5)
    label1 = rownames(label)
    
    colnames(res_tbl)[1] = 'log2FC'
    Significant=ifelse((res_tbl$P.Value < 0.05 & abs(res_tbl$log2FC)> 0.5), ifelse(res_tbl$log2FC > 0.5,"up","down"), "not")
    
    ggplot(res_tbl, aes(log2FC, -log10(P.Value)))+
      geom_point(aes(col=Significant))+
      scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
      labs(title = " ")+
      geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
      geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
      theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
      labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
      theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
      str(res_tbl, max.level = c(-1, 1))+theme_bw()
    
    dev.off()
    
    # PCA
    tExprsData <- t(counts)
    pcaResult <- prcomp(tExprsData)
    summary(pcaResult)
    plot(pcaResult, type="l")
    library(ggplot2)
    pcaData <- pcaResult$x
    pca_df <- as.data.frame(pcaData)
    pca_df$phenotype <- rownames(pca_df)
    ggplot(pca_df, aes(PC1, PC2, color = phenotype)) + geom_point() + labs(x = "PC1", y = "PC2")
    
    # Volcano plot
    volcanoplot(fit,coef=2,highlight=2)
    
    # Mean-difference plot
    plotMD(fit,column=2)
    
    # Q-Q plot of moderated t-statistics
    qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
    abline(0,1)
    
  }

if (opt$tool == "DESeq2") {
  run_deseq2(counts_select, coldata, formula)
} else if (opt$tool == "Limma") {
  run_limma(counts_select, coldata, formula)
}



#---------------------------------------DEGs_analysis------------------------------------

sig_res <- read.csv("sig_genes.csv", row.names=1)
upregulated_genes <- read.csv("upregulated_genes.csv", row.names=1)
downregulated_genes <- read.csv("downregulated_genes.csv", row.names=1)

alldegs <- sig_res$geneid
updegs <- upregulated_genes$geneid
downdegs <- downregulated_genes$geneid

kegg_go <- cbind(sig_res$geneid,upregulated_genes$geneid,downregulated_genes$geneid)
colnames(kegg_go) <- c("all","up","down") 
kegg_go <- na.omit(kegg_go)
kegg_go <- as.data.frame(kegg_go)

#background genes 
if (opt$background == "all") {
  # Use all genes as background
  background_genes <- as.character(background_all$geneid)
  background_genes_entrez <- bitr(as.character(background_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
} else if (opt$background == "no") {
  background_genes_entrez = NULL
} else {
  # Assume opt$background is a file name containing a list of genes
  background_genes <- read_xlsx(opt$background)
  background_genes <- as.character(background_genes[[1]])
  background_genes_entrez <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
}

# Perform enrichment analysis
perform_analysis <- function(gene_list, list_name, plottype) {
  # Convert gene symbols to ENTREZ IDs
  gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
  # Check if any genes were mapped
  if (is.null(gene_list_entrez$ENTREZID)) {
    message("No genes could be mapped in list ", list_name, ". Skipping enrichment analysis.")
    return(NULL)
  }
  
  # Perform GO enrichment analysis
  ontologies <- c("BP", "CC", "MF")
  
  # Initialize an empty list to store the results
  results <- list()
  plots <- list()
  
  # Perform the analysis for each type
  for (ont in ontologies) {
    results[[ont]] <- enrichGO(gene           = gene_list_entrez$ENTREZID,
                               OrgDb          = orgdb,
                               ont            = ont,
                               pAdjustMethod  = "BH",
                               pvalueCutoff   = 0.05,
                               qvalueCutoff   = 0.2,
                               universe       = background_genes_entrez$ENTREZID,
                               readable       = TRUE)
  }
  
  bp_result <- results$BP
  cc_result <- results$CC
  mf_result <- results$MF
  
  # Check if any GO enrichment was found
  if (is.null(results$BP) || nrow(results$BP) == 0) {
    message("No GO enrichment found in list ", list_name, ". Skipping plot.")
    return(NULL)
  }
  
  # To visualize the results, you can use dotplot
  p <- dotplot(bp_result, showCategory=20) + ggtitle("Biological Process")
  ggsave(file=paste0("go_result_bp_", list_name,".",filetype ), plot=p)
  
  
  
  # Add a new column to each result to indicate the ontology
  bp_result <-  as.data.frame(bp_result) 
  bp_result$Ontology <- "BP"
  cc_result <- as.data.frame(cc_result) 
  cc_result$Ontology <- "CC"
  mf_result <- as.data.frame(mf_result)
  mf_result$Ontology <- "MF"

  # Combine the results into a single data frame
  combined_results <- rbind(bp_result, cc_result, mf_result)
  
  # Write the combined results to a CSV file
  write.csv(combined_results, file=paste0("go_result_", list_name, ".csv"), row.names = FALSE)
  
  # Assuming df is your data frame
  combined_results$minus_log_pvalue <- -log10(combined_results$pvalue)
  
 
  # Select the top 10 terms by p-value for each ontology
  combined_results_top10 <- combined_results %>%
    group_by(Ontology) %>%
    top_n(-10, wt = minus_log_pvalue)
  
  # Trim the Term names to a maximum of 15 characters
  combined_results_top10$Description <- strtrim(combined_results_top10$Description, 70)
  
  # Create the bar plot
  q <- ggplot(combined_results_top10, aes(x = reorder(Description, minus_log_pvalue), y = minus_log_pvalue, fill = Ontology)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~Ontology, scales = "free") +
    xlab("Term") +
    ylab("Enrichment score (-log10(pvalue))") +
    ggtitle("Gene Ontology Enrichment")
  
  # Save the plot
  ggsave(file=paste0("go_result_top10_", list_name,".",filetype ), plot=q)
  
  # Perform KEGG enrichment analysis
  kegg_result <- enrichr(genes = gene_list , databases = "KEGG_2019_Mouse")[[1]]
  
  # Check if any KEGG enrichment was found
  if (is.null(kegg_result) || nrow(kegg_result) == 1) {
    message("No KEGG enrichment found in list ", list_name, ". Skipping plot.")
    return(NULL)
  }
  
  # Select the top 30 terms based on Combined Score
  kegg_result_df <- as.data.frame(kegg_result)
  kegg_result_df <- head(kegg_result_df[order(-kegg_result_df$Combined.Score), ], 30)
  
  # Create the bar plot
  bar_plot <- ggplot(kegg_result_df, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +  # To make the bar plot horizontal
    theme_minimal() +
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA)) +
    xlab("KEGG Pathway") +
    ylab("Combined Score")
  
  # Save the bar plot
  ggsave(filename = paste0("kegg_result_bar_", list_name, ".", filetype), plot = bar_plot)
  
  # Create the bubble plot
  bubble_plot <- ggplot(kegg_result_df, aes(x = reorder(Term, -Combined.Score), y = Combined.Score, size = -log10(P.value))) +
    geom_point(alpha = 0.6) +
    scale_size_continuous(range = c(1,10)) +
    coord_flip() +  # To make the bubble plot horizontal
    theme_minimal() +
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA)) +
    xlab("KEGG Pathway") +
    ylab("Combined Score") +
    guides(size=guide_legend(title="Negative log of pvalue"))
  
  # Save the bubble plot
  ggsave(filename = paste0("kegg_result_bubble_", list_name, ".", filetype), plot = bubble_plot)
  
  write.csv(kegg_result, file=paste0("kegg_result_", list_name, ".csv"))
}

string_db <- STRINGdb$new(version = "11.5", species = 10090, score_threshold = 400, input_directory = "")

# Perform network analysis
perform_network_analysis <- function(gene_list, list_name, plottype) {
  
  if (length(gene_list) > 2000) {
    message("Gene list ", list_name, " is too long. Skipping network analysis.")
    return(NULL)
  }
  
  # Convert gene symbols to ENTREZ IDs
  gene_list_entrez <- gene_list %>% bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb, drop = T)
  
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
  ggsave(plot = p, file = paste0("network_", list_name,".",filetype), width = 15, height = 7)
  
  # Determine the degree of each node
  degree_values <- igraph::degree(net)
  
  # Order the nodes by degree in descending order
  nodes_ordered_by_degree <- order(degree_values, decreasing = TRUE)
  
  # Select the top 30 nodes
  top_nodes <- nodes_ordered_by_degree[1:30]
  
  # Create a subgraph with only the top nodes
  subgraph <- igraph::induced_subgraph(net, top_nodes)
  
  # Now plot the subgraph
  p <- ggraph(subgraph, layout = "linear", circular = TRUE)+
    geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
    geom_node_point(aes(size=size), color="orange", alpha=0.7)+
    geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
    scale_edge_width(range = c(0.2,1))+
    scale_size_continuous(range = c(1,10))+
    guides(size=F)+
    theme_graph()
  
  # Save the plot
  ggsave(plot = p, file = paste0("network_new_", list_name,".",filetype), width = 15, height = 7)
  
  
  edgelist <- igraph::get.edgelist(net)
  write.csv(edgelist, file=paste0("network_", list_name, ".csv"))
}


# A helper function to apply perform_analysis to a single column
perform_analysis_on_column <- function(i) {
  
  # Perform enrichment analysis
  perform_analysis(kegg_go[[i]], colnames(kegg_go)[i], filetype)
  
  # Perform protein network analysis
  perform_network_analysis(kegg_go[[i]], colnames(kegg_go)[i], filetype )
  
}

# Apply the function to each column of the data frame
lapply(seq_len(ncol(kegg_go)), perform_analysis_on_column)


