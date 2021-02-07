library(ggplot2)
library(limma)
library(FactoMineR)
library(factoextra)
library(pheatmap)


get_RNA_types <- function(){
    # creates a annotation table
    # Examples:
    # --------------------------------------------------------
    #    ensembl_gene_id      hgnc_symbol       gene_biotype
    # 1   ENSG00000210049       MT-TF            Mt_tRNA
    # 2   ENSG00000211459     MT-RNR1            Mt_rRNA
    # 3   ENSG00000210077       MT-TV            Mt_tRNA

    library(biomaRt)
    
    mart <- useMart('ENSEMBL_MART_ENSEMBL')
    mart <- useDataset('hsapiens_gene_ensembl', mart)
    
    annot <- getBM(mart = mart, 
                   attributes = c("hgnc_symbol",
                                  "ensembl_gene_id",
                                  "ensembl_transcript_id_version",
                                  "refseq_mrna",
                                  "refseq_ncrna",
                                  "gene_biotype"),
                   uniqueRows = T)
    return(annot)
}


perform_pca <- function(expr, group_list, ofname = ""){
    # perform PCA on groups, 
    # to see if there is/are significant difference between them
    # -----------------------------------------------------------------
    # expr: expression matrix
    # group_list: (class: factor) groups that each sample belongs to
    expr_pca <- as.data.frame(t(expr))
    expr_pca <- PCA(expr_pca, graph = FALSE)

    pca_plot <- fviz_pca_ind(expr_pca,
                             geom.ind = "point",
                             col.ind = group_list,
                             palette = c("#00AFBB", "#E7B800"),
                             addEllipses = TRUE,
                             legend.title = "Groups")
    
    if (ofname != "") {
      ggsave(ofname, pca_plot)
    }
    
    return(pca_plot)
}


calc_limma_degs <- function(expr, group_list, logFC_cutoff = 1.5, p_value = 0.05, ofname = ""){
    # calculate differential expressed genes by package "limma";
    # add tags to DEGs based on cutoff values (i.e. logFC and p-value)
    # -----------------------------------------------------------------
    # expr: expression matrix
    # group_list: (class: factor) groups that each sample belongs to

    design <- model.matrix(~group_list)
  
    # calc DEGs
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    deg <- topTable(fit, coef = 2, number = Inf)
    
    # separate down-/up- regulated genes
    down <- (deg$P.Value < p_value) & (deg$logFC < -logFC_cutoff)
    up <- (deg$P.Value < p_value) & (deg$logFC > logFC_cutoff)
    deg$change = ifelse(down,"down",ifelse(up,"up","stable"))
    
    if (ofname != ""){
        write.csv(deg, file = ofname)
    }
    
    return(deg)
}


plot_volcano <- function(deg, labels = "", logFC_cutoff = 1.5, adj_pvalue = 0.05, ofname = ""){
    # generate volcano plot based on previously identified DEGs
    # -----------------------------------------------------------------
    # ofname: output png filename
    # deg: the DEG matrix returned from calc_limma_degs()
  
    # separate down-/up- regulated genes
    down <- (deg$adj.P.Val < adj_pvalue) & (deg$logFC < -logFC_cutoff)
    up <- (deg$adj.P.Val < adj_pvalue) & (deg$logFC > logFC_cutoff)
    deg$change = ifelse(down,"down",ifelse(up,"up","stable"))

    plot <- ggplot(data = deg, mapping = aes(x = logFC, y = -log10(adj.P.Val))) +
            geom_point(alpha = 0.4, size = 3.5, aes(color=change)) +
            geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),lty = 4, col = "black", lwd = 0.8) +
            geom_hline(yintercept = -log10(adj_pvalue), lty = 4, col = "black", lwd = 0.8) +
            ylab("-log10(adjusted P-value)") + 
            scale_color_manual(values = c("blue", "grey","red")) +
            theme_bw()
    
    if (labels != ""){
      labels_df = deg[rownames(deg) %in% labels, ]
        
      plot <- plot + 
              geom_point(size = 3, shape = 1, data = labels_df) +
              ggrepel::geom_label_repel(data = labels_df, aes(label = rownames(labels_df)), color = "black")
      
    }
    
    if (ofname != "") {
        ggsave(ofname, plot)
    }
    
    return(plot)
}


plot_heatmap <- function(expr, transcripts, group_list, ofname = ""){
    # plot heatmap against two groups

    deg_names <- rownames(deg)[deg$change !="stable"]
    diff <- expr[deg_names,]
    
    annotation_col <- data.frame(group=group_list)
    
    rownames(annotation_col) <- colnames(diff)
    
    plot <- pheatmap(diff,
                     annotation_col=annotation_col,
                     scale = "row",
                     show_rownames = F,
                     show_colnames =F,
                     color = colorRampPalette(c("navy", "white", "red"))(50),
                     fontsize = 10,
                     fontsize_row=3,
                     fontsize_col=3)
    return(plot)
}


unicox_survival <- function(covariates, dataset, time = "OS.time", status = "OS", ofname = "", verbose = F){
    # perform univariate Cox regression analysis on a vector of covariates
    # output a csv stating p.value/HR for each covariates
    # ------------------------------------------------------
    # covariates: a vector of covariates of interest, those covariates must present in the dataset
    # dataset: dataset with covariates and time/status on the columns, and observations on the rows
    # verbose: default showing only covariates with p.value < 0.05; set TRUE to get all covariates
    
    univ_formulas <- sapply(covariates, 
                            function(x) as.formula(paste0('Surv(', time, ', ', status, ') ~ ', x)))
  
    univ_models <- lapply(univ_formulas, 
                          function(x){coxph(x, data = dataset)})
  
    univ_results <- lapply(univ_models,
                           function(x){
                               x <- summary(x)
                               p.value <- signif(x$wald["pvalue"], digits=2)
                               HR <- signif(x$coef[2], digits=2);
                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                               HR <- paste0(HR, " (", 
                                            HR.confint.lower, "-", HR.confint.upper, ")")
                               res<-c(p.value,HR)
                               names(res)<-c("p.value","HR (95% CI for HR)")
                               return(res)
                           })
    
    res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
    res$p.value <- as.numeric(as.character(res$p.value))
    
    if (!verbose){
      res <- dplyr::filter(res, p.value < 0.05)
    }
    
    if (ofname != ""){
        write.csv(res, ofname)
    }
    return(res)
}


