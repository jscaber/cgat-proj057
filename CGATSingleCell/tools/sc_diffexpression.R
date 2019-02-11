#' filtering single cell data based on QC metrics
#'
#' WARNING: This script is work-in-progress
#' 
#' Example usage:
#' 
#' cgat sc-counts2counts --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv
#'
#' `feature_counts.tsv` is a table (tab-separated) of ngenes x ncells,
#' that is the genes are in rows and the columns are cells.
#'
#' `phenodata.tsv` is a table (tab-separated) of ncells x nfeatures,
#' that is rows are cells and features are in columns. The table should contain
#' a column called `sample_id` that will match the columns in the table
#' `feature_counts.tsv`.
#'
#' Features can then be selected in the `--factor` option to be
#' plotted.
#'
#' -> todo: parameterize detection of ERCC (pattern?)
#' -> todo: parameterize definition of mitochondrial genes - currently hardcoded for mouse.

## conda dependencies: bioconductor-scater r-cairo

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))
suppressMessages(library(sva))	
suppressMessages(library(zinbwave))
suppressMessages(library(DESeq2))
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
suppressMessages(library(RUVSeq))

source(file.path(Sys.getenv("R_ROOT"), "io.R"))
source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2018.archive.ensembl.org")

getmart_ensembl <- function(values){
    data<- getBM(
        filters= "ensembl_gene_id", 
        attributes= c("ensembl_gene_id", "mgi_symbol", "description"),
        values= values,
        mart= mart)
    data$description <- gsub("\t", "", data$description)
    return(data)
}
getmart_ensembl2 <- function(values){
    data<- getBM(
        filters= "ensembl_gene_id", 
        attributes= c("ensembl_gene_id", "mgi_symbol", "description","entrezgene"),
        values= values,
        mart= mart)
    data$description <- gsub("\t", "", data$description)
    return(data)
}
getmart_symbol <- function(values){
    data<- getBM(
        filters= "mgi_symbol", 
        attributes= c("entrezgene","ensembl_gene_id", "mgi_symbol", "description"),
        values= values,
        mart= mart)
    data$description <- gsub("\t", "", data$description)
    return(data)
}

start_plot <- function(section, height = 6, width = 6, type = "png") {
    file = get_output_filename(paste0(section, ".", type))
    Cairo(file = file,
          type = type,
          width = width,
          height = height,
          units="in",
          dpi = 300,
          bg = "white")
    #opar <- par(lwd=0.5)
}

end_plot <- function() {
    dev.off()
}

makeTPMtable  <- function(genelist, abundance){
  genelist.df <- getmart(genelist[grep("ENS",genelist)])
  
  genelist.df2 <-  as_tibble(t(matrix(rep(genelist[grep("ENS",genelist, invert=TRUE)],each=3), ncol = length(genelist[grep("ENS",genelist, invert=TRUE)]), nrow = 3)))
  colnames(genelist.df2) <- colnames(genelist.df)
  genelist.df <- bind_rows(genelist.df,as_tibble(genelist.df2))
  
  rownames(genelist.df) <- genelist.df$ensembl_gene_id
  genelist.names <- genelist.df[genelist,]$mgi_symbol

  dftemp <- as_tibble(t(abundance[genelist,]), rownames = "sample_id")
  dftemp <- dftemp %>% rename_at(vars(genelist), ~ genelist.names)
  dftemp$group <- colData(sce)[dftemp$sample_id,]$group
  dftemp$total_features <- colData(sce)[dftemp$sample_id,]$total_features
  dftemp
  return(dftemp)
}
plotTPMs <- function(dftemp){
  dftemp %>% 
    gather(key = "var", value="value", -group, -sample_id, -total_features) %>% 
    mutate(var = factor(var, levels=unique(var))) %>%
    ggplot(aes(x = group, y = value, color = total_features)) +
    geom_point(position = position_jitter(w = 0.15, h = 0)) +
    facet_wrap(~ var, scales = "free") + theme_bw() +
    ylab("Log Normalised Counts") + xlab ("Genotype") + scale_colour_gradient(low="red", high="yellow")
}


run <- function(opt) {


    flog.info("Reading in experiment")
    sce <- read_single_cell_experiment_from_rds(opt$rds_filename)
    sce$group = as.factor(sce$group)
    endog_genes <- !rowData(sce)$is_feature_control
    keep <- rowSums(counts(sce) >= 5) >= 25
    if(opt$filter == TRUE){
        flog.info(paste0("After filtering: keeping ",
                         length(keep), " genes out of ", length(rownames(sce))))
        zinb <- sce[keep,]
    }
    else{zinb <- sce}

    flog.info("Running zinbwave...")
    flog.info(paste0("Model: ", opt$zinb_model))
    zinb$condition <- factor(zinb$group)
    # we need to reorganize the assays in the SumExp from splatter
    nms <- c("counts", setdiff(assayNames(zinb), "counts"))
    assays(zinb) <- assays(zinb)[nms]
    # epsilon setting as recommended by the ZINB-WaVE integration paper
    zinb <- zinbwave(zinb, K=1, X=opt$zinb_model, BPPARAM=SerialParam(), epsilon=1e12)
    dds <- DESeqDataSet(zinb, design = opt$deseq_model)
    dds <- DESeq(dds, test="LRT", reduced=~1, minmu=1e-6)

    flog.info("Running SVA and RUVs...")
    # SVA Pacakge script
    dat  <- counts(dds, normalized = TRUE)
    idx  <- rowMeans(dat) > 1
    dat  <- dat[idx, ]
    mod  <- model.matrix(~ group, colData(dds))
    mod0 <- model.matrix(~   1, colData(dds))
    svseq <- svaseq(dat, mod, mod0, n.sv = 2)
    dds$SV1 <- svseq$sv[,1]
    dds$SV2 <- svseq$sv[,2]
    # RUV Package Script
    scIdx <- matrix(-1, ncol = max(table(dds$group)), nrow = length(levels(dds$group)))
    i <- 1
    for(groupitem in levels(dds$group)){    
        tmp <- which(dds$group == groupitem)        
        scIdx[i, 1:length(tmp)] <- tmp
        i <- i + 1
    }
    cIdx <- rownames(dds)
    ruvs <- RUVs(counts(dds), cIdx, k = 2, scIdx = scIdx, isLog = FALSE)
    dds$W_1 <- ruvs$W[,"W_1"]
    dds$W_2 <- ruvs$W[,"W_2"]
    # Plotting versus bias
    p1 <- ggplot(as_tibble(colData(dds)), aes(SV1, total_features)) + geom_point() + ggtitle("SVA Variable 1")
    p2 <- ggplot(as_tibble(colData(dds)), aes(SV2, total_features)) + geom_point() + ggtitle("SVA Variable 2")
    p3 <- ggplot(as_tibble(colData(dds)), aes(W_1, total_features)) + geom_point() + ggtitle("RUV Variable 1")
    p4 <- ggplot(as_tibble(colData(dds)), aes(W_2, total_features)) + geom_point() + ggtitle("RUV Variable 2")
    gridExtra::grid.arrange(p1,p2,p3,p4)
    png('Surrogatevariables_vs_bias.png', width = 6, height = 6, units = 'in', res = 300)
    gridExtra::grid.arrange(p1,p2,p3,p4)
    dev.off()


    flog.info("Running DESeq2...")
    ## Run with RUVSeq batch model
    design(dds) <- opt$deseq_model
    dds <- DESeq(dds, test="Wald", betaPrior = TRUE)
    res <- results(dds,independentFiltering = FALSE)
    flog.info(print(summary(res)))

    flog.info("... plotting dispersion estimates")
    ## Plot dispersion estimates
    keepForDispTrend <- rowSums(counts(dds) >= 10) >= 25
    dds2 <- estimateDispersionsFit(dds[keepForDispTrend,])
    start_plot("Dispersion")
    plotDispEsts(dds2)
    end_plot()

    flog.info("... plotting MA")
    ## MA Plot
    start_plot("MAPlot")
    plotMA(dds, ylim = c(-3,3))
    end_plot()

    flog.info("... saving DE data")
    ## Save DE data
    resSig <- subset(res, padj < 0.1)
    data <- getmart(rownames(resSig))
    resSig$symbol<-data$mgi_symbol[match(rownames(resSig), data$ensembl_gene_id)]
    resSig$desc<-data$description[match(rownames(resSig), data$ensembl_gene_id)]
    write.table(resSig, "results.tsv", sep = "\t")
    write.table(res, "results_full.tsv", sep = "\t")
    resdf <- data.frame(geneid=rownames(res), pvalue=res$pvalue*sign(res$log2FoldChange))
    write.table(resdf, "px_results_pvalue.gene.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    resdf <- data.frame(geneid=rownames(res), l2fc=res$log2FoldChange)
    write.table(resdf, "px_results_l2fc.gene.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    flog.info("... plotting P histogram")
    ## Plot P value Histogram
    start_plot("PHistogram")
    hist(res$pvalue,breaks=50, col='skyblue', xlab="p-value", main="P-value Histogram")
    end_plot()

    flog.info("... performing permutations")
    if(opt$perms >0){    
        designperm <- colData(dds)
        ddsperm = dds
        x = 0
        i = 1
        y= list()
        while(i <= opt$perms) {
            designperm$group = as.factor(sample(c("Ta1neg","Ta1pos"), length(ddsperm$group), replace=TRUE, prob=c(0.5, 0.5)))
            if(sum(designperm$group == "Ta1neg") != table(dds$group)["Ta1neg"]) next
            colData(ddsperm) <- designperm
            ddsperm <- DESeq(ddsperm, test="Wald", betaPrior = TRUE)
            resperm <- results(ddsperm,independentFiltering = FALSE)
            x[i] = length(subset(resperm, padj < 0.1)$padj)
            y[[i]] = ddsperm$group
            i = i+1
        }
        start_plot("Simulations")
        theme_set(theme_gray(base_size = 18))
        sims = qplot(x,
            geom="histogram",
            breaks=seq(0, 20, by = 1),fill=I("grey"), col=I("black"),
            main = "Histogram of DE experiments\n with random group labels", 
            xlab = "Number of differentially expressed genes",
            ylam = "Number of simulations") +
            geom_vline(xintercept = length(rownames(subset(res, padj < 0.1)))) +
            theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=22))
        print(sims)
        end_plot()
        flog.info(paste0("... Permutation p value: ",
                         length(x[x < length(rownames(subset(res, padj < 0.1)))])/length(x)))
        z <- list()
        for(i in 0:length(x)){
          z[i] <- paste( unlist(y[i]), collapse=' ')
        }
        z <- unlist(z)
        df <- data.frame(number=x, combination=z)
        df
        write_tsv(df,"Simulations.tsv")
    }

    flog.info("... plotting downregulated genes")    
    ## Plot Top Downregulated Genes
    genelist <- rownames(res[ order( res$log2FoldChange ), ][0:9,])
    dftemp <- makeTPMtable(genelist, logcounts(sce))
    start_plot("Downregulated")
      plotTPMs(dftemp)
    end_plot()

    flog.info("... plotting downregulated genes")
    genelist <- rownames(res[ order( -res$log2FoldChange ), ][0:9,])
    dftemp <- makeTPMtable(genelist, logcounts(sce))
    start_plot("Upregulated")
      plotTPMs(dftemp)
    dev.off()

    flog.info("... plotting significant genes")
    genelist <- rownames(res[ order( -resSig$padj), ])
    if(length(genelist) > 9){
        genelist <- genelist[0:9]}
    dftemp <- makeTPMtable(genelist, logcounts(sce))
    start_plot("significant")
      plotTPMs(dftemp)
    dev.off()

    ## Plot RRBP1 and its binding partner NDUFA7
    genelist <- c("ENSMUSG00000027422","ENSMUSG00000041881")
    dftemp <- makeTPMtable(genelist, logcounts(sce))
    png('RRbp1 Genes.png', width = 10, height = 6, units = 'in', res = 300)
      plotTPMs(dftemp)
    dev.off()


    resrnk<-res
    data <- getmarte(rownames(resrnk))
    resrnk$symbol <- data$mgi_symbol[match(rownames(resrnk), data$ensembl_gene_id)]
    resrnk$desc <-data$description[match(rownames(resrnk), data$ensembl_gene_id)]
    resrnk$entrezgene <-data$entrezgene[match(rownames(resrnk), data$ensembl_gene_id)]
    rnk.df <- as_tibble(resrnk[,c("entrezgene","log2FoldChange")]) %>% na.omit()
    rnk <- rnk.df$log2FoldChange
    names(rnk) <- rnk.df$entrezgene
    rnk <- rnk[isUnique(names(rnk))]

    for(pathway in opt$pathways){
        pathways <- gmtPathways(pathway)
        fgseaRes <- fgsea(pathways = pathways, 
                          stats = rnk,
                          minSize=15,
                          maxSize=500,
                          nperm=10000)
        fwrite(fgseaRes, file="fgseaResGO.txt", sep="\t", sep2=c("", " ", ""))
        topPathwaysUp <- fgseaRes[ES > 0,][head(order(pval), n=10),]$pathway
        png(paste0('GSEA_GO_up_',file_path_sans_ext(basename(pathway)),'.png'),
            width =15, height = 3, units = 'in', res = 600)
        plotGseaTable(pathways[topPathwaysUp], rnk, fgseaRes, 
                      gseaParam = 0.5, colwidths = c(10,2,1,1,1))
        dev.off()


        topPathwaysDown <- fgseaRes[ES < 0,][head(order(pval), n=10),]$pathway
        plotGseaTable(pathways[topPathwaysDown], rnk, fgseaRes, 
        gseaParam = 0.5)
        png(paste0('GSEA_GO_down_',file_path_sans_ext(basename(pathway)),'.png'),
            width =15, height = 3, units = 'in', res = 600)
        plotGseaTable(pathways[topPathwaysDown], rnk, fgseaRes, 
                      gseaParam = 0.5, colwidths = c(10,2,1,1,1))
        dev.off()
        }
}

main <- function() {

    option_list <- list(
        make_option(
            "--rds-filename",
            dest = "rds_filename",
            type = "character",
            default = "sce.rds",
            help = paste("filename with input data of counts")
        ),
        make_option(
            "--zinbwave-model",
            dest = "zinb_model",
            type = "character",
            default = "~total_features",
            help = paste("model for zinbwave")
        ),
        make_option(
            "--deseq-model",
            dest = "deseq_model",
            type = "character",
            default = "~group",
            help = paste("model for DESeq2")
        ),
        make_option(
            "--filter",
            dest = "filter",
            type = "character",
            default = "TRUE",
            help = paste("model for DESeq2")
        ),
        make_option(
            "--permutations",
            dest = "perms",
            type = "integer",
            default = 0,
            help = paste("number of permutations")
        ),
        make_option(
            "--golibs",
            dest = "golibs",
            type = "character",
            default = "",
            help = paste("number of permutations")
        )
    )
    opt <- experiment_start(option_list = option_list,
                            description = description)

    if (!is.null(opt$golibs)) {
        opt$golibs = unlist(strsplit(opt$golibs, ","))
    }
    run(opt)
    
    experiment_stop()
}

main()
