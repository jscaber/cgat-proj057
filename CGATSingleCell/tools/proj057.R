#' Project Pipeline for CGAT Project 057
#'
#' WARNING: This script is work-in-progress
#' 
#' Example usage:
#' 
#' cgat-singlecell proj057 --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv
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
suppressMessages(library(Cairo))
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(SC3))
suppressMessages(library(Rtsne))
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))	
suppressMessages(library(RUVSeq))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(sva))	
suppressMessages(library(zinbwave))
suppressMessages(library(DESeq))

source(file.path(Sys.getenv("R_ROOT"), "io.R"))
source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


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

# Downsampling funcion (adapted from Hemberg lab)
down_sample_matrix <- function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(unlist(lapply(x, function(y) {
            rbinom(1, y, prob)
        })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}

# Name: plot_normalisation
# Function: Normalises Single Cell dataset using multiple methods and plots the result
# Inputs: single cell dataset, list of endogenous genes (not spike ins, and not mitochondrial	)
# Outputs: png files showing effect of normalisation on PCA
plot_normalisation <- function(sce.raw, endog_genes) {
    sce <- sce.raw

    flog.info("Raw Counts - no normalisation")
    start_plot("rawcounts_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        exprs_values = "counts",
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()

    flog.info("Log count normalisation")
    logcounts(sce) = log(counts(sce)+1)
    start_plot("logcounts_pca")
    print(scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"))
    end_plot()
    start_plot("logcounts_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()

    flog.info("CPM normalisation")
    logcounts(sce) <- log2(calculateCPM(sce) + 1)
    start_plot("cpm_pca")
    print(scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"))
    end_plot()
    start_plot("cpm_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()

    flog.info("Downsampling normalisation")
    logcounts(sce) <- log2(down_sample_matrix(counts(sce)) + 1)
    start_plot("downsampling_pca")
    print(scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"))
    end_plot()
    start_plot("downsampling_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()

    flog.info("SCRAN normalisation with ERCCs")
    sce <- sce.raw
    sce <- computeSpikeFactors(sce, general.use=TRUE)
    sce <- normalize(sce)
    start_plot("scran_ercc_sizefactorcorrelation")
    plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
         xlab="Library size (millions)", ylab="Size factor",
         col=c("red", "black")[sce$group], pch=16)
    legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
           legend=levels(sce$group))
    end_plot()
    start_plot("scran_ercc_pca")
    print(scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"))
    end_plot()
    start_plot("scran_ercc_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()
    fit <- trendVar(sce, parametric=TRUE)
    decomp <- decomposeVar(sce, fit)
    top.hvgs <- order(decomp$bio, decreasing=TRUE)
    start_plot("scran_ercc_fit")
    plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
    o <- order(decomp$mean)
    lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
    points(fit$mean, fit$var, col="red", pch=16)
    end_plot()

    flog.info("SCRAN normalisation without ERCCs")
    sce <- sce.raw
    sce <- computeSumFactors(sce)
    sce <- normalize(sce)
    start_plot("scran_noercc_sizefactorcorrelation")
    plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
         xlab="Library size (millions)", ylab="Size factor",
         col=c("red", "black")[sce$group], pch=16)
    legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
           legend=levels(sce$group))
    end_plot()
    start_plot("scran_noercc_pca")
    print(scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"))
    end_plot()
    start_plot("scran_noercc_explanatory")
    print(plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20))
    end_plot()
    alt.fit <- trendVar(sce, use.spikes=FALSE) 
    alt.decomp <- decomposeVar(sce, alt.fit)
    alt.top.hvgs <- order(alt.decomp$bio, decreasing=TRUE)
    head(alt.decomp[alt.top.hvgs,])
    plot(alt.decomp$mean, alt.decomp$total, xlab="Mean log-expression", ylab="Variance")
    alt.o <- order(alt.decomp$mean)
    lines(alt.decomp$mean[alt.o], alt.decomp$tech[alt.o], col="red", lwd=2)
    start_plot("scran_noercc_fit")
    plot(alt.decomp$mean, alt.decomp$total, xlab="Mean log-expression", ylab="Variance")
    alt.o <- order(alt.decomp$mean)
    lines(alt.decomp$mean[alt.o], alt.decomp$tech[alt.o], col="red", lwd=2)
    end_plot()

}

# Name: check_normalisation
# Function: Checks how well normalisations perform against known ERCC concentrations
# Inputs: single cell dataset, list of endogenous genes (not spike ins, and not mitochondrial	)
# Outputs: png files showing effect of normalisation on PCA
check_normalisation <- function(sce.raw, ERCCconc) {

    #create melted tables for Raw Counts data
    sce <- sce.raw
    ERCC <- as_tibble(counts(sce[rownames(sce)[grep("ER", rownames(sce))],]),rownames = "ID")
    ERCCplot <- inner_join(ERCC, ERCCconc[,1:2], by = "ID")
    dfplot_raw <- melt(ERCCplot,id.vars = c("ID","Mix1"))
    dfplot_nonzero_raw <- subset(dfplot_raw, dfplot_raw$value != 0)
    dfplot_big_raw <- subset(dfplot_raw, log(dfplot_raw$Mix1) >6.25)

    #create melted tables for SCRAN data
    sce <- computeSumFactors(sce)
    sce <- normalize(sce)
    ERCC <- as_tibble(logcounts(sce[rownames(sce)[grep("ER", rownames(sce))],]),rownames = "ID")
    ERCCplot <- inner_join(ERCC, ERCCconc[,1:2], by = "ID")
    dfplot_scran <- melt(ERCCplot,id.vars = c("ID","Mix1"))
    dfplot_nonzero_scran <- subset(dfplot_scran, dfplot_scran$value != 0)
    dfplot_big_scran <- subset(dfplot_scran, log(dfplot_scran$Mix1) >6.25)

    #create melted tables for ERCC normalised SCRAN data
    sce.ercc <- computeSpikeFactors(sce.raw, general.use=TRUE)
    sce.ercc <- normalize(sce.ercc)
    ERCC <- as_tibble(logcounts(sce.ercc[rownames(sce.ercc)[grep("ER", rownames(sce.ercc))],]),
                      rownames = "ID")
    ERCCplot <- inner_join(ERCC, ERCCconc[,1:2], by = "ID")
    dfplot_scranercc <- melt(ERCCplot,id.vars = c("ID","Mix1"))
    dfplot_nonzero_scranercc <- subset(dfplot_scranercc, dfplot_scranercc$value != 0)
    dfplot_big_scranercc <- subset(dfplot_scranercc, log(dfplot_scranercc$Mix1) >6.25)

    #Plotting all counts
    flog.info("Plotting ERCC vs all counts")
    fit1 <- summary(lm(log(dfplot_raw$Mix1)~dfplot_raw$value))
    p1 <- ggplot(dfplot_raw,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) + 
          ggtitle(label = "Raw Counts", paste0("R^2=",as.character(fit1[9]))) + ylab("counts")
    fit2 <- summary(lm(log(dfplot_raw$Mix1)~log(dfplot_raw$value+1)))    
    p2 <- ggplot(dfplot_raw,aes(log(Mix1),log(value+1))) +  geom_point() + geom_smooth(method = lm) +
          ggtitle(label = "Log Counts", paste0("R^2=",as.character(fit2[9]))) + ylab("log(counts)")
    fit3 <- summary(lm(log(dfplot_scran$Mix1)~dfplot_scran$value))    
    p3 <- ggplot(dfplot_scran,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) +
          ggtitle(label = "Scran", paste0("R^2=",as.character(fit3[9]))) + ylab("normalised counts")
    fit4 <- summary(lm(log(dfplot_scranercc$Mix1)~dfplot_scranercc$value))
    p4<- ggplot(dfplot_scranercc,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) +
         ggtitle(label = "Scran ERCC", paste0("R^2=",as.character(fit4[9]))) + ylab("normalised counts")
    start_plot("ERCC_vs_counts_all")
    gridExtra::grid.arrange(p1,p2,p3,p4)
    end_plot()

    # Plotting only high abundance counts (no zeros)
    flog.info("Plotting ERCC vs non-zero counts")
    fit1 <- summary(lm(log(dfplot_big_raw$Mix1)~dfplot_big_raw$value))
    p1 <- ggplot(dfplot_big_raw,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) + 
          ggtitle(label = "Raw Counts", paste0("R^2=",as.character(fit1[9]))) + ylab("counts")
    fit2 <- summary(lm(log(dfplot_big_raw$Mix1)~log(dfplot_big_raw$value+1)))  
    p2 <- ggplot(dfplot_big_raw,aes(log(Mix1),log(value+1))) +  geom_point() + 
          geom_smooth(method = lm) +
          ggtitle(label = "Log Counts", paste0("R^2=",as.character(fit2[9]))) + ylab("log(counts)")
    fit3 <- summary(lm(log(dfplot_big_scran$Mix1)~dfplot_big_scran$value))    
    p3 <- ggplot(dfplot_big_scran,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) +
          ggtitle(label = "Scran", paste0("R^2=",as.character(fit3[9]))) + ylab("normalised counts")
    fit4 <- summary(lm(log(dfplot_big_scranercc$Mix1)~dfplot_big_scranercc$value))
    p4<- ggplot(dfplot_big_scranercc,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) +
         ggtitle(label = "Scran ERCC", paste0("R^2=",as.character(fit4[9]))) + ylab("normalised counts")
    start_plot("ERCC_vs_counts_nozeros")
    gridExtra::grid.arrange(p1,p2,p3,p4)
    end_plot()
}

# Name: runRUV
# Function: Runs RUVSeq on Single Cell Experiment Set
# Inputs: 
# Outputs: 
runRUV <- function(sce, endog_genes, ERCCconc) {

    flog.info("Running RUVs")
    options(stringsAsFactors = FALSE)
    qclust <- quickCluster(sce, min.size = 30)
    sce <- computeSumFactors(sce, sizes = 15, clusters = qclust)
    sce <- normalize(sce)

    # Establishing Group Matrix for RUVs

    scIdx <- matrix(-1, ncol = max(table(sce$group)), nrow = 2)
    tmp <- which(sce$group == "Ta1neg")
    scIdx[1, 1:length(tmp)] <- tmp
    tmp <- which(sce$group == "Ta1pos")
    scIdx[2, 1:length(tmp)] <- tmp
    cIdx <- rownames(sce)

    # Running RUVs
    ruvs <- RUVs(counts(sce), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
    assay(sce, "logcounts") <- log2(
        t(t(ruvs$normalizedCounts) / 
            colSums(ruvs$normalizedCounts) * 1e6) + 1)

    flog.info("Plotting RUVs")
    # PCA
    plot<- scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features"
    )
    start_plot("ruv_pca")
    print(ggplot(plot$data, aes(X, Y, colour = colour_by)) + geom_point() + theme_bw() +
          scale_color_manual(labels = c("Ta1-", "Ta1+"), values=c("red", "green")) + 
          labs(colour = "Genotype") + ylab("Principal Component 2") + 
          xlab("Principal Component 1"))
    end_plot()

    # tSNE
    set.seed(345873945)
    joint_tsne <- Rtsne(t(logcounts(sce)), perplexity = 15)
    cell_type_labels <- factor(c(as.character(colData(sce)$group), as.character(colData(sce)$group)))
    tsnedf <- as.data.frame(joint_tsne$Y)
    tsnedf$group <- sce$group
    start_plot("ruv_tsne")
    print(ggplot(tsnedf, aes(V1, V2, colour = group)) + geom_point() + theme_bw() +
          scale_color_manual(labels = c("Ta1-", "Ta1+"), values=c("red", "green")) +
          labs(colour = "Genotype") + ylab("tSNE Dimension 2") + xlab("tSNE Dimension 1"))
    end_plot()

    # Explanatory
    start_plot("ruv_explanatory")
    print(plotExplanatoryPCs(
          sce[endog_genes, ],
          variable = c("total_features", "total_counts"),
          npcs_to_plot = 20))
    end_plot()

    # Correlation with ERCCs
    ERCC <- as_tibble(logcounts(sce[rownames(sce)[grep("ER", rownames(sce))],]),rownames = "ID")
    ERCCplot <- inner_join(ERCC, ERCCconc[,1:2], by = "ID")
    dfplot_ruv <- melt(ERCCplot,id.vars = c("ID","Mix1"))
    dfplot_nonzero_ruv <- subset(dfplot_ruv, dfplot_ruv$value != 0)
    dfplot_big_ruv <- subset(dfplot_ruv, log(dfplot_ruv$Mix1) >6.25)

    fit1 <- summary(lm(log(dfplot_ruv$Mix1)~dfplot_ruv$value))
    p1 <- ggplot(dfplot_ruv,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) + 
          ggtitle(label = "Log Counts", paste0("R^2=",as.character(fit1[9]))) + ylab("counts")
    fit2 <- summary(lm(log(dfplot_big_ruv$Mix1)~log(dfplot_big_ruv$value+1)))
    p2 <- ggplot(dfplot_big_ruv,aes(log(Mix1),log(value+1))) +  geom_point() + 
          geom_smooth(method = lm) +
          ggtitle(label = "Log Counts", paste0("R^2=",as.character(fit2[9]))) + ylab("log(counts)")
    start_plot("ERCC_vs_counts_RUV")
    gridExtra::grid.arrange(p1,p2, nrow=1)
    end_plot()

    return(sce)    
}

run <- function(opt) {

    options(stringsAsFactors = FALSE)
    set.seed(1234567)

    sce <- read_single_cell_experiment_from_rds(opt$rds_filename)
    sce$group = as.factor(sce$group)
    sce.raw <- sce
    endog_genes <- !rowData(sce)$is_feature_control

    # Check for mouse cell cycle markers
    flog.info("calculating cell cycle markers")
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    assigned <- cyclone(sce, pairs=mm.pairs)
    write.table(table(assigned$phases), "cell_cycle_markers.tsv", 
                quote = FALSE, sep = "\t", row.names = FALSE)

    flog.info("plotting graphs after normalisation")
    plot_normalisation(sce.raw, endog_genes)

    flog.info("Reading in ERCC concentrations")
    ERCCconc <-  read_tsv(opt$ERCCpath,col_names = TRUE )
    ERCCconc <- ERCCconc[,c(2,4,5)]
    colnames(ERCCconc) <- c("ID", "Mix1", "Mix2")
    ERCCconc$ID<- gsub("-",".", ERCCconc$ID)
    
    flog.info("Plotting Normalisations vs ERCCs")
    check_normalisation(sce.raw, ERCCconc)

    flog.info("Running RUVs and plotting vs ERCCs")
    sce.ruv <- runRUV(sce.raw, endog_genes, ERCCconc)

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
            "--ERCC",
            dest = "ERCCpath",
            type = "character",
            default = "ERCC.tsv",
            help = paste("path to file with ERCC concentrations")
        ),
        make_option(
            "--norm",
            dest = "normalisation",
            type = "character",
            default = "group,collection_date",
            help = paste("which normalisation to use.")
        )
    )
    opt <- experiment_start(option_list = option_list,
                            description = description)

    if (!is.null(opt$factors)) {
        opt$factors = unlist(strsplit(opt$factors, ","))
    }
    run(opt)
    
    experiment_stop()
}

main()
