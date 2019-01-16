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
suppressMessages(library(Cairo))
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(SC3))
suppressMessages(library(Rtsne))
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(scRNA.seq.funcs))
suppressMessages(library(RUVSeq))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(kBET))
suppressMessages(library(sva))
suppressMessages(library(edgeR))
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
          bg = "white",
          pointsize = 12)
    opar <- par(lwd=0.5)
}

end_plot <- function() {
    dev.off()
}

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


plot_normalisation <- function(sce.raw, endog_genes) {
    sce <- sce.raw

    flog.info("Raw Counts - no normalisation")
    start_plot("rawcounts_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        exprs_values = "counts",
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
    end_plot()

    flog.info("Log count normalisation")
    logcounts(sce) = log(counts(sce)+1)
    start_plot("logcounts_pca")
    plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features")
    end_plot()
    start_plot("logcounts_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
    end_plot()

    flog.info("CPM normalisation")
    logcounts(sce) <- log2(calculateCPM(sce) + 1)
    start_plot("cpm_pca")
    plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features")
    end_plot()
    start_plot("cpm_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
    end_plot()

    flog.info("Downsampling normalisation")
    logcounts(sce) <- log2(down_sample_matrix(counts(sce)) + 1)
    start_plot("downsampling_pca")
    plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features")
    end_plot()
    start_plot("downsampling_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
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
    plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features")
    end_plot()
    start_plot("scran_ercc_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
    end_plot()
    fit <- trendVar(sce, parametric=TRUE)
    decomp <- decomposeVar(sce, fit)
    top.hvgs <- order(decomp$bio, decreasing=TRUE)]
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
    plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        size_by = "total_features")
    end_plot()
    start_plot("scran_noercc_explanatory")
    plotExplanatoryPCs(
        sce[endog_genes, ],
        variable = c("total_features", "total_counts"),
        npcs_to_plot = 20
    end_plot()
    fit <- trendVar(sce, parametric=TRUE)
    decomp <- decomposeVar(sce, fit)
    top.hvgs <- order(decomp$bio, decreasing=TRUE)]
    start_plot("scran_noercc_fit")
    plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
    o <- order(decomp$mean)
    lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
    points(fit$mean, fit$var, col="red", pch=16)
    end_plot()

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
    write.table(table(assigned$phases), "cell_cycle_markers.tsv")

    flog.info("plotting graphs after normalisation")
    plot_normalisations(sce.raw, endog_genes)

    
    
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
