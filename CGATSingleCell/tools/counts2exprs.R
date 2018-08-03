#' convert single cell counts to expression metrics
#'
#' WARNING: This script is work-in-progress
#' 
#' cgat counts2exprs --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv > exrs.tsv
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

## conda dependencies: bioconductor-scater r-cairo bioconductor-rtracklayer edger

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))
suppressMessages(library(Cairo))
suppressMessages(library(scater))
suppressMessages(library(rtracklayer))
suppressMessages(library(edgeR))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))
source(file.path(Sys.getenv("R_ROOT"), "io.R"))

#' Create a list of GC content
calc_GC_length <- function(x) {
    sum(elementMetadata(x)$widths)
}


start_plot <- function(section, height = 6, width = 10, type = "png") {
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


plotQC <- function(data_sceset, slot, section, opt) {
    endog_genes <- !rowData(data_sceset)$is_feature_control

    for (factor in opt$factors) {
        flog.info(paste("creating plots for factor", factor))

        flog.info("creating RLE plot")
        start_plot(paste("rle_plot", slot, factor, section, sep = "-"))
        print(plotRLE(
            data_sceset[endog_genes, ],
            exprs_mats = list(Raw = "logcounts_raw", Normed = slot),
            exprs_logged = c(TRUE, TRUE),
            colour_by = factor
        ))
        end_plot()
    }
}


run <- function(opt) {

    options(stringsAsFactors = FALSE)
    set.seed(1234567)
    if (is.null(opt$rds_filename)) {
        all_sceset <- read_single_cell_experiment_from_tables(opt$counts_filename, opt$phenotypes_filename)
        flog.info("calculating SCRAN QC metrics using spike-in and mitochondrial controls")
        all_sceset <- scater::calculateQCMetrics(all_sceset,
                                                 feature_controls = list(
                                                     ERCC = isSpike(all_sceset, type="ERCC"),
                                                     Mt = isSpike(all_sceset, type="Mt")))
        all_sceset <- scater::arrange(all_sceset, group)
    } else {
        flog.info(paste("loading single cell data set from", opt$rds_filename))
        all_sceset <- readRDS(opt$rds_filename)
    }
        
    if (is.null(opt$genelengths_filename)) {
        flog.info(paste("reading geneset from", opt$geneset_filename))
        GTF <- rtracklayer::import.gff(opt$geneset_filename,
                                       format="gtf",
                                       genome="GRCm38.88",
                                       feature.type="exon")
        flog.info("building gene lengths")
        grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
        reducedGTF <- unlist(grl, use.names=T)
        elementMetadata(reducedGTF)$gene_id <- rep(names(grl), lengths(grl))
        elementMetadata(reducedGTF)$widths <- width(reducedGTF)
        
        gene_lengths <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length)
        gene_lengths <- as.data.frame(gene_lengths)
        fn = get_output_filename("genelengths.tsv")
        flog.info(paste("writing gene lengths to", fn))

        write.table(gene_lengths,
                file = fn,
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                col.names = NA)
    } else {
        flog.info(paste("reading gene lengths from", opt$geneset_filename))
        gene_lengths <- read.table(opt$genelengths_filename,
                                   sep = "\t",
                                   row.names = 1,
                                   header = TRUE)
    }

    flog.info("normalizing with scater")
    ## calculate log2 transformed expression values
    ## 1. Divide each count by size factor (or scaled library size if no size factor is available)
    ## 2. Add a pseudo-count
    ## 3. Log-transform -> store in log-counts
    ## Find out how to compute custom size factors
    all_sceset <- scater::normalize(all_sceset)

    ## scater:
    ## normalize -> synonym to normalizeSCE
    ## normalizeSCE -> normalized expression values from count data using size factors stored in object
    ## result: normcounts if log=FALSE
    ##         logcounts if log=TRUE
    ## 
    ## normalizeExprs -> normalized expression values (deprecated)
    ## use edgeR::calcNormFactors(), normalize(), limma::removeBatchEffect() directly instead
    
    ## norm_factors <- edgeR::calcNormFactors.default(exprs_mat_for_norm, method=method)
    ## lib_size <- .colSums(exprs_mat_for_norm)
    ##
    ## size_factors <- norm_factors * lib_size
    ## size_factors <- size_factors / mean(size_factors)
    ## sizeFactors(object) <- size_factors
    ## then use cpm, etc using use_size_factors=TRUE
    
    for (method in opt$method) {
        flog.info(paste("building normalized data with method", method))
        if (method == "cpm") {
            SingleCellExperiment::cpm(all_sceset) <- scater::calculateCPM(all_sceset)
            counts_data <- SingleCellExperiment::cpm(all_sceset)
            slot <- "cpm"
        } else if (method == "tmm") {
            all_sce <- normaliseExprs(all_sce, method = "TMM")
        } else if (method == "logcounts") {
            counts_data <- logcounts(all_sceset)
            slot <- "logcounts"
        } else if ((method == "rpkm") || (method == "logcounts-rpkm")) {
            if (method == "rpkm") {
                featurecounts <- counts(all_sceset)
            } else {
                featurecounts <- logcounts(all_sceset)
            }
            shared_genes <- intersect(row.names(gene_lengths),
                                      row.names(featurecounts))
            flog.info(paste("gene names: shared=", length(shared_genes),
                            "out ouf", nrow(gene_lengths), nrow(featurecounts)))
        
            featurecounts <- featurecounts[shared_genes, ]
            gene.lengths <- gene_lengths[shared_genes, ]
        
            counts_data <- edgeR::rpkm(featurecounts,
                                       gene.length = gene.lengths)
            slot <- "logcounts"
        } else {
            flog.warn(paste("unknown method", method))
            next
        }
        fn = get_output_filename(sprintf("%s.tsv", method))
        flog.info(paste("writing data to", fn))
        write.table(counts_data,
                    file = fn,
                    append = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    row.names = TRUE,
                    col.names = NA)

        plotQC(all_sceset, slot, method, opt)
    }
}


main <- function() {

    option_list <- list(
        make_option(
            "--rds-filename",
            dest = "rds_filename",
            type = "character",
            default = "sce.rds",
            help = paste("filename with single cell experiment data")
        ),
        make_option(
            "--counts-filename",
            dest = "counts_filename",
            type = "character",
            default = "featurecounts.tsv",
            help = paste("filename with input data of counts")
        ),
        make_option(
            "--phenotypes-filename",
            dest = "phenotypes_filename",
            type = "character",
            default = "phenodata.tsv",
            help = paste("filename with phenotype data")
        ),
        make_option(
            "--genelengths-filename",
            dest = "genelengths_filename",
            type = "character",
            default = "genelengths.tsv",
            help = paste("filename with geneset data")
        ),
        make_option(
            "--geneset-filename",
            dest = "geneset_filename",
            type = "character",
            default = "geneset.gtf.gz",
            help = paste("filename with geneset data")
        ),
        make_option(
            "--geneset-genome",
            dest = "geneset_genome",
            type = "character",
            default = "GRCm38.88",
            help = paste("geneset filename")
        ),
        make_option(
            "--method",
            dest = "methods",
            type = "character",
            default = "rpkm",
            help = paste("normalization method to apply. Multiple methods can be applied ",
                         "as a ',' separated list")
        ),
        make_option(
            "--factor",
            dest = "factors",
            type = "character",
            default = "group,collection_date",
            help = paste("factors to colour QC plots by.")
        )
    )
    opt <- experiment_start(option_list = option_list,
                            description = description)

    if (!is.null(opt$methods)) {
        opt$methods = unlist(strsplit(opt$methods, ","))
    }
    if (!is.null(opt$factors)) {
        opt$factors = unlist(strsplit(opt$factors, ","))
    }

    run(opt)
    
    experiment_stop()
}

main()
