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
suppressMessages(library(Matrix))
suppressMessages(library(dendextend))
suppressMessages(library(scrattch.hicat))
suppressMessages(library(matrixStats))
suppressMessages(library(tasic2016data))
suppressMessages(library(WGCNA))
suppressMessages(library(biomaRt))
suppressMessages(library(scater))
suppressMessages(library(Rtsne))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))
source(file.path(Sys.getenv("R_ROOT"), "io.R"))

mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2018.archive.ensembl.org")

getmart <- function(values){
    data<- getBM(
        filters= "ensembl_gene_id", 
        attributes= c("ensembl_gene_id", "mgi_symbol", "description", "entrezgene"),
        values= values,
        mart= mart)
    data$description <- gsub("\t", "", data$description)
    return(data)
}
getmartentrez <- function(values){
    data<- getBM(
        filters= "mgi_symbol", 
        attributes= c("ensembl_gene_id", "mgi_symbol", "description", "entrezgene"),
        values= values,
        mart= mart)
    data$description <- gsub("\t", "", data$description)
    return(data)
}



run <- function(opt) {

    flog.info("Reading in Allen design data")
    anno <- read.csv(opt$allen_design)
    anno$cluster <- as.character(anno$cluster)
    if(length(opt$allen_filter) != 0){
        anno <- anno[grep(opt$allen_filter, anno$cluster),]}
    anno$cluster_id = as.numeric(as.factor(anno$cluster))
    ref.cl.df <- as.data.frame(unique(anno[,c("cluster","class","cluster_id")]))
    ref.cl.df <- ref.cl.df[order(ref.cl.df$cluster_id),]
    row.names(ref.cl.df) <- ref.cl.df$cluster_id
    ref.cl <- setNames(factor(anno$cluster_id), anno$sample_name)

    flog.info("Reading in Allen expression data")
    dat <- read.csv(opt$allen_datamatrix,row.names = 1)
    dat.matrix <- as.matrix(dat)
    if(length(opt$allen_filter) != 0){
        dat.matrix <- dat.matrix[,names(ref.cl)]}
    norm.dat <- Matrix(scrattch.hicat::cpm(dat.matrix), sparse = TRUE)
    norm.dat@x <- log2(norm.dat@x+1)

    flog.info("Selecting marker genes")
    marker_results <- select_markers(norm.dat, ref.cl)
    marker_genes <- marker_results$markers

    flog.info("Reading in experiment")
    sce.raw <- readRDS(opt$rds_filename)
    sce <- sce.raw
    cpm(sce) <- calculateCPM(sce)
    data <- getmart(rownames(sce))
    data.unique <- data[match(unique(data$ensembl_gene_id), data$ensembl_gene_id),]
    rownames(data.unique)<-data.unique$ensembl_gene_id
    rownames(sce) <- data.unique[rownames(sce),]$entrezgene
    sc_experiment <- log2(SingleCellExperiment::cpm(sce) + 1)
    marker_genes2 <- marker_genes[marker_genes %in% rownames(sc_experiment)]

    flog.info("Mapping experiment to Allen Data:")
    mapping_results <- map_sampling(train.dat    = norm.dat,
                                    train.cl     = ref.cl,
                                    test.dat     = sc_experiment,
                                    markers      = marker_genes2,
                                    markers.perc = 0.8,
                                    iter         = 100)
    
    colnames(mapping_results$map.freq) <- ref.cl.df[colnames(mapping_results$map.freq),]$cluster
    write.table(as.data.frame.matrix(mapping_results$map.freq), 
                "full_result.tsv", col.names = NA, sep="\t", quote=FALSE)

    map.df <- as_tibble(mapping_results$map.df,rownames = "rowname") %>%
      plyr::mutate(pred.cl = as.numeric(as.character(pred.cl))) %>%
      left_join(ref.cl.df, by = c("pred.cl" = "cluster_id"))
    write_tsv(map.df, "mapping_result.tsv")
    summary_result <- do.call(bind_rows, lapply(as.list(levels(factor(sce$group))),function(level){table(map.df[grep(level,map.df$rowname),]$cluster)}))
    summary_result <- add_column(summary_result,celltype=levels(factor(sce$group)))
    summary_result[is.na(summary_result)] <- 0
    write_tsv(summary_result,"summary_result.tsv")

    if(length(opt$out_filter
) != 0){
        keep <- map.df[grep(opt$out_filter,map.df$cluster),]$rowname
        file = get_output_filename("sce_filtered_hicat.rds")
        flog.info(paste("Saving data matching filter pattern to file", file))
        saveRDS(sce.raw[,keep], file = file)
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
            "--filter",
            dest = "out_filter",
            type = "character",
            default = "",
            help = paste("filename with input data of counts")
        ),
        make_option(
            "--allen-design",
            dest = "allen_design",
            type = "character",
            default = "",
            help = paste("Path to Allen Atlas design file")
        ),
        make_option(
            "--allen-filter",
            dest = "allen_filter",
            type = "character",
            default = "",
            help = paste("regex for filtering Allen data")
        ),
        make_option(
            "--allen-datamatrix",
            dest = "allen_datamatrix",
            type = "character",
            default = "",
            help = paste("Path to Allen Atlas data matrix file")
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
