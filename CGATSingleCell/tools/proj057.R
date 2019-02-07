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



# Downsampling funcion (from Hemberg lab)
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


# Name: normalise_and_plot
# Function: Normalises Single Cell dataset using multiple methods and plots the result
# Inputs: single cell dataset, list of endogenous genes (not spike ins, and not mitochondrial	)
# Outputs: png files showing effect of normalisation on PCA
normalise_and_plot <- function(sce.raw, endog_genes, ERCCconc, method="cpm", subdir=NULL, col_scale=NULL) {

    if(length(subdir)){
        dir.create(subdir)
        if(substr(subdir[1], nchar(subdir), nchar(subdir)) != "/") subdir<-paste0(subdir,"/")
    }
    flog.info("... Raw Counts - no normalisation")
    sce <- sce.raw
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"rawcounts"), assay="counts", col_scale=col_scale)
    
    flog.info("... Log count normalisation")
    logcounts(sce) = log(counts(sce)+1)
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"logcounts"),col_scale=col_scale)
    if(method == "log"){
        sce_out <- sce}

    flog.info("... CPM normalisation")
    logcounts(sce) <- log2(calculateCPM(sce) + 1)
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"cpm"),col_scale=col_scale)
    if(method == "cpm"){
        sce_out <- sce}

    flog.info("... Downsampling normalisation")
    logcounts(sce) <- log2(down_sample_matrix(counts(sce)) + 1)
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"downsampling"),col_scale=col_scale)
    if(method == "downsampling"){
        sce_out <- sce}

    if(!length(grep("ER", rownames(sce)))){
        flog.info("... skipping SCRAN normalisation with ERCCs - no ERCCs found")}
    else{
        flog.info("... SCRAN normalisation with ERCCs")
        sce <- sce.raw
        sce <- computeSpikeFactors(sce, general.use=TRUE)
        sce <- normalize(sce)
        check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"ercc"),col_scale=col_scale)
        fit <- trendVar(sce, parametric=TRUE)
        decomp <- decomposeVar(sce, fit)
        top.hvgs <- order(decomp$bio, decreasing=TRUE)
        start_plot("ercc_fit")
        plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
        o <- order(decomp$mean)
        lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
        points(fit$mean, fit$var, col="red", pch=16)
        end_plot()
        if(method == "ercc"){
            sce_out <- sce}
    }

    flog.info("... SCRAN normalisation without ERCCs")
    sce <- sce.raw
    sce <- computeSumFactors(sce)
    sce <- normalize(sce)
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"scran"),col_scale=col_scale)
    alt.fit <- trendVar(sce, use.spikes=FALSE) 
    alt.decomp <- decomposeVar(sce, alt.fit)
    alt.top.hvgs <- order(alt.decomp$bio, decreasing=TRUE)
    start_plot("scran_fit")
    plot(alt.decomp$mean, alt.decomp$total, xlab="Mean log-expression", ylab="Variance")
    alt.o <- order(alt.decomp$mean)
    lines(alt.decomp$mean[alt.o], alt.decomp$tech[alt.o], col="red", lwd=2)
    end_plot()    
    if(method == "scran"){
        sce_out <- sce}

    flog.info("... RUVs")
    options(stringsAsFactors = FALSE)
    qclust <- quickCluster(sce, min.size = 30)
    sce <- computeSumFactors(sce, sizes = 15, clusters = qclust)
    sce <- normalize(sce)
    # Establishing Group Matrix for RUVs
    scIdx <- matrix(-1, ncol = max(table(sce$group)), nrow = length(levels(sce$group)))
    i <- 1
    write_tsv(as.data.frame(colData(sce)),  "colDataRUVs.tsv")
    write_tsv(as.data.frame(scIdx),  "scIdx.tsv")
    for(groupitem in levels(sce$group)){    
        tmp <- which(sce$group == groupitem)        
        scIdx[i, 1:length(tmp)] <- tmp
        i <- i + 1
    }
    cIdx <- rownames(sce)
    # Running RUVs
    ruvs <- RUVs(counts(sce), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
    assay(sce, "logcounts") <- log2(
        t(t(ruvs$normalizedCounts) / 
            colSums(ruvs$normalizedCounts) * 1e6) + 1)
    check_normalisation(sce, endog_genes, ERCCconc, paste0(subdir,"ruvs"),col_scale=col_scale)
    if(method == "ruvs"){
        sce_out <- sce}

    return(sce_out)
}

# Name: check_normalisation
# Function: Checks how well normalisations perform against known ERCC concentrations
# Inputs: single cell dataset, list of endogenous genes (not spike ins, and not mitochondrial	)
# Outputs: png files showing effect of normalisation on PCA
check_normalisation <- function(sce, endog_genes, ERCCconc, plot.name="no_plot_name_provided",
                                colours=NULL, assay="logcounts", col_scale=NULL) {

    #Sizefactorcorrelation
    if(!length(sizeFactors(sce))){
        start_plot(paste0(plot.name,"_sizefactorcorrelation"))
        plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
             xlab="Library size (millions)", ylab="Size factor",
             col=col_scale[sce$group], pch=16)
        legend("bottomright", col=col_scale, pch=16, cex=1.2,
               legend=levels(sce$group))
        end_plot()
    }
    # PCA
    plot_pca <- scater::plotPCA(
        sce[endog_genes, ],
        colour_by = "group",
        exprs_values = assay,
        size_by = "total_features")
    start_plot(paste0(plot.name,"_pca"))
    print(plot_pca)
    end_plot()
    start_plot(paste0(plot.name,"_pca2"))
    print(ggplot(plot_pca$data, aes(X, Y, colour = colour_by)) + geom_point() + theme_classic() +
          scale_color_manual(labels = levels(sce$group), values=col_scale) + 
          labs(colour = "Genotype") + ylab("Principal Component 2") + 
          xlab("Principal Component 1"))
    end_plot()

    # tSNE
    set.seed(12345678)
    if(assay == "counts") joint_tsne <- Rtsne(t(counts(sce)), perplexity = 15)
    else joint_tsne <- Rtsne(t(logcounts(sce)), perplexity = 15)
    cell_type_labels <- factor(c(as.character(colData(sce)$group), as.character(colData(sce)$group)))
    tsnedf <- as.data.frame(joint_tsne$Y)
    tsnedf$group <- sce$group
    start_plot(paste0(plot.name,"_tsne"))
    print(ggplot(tsnedf, aes(V1, V2, colour = group)) + geom_point() + theme_bw() +
          scale_color_manual(labels = levels(sce$group), values=col_scale) +
          labs(colour = "Genotype") + ylab("tSNE Dimension 2") + xlab("tSNE Dimension 1"))
    end_plot()

    if(packageVersion("scater") <= "1.8.4"){
        plotEx <- plotExplanatoryVariables(sce,
        variable = c("total_features", "total_counts"),
        exprs_values = assay)}
    else{                
        plotEx <- plotExplanatoryPCs(
        sce,
        variable = c("total_features", "total_counts"),
        exprs_values = assay,
        npcs_to_plot = 20)
    }
    start_plot(paste0(plot.name,"_exploratory"))
    print(plotEx)
    end_plot()

    #create melted tables for Counts data
    if(!length(grep("ER", rownames(sce)))){
        # skip correlation with ERCCs if no ERCCs
        return()
    }

    if(assay == "counts"){
        ERCC <- as_tibble(counts(sce[rownames(sce)[grep("ER", rownames(sce))],]),rownames = "ID")
    }
    else {
        ERCC <- as_tibble(logcounts(sce[rownames(sce)[grep("ER", rownames(sce))],]),rownames = "ID")
    }
    ERCCplot <- inner_join(ERCC, ERCCconc[,1:2], by = "ID")
    dfplot <- melt(ERCCplot,id.vars = c("ID","Mix1"))
    dfplot_nonzero <- subset(dfplot, dfplot$value != 0)
    dfplot_big <- subset(dfplot, log(dfplot$Mix1) >6.25)

    #Plotting all counts
    fit <- summary(lm(log(dfplot$Mix1)~dfplot$value))
    p0 <- ggplot(dfplot,aes(log(Mix1),value)) +  geom_point() +
          geom_smooth(method = lm) + 
          ggtitle(label = plot.name, paste0("R^2=",as.character(fit[9]))) + 
          ylab(plot.name) + theme_classic()
    start_plot(paste0("ERCC_vs_", plot.name, "_alldata"))
    print(p0)
    end_plot()

    # Plotting only high abundance counts
    fit1 <- summary(lm(log(dfplot_big$Mix1)~dfplot_big$value))
    p1 <- ggplot(dfplot_big,aes(log(Mix1),value)) +  geom_point() + geom_smooth(method = lm) + 
          ggtitle(label = paste0(plot.name, " - high counts only"), paste0("R^2=",as.character(fit1[9]))) + ylab(plot.name) + theme_classic()
    start_plot(paste0("ERCC_vs_", plot.name, "_highonly"))
    print(p1)
    end_plot()
}


# Name: mergeAllen
# Function: Merges 2018 Allen Dataset
# Inputs: 
# Outputs: 
mergeAllen <- function(list_sces, ERCCconc) {

    sce_counts <- counts(sce)
    design_sce <- colData(sce)[,c("sample_id","group")]
    design_allen.filtered_merge <- design_allen.filtered[,c("sample_name", "cluster")]
    colnames(design_allen.filtered_merge) <- c("sample_id", "group")
    merged_counts <- merge(counts_table3, sce_counts, by="row.names")
    design_merged <- rbind(design_sce, design_allen.filtered_merge)
    rownames(design_merged) <- design_merged$sample_id
    head(merged_counts[,1:10])
    rownames(merged_counts) <- merged_counts$Row.names
    merged_counts <- merged_counts[,rownames(design_merged)]

    table(colnames(as.array.Array(merged_counts)) == rownames(design_merged))

    sce_merged  <- SingleCellExperiment(
        assays = list(counts = as.array.Array(merged_counts)), 
        colData = design_merged
    )

    sce_merged$group <- as.factor(sce_merged$group)
    sce_merged <- calculateQCMetrics(sce_merged)
    endog_genes <- !rowData(sce_merged)$is_feature_control
    sce_merged.raw <-sce_merged
    sum(rownames(sce) %in% rownames(sce_allen))/nrow(sce)
    sum(rownames(sce_allen) %in% rownames(sce))/nrow(sce_allen)

    sce_reduced <- sce[rownames(sce) %in% rownames(sce_allen),]
    sce_allen_reduced <- sce_allen[rownames(sce_allen) %in% rownames(sce),]
    sce_reduced <- sce_reduced[order(rownames(sce_reduced)),]
    sce_allen_reduced <- sce_allen_reduced[order(rownames(sce_allen_reduced)),]


    combined_logcounts <- cbind(logcounts(sce_reduced), logcounts(sce_allen_reduced))
    dataset_labels <- rep(c("ta1", "allen"), times=c(ncol(sce_reduced), ncol(sce_allen_reduced)))
    pheno <- data.frame(Sample_ID = colnames(combined_logcounts),
                    Study_ID=dataset_labels,
                    Celltype=paste(cell_type_labels, dataset_labels, sep="-"))


    var.genes = get_variable_genes(combined_logcounts, pheno)

    corrected <- mnnCorrect(logcounts(sce_reduced), logcounts(sce_allen_reduced), subset.row=var.genes, k=50, sigma=1, pc.approx=TRUE, svd.dim=3 )

    dim(corrected$pairs[[1]]) # sce_reduced
    dim(corrected$pairs[[2]]) # sce_allen_reduced
    head(corrected$pairs[[2]])
    total_pairs <- nrow(corrected$pairs[[2]])
    n_unique_sce_allen_reduced <- length(unique((corrected$pairs[[2]][,1])))
    n_unique_sce_reduced <- length(unique((corrected$pairs[[2]][,1])))
    joint_expression_matrix <- cbind(corrected$corrected[[1]], corrected$corrected[[2]])

    set.seed(345873945)
    joint_tsne <- Rtsne(t(joint_expression_matrix), perplexity = 50)

    cell_type_labels <- factor(c(as.character(colData(sce)$group), as.character(colData(sce_allen)$group)))
    tsnedf <- as.data.frame(joint_tsne$Y)
    tsnedf$group <- cell_type_labels

    ggplot(tsnedf, aes(V1, V2, colour = group)) + geom_point() + 
      scale_color_manual(values=c("azure4", "black", "lightgrey", "red", "green")) +
      ylab("Joint tSNE Dimension 2") + xlab("Joint tSNE Dimension 1")
      
      png('TSNE_mmn_merged.png', width = 6, height = 4, units = 'in', res = 300)
    ggplot(tsnedf, aes(V1, V2, colour = group)) + geom_point() + 
      scale_color_manual(values=c("azure4", "black", "lightgrey", "red", "green")) +
      ylab("Joint tSNE Dimension 2") + xlab("Joint tSNE Dimension 1")
    dev.off()


    require("Seurat")
    set.seed(1234567)

    sce_reduced_seurat <- CreateSeuratObject(raw.data=assays(sce_reduced)[["counts"]]) # raw counts aren't available for sce_reduced
    sce_reduced_seurat@meta.data[, "dataset"] <- 1
    sce_reduced_seurat@meta.data[, "group"] <- colData(sce_reduced)$group

    sce_allen_reduced_seurat <- CreateSeuratObject(raw.data=assays(sce_allen_reduced)[["counts"]])
    sce_allen_reduced_seurat@meta.data[, "dataset"] <- 2
    sce_allen_reduced_seurat@meta.data[, "group"] <- colData(sce_allen_reduced)$group

    sce_reduced_seurat <- NormalizeData(object=sce_reduced_seurat)
    sce_reduced_seurat <- ScaleData(object=sce_reduced_seurat)
    sce_reduced_seurat <- FindVariableGenes(object=sce_reduced_seurat, do.plot=TRUE)

    sce_allen_reduced_seurat <- NormalizeData(object=sce_allen_reduced_seurat)
    sce_allen_reduced_seurat <- ScaleData(object=sce_allen_reduced_seurat)
    sce_allen_reduced_seurat <- FindVariableGenes(object=sce_allen_reduced_seurat, do.plot=TRUE)

    gene.use <- union(rownames(x = head(x = sce_reduced_seurat@hvg.info, n = 2000)),
                      rownames(x = head(x = sce_allen_reduced_seurat@hvg.info, n = 2000)))


    merged_seurat <- RunCCA(object=sce_reduced_seurat, object2=sce_allen_reduced_seurat, genes.use=gene.use, add.cell.id1="m", add.cell.id2="s", num.cc = 5)
    DimPlot(object = merged_seurat, reduction.use = "cca", group.by = "dataset", pt.size = 0.5) # Before correcting 


    merged_seurat <- CalcVarExpRatio(object = merged_seurat, reduction.type = "pca", grouping.var = "dataset", dims.use = 1:5)
    merged.all <- merged_seurat
    merged_seurat <- SubsetData(object=merged_seurat, subset.name="var.ratio.pca", accept.low = 0.5) # CCA > 1/2 as good as PCA
    merged.discard <- SubsetData(object=merged.all, subset.name="var.ratio.pca", accept.high = 0.5)

    summary(factor(merged.discard@meta.data$celltype)) # check the cell-type of the discarded cells.

    merged_seurat <- AlignSubspace(object = merged_seurat, reduction.type = "cca", grouping.var = "dataset", dims.align = 1:5)
    DimPlot(object = merged_seurat, reduction.use = "cca.aligned", group.by = "group", pt.size = 1, cols.use =c(" azure4", "black", "lightgrey","red", "green")) # After aligning subspaces

    png('TSNE_seurat_merged.png', width = 6, height = 6, units = 'in', res = 300)
    DimPlot(object = merged_seurat, reduction.use = "cca.aligned", group.by = "group", pt.size = 1, cols.use = c(" azure4", "black", "lightgrey","red", "green")) # After aligning subspaces

}



run <- function(opt) {

    options(stringsAsFactors = FALSE)
    set.seed(1234567)

    flog.info("Reading in experiment")
    sce <- read_single_cell_experiment_from_rds(opt$rds_filename)
    sce$group = as.factor(sce$group)
    sce.raw <- sce
    endog_genes <- !rowData(sce)$is_feature_control

    flog.info("Reading in ERCC concentrations")
    ERCCconc <-  read_tsv(opt$ERCCpath,col_names = TRUE )
    ERCCconc <- ERCCconc[,c(2,4,5)]
    colnames(ERCCconc) <- c("ID", "Mix1", "Mix2")
    ERCCconc$ID<- gsub("-",".", ERCCconc$ID)

    # Check for mouse cell cycle markers
    flog.info("calculating cell cycle markers")
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    #assigned <- cyclone(sce, pairs=mm.pairs)
    #write.table(table(assigned$phases), "cell_cycle_markers.tsv", 
    #            quote = FALSE, sep = "\t", row.names = FALSE)

    flog.info("Performing normalisations ...")
    sce <-  normalise_and_plot(sce.raw, endog_genes, ERCCconc, 
            method=opt$normalisation, subdir="normalisation_experiment",
            col_scale=c("red","green"))

    flog.info("Plotting SC3")
    rowData(sce)$feature_symbol = rownames(sce)
    sce3 <- sc3(sce, ks = 10, biology = TRUE)
    start_plot("SC3")
    sc3_plot_consensus(sce3, k = 10, show_pdata = "group")
    end_plot()

    flog.info("Importing Allen Atlas Data")
    design_allen <- read_csv(opt$allen_design)
    if(opt$allen_filter == ""){
        design_allen.filtered <- design_allen}
    else{
        design_allen.filtered <- design_allen[grep(opt$allen_filter,design_allen$cluster),]}
    counts_table <- read.csv(opt$allen_datamatrix, row.names = 1 )
    row_table <- read_csv(opt$allen_rowdata)
    counts_table2 <- counts_table[,design_allen.filtered$sample_name]
    rownames(row_table) <- row_table$gene_entrez_id

    flog.info("Converting Allen Atlas Data to ENSEMBL")
    rownames(counts_table2)<- row_table[rownames(counts_table2),]$gene_symbol
    data <- getmart_symbol(rownames(counts_table2))
    data <- data[match(unique(data$mgi_symbol), data$mgi_symbol),]
    rownames(data) <- data$mgi_symbol
    data <- data[match(unique(data$ensembl_gene_id), data$ensembl_gene_id),]
    counts_table3 <- counts_table2[data$mgi_symbol,]
    rownames(data) <- data$mgi_symbol
    rownames(data[!is.na(data[rownames(counts_table3),]$ensembl_gene_id),])
    counts_table3 <- counts_table3[rownames(data[!is.na(data[rownames(counts_table3),]$ensembl_gene_id),]),]
    rownames(counts_table3) <- data[rownames(counts_table3),]$ensembl_gene_id
    sce_allen  <- SingleCellExperiment(
        assays = list(counts = as.array.Array(counts_table3)), 
        colData = design_allen.filtered
    )
    sce_allen$group <- as.factor(sce_allen$cluster)
    sce_allen <- calculateQCMetrics(sce_allen)
    endog_genes <- !rowData(sce_allen)$is_feature_control

    flog.info("Normalising Allen Data ...")
    sce_allen <-  normalise_and_plot(sce_allen, endog_genes, ERCCconc, 
            method=opt$normalisation, subdir="normalisation_allen",
            col_scale=c("azure4", "black", "lightgrey"))

    flog.info("Merging datasets ...")
    
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
            default = "cpm",
            help = paste("which normalisation to use.")
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
        ),
        make_option(
            "--allen-rowdata",
            dest = "allen_rowdata",
            type = "character",
            default = "",
            help = paste("Path to Allen Atlas row data file")
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
