# Make sure to launch the R process from the
# project top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot:
# https://github.com/jennybc/here_here#tldr

pval_cutoff <- .05
lfc_cutoff <- .9

suppressWarnings(suppressPackageStartupMessages(library(docopt)))

"DEA script for GLDS-62 GeneLab entry (raw data).

Usage:
  flocculation.R [-q | --no-color]
  flocculation.R (-h | --help)
  flocculation.R --version
  flocculation.R --qc [-r] [-n] [-t] [--plots] [-q | --no-color]
  flocculation.R --plots [-q | --no-color]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --qc          Produce QC reports.
  --plots       Produce DEA plots.

" -> doc
arguments <- docopt(doc, version = "flocculation 0.1")
qc_selected <- any(arguments$r, arguments$n, arguments$t)

# Import order matters for masking, take care.
# To identify conflicts used by ambiguous function names,
# you can use conflicted by r-lib, loading it into session
# and slowly working your way through the script. See more:
# https://github.com/r-lib/conflicted
suppressMessages(
    suppressWarnings(
        suppressPackageStartupMessages({
            library(dplyr)
            library(affy)
            library(AnnotationDbi)
            library(arrayQualityMetrics)
            library(here)
            library(yeast2.db)
            library(yeast2probe)
            library(limma)
            library(EnhancedVolcano)
            library(ggplot2)
            library(yeast2cdf)
            library(biomaRt)
            library(logger)
            library(tibble)
            library(arrow)
            library(topGO)
            # "suggests" (implicit) dependencies; track carefully.
            library(statmod)
            library(xml2)
            library(gplots)
            library(GO.db)
            library(Rgraphviz)
        })
    )
)

if (arguments$q) {
    log_threshold(SUCCESS)
} else {
    if (!arguments$no_color) {
        log_layout(layout_glue_colors)
    }
}

log_info("Setting up paths...")
# Dataset taken from the GeneLab platform entry:
# https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-62/
# Requires user login/registration.
raw_data_dir <- here(
    "differential-expression-analysis",
    "data",
    "flocculation",
    "raw",
    "GLDS-62_microarray_E-GEOD-64468.raw.1"
)
mask_data_dir <- here(
    "differential-expression-analysis",
    "data",
    "flocculation",
    "mask",
    "s_cerevisiae.msk"
)
annotation_data_dir <- here(
    "differential-expression-analysis",
    "data",
    "flocculation",
    "annotation",
    "Yeast_2.na24.annot.csv"
)

qc_data_dir <- here(
    "differential-expression-analysis",
    "qc",
    "flocculation"
)
helpers_dir <- here(
    "differential-expression-analysis",
    "src",
    "helpers.R"
)

results_dir <- here(
    "differential-expression-analysis",
    "results",
    "flocculation"
)
plots_dir <- here(results_dir, "plots")
gsea_dir <- here(results_dir, "GSEA")
tibbles_dir <- here(
    "differential-expression-analysis",
    "tibbles",
    "flocculation"
)

log_info("Reading in CEL files...")
cel_affybatch <- ReadAffy(filenames = list.celfiles(
    raw_data_dir,
    full.names = TRUE
))

log_info("Reading in maskfile...")
# http://www.affymetrix.com/Auth/support/downloads/maskfiles/scerevisiae.zip
s_cerevisiae_mask <- read.table(
    mask_data_dir,
    skip = 2,
    stringsAsFactors = FALSE
)
probe_filter <- s_cerevisiae_mask[[1]]

log_info("Loading helpers file...")
source(helpers_dir)

log_info("Saving S.cerevisiae-only Yeast 2.0 probe sets tibble...")
# Grab a tibble with only the S. cerevisiae data.
s_cerevisiae_df <- extract_ids(annotation_data_dir, probe_filter)
tibble_path <- "cerevisiae-only-probesets.feather"
arrow::write_feather(s_cerevisiae_df, here::here(tibbles_dir, tibble_path))

cleancdf <- cleancdfname(cel_affybatch@cdfName)
probe_package_name <- "yeast2probe"
# Use the mask to edit the affybatch in-place, removing the pombe probe sets.
# This is needed, since the was created using the Affymetrix Yeast Genome 2.0,
# which also contains probe sets for Schizosaccharomyces pombe:
# https://www.affymetrix.com/support/technical/datasheets/yeast2_datasheet.pdf
log_info("Removing Pombe probe sets...")
remove_probes(probe_filter, cleancdf, probe_package_name)

if (arguments$qc) {
    if (!qc_selected || arguments$r) {
        log_info("Generating QC report for raw data...")
        log_warn(capture.output(suppressMessages(arrayQualityMetrics(
            expressionset = cel_affybatch,
            outdir = here(qc_data_dir, "qc-report-affybatch"),
            force = TRUE,
            do.logtransform = TRUE
        ))))
    }
}

log_info("RMA normalizing...")
eset_rma <- rma(cel_affybatch, verbose = FALSE)

# No need for do.logtransform, since the expression measure from affy::rma
# is already being given in log base 2 scale.
if (arguments$qc) {
    if (!qc_selected || arguments$n) {
        log_info("Generating QC report for normalized data...")
        log_warn(capture.output(suppressMessages(arrayQualityMetrics(
            expressionset = eset_rma,
            outdir = here(
                qc_data_dir,
                "qc-report-rma"
            ),
            force = TRUE
        ))))
    }
}

# Note there is no need to filter low-expressed genes.
# There is no median intensity below 4, a common threshold.

# Adapted from maEndtoEnd
# http://bioconductor.org/packages/devel/workflows/html/maEndToEnd.html
# An end to end workflow for differential gene expression
# using Affymetrix microarrays

log_info("Annotating the expressionset...")
annotated_data <- suppressMessages(AnnotationDbi::select(
    yeast2.db,
    keys = featureNames(eset_rma),
    columns = c("PROBEID", "ENSEMBL", "GENENAME"),
    keytype = "PROBEID"
))

log_info("Removing multiple-mapping probe sets...")
# Grab transcript-cluster identifiers that map to multiple gene identifiers.
annotated_mul_mapping <- annotated_data %>%
    group_by(PROBEID) %>%
    summarize(no_of_matches = n_distinct(ENSEMBL)) %>%
    dplyr::filter(no_of_matches > 1)

mul_mapping_ids <- (featureNames(eset_rma) %in% annotated_mul_mapping$PROBEID)
eset_final <- subset(eset_rma, !mul_mapping_ids)

# Get the feature data only for the filtered probes.
fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- suppressMessages(
    left_join(
        fData(eset_final),
        annotated_data
    )
)
# Restore the rownames after performing the left join.
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID

log_info("Removing control probe sets...")
control_affymetrix <- grep("AFFX", featureNames(eset_final))
eset_final <- eset_final[-control_affymetrix, ]

control_reporter_genes <- grep("RPTR", featureNames(eset_final))
eset_final <- eset_final[-control_reporter_genes, ]

# yeast2.db doesn't have SYMBOL.
# However, ENSEMBL also maps to ENTREZ.
# Remove the probes that do not map to an ENSEMBL/ENTREZ ID.
log_info("Removing probe sets without ENSEMBL ID...")
no_ensembl_ids <- is.na(fData(eset_final)$ENSEMBL)
eset_final <- eset_final[!no_ensembl_ids, ]

# Group membership for all samples.
# Ground vs microgravity.
group_membership_ground <- "00011100011000111"
sml <- strsplit(group_membership_ground, split = "")[[1]]

gs <- factor(sml)

groups <- make.names(c("onground", "micro"))
levels(gs) <- groups
eset_final$group <- gs

log_info("Creating design matrix...")
# Create an independent t-test design matrix.
# Linear model equation:
# y = mean(on ground) + mean(micro)
design_matrix <- model.matrix(~ group + 0, eset_final)
colnames(design_matrix) <- levels(gs)

# Fit multiple linear models by generalized least squares.
# Least squares is used instead of robust regression,
# because only three replicates are available.
# Due to the low number of replicates (common in microarray datasets),
# using the robust estimation could remove real variation,
# resulting to more false-positives.
log_info("Fitting linear models...")
fit <- lmFit(eset_final, design_matrix)

# Set up contrasts of interest and recalculate model coefficients.
contrast <- paste(groups[1], groups[2], sep = "-")
log_info("Creating contrast matrix...")
contrast_matrix <- makeContrasts(contrasts = contrast, levels = design_matrix)
# Re-orientate the fitted model from the coefficients of the
# design matrix to the set of contrasts of the original coefficients.
log_info("Re-orientating fitted model to the set of contrasts...")
fit2 <- contrasts.fit(fit, contrast_matrix)

log_info("Computing statistics and metrics by empirical Bayes...")
fit_eb <- eBayes(fit2, robust = TRUE)

log_info("Selecting DE genes...")
# Grab DE genes with a FDR (Benjamini-Hochberg adjusted p-values),
# as well as with a absolute log2 fold-change cutoff.
de_genes <- topTable(fit_eb,
    adjust.method = "BH",
    sort.by = "B",
    p.value = pval_cutoff,
    lfc = lfc_cutoff
)
de_genes <- subset(de_genes,
    select = c(
        "PROBEID",
        "adj.P.Val",
        "P.Value",
        "t",
        "B",
        "logFC",
        "ENSEMBL",
        "GENENAME"
    )
)

log_info("Adding ENTREZ IDs...")
# Select the Affymetrix Yeast Genome 2.0 database from ensembl,
# and submit the query for the ENTREZ IDs.
ensembl <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
de_genes$ENTREZ <- getBM(
    attributes = c("affy_yeast_2", "entrezgene_id"),
    filters = "affy_yeast_2",
    values = de_genes$PROBEID,
    mart = ensembl
)[, 2]

de_genes <- as_tibble(de_genes)
log_info("Saving DE genes tibble...")
tibble_path <- "de-genes.feather"
arrow::write_feather(de_genes, here::here(tibbles_dir, tibble_path))

log_info("Running decideTests...")
# Identify the significantly differentially expressed genes for each contrast.
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = pval_cutoff,
    lfc = lfc_cutoff
)

if (arguments$qc || arguments$plots) {
    # Venn diagram needs TestResults-class.
    results_dT <- results
}
# tibble gets all confused with the TestResults-class...
results <- as.data.frame(results)
# tibble also won't preserve rownames.
# https://tibble.tidyverse.org/reference/rownames.html
results <- rownames_to_column(results, var = "probe_id")
results <- as_tibble(results)

log_info("Saving decideTests results tibble...")
tibble_path <- "decideTests.feather"
arrow::write_feather(results, here::here(tibbles_dir, tibble_path))

full_tt <- topTable(fit_eb,
    adjust.method = "BH",
    sort.by = "B",
    number = Inf
)

full_tt <- as_tibble(full_tt)
log_info("Saving full gene list from LMF tibble...")
tibble_path <- "full-de-genes.feather"
arrow::write_feather(full_tt, here::here(tibbles_dir, tibble_path))

log_info("Generating gene universe...")
all_genes <- full_tt$adj.P.Val
names(all_genes) <- full_tt$PROBEID

log_info("Generating BP ID2GO mapping...")
id_to_go_bp <- get_id_to_go_ontology("BP")

log_info("Generating CC ID2GO mapping...")
id_to_go_cc <- get_id_to_go_ontology("CC")

log_info("Generating MF ID2GO mapping...")
id_to_go_mf <- get_id_to_go_ontology("MF")

log_info("Generating the BP topGOdata object...")
GOdata_bp <- generate_topGOdata("BP", all_genes, id_to_go_bp, pval_cutoff)

log_info("Generating the CC topGOdata object...")
GOdata_cc <- generate_topGOdata("CC", all_genes, id_to_go_cc, pval_cutoff)

log_info("Generating the MF topGOdata object...")
GOdata_mf <- generate_topGOdata("MF", all_genes, id_to_go_mf, pval_cutoff)

# NOT adjusted p-values (topGO.pdf, 3.3 Analysis of Results)
# "The p-values computed by the runTest function are unadjusted for multiple
# testing. We do not advocate against adjusting the p-values of the tested
# groups, however in many cases adjusted p-values might be misleading."
log_info("Running BP Fisher GSEA test...")
fisher_bp <- suppressMessages(
    runTest(GOdata_bp, algorithm = "classic", statistic = "fisher")
)
log_info("Running BP Kolmogorov-Smirnov GSEA tests...")
ks_bp <- suppressMessages(
    runTest(GOdata_bp, algorithm = "classic", statistic = "ks")
)
elim_ks_bp <- suppressMessages(
    runTest(GOdata_bp, algorithm = "elim", statistic = "ks")
)

log_info("Running CC Fisher GSEA test...")
fisher_cc <- suppressMessages(
    runTest(GOdata_cc, algorithm = "classic", statistic = "fisher")
)
log_info("Running CC Kolmogorov-Smirnov GSEA tests...")
ks_cc <- suppressMessages(
    runTest(GOdata_cc, algorithm = "classic", statistic = "ks")
)
elim_ks_cc <- suppressMessages(
    runTest(GOdata_cc, algorithm = "elim", statistic = "ks")
)

log_info("Running MF Fisher GSEA test...")
fisher_mf <- suppressMessages(
    runTest(GOdata_mf, algorithm = "classic", statistic = "fisher")
)
log_info("Running MF Kolmogorov-Smirnov GSEA tests...")
ks_mf <- suppressMessages(
    runTest(GOdata_mf, algorithm = "classic", statistic = "ks")
)
elim_ks_mf <- suppressMessages(
    runTest(GOdata_mf, algorithm = "elim", statistic = "ks")
)

log_info("Gathering BP GSEA results...")
bp_results <- GenTable(GOdata_bp,
    classicFisher = fisher_bp,
    classicKS = ks_bp,
    elimKS = elim_ks_bp,
    orderBy = "elimKS",
    ranksOf = "classicFisher",
    topNodes = 10
)

log_info("Gathering CC GSEA results...")
cc_results <- GenTable(GOdata_cc,
    classicFisher = fisher_cc,
    classicKS = ks_cc,
    elimKS = elim_ks_cc,
    orderBy = "elimKS",
    ranksOf = "classicFisher",
    topNodes = 10
)

log_info("Gathering MF GSEA results...")
mf_results <- GenTable(GOdata_mf,
    classicFisher = fisher_mf,
    classicKS = ks_mf,
    elimKS = elim_ks_mf,
    orderBy = "elimKS",
    ranksOf = "classicFisher",
    topNodes = 10
)

bp_results <- as_tibble(bp_results)
log_info("Saving GSEA BP tibble...")
tibble_path <- "gsea_bp.feather"
arrow::write_feather(bp_results, here::here(tibbles_dir, tibble_path))

cc_results <- as_tibble(cc_results)
log_info("Saving GSEA CC tibble...")
tibble_path <- "gsea_cc.feather"
arrow::write_feather(cc_results, here::here(tibbles_dir, tibble_path))

mf_results <- as_tibble(mf_results)
log_info("Saving GSEA MF tibble...")
tibble_path <- "gsea_mf.feather"
arrow::write_feather(mf_results, here::here(tibbles_dir, tibble_path))

if (arguments$plots) {
    log_info("Generating Fisher BP GSEA graph...")
    plot_gsea(
        GOdata_bp,
        score(fisher_bp),
        5,
        gsea_dir,
        "fisher_bp.pdf"
    )

    log_info("Generating KS BP GSEA graph...")
    plot_gsea(
        GOdata_bp,
        score(ks_bp),
        5,
        gsea_dir,
        "ks_bp.pdf"
    )

    log_info("Generating Elimination KS BP GSEA graph...")
    plot_gsea(
        GOdata_bp,
        score(elim_ks_bp),
        5,
        gsea_dir,
        "elim_ks_bp.pdf"
    )

    log_info("Generating Fisher CC GSEA graph...")
    plot_gsea(
        GOdata_cc,
        score(fisher_cc),
        5,
        gsea_dir,
        "fisher_cc.pdf"
    )

    log_info("Generating KS CC GSEA graph...")
    plot_gsea(
        GOdata_cc,
        score(ks_cc),
        5,
        gsea_dir,
        "ks_cc.pdf"
    )

    log_info("Generating Elimination KS CC GSEA graph...")
    plot_gsea(
        GOdata_cc,
        score(elim_ks_cc),
        5,
        gsea_dir,
        "elim_ks_cc.pdf"
    )

    log_info("Generating Fisher MF GSEA graph...")
    plot_gsea(
        GOdata_mf,
        score(fisher_mf),
        5,
        gsea_dir,
        "fisher_mf.pdf"
    )

    log_info("Generating KS MF GSEA graph...")
    plot_gsea(
        GOdata_mf,
        score(ks_mf),
        5,
        gsea_dir,
        "ks_mf.pdf"
    )

    log_info("Generating Elimination KS MF GSEA graph...")
    plot_gsea(
        GOdata_mf,
        score(elim_ks_mf),
        5,
        gsea_dir,
        "elim_ks_mf.pdf"
    )
}

# https://support.bioconductor.org/p/123219/#123228
log_info("Linking genes to GO IDs...")
gene_to_go_ids <- toTable(yeast2GO2ALLPROBES)

log_info("Linking GO IDs to GO terms...")
go_ids_to_terms <- toTable(GOTERM)
go_ids_to_terms <- go_ids_to_terms[, c("go_id", "Term")]

log_info("Performing over-representation GO analysis...")
go_analysis <- kegga(fit_eb,
    gene.pathway = gene_to_go_ids,
    pathway.names = go_ids_to_terms
)
log_info("Extracting top GO terms...")
top_go_terms <- topKEGG(go_analysis)
colnames(top_go_terms)[1] <- "GO Term"

top_go_terms <- rownames_to_column(top_go_terms, var = "GO ID")
top_go_terms <- as_tibble(top_go_terms)

log_info("Saving kegga/goana GO terms tibble...")
tibble_path <- "kegga-GOTerms.feather"
arrow::write_feather(top_go_terms, here::here(tibbles_dir, tibble_path))

# Histogram of the adjusted p-value distribution
# meant for QC of the test results.
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6164648/ or
# http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

# Note: normal test assumption (limma::eBayes proportion)
# is most genes are not differentially expressed.
if (arguments$qc) {
    if (!qc_selected || arguments$t) {
        log_info("Generating p-value distribution histogram QC...")
        pdf(file = here(qc_data_dir, "p-val-hist.pdf"))
        hist(full_tt$P.Value,
            breaks = "Scott",
            col = "grey",
            border = "white",
            xlab = "p-val",
            ylab = "Number of genes",
            main = "p-value distribution"
        )
        invisible(dev.off())

        log_info("Generating adj. p-value distribution histogram QC...")
        pdf(file = here(qc_data_dir, "adj-p-val-hist.pdf"))
        hist(full_tt$adj.P.Val,
            breaks = "Scott",
            col = "grey",
            border = "white",
            xlab = "p-adj",
            ylab = "Number of genes",
            main = "Adjusted p-value distribution"
        )
        invisible(dev.off())

        good_test_probes <- which(!is.na(fit_eb$F))
        log_info("Generating Q-Q for the mod t-statistics...")
        pdf(file = here(
            qc_data_dir,
            "q-q-mod-t-stat.pdf"
        ))
        qqt(fit_eb$t[good_test_probes],
            fit_eb$df.total[good_test_probes],
            main = "Moderated t-statistics"
        )
        invisible(dev.off())
    }
}

if (arguments$plots) {
    ct <- 1
    log_info("Generating volcanoplot...")
    pdf(file = here(plots_dir, "volcano.pdf"))
    volcanoplot(fit_eb,
        coef = ct,
        main = colnames(fit_eb)[ct],
        pch = 20,
        highlight = length(which(results_dT[, ct] != 0)),
        names = rep("+", nrow(fit_eb))
    )
    invisible(dev.off())

    # Generate MD plot (log2 fold-change vs mean log2 expression).
    # Highlight statistically significant (p-adj < .05) probes.
    log_info("Generating MD plot...")
    pdf(file = here(plots_dir, "MD.pdf"))
    plotMD(fit_eb,
        column = ct,
        status = results_dT[, ct],
        legend = F,
        pch = 20,
        cex = 1
    )
    abline(h = 0)
    invisible(dev.off())

    log_info("Generating Venn diagram...")
    pdf(file = here(plots_dir, "venn.pdf"))
    vennDiagram(results_dT, circle.col = palette())
    invisible(dev.off())

    # Generate heatmap where genes are clustered
    # by relative changes in expression.
    # To cluster together genes with similar DE patterns,
    # the genes are clustered by Pearson correlation,
    # and the log-expression values are mean-corrected by rows for the plot.
    log_info("Generating heatmap...")
    pdf(file = here(plots_dir, "heatmap.pdf"))
    coolmap(eset_final[de_genes$PROBEID],
        labRow = de_genes$GENENAME
    )
    invisible(dev.off())

    log_info("Generating EnhancedVolcano...")
    pdf(file = here(plots_dir, "p-val-volcano.pdf"))
    ev <- EnhancedVolcano(full_tt,
        lab = full_tt$GENENAME,
        x = "logFC",
        y = "P.Value",
        ylab = bquote(~ -Log[10] ~ "  P-value"),
        ylim = c(0, 8),
        selectLab = de_genes$GENENAME,
        title = "On-ground vs Microgravity",
        subtitle = "Differential Expression",
        caption = bquote(~ Log[2] ~
        paste("fold change cutoff, ",
            lfc_cutoff,
            "; p-value cutoff, ",
            pval_cutoff,
            sep = ""
        )),
        pCutoff = pval_cutoff,
        FCcutoff = lfc_cutoff,
        pointSize = 1.2,
        labSize = 5.0,
        colAlpha = .8,
        boxedLabels = TRUE,
        labCol = "black",
        labFace = "bold",
        cutoffLineType = "blank",
        cutoffLineCol = "black",
        cutoffLineWidth = .8,
        hline = c(.05, .005, .0005, .00005),
        hlineCol = c("black"),
        hlineType = c("dotted"),
        gridlines.major = TRUE,
        gridlines.minor = FALSE,
        legendLabels = c(
            "Non-significant",
            "FC",
            "P-value",
            "P-value & FC"
        ),
        legendPosition = "bottom",
        legendLabSize = 14,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = .85,
        colConnectors = "black"
    )
    suppressMessages(ev +
        coord_cartesian(xlim = c(-2.5, 2.5)) +
        scale_x_continuous(
            breaks = seq(-2.5, 2.5, .25)
        ))
    # Force EnhancedVolcano to return a plot object instead of a ggplot one.
    # Idea from https://github.com/kevinblighe/EnhancedVolcano/issues/3
    plot(ev)
    invisible(dev.off())
}

log_success("The script finished running successfully!")