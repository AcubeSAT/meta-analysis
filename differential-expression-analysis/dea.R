# Make sure to launch the R process from the
# project top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot:
# https://github.com/jennybc/here_here#tldr

suppressWarnings(suppressPackageStartupMessages(library(docopt)))

"DEA script for GLDS-62 GeneLab entry.

Usage:
  dea.R [-q | --no-color]
  dea.R (-h | --help)
  dea.R --version
  dea.R --qc [-r] [-n] [-t] [--plots] [-q | --no-color]
  dea.R --plots [-q | --no-color]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --qc          Produce QC reports.
  --plots       Produce DEA plots.

" -> doc
arguments <- docopt(doc, version = "DEA 0.1")
qc_selected <- any(arguments$r, arguments$n, arguments$t)

# Import order matters for masking, take care.
# To identify conflicts used by ambiguous function names,
# you can use conflicted by r-lib, loading it into session
# and slowly working your way through the script. See more:
# https://github.com/r-lib/conflicted
suppressWarnings(suppressPackageStartupMessages({
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
    # "suggests" dependencies; track carefully.
    library(statmod)
    library(xml2)
    library(gplots)
    library(GO.db)
}))

if(arguments$q) {
    log_threshold(SUCCESS)
} else {
    if(!arguments$no_color) {
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
    "raw",
    "GLDS-62_microarray_E-GEOD-64468.raw.1"
)
mask_data_dir <- here(
    "differential-expression-analysis",
    "data",
    "mask",
    "s_cerevisiae.msk"
)

qc_data_dir <- here("differential-expression-analysis", "qc")
helpers_dir <- here("differential-expression-analysis", "helpers.R")

results_dir <- here("differential-expression-analysis", "results")
plots_dir <- here(results_dir, "plots")

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

log_info("Removing Pombe probe sets...")
# Grab a dataframe with only the S. cerevisiae data.
s_cerevisiae_df <- extract_ids(probe_filter)

cleancdf <- cleancdfname(cel_affybatch@cdfName)
probe_package_name <- "yeast2probe"
# Use the mask to edit the affybatch in-place, removing the pombe probe sets.
# This is needed, since the was created using the Affymetrix Yeast Genome 2.0,
# which also contains probe sets for Schizosaccharomyces pombe:
# https://www.affymetrix.com/support/technical/datasheets/yeast2_datasheet.pdf
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
fData(eset_final) <- suppressMessages(left_join(fData(eset_final), annotated_data))
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
    p.value = .05,
    lfc = .9
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
ensembl = useMart("ensembl",dataset="scerevisiae_gene_ensembl")
de_genes$ENTREZ <- getBM(
    attributes = c("affy_yeast_2", "entrezgene_id"),
    filters = "affy_yeast_2",
    values = de_genes$PROBEID,
    mart = ensembl
)[, 2]

log_info("Running decideTests...")
# Identify the significantly differentially expressed genes for each contrast.
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = .05,
    lfc = .9
)

if (arguments$qc || arguments$plots) {
    full_tt <- topTable(fit_eb,
        adjust.method = "BH",
        sort.by = "B",
        number = Inf
    )
}

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
        highlight = length(which(results[, ct] != 0)),
        names = rep("+", nrow(fit_eb))
    )
    invisible(dev.off())

    # Generate MD plot (log2 fold-change vs mean log2 expression).
    # Highlight statistically significant (p-adj < .05) probes.
    log_info("Generating MD plot...")
    pdf(file = here(plots_dir, "MD.pdf"))
    plotMD(fit_eb,
        column = ct,
        status = results[, ct],
        legend = F,
        pch = 20,
        cex = 1
    )
    abline(h = 0)
    invisible(dev.off())

    log_info("Generating Venn diagram...")
    pdf(file = here(plots_dir, "venn.pdf"))
    vennDiagram(results, circle.col = palette())
    invisible(dev.off())

    # Generate heatmap where genes are clustered by relative changes in expression.
    # To cluster together genes with similar DE patterns,
    # the genes are clustered by Pearson correlation,
    # and the log-expression values are mean-corrected by rows for the plot.
    log_info("Generating heatmap...")
    pdf(file = here(plots_dir, "heatmap.pdf"))
    coolmap(eset_final[rownames(de_genes), ],
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
        "fold change cutoff, 0.9; p-value cutoff, 0.05"),
        pCutoff = .05,
        FCcutoff = .9,
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