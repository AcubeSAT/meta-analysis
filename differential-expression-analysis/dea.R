# Make sure to launch the R process from the
# project top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot:
# https://github.com/jennybc/here_here#tldr

suppressPackageStartupMessages(library(docopt))

"DEA script for GLDS-62 GeneLab entry.

Usage:
  dea.R
  dea.R (-h | --help)
  dea.R --version
  dea.R --qc [-r] [-n] [-t] [--plots]
  dea.R --plots

Options:
  -h --help     Show this screen.
  --version     Show version.
  --qc          Produce QC reports.
  --plots       Produce DEA plots.

" -> doc
arguments <- docopt(doc, version = "DEA 0.1")
qc_selected <- any(arguments$r, arguments$n, arguments$t)

# Import order matters for masking, take care.
suppressPackageStartupMessages({
    library(dplyr)
    library(affy)
    library(AnnotationDbi)
    library(arrayQualityMetrics)
    library(here)
    library(yeast2.db)
    library(yeast2probe)
    library(limma)
})

# Set up related paths.

# Dataset taken from the GeneLab platform entry:
# https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-62/
# Requires user login/registration.
raw_data_dir <- here(
    "data",
    "raw",
    "GLDS-62_microarray_E-GEOD-64468.raw.1"
)
mask_data_dir <- here(
    "data",
    "mask",
    "s_cerevisiae.msk"
)

qc_data_dir <- here("qc")
helpers_dir <- here("helpers.R")

results_dir <- here("results")
plots_dir <- here(results_dir, "plots")

# Read in the .CEL files generated from the Affymetrix.
cel_affybatch <- ReadAffy(filenames = list.celfiles(
    raw_data_dir,
    full.names = TRUE
))

# Read in the S. cerevisiae mask provided by Affymetrix.
# http://www.affymetrix.com/Auth/support/downloads/maskfiles/scerevisiae.zip
s_cerevisiae_mask <- read.table(
    mask_data_dir,
    skip = 2,
    stringsAsFactors = FALSE
)
probe_filter <- s_cerevisiae_mask[[1]]

source(helpers_dir)
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
        arrayQualityMetrics(
            expressionset = cel_affybatch,
            outdir = here(qc_data_dir, "qc-report-affybatch"),
            force = TRUE,
            do.logtransform = TRUE
        )
    }
}

eset_rma <- rma(cel_affybatch)

# No need for do.logtransform, since the expression measure from affy::rma
# is already being given in log base 2 scale.
if (arguments$qc) {
    if (!qc_selected || arguments$n) {
        arrayQualityMetrics(
            expressionset = eset_rma,
            outdir = here(
                qc_data_dir,
                "qc-report-rma"
            ),
            force = TRUE
        )
    }
}

# Note there is no need to filter low-expressed genes.
# There is no median intensity below 4, a common threshold.

# Adapted from maEndtoEnd
# http://bioconductor.org/packages/devel/workflows/html/maEndToEnd.html
# An end to end workflow for differential gene expression
# using Affymetrix microarrays
annotated_data <- select(
    yeast2.db,
    keys = featureNames(eset_rma),
    columns = c("PROBEID", "ENSEMBL", "GENENAME"),
    keytype = "PROBEID"
)

# Grab transcript-cluster identifiers that map to multiple gene identifiers.
annotated_mul_mapping <- annotated_data %>%
    group_by(PROBEID) %>%
    summarize(no_of_matches = n_distinct(ENSEMBL)) %>%
    filter(no_of_matches > 1)

# Generate an expression set without the probes with multiple mappings.
mul_mapping_ids <- (featureNames(eset_rma) %in% annotated_mul_mapping$PROBEID)
eset_final <- subset(eset_rma, !mul_mapping_ids)

# Get the feature data only for the filtered probes.
fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- left_join(fData(eset_final), annotated_data)
# Restore the rownames after performing the left join.
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID

# Remove control probe sets prior to the DEA.
control_affymetrix <- grep("AFFX", featureNames(eset_final))
eset_final <- eset_final[-control_affymetrix, ]

control_reporter_genes <- grep("RPTR", featureNames(eset_final))
eset_final <- eset_final[-control_reporter_genes, ]

# yeast2.db doesn't have SYMBOL.
# However, ENSEMBL also maps to ENTREZ.
# Remove the probes that do not map to an ENSEMBL/ENTREZ ID.
no_ensembl_ids <- is.na(fData(eset_final)$ENSEMBL)
eset_final <- eset_final[!no_ensembl_ids, ]

# Group membership for all samples.
# Ground vs microgravity.
group_membership_ground <- "00011100011000111"
sml <- strsplit(group_membership_ground, split = "")[[1]]

# Factorize.
gs <- factor(sml)

# Add levels to the factors.
groups <- make.names(c("onground", "micro"))
levels(gs) <- groups
eset_final$group <- gs

# Create an independent t-test design matrix.
# Linear model equation:
# y = mean(on ground) + mean(micro)
design_matrix <- model.matrix(~ group + 0, eset_final)
# Overwrite the default generated column names.
colnames(design_matrix) <- levels(gs)

# Fit multiple linear models by generalized least squares.
# Least squares is used instead of robust regression,
# because only three replicates are available.
# Due to the low number of replicates (common in microarray datasets),
# using the robust estimation could remove real variation,
# resulting to more false-positives.
fit <- lmFit(eset_final, design_matrix)

# Set up contrasts of interest and recalculate model coefficients.
contrast <- paste(groups[1], groups[2], sep = "-")
contrast_matrix <- makeContrasts(contrasts = contrast, levels = design_matrix)
# Re-orientate the fitted model from the coefficients of the
# design matrix to the set of contrasts of the original coefficients.
fit2 <- contrasts.fit(fit, contrast_matrix)

# From a general linear model fit,
# compute moderated t-statistics, moderated F-statistic and
# log-odds of differential expression by empirical Bayes moderation
# of the standard errors towards a common value.

# Based on: Smyth, G. K. (2004).
# Linear models and empirical bayes methods for
# assessing differential expression in microarray experiments.
# Statistical applications in genetics and molecular biology, 3(1).
fit_eb <- eBayes(fit2, robust = TRUE)

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

# Identify the significantly differentially expressed genes
# for each contrast from the fit object with the p-values etc.
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = .05,
    lfc = .9
)

# Histogram of the adjusted p-value distribution
# meant for QC of the test results.
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6164648/

# Normal test assumption (limma::eBayes proportion)
# is most genes are not differentially expressed.
if (arguments$qc) {
    if (!qc_selected || arguments$t) {
        qc_tt <- topTable(fit_eb,
            adjust.method = "BH",
            sort.by = "B",
            number = Inf
        )
        pdf(file = here(qc_data_dir, "adj-p-val-hist.pdf"))
        hist(qc_tt$adj.P.Val,
            breaks = "Scott",
            col = "grey",
            border = "white",
            xlab = "P-adj",
            ylab = "Number of genes",
            main = "Adjusted p-value distribution"
        )
        invisible(dev.off())

        # Filter out probes with invalid moderated F-statistics.
        good_test_probes <- which(!is.na(fit_eb$F))
        # Quantile-Quantile plot for the moderated t-statistics.
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
    # Generate volcanoplot with the DE genes as selected from the
    # adjusted p-value and log2 fold-change cutoff.
    ct <- 1
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

    pdf(file = here(plots_dir, "venn.pdf"))
    vennDiagram(results, circle.col = palette())
    invisible(dev.off())
}