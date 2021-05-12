# Import order matters for masking, take care.
library(dplyr)
library(affy)
library(AnnotationDbi)
library(arrayQualityMetrics)
library(here)
library(yeast2.db)
library(yeast2probe)
library(limma)


# Make sure to launch the R process from the
# project top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot:
# https://github.com/jennybc/here_here#tldr

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

arrayQualityMetrics(
    expressionset = cel_affybatch,
    outdir = here(qc_data_dir, "qc-report-affybatch"),
    force = TRUE,
    do.logtransform = TRUE
)

eset_rma <- rma(cel_affybatch)
# No need for do.logtransform, since the expression measure from affy::rma
# is already being given in log base 2 scale.
arrayQualityMetrics(
    expressionset = eset_rma,
    outdir = here(
        qc_data_dir,
        "qc-report-rma"
    ),
    force = TRUE
)

# Note there is no need to filter low-expressed genes.
# There is no median intensity below 4, a common threshold.
eset_medians <- rowMedians(exprs(eset_rma))
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
# Ground vs microgravity, WT vs FLO1 vs FLO8.
group_membership_ground <- "00011100011000111"
sml <- strsplit(group_membership_ground, split = "")[[1]]

group_membership_wt <- "00000011111222222"
sml2 <- strsplit(group_membership_wt, split = "")[[1]]

# Factorize both grouping variables.
gs <- factor(sml)
gs2 <- factor(sml2)

# Add levels to the factors.
groups <- make.names(c("onground", "micro"))
levels(gs) <- groups
groups2 <- make.names(c("WT", "FLO1", "FLO8"))
levels(gs2) <- groups2

# Create a paired-test design matrix.
# Linear model equation:
# onground * WT mean + micro * offset + difference(FLO8 - FLO1- WT)
paired_design_matrix <- model.matrix(~ gs + gs2)
colnames(paired_design_matrix) <- c(
    "intercept",
    "FLO1vsWT",
    "FLO8vsWT",
    "micro vs ground"
)

# Fit multiple linear models by generalized least squares.
# Least squares is used instead of robust regression,
# because only three replicates are available.
# Due to the low number of replicates (common in microarray datasets),
# using the robust estimation could remove real variation,
# resulting to more false-positives.
data_fit <- lmFit(eset_final, paired_design_matrix)

# From a general linear model fit,
# compute moderated t-statistics, moderated F-statistic and
# log-odds of differential expression by empirical Bayes moderation
# of the standard errors towards a common value.

# Based on: Smyth, G. K. (2004).
# Linear models and empirical bayes methods for
# assessing differential expression in microarray experiments.
# Statistical applications in genetics and molecular biology, 3(1).
data_fit_eb <- eBayes(data_fit, robust = TRUE)

# Generate the volcano plot for the "micro vs ground" linear model coefficient.
volcanoplot(data_fit_eb, coef = 4, highlight = 10)

# Identify the significantly differentially expressed genes for each contrast,
# from the fit object with the p-values etc.
results <- decideTests(data_fit_eb)
summary(results)

# Grab DE genes with a FDR (Benjamini-Hochberg adjusted p-values),
# as well as with a absolute log2 fold-change cutoff.
de_genes <- topTable(data_fit_eb,
    coef = "micro vs ground",
    adjust.method = "BH",
    p.value = .05,
    lfc = .9
)

upregulated_cutoff <- de_genes[de_genes[, "logFC"] > .9, ]
downregulated_cutoff <- de_genes[de_genes[, "logFC"] < -.9, ]

# Grab the DE genes as produced by limma::decideTests()
upregulated_dtests <- subset(results, results[, 4] == 1)
downregulated_dtests <- subset(results, results[, 4] == -1)