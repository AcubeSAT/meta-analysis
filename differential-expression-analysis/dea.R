# Import order matters for masking, take care.
library(dplyr)
library(affy)
library(AnnotationDbi)
library(arrayQualityMetrics)
library(here)
library(yeast2.db)
library(yeast2probe)


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
# Filter out probes that do not map to a gene.
# yeast2.db doesn't have SYMBOL.
# However, ENSEMBL also maps to ENTREZ.
annotated_data <- subset(annotated_data, !is.na(ENSEMBL))

# Grab transcript-cluster identifiers that map to multiple gene identifiers.
annotated_grouped <- group_by(annotated_data, PROBEID)
annotated_summarized <-
    dplyr::summarize(annotated_grouped, no_of_matches = n_distinct(ENSEMBL))
annotated_filtered <- filter(annotated_summarized, no_of_matches > 1)
probe_stats <- annotated_filtered

# Generate an expression set without the probes with multiple mappings.
ids_to_exclude <- (featureNames(eset_rma) %in% probe_stats$PROBEID)
eset_final <- subset(eset_rma, !ids_to_exclude)

# Get the feature data only for the filtered probes.
fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- left_join(fData(eset_final), annotated_data)
# Restore the rownames after performing the left join.
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID