library(affy)
library(arrayQualityMetrics)
library(here)
library(yeast2probe)


# Make sure to launch the R process from the
# project’s top-level directory for here() to work.
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