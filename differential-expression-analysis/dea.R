library(affy)
library(arrayQualityMetrics)
library(here)
library(gcrma)
library(yeast2probe)


# Make sure to launch the R process from the
# project’s top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot.
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

cel_affybatch <- ReadAffy(filenames = list.celfiles(
    raw_data_dir,
    full.names = TRUE
))

s_cerevisiae_mask <- read.table(
    mask_data_dir,
    skip = 2,
    stringsAsFactors = FALSE
)
probe_filter <- s_cerevisiae_mask[[1]]

source(helpers_dir)
s_cerevisiae_df <- extract_ids(probe_filter)

cleancdf <- cleancdfname(cel_affybatch@cdfName)
probe_package_name <- "yeast2probe"
remove_probes(probe_filter, cleancdf, probe_package_name)

arrayQualityMetrics(
    expressionset = cel_affybatch,
    outdir = here(qc_data_dir, "qc-report-affybatch"),
    force = TRUE,
    do.logtransform = TRUE
)

eset_rma <- rma(cel_affybatch)
arrayQualityMetrics(
    expressionset = eset_rma,
    outdir = here(
        qc_data_dir,
        "qc-report-rma"
    ),
    force = TRUE
)

eset_gcrma <- gcrma(cel_affybatch, type = "fullmodel")
arrayQualityMetrics(
    expressionset = eset_gcrma,
    outdir = here(qc_data_dir, "qc-report-gcrma"),
    force = TRUE
)