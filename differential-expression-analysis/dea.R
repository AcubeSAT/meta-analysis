library(affy)
library(arrayQualityMetrics)
library(here)
library(gcrma)


raw_data_dir <- here(
    "data",
    "raw",
    "GLDS-62_microarray_E-GEOD-64468.raw.1"
)
qc_data_dir <- here("qc")

cel_affybatch <- ReadAffy(filenames = list.celfiles(
    raw_data_dir,
    full.names = TRUE
))

arrayQualityMetrics(
    expressionset = cel_affybatch,
    outdir = here(qc_data_dir, "qc-report-affybatch"),
    force = TRUE,
    do.logtransform = TRUE
)

eset_gcrma <- gcrma(cel_affybatch, type = "fullmodel")
arrayQualityMetrics(
    expressionset = eset_gcrma,
    outdir = here(qc_data_dir, "qc-report-gcrma"),
    force = TRUE
)