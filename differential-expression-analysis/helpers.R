# Taken from https://github.com/csgillespie/bmc-microarray
# As mentioned in:
# Gillespie, C. S., Lei, G., Boys,
# R. J., Greenall, A., & Wilkinson, D. J. (2010).
# Analysing time course microarray data using Bioconductor:
# a case study using yeast2 Affymetrix arrays.
# BMC research notes, 3(1), 1-10.

## @knitr extract_ids
extract_ids <- function(probe_filter) {
    # probe_filter: a vector of S. cerevisiae genes
    # Get both S. pombe & S. cerevisiae ids from yeast2GENENAME library.
    require(yeast2.db)
    genenames <- as.list(yeast2GENENAME)
    probes <- names(genenames)

    # Get all transcript ids from yeast2annotation.csv.
    # http://www.affymetrix.com/Auth/analysis/downloads/na24/ivt/Yeast2.na24.annot.csv.zip
    annotation_data_dir <- here(
        "data",
        "annotation",
        "Yeast_2.na24.annot.csv"
    )
    annotations <- read.csv(
        file = annotation_data_dir, header = TRUE,
        stringsAsFactors = FALSE
    )
    transcript_id <- annotations[, 3]
    probeset_id <- annotations[, 1]

    # Reorder the transcript_id to match probes.
    transcript_id <- transcript_id[match(probes, probeset_id)]

    # Retrieve the probeset and transcript ids for S. cerevisiae.
    c_probe_id <- probes[-match(probe_filter, probes)]
    c_transcript_id <- transcript_id[-match(probe_filter, probes)]

    # We need the TranscriptID if the gene name is NA.
    yeast_genenames <- transcript_id
    for (i in seq(along = probeset_id)) {
        gname <- genenames[i][[1]]
        if (!is.na(gname)) {
            yeast_genenames[i] <- gname
        }
    }

    # Set the gene name.
    c_genename <- yeast_genenames[-match(probe_filter, probes)]
    df <- data.frame(
        probe = c_probe_id, transcript = c_transcript_id,
        genename = c_genename, stringsAsFactors = FALSE
    )
    return(df)
}

## @knitr remove_probes
remove_probes <- function(list_out_probe_sets,
                          cdfpackagename,
                          probepackagename) {
    # list_out_probe_sets: Probes sets that are removed.
    # cdfpackagename: The cdf package name.
    # probepackagename: The probe package name.
    require(cdfpackagename, character.only = TRUE)
    require(probepackagename, character.only = TRUE)

    probe_env_org <- get(probepackagename)

    # Remove probesets from the CDF environment.
    rm(list = list_out_probe_sets, envir = get(cdfpackagename))

    # Set the PROBE env accordingly.
    # (idea originally from gcrma compute.affinities.R)
    tmp <- get("xy2indices", paste("package:", cdfpackagename, sep = ""))

    new_affybatch <- new("AffyBatch", cdfName = cleancdf)
    pm_index <- unlist(indexProbes(new_affybatch, "pm"))
    sub_index <- match(tmp(probe_env_org$x, probe_env_org$y,
        cdf = cdfpackagename
    ), pm_index)

    i_na <- which(is.na(sub_index))

    # Need to unlock the environment binding to alter the probes.
    ipos <- grep(probepackagename, search())
    env <- as.environment(search()[ipos])

    unlockBinding(probepackagename, env)
    assign(probepackagename, probe_env_org[-i_na, ], env = env)
    lockBinding(probepackagename, env)
}