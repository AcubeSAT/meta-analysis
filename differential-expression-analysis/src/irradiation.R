library(logger)

log_info("Importing libraries...")

library(GEOquery)
library(limma)
library(umap)
# implicit dependencies
library(statmod)

log_info("Loading data from GEO...")
gset <- getGEO("GSE4136", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXX111XXX000"
sml <- strsplit(gsms, split = "")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

log_info("log2 transforming...")
# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
log_c <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
if (log_c) {
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex)
}

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test", "control"))
levels(gs) <- groups
gset$group <- gs

log_info("Creating design matrix...")
design_matrix <- model.matrix(~ group + 0, gset)
colnames(design_matrix) <- levels(gs)

log_info("Fitting linear models...")
# fit linear model
fit <- lmFit(gset, design_matrix)

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep = "-")
log_info("Creating contrast matrix...")
cont_matrix <- makeContrasts(contrasts = cts, levels = design_matrix)

log_info("Re-orientating fitteed model to the set of contrasts...")
fit2 <- contrasts.fit(fit, cont_matrix)

log_info("Computing statistics and metrics by empirical Bayes...")
fit_eb <- eBayes(fit2, robust = TRUE)

log_info("Selecting DE genes...")
tt <- topTable(fit_eb,
    adjust.method = "BH",
    sort.by = "B",
    number = Inf,
    p.value = .05,
    lfc = 2
)

log_info("Running decideTests...")
# Identify the significantly differentially expressed genes for each contrast.
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = .05,
    lfc = 2
)

log_success("The script finished running successfully!")
# ct <- 1
# volcanoplot(fit_eb,
#     coef = ct,
#     main = colnames(fit_eb)[ct],
#     pch = 20,
#     highlight = length(which(results[, ct] != 0)),
#     names = rep("+", nrow(fit_eb))
# )