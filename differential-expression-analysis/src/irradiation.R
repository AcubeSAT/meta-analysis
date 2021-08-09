library(GEOquery)
library(limma)
library(umap)
# implicit dependencies
library(statmod)

# load series and platform data from GEO

gset <- getGEO("GSE4136", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL2529", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXX111XXX000"
sml <- strsplit(gsms, split = "")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

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

design <- model.matrix(~ group + 0, gset)
colnames(design) <- levels(gs)

# fit linear model
fit <- lmFit(gset, design)

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep = "-")
cont_matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)

# compute statistics and table of top significant genes
fit_eb <- eBayes(fit2, robust = TRUE)
tt <- topTable(fit_eb,
    adjust.method = "BH",
    sort.by = "B",
    number = Inf,
    p.value = .05,
    lfc = 2
)

# summarize test results as "up", "down" or "not expressed"
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = .05,
    lfc = 2
)

# ct <- 1
# volcanoplot(fit_eb,
#     coef = ct,
#     main = colnames(fit_eb)[ct],
#     pch = 20,
#     highlight = length(which(results[, ct] != 0)),
#     names = rep("+", nrow(fit_eb))
# )