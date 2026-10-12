#!/usr/bin/env Rscript
packages <- c("data.table", "Matrix", "Seurat", "coloc", "ggplot2", "gridExtra", "png", "scales", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
versions <- data.frame(package = packages, version = vapply(packages, function(x) as.character(packageVersion(x)), character(1)))
print(versions, row.names = FALSE)
