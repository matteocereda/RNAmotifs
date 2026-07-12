# Install R dependencies for rnamotifs-mars pipeline
# Run: Rscript install_mars_deps.R

cran_pkgs <- c("circlize", "viridis", "BAMMtools", "snow", "lsa",
               "optparse", "data.table", "stringr", "rlang", "dplyr",
               "plyr", "tidyr", "lattice", "latticeExtra", "bootstrap",
               "reshape2")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
        cat(pkg, "already installed.\n")
    }
}

# ComplexHeatmap is from Bioconductor
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    cat("Installing ComplexHeatmap from Bioconductor...\n")
    BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
} else {
    cat("ComplexHeatmap already installed.\n")
}

cat("\nAll dependencies checked.\n")
