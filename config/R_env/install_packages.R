## This script installs the necessary packages to run the R scripts in the Piol,
## et al. manuscript. This script assumes that R v4.2.0 has already been installed

# To install the packages using Bioconductor's BiocManager::install() function
# with specified versions, you can run this R script in the command line, or
# copy the code below and run within RStudio:

# Example for Windows
# 
# Open Command Prompt or PowerShell.
# 
# Navigate to the directory containing your script:
#   
#   sh
# 
# cd C:\path\to\your\script
# 
# Run the script:
#   
#   sh
# 
# Rscript install_packages.R
# 
# Example for macOS or Linux
# 
# Open Terminal.
# 
# Navigate to the directory containing your script:
#   
#   sh
# 
# cd /path/to/your/script
# 
# Run the script:
#   
#   sh
# 
# Rscript install_packages.R

# Load BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install packages with specific versions
BiocManager::install(c(
  "AnnotationDbi@1.66.0",
  "Biobase@2.64.0",
  "BiocFileCache@2.12.0",
  "BiocGenerics@0.50.0",
  "BiocIO@1.14.0",
  "BiocManager@1.30.23",
  "BiocParallel@1.38.0",
  "BiocStyle@2.32.1",
  "Biostrings@2.72.0",
  "DBI@1.2.2",
  "DESeq2@1.44.0",
  "DT@0.33",
  "DelayedArray@0.30.0",
  "EnvStats@2.8.1",
  "GGally@2.2.1",
  "GO.db@3.19.1",
  "GenomeInfoDbData@1.2.12",
  "GenomeInfoDb@1.40.1",
  "GenomicAlignments@1.40.0",
  "GenomicRanges@1.56.1",
  "GeoMxWorkflows@1.10.0",
  "GeomxTools@3.8.0",
  "IRanges@2.38.0",
  "KEGGREST@1.44.0",
  "MASS@7.3-61",
  "Matrix@1.7-0",
  "MatrixGenerics@1.16.0",
  "NanoStringNCTools@1.12.0",
  "R6@2.5.1",
  "RColorBrewer@1.1-3",
  "RCurl@1.98-1.16",
  "RSQLite@2.3.7",
  "RSpectra@0.16-2",
  "Rcpp@1.0.12",
  "Rsamtools@2.20.0",
  "Rtsne@0.17",
  "S4Arrays@1.4.0",
  "S4Vectors@0.42.0",
  "SeuratObject@5.0.2",
  "SparseArray@1.4.1",
  "SummarizedExperiment@1.34.0",
  "UCSC.utils@1.0.0",
  "XML@3.99-0.17",
  "XVector@0.44.0",
  "abind@1.4-5",
  "askpass@1.2.0",
  "babelgene@22.9",
  "beeswarm@0.4.0",
  "biomaRt@2.60.0",
  "bit64@4.0.5",
  "bit@4.0.5",
  "bitops@1.0-8",
  "blob@1.2.4",
  "boot@1.3-30",
  "bslib@0.7.0",
  "cachem@1.1.0",
  "cellranger@1.1.0",
  "cli@3.6.2",
  "codetools@0.2-20",
  "colorspace@2.1-0",
  "compiler@4.4.0",
  "cowplot@1.1.3",
  "crayon@1.5.2",
  "crosstalk@1.2.1",
  "curl@5.2.1",
  "data.table@1.15.4",
  "dbplyr@2.5.0",
  "digest@0.6.35",
  "dotCall64@1.1-1",
  "dplyr@1.1.4",
  "evaluate@0.23",
  "fansi@1.0.6",
  "farver@2.1.2",
  "fastmap@1.2.0",
  "fastmatch@1.1-4",
  "fgsea@1.30.0",
  "filelock@1.0.3",
  "future@1.34.0",
  "future.apply@1.11.2",
  "generics@0.1.3",
  "ggbeeswarm@0.7.2",
  "ggforce@0.4.2",
  "ggiraph@0.8.10",
  "ggplot2@3.5.1",
  "ggrepel@0.9.5",
  "ggstats@0.6.0",
  "ggthemes@5.1.0",
  "globals@0.16.3",
  "glue@1.7.0",
  "gridExtra@2.3",
  "grid@4.4.0",
  "gtable@0.3.5",
  "highr@0.11",
  "hms@1.1.3",
  "htmltools@0.5.8.1",
  "htmlwidgets@1.6.4",
  "httpuv@1.6.15",
  "httr2@1.0.1",
  "httr@1.4.7",
  "hypeR@2.2.0",
  "igraph@2.0.3",
  "janitor@2.2.0",
  "jquerylib@0.1.4",
  "jsonlite@1.8.8",
  "kableExtra@1.4.0",
  "knitr@1.46",
  "labeling@0.4.3",
  "later@1.3.2",
  "lattice@0.22-6",
  "lifecycle@1.0.4",
  "listenv@0.9.1",
  "lme4@1.1-35.5",
  "lmerTest@3.1-3",
  "locfit@1.5-9.9",
  "lubridate@1.9.3",
  "magrittr@2.0.3",
  "matrixStats@1.3.0",
  "memoise@2.0.1",
  "mgcv@1.9-1",
  "mime@0.12",
  "minqa@1.2.7",
  "msigdbr@7.5.1",
  "munsell@0.5.1",
  "networkD3@0.4",
  "nlme@3.1-164",
  "nloptr@2.1.1",
  "numDeriv@2016.8-1.1",
  "openssl@2.2.0",
  "openxlsx@4.2.6.1",
  "org.Mm.eg.db@3.19.1",
  "outliers@0.15",
  "parallel@4.4.0",
  "parallelly@1.38.0",
  "pheatmap@1.0.12",
  "pillar@1.9.0",
  "pkgconfig@2.0.3",
  "plyr@1.8.9",
  "plyranges@1.24.0",
  "png@0.1-8",
  "polyclip@1.10-7",
  "prettyunits@1.2.0",
  "progress@1.2.3",
  "progressr@0.14.0",
  "promises@1.3.0",
  "purrr@1.0.2",
  "rappdirs@0.3.3",
  "reactable@0.4.4",
  "readr@2.1.5",
  "readxl@1.4.3",
  "reshape2@1.4.4",
  "restfulr@0.0.15",
  "reticulate@1.38.0",
  "rjson@0.2.21",
  "rlang@1.1.3",
  "rmarkdown@2.27",
  "rstudioapi@0.16.0",
  "rtracklayer@1.64.0",
  "sass@0.4.9",
  "scales@1.3.0",
  "sessioninfo@1.2.2",
  "shiny@1.9.1",
  "snakecase@0.11.1",
  "sp@2.1-4",
  "spam@2.10-0",
  "splines@4.4.0",
  "stringi@1.8.4",
  "stringr@1.5.1",
  "svglite@2.1.3",
  "systemfonts@1.1.0",
  "tibble@3.2.1",
  "tidyr@1.3.1",
  "tidyselect@1.2.1",
  "timechange@0.3.0",
  "tools@4.4.0",
  "units@0.7-2",
  "utf8@1.2.3",
  "uuid@1.1-2",
  "vctrs@1.6.0",
  "viridis@0.6.3",
  "viridisLite@0.4.2",
  "webshot@0.5.4",
  "withr@2.5.1",
  "xfun@0.35",
  "xml2@1.3.4",
  "yaml@2.3.7",
  "zoo@1.9-6"
))
