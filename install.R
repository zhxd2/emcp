if (!requireNamespace("jsonlite", quietly = TRUE)) {install.packages("jsonlite")}
if (!requireNamespace("tidyverse", quietly = TRUE)) {install.packages("tidyverse")}
if (!requireNamespace("argparser", quietly = TRUE)) {install.packages("argparser")}
if (!requireNamespace("seqinr", quietly = TRUE)) {install.packages("seqinr")}
if (!requireNamespace("formattable", quietly = TRUE)) {install.packages("formattable")}

if (!requireNamespace("AnnotationForge", quietly = TRUE)) {BiocManager::install("AnnotationForge")}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {BiocManager::install("clusterProfiler")}
if (!requireNamespace("pathview", quietly = TRUE)) {BiocManager::install("pathview")}
