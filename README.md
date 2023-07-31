# Transcriptome-Explorer
A comprehensive RNA-Seq data analysis pipeline

This project encompasses a comprehensive pipeline for analyzing RNA-Seq data, focusing on the identification of differentially expressed genes using DESeq2 and Limma and subsequent functional enrichment analysis. The script takes raw count data as input, normalizes it, and applies rigorous statistical methods to detect differentially expressed genes. Subsequently, it leverages various R packages to conduct functional enrichment, aiming to unravel the biological implications of the differential expression. The script also includes visualization tools to produce clear, interpretable plots of the data and results.


# Installation Guide

## Prerequisites

This script requires R (version 3.6 or later) and several R packages. Make sure you have R installed on your system. If you do not, you can download it from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

## Cloning the Repository

To get the script on your local machine, clone the repository from GitHub. Open your terminal and run the following command (replace "your-repository-url" with the URL of your GitHub repository):

```
git clone your-repository-url
```

This will create a copy of the repository in your local machine.

## Installing R Packages

The script requires several R packages. Here's how to install them:

Open R or RStudio, and run the following commands:

```R
packages <- c("DESeq2", "pheatmap", "tidyverse", "cowplot", "Matrix.utils", "edgeR", "dplyr", "magrittr", 
"Matrix", "purrr", "reshape2", "S4Vectors", "tibble", "apeglm", "png", "RColorBrewer", "optparse", "limma", 
"edgeR", "openxlsx", "ggplot2", "ggrepel", "clusterProfiler", "enrichplot", "gprofiler2", "ggplot2", "enrichR", 
"tidyverse", "STRINGdb", "igraph")

install.packages(setdiff(packages, rownames(installed.packages())))
```

The script will automatically install any package in the list that is not already installed.

## Running the Script

To run the script, navigate to the directory containing the script file (`Transcriptome-Explorer.R`) in the bash terminal, then source the script:

```
"Rstript Transcriptome-Explorer.R -h"
```

## Script Options

The script supports the following command-line options:

| Option | Description |
| --- | --- |
| `-c CHARACTER, --counts=CHARACTER` | Counts matrix file |
| `-f CHARACTER, --samplefile=CHARACTER` | Sample file |
| `-b CHARACTER, --background=CHARACTER` | Background gene list file (e.g. all/no/a list of genes) |
| `-t CHARACTER, --tool=CHARACTER` | Tool selection (DESeq2/Limma), default=DESeq2 |
| `-d CHARACTER, --datatype=CHARACTER` | Data type (raw/norm), default=raw |
| `-g CHARACTER, --groupnames=CHARACTER` | Group names (e.g., condition,wt,ko) |
| `-p NUMERIC, --pvalue=NUMERIC` | P-value cutoff, default=0.05 |
| `-q NUMERIC, --fdr=NUMERIC` | FDR cutoff, default=0.05 |
| `-F NUMERIC, --foldchange=NUMERIC` | Fold change cutoff, default=0.58 |
| `-T INTEGER, --top=INTEGER` | Number of top DEGs, default=20 |
| `-D CHARACTER, --formula=CHARACTER` | Design formula, e.g. ~sex + age + treatment |
| `-P CHARACTER, --plottype=CHARACTER` | Plot type (png/pdf/svg) |
| `-r REFERENCE, --reference=REFERENCE` | Reference genome, can be 'mouse' or 'human' |
| `-m MISSING_OPTION, --missing_option=MISSING_OPTION` | |
| `-h, --help` | Show help message and exit |

Please replace the placeholders (e.g., `CHARACTER`, `NUMERIC`, `INTEGER`, `REFERENCE`, `MISSING_OPTION`) with your actual values when running the script.

example:
```
Rscript Transcriptome-Explorer.R -c NormCounts.xlsx -f samplefile.xlsx -b no  -g "condition,KO,WT" -D  "~ age + condition" -P png -F 0.2 -q 0.6 -T 10 -d norm
```

Follow the instructions in the script to input your data and run the analysis.

## Getting Help

If you encounter any issues or have any questions, please check the [GitHub Issues page](https://github.com/Kai6662/Transcriptome-Explorer/issues) for this repository. If you don't see your issue, feel free to open a new one. 

