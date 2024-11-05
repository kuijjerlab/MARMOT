# MARMOT

## Introduction

There have been many tools for joint dimensionality reduction (JDR) of multi-omics data, many of which have been benchmarked^1^, however the suitability of various tools seems to be largely dependent on the data and the downstream analysis being performed.

Here we present MARMOT (<u>M</u>odel comp<u>AR</u>ison for <u>M</u>ulti-<u>O</u>mics <u>T</u>ools), an R tool for comparing JDR models in different conditions or with different inputs.

## Installation guide
MARMOT is available as a GitHub package and can be installed using `devtools`. 

```
install.packages("devtools") # if not already installed
library(devtools)
install_github("rtpop/MARMOT")
```

## References
1. Cantini, L., Zakeri, P., Hernandez, C. et al. Benchmarking joint multi-omics dimensionality reduction approaches for the study of cancer. Nat Commun 12, 124 (2021). https://doi.org/10.1038/s41467-020-20430-7
