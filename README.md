SCATE: Single-cell ATAC-seq Signal Extraction and Enhancement
====

## Introductions
Single-cell sequencing assay for transposase-accessible chromatin (scATAC-seq) is
the state-of-the-art technology for analyzing genome-wide regulatory landscape in
single cells. Due to data sparsity and discreteness, analyzing scATAC-seq data is
challenging. Existing computational methods cannot accurately reconstruct
activities of individual cis-regulatory elements (CREs) in individual cells. We
present a new statistical framework, SCATE, that adaptively integrates
information from co-activated CREs, similar cells, and publicly available regulome
data to substantially increase the accuracy for estimating individual CRE
activities in single cell and rare cell subpopulations. 

## SCATE Installation

SCATE software can be installed via Github.
Users should have R installed on their computer before installing SCATE. R version needs to be at least 3.5.x or higher. R can be downloaded here: http://www.r-project.org/.

For Windows users, Rtools is also required to be installed. Rtools can be downloaded here: (https://cloud.r-project.org/bin/windows/Rtools/). For R version 3.5.x, Rtools35.exe is recommended. Use default settings to perform the installation.

For mac users, if there is any problem with installation problem, please try download and install clang-8.0.0.pkg from the following URL: https://cloud.r-project.org/bin/macosx/tools/clang-8.0.0.pkg

To install the latest version of SCATE package via Github, run following commands in R:
```{r }
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicAlignments","preprocessCore"))
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/SCATE")
```
<!---
In rare case where devtools cannot be installed, SCATE software package can also be downloaded from this link:

http://jilab.biostat.jhsph.edu/projects/scate/SCATE_1.0.tar.gz

After downloading the package, open R, type in the following commands:
```{r }
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicAlignments")
```
After the installation process, go to 'Packages & Data' menu, go to 'Package Installer', choose 'Local Source Package' from the pull-down menu the top, click 'Install' button on the bottom, choose the path to the software package 'SCATE_1.0.tar.gz' just downloaded, and click 'open'. The installation process will be completed soon. 
--->

If there is any problem with the installation process, please make sure you have R version at least 3.5.x and you have installed Rtools (Windows users) or clang (mac users). If the problem still occurs, please contact the author (see below)

## User Manual
Check the following page for PDF version of the user manual:
https://github.com/zji90/SCATE/raw/master/inst/doc/SCATE.pdf

Check below link for the R code of the user manual.

https://github.com/zji90/SCATE/blob/master/inst/doc/SCATE.R

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
