# Microbiome Tutorial 2022 -- 16S Data Analysis

# Installing Required Packages
Generally speaking, most packages for microbial ecology exist in 3 places: CRAN, Bioconductor, and Github. CRAN is the main R repository from which packages can be installed using install.packages(). Bioconductor is a repository specializing in bioinformatics. Bioconductor packages are installed via its own package (available from CRAN) called BiocManager which has a function called BiocManager::install(). Finally, the newest versions of packages or in-development packages are often found on github which can be installed using a package called devtools which offers devtools::install_github(). To install the required packages, copy and past the code below line by line into your R session. If asked to compile type yes and if asked to update packages type no.

```{r setup, include=FALSE}
install.packages("tidyverse") #A group of related packages to facilitate a wide variety of useful tasks including plotting
install.packages("vegan") #A commonly used ecology package offering diversity calculations and statistical tests
install.packages("ape") #A package for phylogenetic analysis offering several useful functions for microbial ecology
install.packages("lme4") #A package for advanced statistical testing
install.packages("lmerTest") #A package for advanced statistical testing

install.packages("BiocManager") #The package manager for Bioconductor
BiocManager::install("dada2") #A one-stop-shop for processing amplicon data
BiocManager::install("phyloseq") #A R package for microbiome analysis

install.packages("devtools") #A multi-purpose tool for developing packages
devtools::install_github("jbisanz/qiime2R") #A multi-purpose microbiome import/processing package written by Jordan Bisanz
```

```
dir.create("~/Desktop/16S_tutorial") #Note: windows users would use setwd("C:\Users\YourUserName\Desktop\16S_tutorial")
```
