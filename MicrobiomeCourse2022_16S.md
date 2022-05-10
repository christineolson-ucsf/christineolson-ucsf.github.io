# Microbiome Tutorial 2022 -- 16S Data Analysis

# Installing R and R Studio
We will be making use of the R programming language which is commonly used across many fields for statistical analysis. To begin, please install and/or update to version of R (4.2.0 by downloading from the following [link](https://ftp.osuosl.org/pub/cran/). R is commonly accessed through an easy-to-use interface called R Studio. This can be downloaded [here](https://www.rstudio.com/products/rstudio/download/). Microbiome-analysis software is frequently updated and is often not back-compatible so please ensure you have the latest version of R. You can find how to update R and R Studio [here](https://bootstrappers.umassmed.edu/bootstrappers-courses/courses/rCourse/Additional_Resources/Updating_R.html#updating-on-mac-and-ubuntu). 

If you are not familiar with R, it is an extremely useful tool for analyzing and visualizing data. I would recommend these videos:

- [Introduction to R Programming with Data Camp](https://www.youtube.com/watch?v=HkNFn6eosaU)
- [Getting started with R and RStudio](https://www.youtube.com/watch?v=lVKMsaWju8w&t=3s)
- If you are interested in learning more R, I recommend this [book](https://www.amazon.com/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=sr_1_9?dchild=1&keywords=tidyverse&qid=1620329413&s=books&sr=1-9).
To ensure you have successfully installed the software, please open R studio and type the following code into your console:

```
R.Version()$version.string
```
The console should print your version of R (4.2.0 if you downloaded as followed instructions above). 

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

After you have installed all of the above packages, it is a good idea to try loading them one by one as below. If a package fails to load, read the error message and then try to reinstall. Please contact me if there are issues which you can't resolve.
```
library(tidyverse)
library(dada2)
library(phyloseq)
library(ape)
library(vegan)
library(qiime2R)
library(lme4)
```

```
dir.create("~/Desktop/16S_tutorial") #Note: windows users would use setwd("C:\Users\YourUserName\Desktop\16S_tutorial")
```

Then we can specifically set our working directory to be inside this folder:

```
setwd("~/Desktop/16S_tutorial")
```

# Downloading Data
Please download the "reads.tar.gz" file from this [Google Drive folder](https://drive.google.com/drive/u/1/folders/1yINuUX0RExD1wWC-KERXmcH_uiB5-7XB) and move this file to the "16S_tutorial" folder on your Desktop.
