# Microbiome Tutorial 2022 -- 16S Data Analysis

# Prior to Class:
# Installing R and R Studio
We will be making use of the R programming language which is commonly used across many fields for statistical analysis. To begin, please install and/or update to version of R (4.2.0 by downloading from the following [link](https://ftp.osuosl.org/pub/cran/). R is commonly accessed through an easy-to-use interface called R Studio. This can be downloaded [here](https://www.rstudio.com/products/rstudio/download/). Microbiome-analysis software is frequently updated and is often not back-compatible so please ensure you have the latest version of R. You can find how to update R and R Studio [here](https://bootstrappers.umassmed.edu/bootstrappers-courses/courses/rCourse/Additional_Resources/Updating_R.html#updating-on-mac-and-ubuntu). 

If you are not familiar with R, it is an extremely useful tool for analyzing and visualizing data. I would recommend these resources:

- [Data Carpentry: Data Analysis for Ecologists](https://datacarpentry.org/R-ecology-lesson/)
- [Getting started with R and RStudio](https://www.youtube.com/watch?v=lVKMsaWju8w&t=3s)
- If you are interested in learning more R, I recommend this [book](https://r4ds.had.co.nz/).
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
BiocManager::install("ggtree") #a tool for making phylogenetic trees

install.packages("devtools") #A multi-purpose tool for developing packages
devtools::install_github("jbisanz/qiime2R") #A multi-purpose microbiome import/processing package written by Jordan Bisanz
devtools::install_github("blekhmanlab/biomehorizon") #a tool for looking at microbiome time-series data
```

After you have installed all of the above packages, it is a good idea to try loading them one by one as below. If a package fails to load, read the error message and then try to reinstall. Please contact me (christine.olson@ucsf.edu) if there are issues which you can't resolve.
```
library(tidyverse)
library(dada2)
library(phyloseq)
library(ape)
library(vegan)
library(qiime2R)
library(lme4)
library(biomehorizon)
library(ggtree)
```

```
dir.create("~/Desktop/16S_tutorial") #Note: windows users would use setwd("C:\Users\YourUserName\Desktop\16S_tutorial")
```

Then we can specifically set our working directory to be inside this folder:

```
setwd("~/Desktop/16S_tutorial")
```

# Downloading Data
Please download the "reads.zip" file from this [Google Drive folder](https://drive.google.com/drive/u/1/folders/1yINuUX0RExD1wWC-KERXmcH_uiB5-7XB) and move this file to the "16S_tutorial" folder on your Desktop.

Now that we have downloaded the reads, we can unpack the folder as below:

```
unzip("~/Desktop/16S_tutorial/reads.zip")
```

You will now see that all sequencing files are in the 16S_tutorial folder.


# Metadata
Now that we have the reads, we need to know which samples they belong to! This is called metadata, and a spreadsheet (in tab separated format) is available to download as below. Note that you could just as easily use excel spread sheets; however, there are some special considerations for reading/writing excel files and they love to auto-convert text to dates which can cause issues.


```
download.file("https://raw.githubusercontent.com/christineolson-ucsf/MicrobiomeTutorialData/main/metadata.tsv", "~/Desktop/16S_tutorial/metadata.tsv")
```

# Taxonomy Database
Finally, we will download a database for taxonomic assignment. We will be using the RDP database; however, we will discuss alternate (and perhaps better) choices during the tutorial.

```
download.file("https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz", "~/Desktop/16S_tutorial/rdp_species_assignment_18.fa.gz")
download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz", "~/Desktop/16S_tutorial/rdp_train_set_18.fa.gz")
```
