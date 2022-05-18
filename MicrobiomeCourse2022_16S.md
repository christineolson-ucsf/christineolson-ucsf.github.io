# Microbiome Tutorial 2022 -- 16S Data Analysis

## Prior to Class:
#### Installing R and R Studio
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

#### Installing Required Packages
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

#### Downloading Data
Please download the "reads.zip" file from this [Google Drive folder](https://drive.google.com/drive/u/1/folders/1yINuUX0RExD1wWC-KERXmcH_uiB5-7XB) and move this file to the "16S_tutorial" folder on your Desktop.

Now that we have downloaded the reads, we can unpack the folder as below:

```
unzip("~/Desktop/16S_tutorial/reads.zip")
```

You will now see that all sequencing files are in the 16S_tutorial folder.


#### Metadata
Now that we have the reads, we need to know which samples they belong to! This is called metadata, and a spreadsheet (in tab separated format) is available to download as below. Note that you could just as easily use excel spread sheets; however, there are some special considerations for reading/writing excel files and they love to auto-convert text to dates which can cause issues.


```
download.file("https://raw.githubusercontent.com/christineolson-ucsf/MicrobiomeTutorialData/main/metadata.tsv", "~/Desktop/16S_tutorial/metadata.tsv")
```

#### Taxonomy Database
Finally, we will download a database for taxonomic assignment. We will be using the RDP database; however, we will discuss alternate (and perhaps better) choices during the tutorial.

```
download.file("https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz", "~/Desktop/16S_tutorial/rdp_species_assignment_18.fa.gz")
download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz", "~/Desktop/16S_tutorial/rdp_train_set_18.fa.gz")
```

## Class:

First, set working directory to where the data and metadata are:
```
setwd("~/Desktop/16S_tutorial")
```
Then, load the proper packages:
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

Next we might want to set up some directories for our outputs as below:
```
dir.create("figures")
dir.create("tables")
```
Finally, it is always helpful to keep a record of your work environment. This can easily be accomplished with SessionInfo():

```
sessionInfo()
```

Before going forward, we will read in our metadata table with the help of the read_tsv() function.
```
metadata<-read_tsv("metadata.tsv")
metadata <- metadata[!is.na(metadata$every_sample),]
```
The most standard format for sequencing data is called a FastQ. It is composed of a repeating pattern of 4 lines as below:

@M01869:152:000000000-AMWM8:1:1101:21609:1687 1:N:0:0
TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCATGGCCG...
+
BBCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGC,-...
The first line contains a unique identifier for each read starting with an @ sign. This data which comes from an Illumina MiSeq shows the machine's serial number (M01869), the flow cell ID (AMWM8; unique to the sequencing run), and the read's physical location (1:1101:21609:1687) on the flow cell. The next line stores the actual DNA sequence for the read which is followed by a + sign. The last line of the pattern stores quality information wherein the probability of an incorrect base call is represented by a letter. 


Given that it is hard to quite hard to read the quality Profile by eye, we can instead plot it using the plotQualityProfile() function.
```
plotQualityProfile(metadata$read1[1:3]) + ggtitle("Forward Read Quality (R1)")
```

The forward reads look of high quality...let's look at the reverse reads. 
```
plotQualityProfile(metadata$read2[1:3]) + ggtitle("Reverse Read Quality (R2)")
```

It is common for the reverse read to have lower quality than the forward read. 

We should do some QC which entails trimming back the ends of the reads and removing those which are suspected to be highly error prone.

```
dir.create("filtered_reads")

metadata<-
  metadata %>%
  mutate(read1_filtered=paste0("filtered_", read1)) %>%
  mutate(read2_filtered=paste0("filtered_", read2))
  

trimlog <- filterAndTrim(
                      fwd=metadata$read1,
                      rev=metadata$read2,
                      filt=metadata$read1_filtered,
                      filt.rev=metadata$read2_filtered,
                      compress=TRUE,
                      truncLen=c(240,140), # trim back the reads on the 3' end
                      trimLeft=5, # remove 5 bases from the beginning as these q scores may not be reliable
                      maxEE=2, # no more than 2 expected errors per reads
                      truncQ=2, # trim the read after a q score of 2
                      multithread=TRUE, # use multiple processors, you may need to set this to FALSE on PCs
                      rm.phix=TRUE, # Remove reads mapping to phiX
                      verbose=TRUE
                   )
```
QUESTION: What is phiX and why might we remove it?

Now we can check how many reads have been removed:

```
trimlog %>%
  as.data.frame() %>%
  mutate(Percent_Passed=reads.out/reads.in*100) %>%
  interactive_table()
```

Amplicon sequencing tends to be very noisy. To deal with this, operational taxonomic units (OTUs) have been heavily used for analysis. In this type of analysis, sequences are clustered based on a % identity threshold, often 97%. This has the advantage of negating sequencing errors but causes a lack of resolution. 97% is generally based on the notion that the 16S rRNA of two members of the same species is >=97%; however, this is based on the whole 16S rRNA gene, and not a small fragment as we are commonly analyzing so this logic is not entirely sound.

As a more modern approach, denoising is commonly used to identify sequencing errors and correct them. As for how it works, a simple explanation is that errors are relatively rare and as such can be identified when the error is substitution from a more abundant sequence. For a more complete explanation and alternate approaches see the manuscripts describing Dada2, or any of the other popular denoising algorithms available:

dada2 https://www.nature.com/articles/nmeth.3869
deblur https://msystems.asm.org/content/2/2/e00191-16
unoise https://www.biorxiv.org/content/10.1101/081257v1
Comparison of all three https://peerj.com/articles/5364/

The first step is to learn the error profile for the sequencing run. This is the longest part of the analysis and is based on a sampling of the data.

```
forward_error<-learnErrors(metadata$read1_filtered, multithread = TRUE)
```

```
reverse_error<-learnErrors(metadata$read2_filtered, multithread = TRUE)
```

Next we can actually correct the data:

```
forward_denoised<- dada(metadata$read1_filtered, err=forward_error, multithread=TRUE)
```

```
reverse_denoised<- dada(metadata$read2_filtered, err=reverse_error, multithread=TRUE)
```

Up until this point, we have been treating our forward and reverse reads separately. We can now merge them together to represent the complete V4 sequence using merge pairs.

```
merged_reads <- mergePairs(
                          dadaF=forward_denoised, 
                          derepF=metadata$read1_filtered, 
                          dadaR=reverse_denoised, 
                          derepR=metadata$read2_filtered, 
                          verbose=TRUE
                          )
```

Now we can get can get the output we really want: a table with each sample and the number of times each sequence is observed.
```
merged_table <- makeSequenceTable(merged_reads)
```

When dealing with amplicon data, it is incredibly important to remove chimeric reads. These sequences are generated through a number of mechanisms but result in a read that is derived from multiple original pieces of template DNA. Certain protocols for generating 16S rRNA amplicons, including the methods used to generate this data, are highly prone to chimera formation. For more information on how to prevent their formation, see the manuscript by [Gohl et al. 2016.](https://www.nature.com/articles/nbt.3601)

```
asv_table<-removeBimeraDenovo(merged_table, method="pooled", multithread=TRUE, verbose=TRUE)
```
Question: How many ASVs did we remove? Do they represent a lot of our data?


Frustratingly, some people like to have samples as rows and features (ASVs) as columns, others like the reverse. Different tools have different expectations and you must always be aware which orientation is expected. I, unlike the authors of Dada2, prefer samples as columns. The good news is that we can easily transpose the table back and forth using the t() function.

```
asv_table<-t(asv_table)
```

If you look at the sample names in your table, it is the currently the name of the fastq file rather than the sample ID, we can switch them out as below.
```
asv_table<-asv_table[,basename(metadata$read1_filtered)] #put columns in the same order as the metadata
names(metadata)[names(metadata) == "sample-id"] <- "sample_id" #change sample name column due to issue with dashes in r variable names
colnames(asv_table)<-metadata$sample_id #overwrite the column names
```
You can also see in our table that the ASVs are being stored as their literal sequence. We can swap this out with something a little bit more readable as below while storing a copy of the sequences:
```
sequences<-tibble(ASV=paste0("ASV_",1:nrow(asv_table)), Sequence=rownames(asv_table))
rownames(asv_table)<-sequences$ASV
```

One of the most obvious to describe a community is in terms of its diversity. While we probably all intuitively understand the concept, how might we actually define it quantitatively?

Question: how even is our sequencing depth across samples?

In this next code chunk our goal is to calculate our diversity metrics and then put them in a table with our metadata.

```
diversity_table<-
  data.frame(
    Shannon=diversity(asv_table, index="shannon", MARGIN=2),
    Obs_ASVs=specnumber(asv_table, MARGIN=2)
  ) %>%
  rownames_to_column("sample_id") %>%
  full_join(metadata) %>%
  select(sample_id, every_sample, age, years_since_ldopastartage_T0, sample_timepoint, Shannon, Obs_ASVs, patient_number, gender, age_5y_bins, T0_patientID)
```

Now that we have our table made, we can start making our plot.

```
diversity_table %>%
  pivot_longer(c(Shannon, Obs_ASVs), names_to="Metric", values_to="Diversity") %>%
  ggplot(aes(x=sample_timepoint, y=Diversity, fill=age)) +
  stat_summary(geom="ribbon", alpha=0.5) +
  geom_jitter(width=1, height=0, shape=21) +
  facet_wrap(~Metric, scales="free") +
  xlab("Sample Timepoint") +
  ylab("Alpha Diversity") 
```
Question: Do you see any trends in the data?


```
ggsave("figures/diversity.pdf", height=3, width=5)
```
Now we can do a statistical test accounting for the nature of the data:
```
fit<-glmer.nb(Obs_ASVs~age_5y_bins*years_since_ldopastartage_T0+(1|sample_id), data=diversity_table)
summary(fit)
```
Question: What trends are significant?

Question: What conclusions do you draw from this diversity analysis? 

Question: Did you get an error message about a singular fit? Should this be ignored?

Beta diversity refers to a way of measuring similarity between samples. These are often called distances or dissimilarities. Conceptually it is easy to measure the distance between two physical points, but how do we measure the distance between two compositions? There is an entire universe of metrics used for this purpose, but today we will use a newer metric called Aitchison distance, also called CLR Euclidean. This metric comes from the field of compositional data analysis (CoDA) and deals with a number of shortcomings of more traditional metrics were developed for macro-ecology. Read more [here][CoDA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5695134/).


```
clr_euc<-asv_table %>% make_clr() %>% t() %>% vegdist(method="euclidean")
```

Now we can look at our distances:
```
as.matrix(clr_euc) %>% corner()
```

We have calculated our distances; however, now we need a human-readable way to interpret them. We need a method for reducing the number of dimensions to something that can be plotted in 2 or 3 dimensions. There are a number of methods, of which you may be familiar with including PCA, PCoA, tSNE, NMDS, etc. Today we will use PCoA (principal coordinates analysis) which is fairly versatile and can work with a number of microbiome-relevant distance metrics.

```
pc<-pcoa(clr_euc)
```

We can now visualize the plot as below, with lines connecting each timepoint per patient:
```
pc$vectors %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata) %>%
  arrange(age) %>%
  ggplot(aes(x=Axis.1, y=Axis.2, group=patient_number, label=sample_timepoint, color=patient_number)) +
  scale_colour_gradient2() +
  geom_path(arrow=arrow(type="closed", length=unit(0.1,"cm"))) +
  #geom_text() +
  theme_q2r() +
  xlab(paste("PC1:",round(pc$values$Relative_eig[1]*100,2),"%")) +
  ylab(paste("PC2:",round(pc$values$Relative_eig[2]*100,2),"%"))

```

Based on the beta diversity distances from the starting point, which patients do you think may have received antibiotics? What are other reasons there could be for a change in beta diversity within patient over time?

```
ggsave("figures/pcoa.pdf", height=2, width=3)  
```

Taxonomic Analysis

So we have seen that perhaps patients of different ages have different microbiotas... but we don't actually know what makes them different yet. In this section we will assign the ASVs to known taxa and create some plots to help us conduct some exploratory analysis.


There are multiple choices of database for this tool. For the purposes of today, we will use the Ribosomal Database Project as it is considerably smaller in size and will be fairly fast for our purposes today. I would recommend the SILVA database which is far more inclusive and has up to date taxonomy. There are two steps:

```
taxonomy <- assignTaxonomy(sequences$Sequence, "rdp_train_set_18.fa.gz", multithread=TRUE)
taxonomy <- addSpecies(taxonomy, "rdp_species_assignment_18.fa.gz", allowMultiple = TRUE)

taxonomy<-
  taxonomy %>%
  as.data.frame() %>%
  rownames_to_column("Sequence") %>%
  left_join(sequences) %>%
  select(ASV, everything()) %>%
  select(-Sequence) %>%
  column_to_rownames("ASV")
```

Question: How many ASVs were assigned to a species?

Question: Why would we allow multiple species assignments?

The stacked bar plot is a staple of microbiome manuscripts. We can easily create one as below:
```
taxa_sums<-summarize_taxa(asv_table, taxonomy)

taxa_barplot(features = taxa_sums$Phylum, metadata = metadata, category = "patient_number")
```
```
ggsave("figures/barplot.pdf", height=3, width=8)
```
Question: What conclusion can you draw from this plot?

Question: What are the limitations of this plot?

Question: How might we make a better plot to compare specific phyla?


Depending on the nature of the data, trends can sometimes be spotted in heat maps. Let's create one.
```
taxa_heatmap(
              features = taxa_sums$Genus,
              metadata = metadata,
              category = "patient_number"
            )
```

```
ggsave("figures/heatmap.pdf", height=4, width=6)
```

Question: What conclusion can you draw from this plot?

Question: What are the limitations of this plot?

Ultimately, what people often care most about is individual features/strains. After all, scientists tend to take a reductionist approach and want to identify novel pathogens or new probiotics in microbiome data. While fundamentally, 16S rRNA gene sequencing is a tool designed to look at communities as a whole, many approaches exist for looking at individual features. There are also a wide range of approaches and ethos for these. Generally speaking, it is a good idea to use multivariate statistics such as the beta diversity analysis above before looking for differences in individual taxa.

While there is not necessarily a single correct way to find significant taxa, there is certainly a wrong way: Do not take raw sequence counts and do a series of uncorrected t-tests. Remember that data is a technical sampling (if subsampled) of a technical sampling (which DNA binds to sequencing flow cell), of a technical sampling (which DNA was amplified by PCR primers), of a technical sampling (which portion of a sample was extracted for DNA), of a biological sampling (one arbitrary sample collection from an individual). As such it is important to realize that if you were to analyze replicate samples at every step of the process you would end up with slightly different numbers (but hopefully the same trends). Even the same sequencing library sequenced multiple times will yield slightly different results. Some tools are more false positive prone then others, and some create better visualization than others. 

It is also important to consider the nature of your dataset. In our case, we have a repeated sampling design over time, and as such it would not be kosher to treat every sample as an independent unit of measurement. We will use a linear mixed effects model after applying a centered log ratio transformation which helps normalize the data to be appropriate for this analysis. An alternate approach which is conceptually similar is available in ANCOMII.

As a first step, we can remove features that are very rare as they won't be helpful and every additional test performed reduces our power to detect significant changes when we consider multiple testing corrects.

```
asv_table_filtered<-filter_features(asv_table, minsamples = 3, minreads = 10)
```

Question: Does this seem like an appropriate filter? What happens if you change the filtering criteria?

Now we can generate our linear mixed effects model on a per-feature basis:
```
results<-
asv_table_filtered %>%
  make_clr() %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(!ASV, names_to="sample_id", values_to="CLR") %>%
  left_join(metadata) %>%
  group_by(ASV) %>%
  do(
    lmerTest::lmer(CLR~patient_number*age_5y_bins+(1|patient_number), data=.) %>%
    summary() %>%
    .$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("Term")
  )
```
Now we can filter for organisms which are significantly different between patients aged 70-74 and other patients.

```
results2<-
  results %>%
  filter(Term=="age_5y_bins70_74") %>%
  ungroup() %>%
  mutate(FDR=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  mutate(Significant=if_else(FDR<0.1 & abs(Estimate)>1, "Significant","ns")) %>%
  left_join(taxonomy %>% rownames_to_column("ASV"))
```

Question: How many ASVs are signficantly higher in the samples from patients aged 70-74?

Question: What term would we look at if we were interested in the organisms which change over time differentially between groups?

One way to represent the this data for exploring the number of significant results is with a volcano plot where we plot fold change against P-value as below:

```
results2 %>%
  ggplot(aes(x=Estimate, y=-log10(`Pr(>|t|)`), color=Significant)) +
  geom_point(alpha=0.5, shape=16) +
  xlab("log2(fold change)") +
  ylab("-log10(P-value)") +
  theme_q2r() +
  scale_color_manual(values=c("grey50","indianred")) +
  theme(legend.position="none")
```

Question: How many ASVs are signficantly higher in the age 70-74 group?

Question: Which ASVs do you think are most important?

Question: What would your next analysis or experiment be?
