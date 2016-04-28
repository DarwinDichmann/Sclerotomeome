### TITLE: Sclerotomeome analysis including second sequencing of peat mutants.
### AUTHOR: Darwin Sorento Dichmann, (C) 2016

#################################################
##########        LOAD PACKAGES        ##########
#################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("vsn") 
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(cluster)
library(vsn)
# TODO: tests for packages and use require() to install only if necessary.

# sessionInfo()

#################################################
#######    EXPERIMENTAL DESIGN TABLES      ######
#################################################

### Use pattern matching to get HT-Seq count files.
count_dir <- file.path("countFiles/")
count_files <- list.files(path = count_dir, 
                          pattern = c("(WT*)|(MUT*)|",
                                      "(SOMA*)|(SOMB*)|(SOMC*)|(SOMD*)|",
                                      "(CTRL*)|(DIR*)|", 
                                      "(HEAD*)|(TRUNK*)"))

### Create experimental design list corresponding to file order.
sample_IDs <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", 
                "DIR1", "DIR2", "DIR3", "DIR4", 
                "HEAD",
                "MUT1", "MUT2", "MUT3",
                "SOMA1", "SOMA3", "SOMA4", # SOMA2 removed.
                "SOMB1", "SOMB2", "SOMB3", "SOMB4", 
                "SOMC1", "SOMC2", "SOMC3", "SOMC4", 
                "SOMD1", "SOMD2", "SOMD3", "SOMD4",
                "TRUNK",
                "WT1", "WT2", "WT3")

exp_groups <- c(rep("Control", 4), 
                rep("Directed", 4),
                "Head",
                rep("Mutant", 3),
                rep("SomA", 3),
                rep("SomB", 4),
                rep("SomC", 4),
                rep("SomD", 4),
                "Trunk",
                rep("Wild-type", 3))

exp_seq <- c(rep("Paired-end", 9),
             rep("Single-read", 3),
             rep("Paired-end", 16),
             rep("Single-read", 3))

tissue_type <- c(rep("Sorted-cells", 8),
                 rep("Dissected-tissue", 4),
                 rep("Sorted-cells", 15),
                 rep("Dissected-tissue", 4))

design_full <- data.frame(sampleName = sample_IDs, 
                         fileName = count_files,
                         lib_type = exp_seq,
                         tissue = tissue_type,
                         condition = exp_groups)

# Separate analysis of peat and somite data need subsets of full design. 
design_peat <- data.frame(droplevels(design_full[c(10:12, 29:31), ]))
design_somites <- data.frame(droplevels(design_full[c(1:9, 13:28), ]))
# design_full; design_somites; design_peat  # OK
rm(count_files, sample_IDs, exp_seq, exp_groups, tissue_type)  # Clean up.


#################################################
##########        DESeq2 Results        #########
#################################################

# DESeq analysis and results are performed here. 
# Resulting objects are saved for loading in other scrips that detail either 
# the somite or peat section of the analysis.

# Full data analysis.
dds_full <- DESeqDataSetFromHTSeqCount(sampleTable = design_full, 
                                      directory = count_dir, 
                                      design = ~ condition)
dds_full <- DESeq(dds_full)
res_full <- results(dds_full)

# Somite data only.
dds_som <- DESeqDataSetFromHTSeqCount(sampleTable = design_somites, 
                                  directory = count_dir, 
                                  design = ~ condition)
dds_som <- DESeq(dds_som)
res_som <- results(dds_som)


# Peat data only.
dds_peat <- DESeqDataSetFromHTSeqCount(sampleTable = design_peat, 
                                      directory = count_dir, 
                                      design = ~ condition)
dds_peat <- DESeq(dds_peat)
res_peat <- results(dds_peat)


#################################################
######        DATA TRANSFORMATIONS         ######
#################################################

# Use regularized log and variance stabilization transformations to 
# determine the lowest SD-mean dependency.

rld_full <- rlog(dds_full, blind = TRUE)
vsd_full <- varianceStabilizingTransformation(dds_full, blind = TRUE)

rld_som <- rlog(dds_som, blind = TRUE)
vsd_som <- varianceStabilizingTransformation(dds_som, blind = TRUE)

rld_peat <- rlog(dds_peat, blind = TRUE)
vsd_peat <- varianceStabilizingTransformation(dds_peat, blind = TRUE)


# Mean-sd plots.

# Full data set.
not_all_zero <- rowSums(counts(dds_full)) > 0  # Exlude genes without counts.
meanSdPlot(log2(counts(dds_full, normalized = TRUE)[not_all_zero, ] + 1))
meanSdPlot(assay(rld_full[not_all_zero, ]))  # Better.
meanSdPlot(assay(vsd_full[not_all_zero, ]))  # Not much different from rld.

# MeanSD plots for somites
not_all_zero <- rowSums(counts(dds_som)) > 0  # Exlude genes without counts.
meanSdPlot(log2(counts(dds_som, normalized = TRUE)[not_all_zero, ] + 1))
meanSdPlot(assay(rld_som[not_all_zero, ]))  # Better than vsd.
meanSdPlot(assay(vsd_som[not_all_zero, ]))  # Not much different from rld.

# MeanSD plots for peat
not_all_zero <- rowSums(counts(dds_peat)) > 0  # Exlude genes without counts.
meanSdPlot(log2(counts(dds_peat, normalized = TRUE)[not_all_zero, ] + 1))
meanSdPlot(assay(rld_peat[not_all_zero, ]))  # Better than vsd.
meanSdPlot(assay(vsd_peat[not_all_zero, ]))  # Not much different from rld.

# TODO: facet plots and save them.
# In all data sets rld transformation performs marginally better than vsd.
# All further analysis will use rld transformed data.
rm(not_all_zero, list = ls(pattern = "vsd"))  # Clean.

# Save data objects for use in proceeding analysis.
dir.create("DESeq-objects")
save(dds_full, res_full, rld_full, file = "DESeq-objects/DDS_full.rda")
save(dds_som, res_som, rld_som, file = "DESeq-objects/DDS_somites.rda")
save(dds_peat, res_peat, rld_peat, file = "DESeq-objects/DDS_peat.rda")


#################################################
######     EXPLORATORY SAMPLE ANALYSIS    #######
#################################################












# TODO: EDA plots
# 