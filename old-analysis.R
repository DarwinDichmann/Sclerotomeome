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


#################################################

### DEFINE COLOR PALETTES
### TODO: Move to relevant sections.
### Heatmap colors
# hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu"))) (255)

#################################################
######        FUNCTION DEFINITIONS        #######
#################################################
writeDESeqResults <- function(dds, num, denom, 
                              condition = "condition", 
                              log2_cut = 1,  # 
                              padj_cut = 0.1,
                              return_data = TRUE,
                              filename = NULL) {
        # Convenience function that extracts results from a DESeqResult 
        # object and returns a dataframe or writes a file if a filename is
        # provided.
        #  
        # Args:
        #       dds:            A DESeqDataSet.
        #       num:            Group used as numerator in contrast.
        #       denom:          Group used denominator in contrast.
        #       condition:      Factor be used for contrast. Default: "condition".
        #       log2_cut:       Log2FC threshold for excluding genes. Default: 1.
        #       padj_cut:       Adjusted p-value threshold. Default: 0.1. 
        #       return_data:    If TRUE, returns results as dataframe.
        #       filename:       Results will be written to file if provided.
        #                       Can be combined with return_data = TRUE.
        condition <- as.character(condition)
        num <- as.character(num)
        denom <- as.character(denom)
        deg <- results(dds, c(condition, num, denom))
        deg <- as.data.frame(deg)
        deg <- subset(deg, padj < padj_cut & abs(log2FoldChange) > log2_cut)
        deg <- deg[order(- deg$log2FoldChange), ]
        deg <- cbind(Row.Names = rownames(deg), deg)
        colnames(deg) [1] <- "Gene"
        if (! is.null(filename)) {
                write.table(x = deg, 
                            file = filename, 
                            quote = FALSE, 
                            sep = "\t")
        }
        if (return_data == TRUE) {
                return(deg)
        }
}


#################################################
#######    EXPERIMENTAL DESIGN TABLE      #######
#################################################

### Path to HT-Seq count files.
count_dir <- file.path("countFiles/")
count_files <- list.files(path = count_dir, 
                          pattern = "(WT*)|(MUT*)|(SOMA*)|(SOMB*)|(SOMC*)|(SOMD*)|(CTRL*)|(DIR*)|(HEAD*)|(TRUNK*)" ) # No space between patterns

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

exp_design <- data.frame(sampleName = sample_IDs, 
                         fileName = count_files,
                         lib_type = exp_seq,
                         tissue = tissue_type,
                         condition = exp_groups)

rm(count_files, sample_IDs, exp_seq, exp_groups, tissue_type)  # Clean up.
# exp_design # OK


#################################################
##########           DESeq2            ##########
#################################################
dds <- DESeqDataSetFromHTSeqCount(sampleTable = exp_design, 
                                  directory = count_dir, 
                                  design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)


#################################################
######        DATA TRANSFORMATIONS         ######
#################################################
### Regularized log transformation
rld <- rlog(dds, blind = TRUE)

### Variance stabilized transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)


#################################################
######     EXPLORATORY SAMPLE ANALYSIS    #######
#################################################
### Determine best transformation
not_all_zero <- rowSums(counts(dds)) > 0 
meanSdPlot(log2(counts(dds, normalized = TRUE)[not_all_zero, ] + 1))
meanSdPlot(assay(rld[not_all_zero, ])) # Better
meanSdPlot(assay(vsd[not_all_zero, ]))
# rm(not_all_zero)  # Might need this later
## TODO: facet plots and save them
#################################################


#################################################
### PCA plot
### Color palette: Qualitiative, 10-class paired from colorbrewer2.org
### A reference in a journal article like this:
### Brewer, Cynthia A., 200x. http://www.ColorBrewer2.org, accessed 3/27/2016.
dir.create("EDA-plots/")
pc_data <- plotPCA(rld, intgroup = 'condition', returnData = TRUE)
# PC1:49% variance, PC2: 23% variance
pc_data$condition <- factor(x = pc_data$condition, 
                            levels = c("Control", "Directed",
                                       "SomA", "SomB", "SomC", "SomD",
                                       "Wild-type", "Mutant",
                                       "Head", "Trunk"))

# Color palette: Qualitiative, 10-class paired from colorbrewer2.org
# A reference in a journal article like this:
# Brewer, Cynthia A., 200x. http://www.ColorBrewer2.org, accessed 3/27/2016.
pc_col <- c('#a6cee3','#1f78b4',
            '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', 
            '#b2df8a', '#33a02c', 
            '#cab2d6', '#6a3d9a')

pc_plot <- ggplot(data = pc_data, aes(pc_data[, 1], pc_data[, 2])) 
pc_plot <- pc_plot + theme_bw()
pc_plot <- pc_plot + geom_point(size = 4, alpha = 0.7, aes(colour = condition))
pc_plot <- pc_plot + scale_colour_manual(values = pc_col, guide_legend( "" )) 
pc_plot <- pc_plot + ggtitle("Principal Component Analysis\nof Somite Data")
pc_plot <- pc_plot + labs(x = "PC1: 49% variance", y = "PC2: 23% variance")
pc_plot <- pc_plot + theme(legend.position = c(0.85, 0.2))
# pc_plot <- pc_plot + geom_text(aes(label = name, color = group))  # Too messy.
pc_plot
ggsave("EDA_plots/PCA.pdf", width = 8, height = 8 )
rm(pc_plot, pc_data) # Clean up.
# TODO: Fix title

#################################################
#########       OTHER EDA PLOTS        ##########
#################################################
plotDispEsts(dds) # OK

### Correlation matrix plots
rld_dist <- as.matrix( dist( t( assay( rld ) ) ) )
distcol <- colorRampPalette( rev( brewer.pal(9, "GnBu") ) ) ( 20 )
heatmap.2 (rld_dist, trace= "none", col= distcol, Colv= FALSE, Rowv= FALSE, dendrogram= "none" ,
           main= "rld",
           key.title= NA, key.ylab=NA, key.ytick=NA, density.info="none",
           key.xlab= "Difference",
           keysize= .8 )


#### OK UNTIL HERE 3/27/2016







### GET DESeq2 RESULTS

# Use writeDESeqResults() instead.
# # Create convenience function for getting DESeq2 results
# getDDSResults <- function( dds, 
#                               condition, numerator, denominator, 
#                               log2FC = 1, padj_cutoff = 0.1 ) {
#         condition <- as.character(condition)
#         numerator <- as.character(numerator)
#         denominator <- as.character(denominator)
#         x <- results(dds, c(condition, numerator, denominator))
#         x <- as.data.frame( x )
#         x <- subset(x, padj < padj_cutoff & abs( log2FoldChange ) > log2FC )
#         x <- x[ order( - x$log2FoldChange ), ]
#         x <- cbind( Row.Names = rownames( x ), x )
#         colnames( x ) [ 1 ] <- "Gene"
#         x
# }



# Add contrasts to list of dataframes.
# TODO: Add contrasts to list of dataframes. Run all with apply()?
# resultsList <- list()
# resultsList[[1]] <- getDDSResults(dds = dds, 
#                                      condition = "condition", 
#                                      numerator = "Directed", 
#                                      denominator = "Control", 
#                                      log2FC = 1, padj_cutoff = 0.1 )
# 
# resultsList[[2]] <- getDDSResults(dds = dds, 
#                                      condition = "condition", 
#                                      numerator = "Mutant", 
#                                      denominator = "Wild-type", 
#                                      log2FC = 0, padj_cutoff = 0.1)


# Write Scl and peat DEGs to files.
# TODO: Modify function so that it warns if file is already present.
dir.create("DEG-lists")


deg_scl <- writeDESeqResults(dds = dds, num = "Directed", denom = "Control",return_data = TRUE, filename = "DEG-lists/deg_scl.txt")

deg_peat <- writeDESeqResults(dds = dds, num = "Mutant", denom = "Wild-type", log2_cut = 0, return_data = TRUE, filename = "DEG-lists/deg_peat.txt")

deg_peat_2fc <- writeDESeqResults(dds = dds, num = "Mutant", denom = "Wild-type", return_data = TRUE, filename = "DEG-lists/deg_peat_2fc.txt")
