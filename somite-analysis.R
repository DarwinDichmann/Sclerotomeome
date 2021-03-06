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
library(cowplot)
library(plyr)
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
                              padj_cut = 0.05,
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
# TODO: Modify function so that it warns if file is already present.

genePlot <- function(gene_2_plot, dds = dds) {
        # TODO: Add function information
        # gene_2_plot <- as.character(gene2plot) # Probably redundant.
        # dds <- dds
        # plotCount() from DESeq2 returns a dataframe of normalized counts.
        gene_counts <- plotCounts(dds, gene_2_plot, returnData = TRUE) 
        
        ## Summarize and calculate meand and variation.
        gene_counts_sum <- ddply(.data = gene_counts, .(condition), summarise, 
                                 N = length(condition), 
                                 meanExp = mean(count), 
                                 sdExp = sd(count), 
                                 seExp = sdExp / sqrt(N))
        
        ## Set the limits for the error-bars
        limits <- aes(ymax = meanExp + seExp, ymin = meanExp - seExp)
        
        ## Create plot with individual samples, mean, and SEM errorbars.
        pl <-ggplot(gene_counts_sum, aes(x = condition, y = meanExp)) 
        pl <- pl + geom_point(size = 2, alpha = 0.5)  # Do I need this?
        pl <- pl + theme_bw() + ggtitle(gene_2_plot)
        pl <- pl + xlab("") 
        pl <- pl + ylab("Normalized counts (mean + SEM)") 
        pl <- pl + geom_errorbar(limits, 
                                 width = 0.25, 
                                 size = 1, 
                                 alpha = 0.5)
        pl <- pl + geom_line(aes(group = "meanExp"), 
                             size = 1, 
                             alpha = 0.5) 
        pl <- pl + geom_point(data = gene_counts, 
                              aes(condition, count), 
                              position = position_jitter(0.1, 0), 
                              color = "dodgerblue2",
                              size = 4, 
                              alpha = 0.7)
        pl  # TODO: I could add a ggsave option or use cowplot.
}


#### Load DESeq2 objects. ####
load("DESeq-objects/DDS_somites.rda")



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
dir.create("DEG-lists")

deg_scl <- writeDESeqResults(dds = dds_som, 
                             num = "Directed", 
                             denom = "Control",
                             return_data = TRUE, 
                             filename = "DEG-lists/deg_scl.txt")

deg_scl_all <- writeDESeqResults(dds = dds_som, 
                             num = "Directed", 
                             denom = "Control", 
                             log2_cut = 0,
                             return_data = TRUE, 
                             filename = "DEG-lists/deg_scl_all.txt")

save(deg_scl, deg_scl_all, file = "DESeq-objects/deg_scl.rda")


#### Scatter plot of DE genes. ####

# Extract mean normalized counts from scl data for scatter plot.
scl_counts <- as.data.frame(counts(dds_som[, 1:8], normalized = TRUE))
scl_means <- summarise(.data = scl_counts,
                       meanCTR = rowMeans(scl_counts[, 1:4]), 
                       meanDIR = rowMeans(scl_counts[, 5:8]))
row.names(scl_means) <- rownames(scl_counts)

# Make subsets with DEG to highlight on plot.
deg_all_means <- subset(scl_means, 
                        row.names(scl_means) %in% row.names(deg_scl_all))
deg_scl_means <- subset(scl_means, 
                        row.names(scl_means) %in% row.names(deg_scl))
# scl_means[c("Pax1", "Pax9", "AI646519"), ]

# Fit a linear model for plot.
lm_fit <- lm(meanDIR ~ meanCTR, scl_means)

# Extract means of interesting genes to plot.
deg_scl_up <- c("Pax1", "Pax9", "AI646519", "Nkx3-1", "Heyl", "Bmp3", "Nbl")
deg_scl_dn <- c("Id1", "Id2", "Id3", "Myf5", "Pax8")
up_means <- subset(deg_scl_means, rownames(deg_scl_means) %in% deg_scl_up )
dn_means <- subset(deg_scl_means, rownames(deg_scl_means) %in% deg_scl_dn)

# Plot.
de_plot <- ggplot(scl_means, aes(log2(meanCTR), log2(meanDIR))) 
de_plot <- de_plot + theme_cowplot()
de_plot <- de_plot + xlab("Control\nlog2(Normalized Counts)")
de_plot <- de_plot + ylab("Directed\nlog2(Normalized counts)")
de_plot <- de_plot + geom_point(colour="grey50", size= 2, alpha= 0.3)
# Add regression line.
de_plot <- de_plot + geom_abline(data = lm_fit, colour = "grey20")
de_plot <- de_plot + annotate("text", colour = "grey20", 
                              label = "slope = 1.095",
                              x = 14, y = 18, size = 4)
# Add all DE genes.
de_plot <- de_plot + geom_point(data = deg_all_means,
                                colour = '#a6cee3',
                                size = 3,
                                alpha = 0.7)

# Add DE +2FC genes
de_plot <- de_plot + geom_point(data = deg_scl_means, 
                                colour = "#1f78b4", 
                                size = 3, 
                                alpha = 0.7)
# Add interesting up genes.
de_plot <- de_plot + geom_point(data = up_means, 
                                colour = "#e31a1c",
                                size = 3)
de_plot <- de_plot + geom_text(data = up_means, 
                               aes(label = rownames(up_means)), 
                               colour = "#e31a1c",
                               size = 4,
                               hjust = 1, vjust = -1)
# Add interesting down genes.
de_plot <- de_plot + geom_point(data = dn_means, 
                                colour= "#ff7f00", 
                                size= 3) 
de_plot <- de_plot + geom_text(data = dn_means, 
                               aes(label = rownames(dn_means) ), 
                               colour = "#ff7f00", 
                               size = 4,
                               hjust = 0, vjust= 1.5)
de_plot
# ggsave(de_plot, file = "Figure-plots/Fig2_wAll.pdf", width = 6, height = 6)
rm(deg_all_means, deg_scl_means, deg_scl_dn, deg_scl_up, dn_means, up_means,
   scl_means, scl_counts, de_plot, lm_fit)  # Do you really need all this?

#### Gene Ontology analysis of scl   #####
# TODO: In other scripts; combine GO-tmp.R and scl_GO.R

#### Write all DEG lists for stages/dissected somites ####

# ##########################################################
# Write DEG lists with default p < 0.05 and any fold change.
# ##########################################################

# Contrast SomA vs SomB
deg_soma_b <- writeDESeqResults(dds_som,
                                num = "SomB",
                                denom = "SomA",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_somab.txt")
deg_soma_b$Contrast <- rep("SomA_vs_SomB", nrow(deg_soma_b))
# Contrast SomA vs SomC
deg_soma_c <- writeDESeqResults(dds_som,
                                num = "SomC",
                                denom = "SomA",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_somac.txt")
deg_soma_c$Contrast <- rep("SomA_vs_SomC", nrow(deg_soma_c))
# Contrast SomA vs SomD
deg_soma_d <- writeDESeqResults(dds_som,
                                num = "SomD",
                                denom = "SomA",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_somad.txt")
deg_soma_d$Contrast <- rep("SomA_vs_SomD", nrow(deg_soma_d))
# Contrast SomB vs SomC
deg_somb_c <- writeDESeqResults(dds_som,
                                num = "SomC",
                                denom = "SomB",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_sombc.txt")
deg_somb_c$Contrast <- rep("SomB_vs_SomC", nrow(deg_somb_c))
# Contrast SomB vs SomD
deg_somb_d <- writeDESeqResults(dds_som,
                                num = "SomD",
                                denom = "SomB",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_sombd.txt")
deg_somb_d$Contrast <- rep("SomB_vs_SomD", nrow(deg_somb_d))
# Contrast SomC vs SomD
deg_somc_d <- writeDESeqResults(dds_som,
                                num = "SomD",
                                denom = "SomC",
                                log2_cut = 0,
                                filename = "DEG-lists/deg_somcd.txt")
deg_somc_d$Contrast <- rep("SomC_vs_SomD", nrow(deg_somc_d))



