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
#### LOAD GO LIBRARIES
library(clusterProfiler)
library(DOSE)

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

egoFunc <- function(genes, universe, ont) {
        # Convenience function to run GO enrichment.
        # Use default values and provide gene lists and ontogony.
        ego <- enrichGO(gene = genes,
                        universe = universe,
                        ont = ont,
                        organism = 'mouse',
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.01,
                        readable = TRUE)
        return(ego)
}


#################################################
##### Load DESeq2 objects and get DE genes. #####
load("DESeq-objects/DDS_peat.rda")

deg_peat <- writeDESeqResults(dds = dds_peat, 
                              num = "Mutant", 
                              denom = "Wild-type", 
                              log2_cut = 0, 
                              return_data = TRUE,
                              filename = "DEG-lists/deg_peat.txt")

#### Get peat DEG Ids for GO enrichment ####
deg_peat_ids <- bitr(rownames(deg_peat), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID",
                     annoDb = "org.Mm.eg.db")
# # 28 peat genes didn't map. Which?
# rownames(deg_peat)[ ! rownames(deg_peat) %in% deg_peat_ids$SYMBOL]  # Who cares.

### Explore number of ribosomal protein genes.
# universes were first coded for scl in GO-tmp.R



##########################
#### Define universes ####
peat_uni0 <- rowSums(counts(dds_peat, normalized = TRUE) > 0)
peat_uni10 <- rowSums(counts(dds_peat, normalized = TRUE) > 10) 
peat_uni100 <- rowSums(counts(dds_peat, normalized = TRUE) > 100) 
# save(peat_uni0, peat_uni10, peat_uni100, 
#      file = "DESeq-objects/peat_universes.Rda")
# Ids from universe with min 10 reads.
uni_peat <- names(peat_uni10[peat_uni10 >= 3])  # 14485 genes
uni_peat_ids <- bitr(uni_peat, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID",
                     annoDb = "org.Mm.eg.db")


# Histogram of universes.
plot_reads0 <- qplot(peat_uni0, binwidth = 1, 
                     fill = I('lightblue'), color = I('grey40'), 
                     ylim = c(0, 15000),
                     xlab = 'Number of samples', 
                     ylab = 'Number of genes',
                     main = "Cut-off 0 reads") + theme_cowplot()
plot_reads10 <- qplot(peat_uni10, binwidth = 1, 
                      fill = I('pink'), color = I('grey40'),
                      ylim = c(0, 15000),
                      xlab = 'Number of samples', 
                      ylab = 'Number of genes',
                      main = "Cut-off 10 reads") + theme_cowplot()
plot_reads100 <- qplot(peat_uni100, binwidth = 1, 
                       fill = I('palegreen2'), color = I('grey40'), 
                       ylim = c(0, 15000),
                       xlab = 'Number of samples', 
                       ylab = 'Number of genes',
                       main = "Cut-off 100 reads") + theme_cowplot()

plot_peat_reads <- plot_grid(plot_reads0, plot_reads10, plot_reads100, ncol = 1,
                             labels = c('A', 'B', 'C'))
# plot_peat_reads  # B looks reasonable.
ggsave(plot_peat_reads, filename = "tmp-plots/peat-readsample-dist.pdf", 
       height = 8, width = 4)
rm(list = ls(pattern = 'plot'))


###########################################
#### Biological Process GO enrichment. ####
ego_peat_bp <- egoFunc(genes = deg_peat_ids$ENTREZID, 
                       universe = uni_peat_ids$ENTREZID, 
                       ont = "BP")
# barplot(ego_peat_bp)  # translation number 1
# dotplot(ego_peat_bp)  # ribosomal not the most prominent
# Get the results in df and genes in translation.
ego_peat_bp_df <- as.data.frame(ego_peat_bp@result)
trns_genes <- ego_peat_bp_df[ego_peat_bp_df$Description == 'translation', 
                             c(1, 2, 8)]
trns_genes <- strsplit(trns_genes$geneID, '/')[[1]]

# Combine patterns for Rps, Rpl, Mrps, Mprl; should find 73 genes.
deg_rp <- deg_peat[grep('(Mr|R)p(s|l)', deg_peat$Gene), ] 

# which translation genes are not in deg_rp and what are their FC?
# Not all ribosomal proteins are annotated so
# merge results and remove duplicates.
deg_trns <- deg_peat[deg_peat$Gene %in% trns_genes, ]
deg_trns_rp <- rbind(deg_trns, deg_rp)
deg_trns_rp <- unique(deg_trns_rp)
qplot(deg_trns_rp$log2FoldChange)

# plot DEG of different classes.
plot_a <- qplot(deg_peat$log2FoldChange, 
                bins = 27,
                xlim = c(-0.6, 1.1),
                ylim = c(0, 125),
                main = "All DEG (700)\nin peat mutants",
                fill = I('lightblue'), color = I('grey40'),
                xlab = "log2(Fold Change)") + theme_cowplot()
plot_b <- qplot(deg_rp$log2FoldChange, 
                bins = 27,
                xlim = c(-0.6, 1.1),
                ylim = c(0, 25),
                main = "Ribosomal DEG (73)\nin peat mutants",
                fill = I('pink'), color = I('grey40'),
                xlab = "log2(Fold Change)") + theme_cowplot()
plot_c <- qplot(deg_trns_rp$log2FoldChange, 
                bins = 27,
                xlim = c(-0.6, 1.1),
                ylim = c(0, 25),
                main = "Rp and translation\nDEG (102) in peat mutants",
                fill = I('palegreen'), color = I('grey40'),
                xlab = "log2(Fold Change)") + theme_cowplot()
plot_peat_hist <- plot_grid(plot_a, plot_b, plot_c,
                            ncol = 1,
                            labels = c('A', 'B', 'C'))
# plot_peat_hist
ggsave(plot_peat_hist, filename = "tmp-plots/peat_DEG-Rp_hist.pdf", 
       height = 8, width = 4)
rm(list = ls(pattern = 'plot'))


###########################################
#### Molecular Function GO enrichment. ####
ego_peat_mf <- egoFunc(genes = deg_peat_ids$ENTREZID, 
                       universe = uni_peat_ids$ENTREZID, 
                       ont = "MF")
# barplot(ego_peat_mf)  # GO:0003735: structural constituent of ribosome
# dotplot(ego_peat_mf)  # junk
# Get the results in df and genes in ribosome
ego_peat_mf_df <- as.data.frame(ego_peat_mf@result)
ribo_genes <- ego_peat_mf_df[ego_peat_mf_df$ID == 'GO:0003735', 
                             c(1, 2, 8)]
ribo_genes <- strsplit(ribo_genes$geneID, '/')[[1]]

# Plot the ribo genes
# Nice but doesn't add anything.
deg_ribo <- deg_peat[deg_peat$Gene %in% ribo_genes, ]
plot_ribo <- qplot(deg_ribo$log2FoldChange, 
                   bins = 27,
                   xlim = c(-0.6, 1.1),
                   ylim = c(0, 15),
                   main = "Ribosomal (MF GO:0003735) (60)\nin peat mutants",
                   fill = I('grey80'), color = I('grey40'),
                   xlab = "log2(Fold Change)") + theme_cowplot()
# plot_ribo
ggsave(plot_ribo, filename = "tmp-plots/peat_DEG-MF-ribo_hist.pdf", 
       height = 3, width = 4)


###########################################
#### Cellular Component GO enrichment. ####
ego_peat_cc <- egoFunc(genes = deg_peat_ids$ENTREZID, 
                       universe = uni_peat_ids$ENTREZID, 
                       ont = "CC")
barplot(ego_peat_cc)  # GO:0030529 ribonucleoprotein complex
dotplot(ego_peat_cc)  # junk
# Get the results in df and genes in ribosome
ego_peat_cc_df <- as.data.frame(ego_peat_cc@result)
ribo_cc_genes <- ego_peat_cc_df[ego_peat_cc_df$ID == 'GO:0030529', 
                                c(1, 2, 8)]
ribo_cc_genes <- strsplit(ribo_cc_genes$geneID, '/')[[1]]

# Plot the CC RNP complex genes
# Nice, and perhaps more comprehensive.
deg_cc_ribo <- deg_peat[deg_peat$Gene %in% ribo_cc_genes, ]
plot_ribo_cc <- qplot(deg_cc_ribo$log2FoldChange, 
                      bins = 27,
                      xlim = c(-0.6, 1.1),
                      ylim = c(0, 25),
                      main = "RNP Complex (CC GO:0030529) (123)\nin peat mutants",
                      fill = I('slategray2'), color = I('grey40'),
                      xlab = "log2(Fold Change)") + theme_cowplot()
plot_ribo_cc
ggsave(plot_ribo_cc, filename = "tmp-plots/peat_DEG-CC-RNP_hist.pdf", 
       height = 3, width = 4)
