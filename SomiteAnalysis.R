### TITLE: Scleome analysis including second sequencing of peat mutants.
### AUTHOR: Darwin Sorento Dichmann, (C) 2016

#################################################
##########        LOAD PACKAGES        ##########
#################################################
# source( "http://bioconductor.org/biocLite.R" )
# biocLite( "vsn" ) 
library( vsn )
library( DESeq2 )
library( gplots )
library( RColorBrewer )
library( ggplot2 )
library( cluster )
# TODO: tests for packages and install only if necessary.
#################################################
### DIRECTORY PATHS.

### TODO: Create other dirs here?

### DEFINE COLOR PALETTES
### TODO: Move to relevant sections.
### Heatmap colors
hmcol <- colorRampPalette( rev( brewer.pal( 9, "RdBu" ) ) ) ( 255 )
### Distance matrix colors
distcol <- colorRampPalette( rev( brewer.pal(9, "GnBu") ) ) ( 20 )


#################################################
#######    EXPERIMENTAL DESIGN TABLE      #######
#################################################

### Path to HT-Seq count files.
inDir <- file.path( "countFiles/" )
countFiles <- list.files( path = inDir, 
                          pattern = "(WT*)|(MUT*)|(SOMA*)|(SOMB*)|(SOMC*)|(SOMD*)|(CTRL*)|(DIR*)|(HEAD*)|(TRUNK*)" ) # No space between patterns

### Create experimental design list corresponding to file order.
sampleIDs <- c( "CTRL1", "CTRL2", "CTRL3", "CTRL4", 
                "DIR1", "DIR2", "DIR3", "DIR4", 
                "HEAD",
                "MUT1", "MUT2", "MUT3",
                "SOMA1", "SOMA3", "SOMA4", # SOMA2 removed.
                "SOMB1", "SOMB2", "SOMB3", "SOMB4", 
                "SOMC1", "SOMC2", "SOMC3", "SOMC4", 
                "SOMD1", "SOMD2", "SOMD3", "SOMD4",
                "TRUNK",
                "WT1", "WT2", "WT3" )

expConditions <- c( rep( "Control", 4 ), 
                    rep( "Directed", 4 ),
                    "Head",
                    rep( "Mutant", 3 ),
                    rep( "SomA", 3 ),
                    rep( "SomB", 4 ),
                    rep( "SomC", 4 ),
                    rep( "SomD", 4 ),
                    "Trunk",
                    rep( "Wild-type", 3 ) )

expType <- c( rep( "Paired-end", 9 ),
              rep( "Single-read", 3 ),
              rep( "Paired-end", 16 ),
              rep( "Single-read", 3) )

tissueType <- c( rep( "Sorted-cells", 8),
                 rep( "Dissected-tissue", 4),
                 rep( "Sorted-cells", 15),
                 rep( "Dissected-tissue", 4) )

expDesign <- data.frame( sampleName = sampleIDs, 
                         fileName = countFiles,
                         type = expType,
                         tissue = tissueType,
                         condition = expConditions )
# Clean up.
rm( countFiles, sampleIDs, expConditions, expType, tissueType )
expDesign # OK

#################################################
##########           DESeq2            ##########
#################################################
dds <- DESeqDataSetFromHTSeqCount(sampleTable = expDesign, 
                                  directory = inDir, 
                                  design = ~condition)
dds <- DESeq( dds )

#################################################
######        DATA TRANSFORMATIONS         ######
#################################################
### Regularized log transformation
rld <- rlog( dds, blind = TRUE )

### Variance stabilized transformation
vsd <- varianceStabilizingTransformation( dds, blind = TRUE )


#################################################
######     EXPLORATORY SAMPLE ANALYSIS    #######
#################################################
### Determine best transformation
notAllZero <- rowSums(counts(dds)) > 0 
meanSdPlot(log2(counts(dds, normalized = TRUE)[notAllZero,]+1))
meanSdPlot(assay(rld[notAllZero,])) # Better
meanSdPlot(assay(vsd[notAllZero,]))
rm(notAllZero)
## TODO: facet plots and save them
#################################################


#################################################
### PCA plot
### Color palette: Qualitiative, 10-class paired from colorbrewer2.org
### A reference in a journal article like this:
### Brewer, Cynthia A., 200x. http://www.ColorBrewer2.org, accessed 3/27/2016.
dir.create("EDA_plots/")
pcaData <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
# PC1:49% variance, PC2: 23% variance
pcaData$condition <- factor(x = pcaData$condition, 
                            levels = c("Control", "Directed",
                                       "SomA", "SomB", "SomC", "SomD",
                                       "Wild-type", "Mutant",
                                       "Head", "Trunk" ))
pccol <- c( '#a6cee3','#1f78b4',
            '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', 
            '#b2df8a', '#33a02c', 
            '#cab2d6', '#6a3d9a' )

pca <- ggplot( data = pcaData, aes( pcaData[,1], pcaData[,2] ) ) + theme_bw()
pca <- pca + geom_point( size = 4, alpha = 0.8, aes( colour = condition ) )
pca <- pca + scale_colour_manual( values = pccol, guide_legend( "" ) ) 
pca <- pca + ggtitle( "PCA plot\nrld transformed values" )
pca <- pca + labs( x = "PC1: 49% variance", y = "PC2: 23% variance" )
pca <- pca + theme( legend.position = c( 0.85, 0.2 ) )
# pca <- pca + geom_text( aes( label = name, color = group ) )
pca
ggsave( "EDA_plots/PCA.pdf", width = 8, height = 8 )




#### OK UNTIL HERE WED 2/10/2016
res1 <- results(dds)
res1

plot( attr( res1, "filterNumRej" ), type = 'b', ylab= "number of rejections")

#plotDispEsts(dds)


# Correlation matrix plots
# rld transformed
# rld_dist <- as.matrix( dist( t( assay( rld ) ) ) )
# heatmap.2 (rld_dist, trace= "none", col= distcol, Colv= FALSE, Rowv= FALSE, dendrogram= "none" ,
#            main= "rld", 
#            key.title= NA, key.ylab=NA, key.ytick=NA, density.info="none", 
#            key.xlab= "Similarity",
#            keysize= .8 )
# # vsd transformed
# vsd_dist <- as.matrix( dist( t( assay(vsd ) ) ) )
# rld_cor <- heatmap.2 (vsd_dist, trace= "none", col= distcol, Colv= FALSE, Rowv= FALSE, dendrogram= "none" ,
#                       main= "vsd", 
#                       key.title= NA, key.ylab=NA, key.ytick=NA, density.info="none", 
#                       key.xlab= "Similarity",
#                       keysize= .8 )


### GET DESeq2 RESULTS

# Create convenience function for getting DESeq2 results
getDDSResults <- function( dds, 
                              condition, numerator, denominator, 
                              log2FC = 1, padj_cutoff = 0.1 ) {
        condition <- as.character(condition)
        numerator <- as.character(numerator)
        denominator <- as.character(denominator)
        x <- results(dds, c(condition, numerator, denominator))
        x <- as.data.frame( x )
        x <- subset(x, padj < padj_cutoff & abs( log2FoldChange ) > log2FC )
        x <- x[ order( - x$log2FoldChange ), ]
        x <- cbind( Row.Names = rownames( x ), x )
        colnames( x ) [ 1 ] <- "Gene"
        x
}



# Add contrasts to list of dataframes.
resultsList <- list()
resultsList[[1]] <- getDDSResults(dds = dds, 
                                     condition = "condition", 
                                     numerator = "Directed", 
                                     denominator = "Control", 
                                     log2FC = 1, padj_cutoff = 0.1 )

resultsList[[2]] <- getDDSResults(dds = dds, 
                                     condition = "condition", 
                                     numerator = "Mutant", 
                                     denominator = "Wild-type", 
                                     log2FC = 0, padj_cutoff = 0.1)



