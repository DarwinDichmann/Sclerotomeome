### TITLE: Scleome analysis including second sequencing of peat mutants.


### LOAD PACKAGES
# TODO: add a section that tests for install packages and install if necessary.
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(cluster)


### DIRECTORY PATHS.
### Your workdir should contain the folder with the countfiles.
setwd("./")
inDir <- file.path( "countFiles/" )
# inDir <- file.path( paste( getwd(), "count_files/", sep = "/" ) )


### DEFINE COLOR PALETTES
# PCA plots; adjust number to conditions. See other palettes PCA section.
# pccol <- colorRampPalette( brewer.pal( 10, "Spectral" ) ) ( 10 )

# Heatmap colors
hmcol <- colorRampPalette( rev( brewer.pal( 9, "RdBu" ) ) ) ( 255 )
# Distance matrix colors
distcol <- colorRampPalette( rev( brewer.pal(9, "GnBu") ) ) ( 20 )


### CREATE EXPERIMENTAL DESIGN TABLE
countFiles <- list.files( path = inDir, 
                          pattern = "(WT*)|(MUT*)|(SOMA*)|(SOMB*)|(SOMC*)|(SOMD*)|(CTRL*)|(DIR*)|(HEAD*)|(TRUNK*)" ) # No space between patterns
# Create experimental design list corresponding to file order.
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
# Clean up temp vars.
rm( countFiles, sampleIDs, expConditions, expType, tissueType )
expDesign # OK


### DESeq2 
## TODO Investigate MultiFactorial design
dds <- DESeqDataSetFromHTSeqCount(sampleTable = expDesign, 
                                  directory = inDir, 
                                  design = ~condition)
dds <- DESeq( dds )


### INITIAL TRANSFORMATIONS
# Regularized log transformation
rld <- rlog( dds, blind = TRUE )
# Variance stabilized transformation
vsd <- varianceStabilizingTransformation( dds, blind = TRUE )


### EXPLORATORY SAMPLE ANALYSIS
# PCA
# Better colors for PCA plot.
pccol <- c( '#a6cee3','#1f78b4','#6a3d9a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#b2df8a')
#6a3d9a previous pccol2[10]
#b2df8a previous pccol[3]
# rld_pca <- plotPCA( rld, intgroup = "condition", returnData = FALSE ) +
#         theme_bw() + ggtitle("rld transformed") +
#         scale_colour_manual( values = pccol ) +
#         geom_text( aes( label = name ))
# rld_pca

vsd_pca <- plotPCA( vsd, intgroup = "condition", returnData = FALSE ) +
        geom_point( size = 4, alpha = 0.5 ) +
        theme_bw() + ggtitle("vsd transformed") +
        scale_colour_manual( values = pccol ) +
        geom_text( aes( label = name ) )
# TODO: fix text position adn size. (Could add a 'series' component?)
vsd_pca

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



