##----Set home directory & create the main folder----
homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(ggplot2, DESeq2, AnnotationDbi, org.Ce.eg.db, pheatmap, genefilter, RColorBrewer, rafalib)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Load saved featureCounts----
load(paste0(rnaLocal, "Rdata/featurecounts.RData"))

countMatrix <- featurecounts$counts # combine counts as a matrix
colnames(countMatrix) <- gsub("_R1.fastq.gz.bam$", "", colnames(countMatrix)) # rename columns
coldata <- read.csv("coldata.csv", header = TRUE) # column information
countMatrix <- countMatrix[, coldata$X]

##----Differential Expression Analysis----
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = coldata,
                              design = ~ condition) # construct DESeq2 matrix

##----Pre-filtering----
keep <- rowSums(counts(dds) >= 10) >= 5 # remove low counts less than 10 for at least 5 biological replicates
dds <- dds[keep,]

dds <- DESeq(dds)
colnames(dds) <- coldata$names # rename matrix columns

##----treatment vs control----
res <- results(dds, contrast = c("condition", "treatment1", "control"),
               lfcThreshold = 0.5, # only +/- 0.5 fold change and above
               alpha = 0.05, # pvalue
               altHypothesis = "greaterAbs") # alternative hypothesis

##----Annotation----
res$symbol <- mapIds(org.Ce.eg.db, # annotation database for organism
                     keys=row.names(res), #gene names as reference
                     column="SYMBOL", # to add symbols 
                     keytype="WORMBASE", # from this database
                     multiVals="first")
res$name <- mapIds(org.Ce.eg.db,
                    keys=row.names(res),
                    column="GENENAME", # to add gene functions
                    keytype="WORMBASE",
                    multiVals="first")
res$entrez <- mapIds(org.Ce.eg.db,
                     keys=row.names(res),
                     column="ENTREZID", # to add entrez ids
                     keytype="WORMBASE",
                     multiVals="first")
res$GO <- mapIds(org.Ce.eg.db,
                 keys=row.names(res),
                 column="GO", #to add GO terms
                 keytype="WORMBASE",
                 multiVals="first")

##----Filter results----

res <- res[order(res$pvalue), ] # sort based on p values, ascending
write.csv(as.data.frame(res), file = "treatment1_vs_control.csv") #save raw results as .CSV
summary(res) # results at a glance

report <- res[which(res$name != "hypothetical protein" & res$pvalue <= 0.05), ] # if required, filter unknown proteins and those insignificant
write.csv(as.data.frame(report), file="filtered_treatment1_vs_control.csv") # save filtered results as .CSV

##----VSD normalisation----
# if less than 100 samples:
# vsd <- varianceStabilizingTransformation(dds)
vsd <- vst(dds, blind=FALSE, nsub = 10000) # top 10000 significant differentially expressed genes

##----Heatmap----
plotsLocal <- paste0(rnaLocal, "plots") # dir for plots
dir.create(plotsLocal, recursive = TRUE, showWarnings = FALSE) # create dir for plots

topVarGenes <- head(order(-rowVars(assay(vsd))), 10000) # top 10000 variance genes, ordered
mat <- assay(vsd)[topVarGenes, ] # matrix for heat map
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[,"condition"]) # columns for labelling heatmap
colnames(df) <- "Diet" # Treatment name
rownames(df) <- colnames(vsd) # gene names for heatmap
ann_colors <- list(Diet = c(nd = "#808080", hgd1 = "#88419D", hgd5 = "#E6550D")) # colours for treatments

mapColours <- brewer.pal(n = 11, name = "PiYG") #colours for heat map gradient
mapColours <- colorRampPalette(mapColours)(100)

pheatmap(mat, cluster_rows = TRUE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df, annotation_colors = ann_colors,
         angle_col = 45, color = mapColours) # plot heatmap with gene names, sample names angled

##----Principle Component Analysis----
plotPCA(vsd, intgroup = c("condition")) # PCA analysis

##----Distance Matrix----
sampleDists <- dist(t(assay(vsd))) # transpose vsd matrix and calculate distances
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, col = mapColours, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45,
         annotation_col = df, annotation_colors = ann_colors) # plot sample distance matrix as a heat map

##----Volcano plots----
mypar(1, 1) # to show plots side by side (row, col)
# with to plot two overlapping plots
with(res, plot(log2FoldChange, # fold changes to plot
               -log10(pvalue), # log10 p values for Z-score
               pch = 20, # type of data points
               main = "Treatment1 vs Control", # title for plot
               xlim = c(-6,6), # x axis limits
               ylim = c(0,22), # y  axis limits
               col="#808080")) # colour for Control
with(subset(res, padj < 0.1 ), # adjusted p values less than 0.1
     points(log2FoldChange,
            -log10(pvalue),
            pch=20,
            col="#88419D")) # colour for Treatment1

##----MA plot and histograms----
# MA plot
plotMA(res,
       colNonSig = "#808080", # colour for Control
       colSig = "#88419D", # colour for Treatment1
       main = "Treatment 1",
       colLine = "black") # title for Treatment1

# Histogram
hist(res$pvalue[res$baseMean > 1], 
     breaks = 20, # bin size
     col = "#88419D", # colour for Treatment1
     border = NULL, 
     main = "Treatment 1", # title 
     xlab = "p-value") # x axis label

##----Plot counts for individual genes----
col <- c("#808080", "#88419D", "#E6550D") # colours for treatments
col.dots <- rep(c("#808080", "#88419D", "#E6550D"), each = 6) # colours for data points

#Gene1
Gene1 <- plotCounts(dds,
                    gene = "WBGene00016097", #ENSEMBL id for a gene of interest
                    intgroup = "condition", # compare by
                    returnData = TRUE)

Gene1 <- ggplot(d, 
                aes(x = forcats::fct_relevel(condition, c("Control", "Treatment1","Treatment2")), # by "condition", reordered
                    y = count, # gene counts
                    fill = col.dots)) + # colours for data points
        labs(x = "", y = "Counts", title = "GENE 1") +
        theme(panel.background = element_rect(colour = "black", fill = "transparent"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.position = "none") +
        geom_violin(trim = FALSE) + # violin plot
        geom_boxplot(width = 0.1, fill = "white") + # box plot in a violin plot
        scale_fill_manual(values = col) + # colours for box plots
        scale_x_discrete(labels=c("Control" = "Control", "Treatment1" = "Treatment 1", "Treatment 2" = "Treatment 2")) # relabel x axis

##----Session info----
sessionInfo()