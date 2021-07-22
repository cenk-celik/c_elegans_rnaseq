##----Set home directory & create the main folder----
homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(DESeq2, AnnotationDbi, org.Ce.eg.db, pathview, gage, gageData, dplyr)
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

##----KEGG pathway analysis----
# Set up KEGG database
kg.cel <- kegg.gsets("cel", id.type = 'entrez') # cel = C elegans
kegg.gs <- kg.cel$kg.sets[kg.cel$sigmet.idx]
head(kegg.gs, 3) # top three

foldchanges <- res$log2FoldChange # get fold changes
names(foldchanges) <- res$entrez # assign gene names to fold changes
keggres <- gage(foldchanges, gsets = kegg.gs, same.dir = TRUE) # find pathways with fold changes
lapply(keggres, head)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
  .$id %>%
  as.character()

keggresids <- substr(keggrespathways, start = 1, stop = 8) # substring found pathway names

# Define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data = foldchanges, # fold changes
                                       pathway.id = pid, # pathway ids
                                       species = "cel", # organism
                                       new.signature = FALSE)

detach("package:dplyr", unload=T) # Unload dplyr since it conflicts with the next line
keggLocal <- paste0(rnaLocal,"kegg/treatment1_vs_control") # directory for KEGG plots
dir.create(keggLocal, recursive = TRUE, showWarnings = FALSE) # create dir
setwd(keggLocal) # set dir to KEGG pathways
tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges, 
                                                 pathway.id = pid, 
                                                 species = "cel")) # plot found pathways
setwd(rnaLocal) # set back to main dir

##----Session info----
sessionInfo()