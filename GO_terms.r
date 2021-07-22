##----Set home directory & create the main folder----

homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(clusterProfiler, DOSE, ggplot2, enrichplot)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Create gene list----
d <- read.csv("treatment1_vs_control.csv") # read raw data
genelist <- d[,3] # get log2FoldChange column into a geneList matrix
names(genelist) <- as.character(d[,10]) # get ENSEMBL names for log2FoldChange
genelist <- na.omit(genelist) # omit NA values
genelist <- sort(genelist, decreasing = TRUE) # sort descending

#Gene Set Enrichment Analysis of Gene Ontology
gse <- gseGO(geneList = genelist,
              ont = "CC", #BP: biological pathways, MF: molecular function, CC: cell cycle, ALL
              minGSSize = 3, # minimal size of each gene set for analysing
              maxGSSize = 800, # maximal size of each gene set for analysing
              pvalueCutoff = 0.05,
              verbose = TRUE, # progress
              OrgDb = org.Ce.eg.db, #organism database
              pAdjustMethod = "none",
              eps = 0) # boundary for calculating the p value

##----Dot plot for enrichment result----
dp <- dotplot(gse, showCategory = 10, split = ".sign", orderBy = "x") + 
  facet_grid(.~.sign) + labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Ridge plot for enrichment result----
rp <- ridgeplot(gse, showCategory = 10) + 
  labs(x = "Enrichment distribution for Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Enrichment analysis----
genes <- names(genelist)
ego <- enrichGO(gene = genes,
                OrgDb = org.Ce.eg.db,
                ont = "ALL", #BP, MF or CC
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

##----Bar plot----
bp <- barplot(ego, 
              drop = TRUE, 
              showCategory = 10, 
              title = "GO All Pathways", 
              font.size = 8)

##----Dot plot----
dp2 <- dotplot(ego,  showCategory = 10) + 
  labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Emap Plots----
# get similarity matrix
pt <- pairwise_termsim(ego,
                       method = "JC", # Jaccard similarity coefficient
                       semData = NULL, # GOSemSimDATA object
                       showCategory = 200) # number of enriched terms to display

ep <- emapplot(pt) + 
  labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Cnet plots----
cp <- cnetplot(ego, 
               showCategory = 9, 
               categorySize = "pvalue", 
               foldChange = genelist, 
               layout = "kk", 
               colorEdge = TRUE) +
  labs(title = "HGD-1") + 
  scale_color_gradient(low = "#808080", high = "#88419d")

##----Session info----
sessionInfo()