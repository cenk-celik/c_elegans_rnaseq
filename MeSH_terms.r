##----Set home directory & create the main folder----
homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(meshes, MeSH.Cel.eg.db)
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

##----Create MeSH object----
mesh <- gseMeSH(genelist, 
                MeSHDb = "MeSH.Cel.eg.db", 
                database = 'gene2pubmed', 
                category = "C", 
                eps = 0)

gseaplot <- gseaplot2(mesh, 
                      geneSetID = 1, # geneSetID could be a vector, 1:x
                      pvalue_table = TRUE) 

rp <- ridgeplot(mesh) + 
  labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Gene-Concept network plot----
cp <- cnetplot(mesh, 
               categorySize = "pvalue", 
               foldChange = genelist, 
               node_label = "all") + 
  labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----Heat map-like plot for functional classification----
hp <- heatplot(mesh, 
               foldChange = genelist, 
               showCategory = 10) + 
  labs(title = "Treatment 1") + 
  scale_color_gradient(low = "#56B1F7", high = "#132B43")

##----UpSet plot----
up <- upsetplot(mesh, 11) + # n, number of categories to be plotted
  labs(title = "Treatment 1")

##----Session info----
sessionInfo()