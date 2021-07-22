##----Set home directory & create the main folder----
homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(Rsubread)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Set the directory for reference genome & annotation file
refLocal <- paste0(rnaLocal,"ref") # directory for reference genome and annotation
dir.create(refLocal, recursive = TRUE, showWarnings = FALSE) # create directory for reference genome and annotation

##----Download genome and anootation files----
faFTP <- "ftp://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/dna/" # link for reference genome
gtfFTP <- "ftp://ftp.ensembl.org/pub/release-104/gtf/caenorhabditis_elegans/" # link for annotation file

faFile <- "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz" # name of reference genome
gtfFile <- "Caenorhabditis_elegans.WBcel235.104.gtf.gz" # name of annotation file

download.file(paste0(faFTP, faFile), paste0(refLocal, faFile)) # download reference genome
download.file(paste0(gtfFTP, gtfFile), paste0(refLocal, gtfFile)) # download annotation

##----Build the index----
idxLocal <- paste0(rnaLocal,"index") # create index using reference genome
dir.create(idxLocal, recursive= TRUE, showWarnings = FALSE) # directory for index

system(paste0("gunzip ",paste0(refLocal, faFile))) # unzip reference genome archive
refFile <- list.files(refLocal, pattern = "fa$", full = TRUE) # find reference genome file
buildindex(basename = paste0(idxLocal,"Ce104"), reference = refFile) # unzip with the name

##----Alignment----
system(paste0("tar -zxvf ", rnaLocal, "FILENAME.tar.gz -C ",rnaLocal)) # untar the data

mcIndex <- paste0(idxLocal,"ce104") # name of index
pair1Files <- list.files(paste0(rnaLocal, "data"), pattern=".*_R1.*", full.names = TRUE) # raw data
bamLocal <- paste0(rnaLocal,"bam") # directory to save aligned samples
dir.create(bamLocal, recursive = TRUE, showWarnings = FALSE) # create directory for .BAM

lapply(pair1Files, function(fn1){ # function to align all raw data
  fileStub <- gsub("_lib.*", "", basename(fn1)) 
  align(index=mcIndex, # reference
        readfile1=fn1, # 5' read
        readfile2=gsub("_R1", "_R2", fn1), # 3' read
        output_file=paste0(bamLocal,fileStub,".bam")) # rename the .BAM files
})

dir(bamLocal)

##----Feature Counts----
system(paste0("gunzip ", paste0(refLocal, gtfFile)))
gtf <- list.files(refLocal, pattern=".*gtf$", full.names = TRUE)

bamFiles <- list.files(bamLocal, pattern = ".*bam$")
curdir <- getwd() # remember current directory
setwd(bamLocal) #temporarily set directory to .BAM folder

featurecounts <- featureCounts(files = bamFiles, GTF.featureType = "exon",
                               GTF.attrType = "gene_id", annot.ext = gtf,
                               isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)

setwd(curdir) # set back to main directory

knitr::kable(featurecounts$stat) # have a look at featureCounts stats

##----Save featureCounts table----
dir.create(paste0(rnaLocal, "Rdata"), recursive = FALSE, showWarnings = FALSE)
save(featurecounts, file = paste0(rnaLocal, "Rdata/featurecounts.RData"))

##----Session info----
sessionInfo()