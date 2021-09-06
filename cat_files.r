homeDir <- system("echo ${HOME}", intern=T)
rnaLocal <- paste0(homeDir,"/glu1920/")
dir.create(rnaLocal, recursive=T, showWarnings = F)

dataLocal <- paste0(rnaLocal,"data/")
dir.create(dataLocal, recursive=T, showWarnings = F)

rawFiles <- list.files(paste0(rnaLocal, "data"), pattern=".*RCW*", full.names = T) # raw data

untarLocal <- paste0(rnaLocal,"untar")
dir.create(untarLocal, recursive=TRUE, showWarnings = FALSE)

library(plyr)
ldply(.data = rawFiles, .fun = untar, exdir = untarLocal)

folderList <- list.dirs(path = untarLocal, full.names = TRUE, recursive = FALSE)

lapply(folderList, function(x){
  setwd(x)
  system(paste0("cat *_R1_unaligned_00*.fastq.gz > ", paste0(basename(x), "_R1_cat.fastq.gz")))
  system(paste0("cat *_R2_unaligned_00*.fastq.gz > ", paste0(basename(x), "_R2_cat.fastq.gz")))
  system('for x in *_R*_cat.fastq.gz; do gunzip -c "$x" | wc -l && echo $x; done')
  setwd(untarLocal)
})






