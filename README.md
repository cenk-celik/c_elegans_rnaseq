# RNA-seq Data Analysis

## Before alignment
Untar raw data in individual folders, concatenate if the samples are run in two lanes, then count number of lines in 5' and 3' direction. [cat_files.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/cat_files.r)

## Alignment
Using _Rsubread_ package, this piece of code will help analyse from constructing reference genome to annotation and from alignment to .BAM files. Finally, resulting .BAM files can be used to build _featureCounts_ and combined in a matrix. [alignment.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/alignment.r)

## Data pre-processing and data exploration
Using _DESeq2_ package, the _featureCounts_ matrix will be used to analyse differentially-expressed genes. A general overview of results can be seen by creating histograms, volcano plots, heat maps, distance matrices and _Principle Component Analysis_. [data_exploration.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/data_exploration.r)

## KEGG pathway analysis
Results from _DESeq2_ analysis can be constructed into a new matrix with the most differentially expressed genes to identify activated KEGG pathways. [kegg_pathways.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/kegg_pathways.r)

## Gene Ontology terms
Again, using the same results, a gene list can be used to find out "Gene Ontology" terms to visualise enriched pathways. [GO_terms.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/GO_terms.r)

## Violin Superplots
Finally, to display genes associated with certain metabolic pathways, violin superplots can be useful. [violin_superplots.r](https://github.com/cenk-celik/complete-rna-seq/blob/main/violin_superplots.r)
