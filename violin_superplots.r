##----Set home directory & create the main folder----
homeDir <- system("echo ${HOME}", intern = TRUE) # set main directory
rnaLocal <- paste0(homeDir,"MAIN_DIRECTORY/") # name of main directory
dir.create(rnaLocal, recursive = TRUE, showWarnings = FALSE) # create main directory

##----Install libraries----
libraries <- rlang::quos(ggplot2, ggpubr, dplyr, ggbeeswarm)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Create a dir for superplots
superPlots <- paste0(homeDir,"superPlots/")
dir.create(superPlots, recursive = TRUE, showWarnings = FALSE)
setwd("superPlots") # set working dir

##----Set colours for treatments
col <- c("#88419D", "#E6550D", "#808080")

##----Read data----
# data set should have three columns
# Genes   Treatment   Counts
combined <- read.csv("SUPERPLOT_DATASET.csv")
ReplicateAverages <- combined %>% 
  group_by(treatment, genes) %>% 
  summarise_each(list(mean))

title <- expression(paste(italic("gene-1"), " text text")) # title for the plot with partially italicised text

ggplot(combined,
       aes(x = forcats::fct_relevel(treatment, c("Control", "Treatment 1","Treatment 5")), y = counts, color = factor(genes))) + 
  labs(x = "", y = "Log2FC", title = title, colour = "Genes") +
  stat_ydensity(geom = "violin", # violin layout
                aes(x = forcats::fct_relevel(treatment, c("Control", "Treatment 1","Treatment 5")), y = counts, group = treatment),
                show.legend = FALSE) +
  geom_beeswarm(cex = 1.5, pch = 20, priority = "density") + # add data points
  scale_colour_brewer(palette = "Set3", aesthetics = "fill") + # colours for genes
  scale_fill_manual(values = col) + 
  geom_beeswarm(data = ReplicateAverages, # add averages for genes
                size = 2, 
                show.legend = TRUE,
                pch = 20, 
                stroke = 1) + 
  stat_compare_means(data = ReplicateAverages, # add significance levels in the plot
                     comparisons = list(c("Control", "Treatment 1"), c("Control", "Treatment 2")),
                     method = "t.test", #label = "p.signif",
                     paired = TRUE) + 
  theme(panel.background = element_rect(colour = "black", fill = "transparent"), # aesthetics for plot canvas
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##----Session info----
sessionInfo()