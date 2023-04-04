#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig3
# Date: 08/26/2022
# Author: Chenyang Skylar Li
# Note: 
# Use the codes in the example folder
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list = ls())
library(TimiGP)

#[1]Bindea2013_Cancer(Related to Figure 3A-D, Talbe S3) ########################
# The tutorial is example/example01_Bindea2013_Cancer.R
# You can learn how to generate data("Bindea2013c_enrich") there
rm(list=ls())
outdir <- "~/Mypackage/MSofTimiGP/Fig3/Bindea2013_Cancer/"
data("Bindea2013c_enrich")
res <- Bindea2013c_enrich 

# a) Cell network files (related to Figure 3C) ---------------------------------
NET <- TimiCellNetwork(resdata = res,dataset = "Bindea2013_Cancer",
                       export =TRUE, path = outdir)
# The files are used to visualize network with Cytoscape
# b) Dot plot  (related to Figure 3A) ------------------------------------------
p <- TimiDotplot(resdata = res,select = c(1:10))
pdf(paste0(outdir,"dotplot.pdf"),width = 8,height = 6)
print(p)
dev.off()

# c)Cell Chord Diagram (related to Figure 3B) ----------------------------------
pdf(paste0(outdir,"circle.pdf"),width = 8,height = 9,)
TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer")
dev.off()

# d) Favorability score (related to Figure 3D) ---------------------------------
# Calculate
score <- TimiFS(res)
head(score)
# barplot
p <- TimiFSBar(score,select = c(1:5,(nrow(score)-2):nrow(score)))
pdf(paste0(outdir,"score_select.pdf"),width = 5,height = 5)
print(p)
dev.off()

