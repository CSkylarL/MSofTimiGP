#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig6
# Date: 8/27/2022
# Author: Chenyang Skylar Li
# Note: 
# Use the codes in the example folder
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list = ls())
library(TimiGP)

#[1] Tirosh2016_melanoma_TME(Related to Figure 6B,C)############################
# The tutorial is example/example04_Tirosh2016_melanoma_TME.R
# You can learn how to generate data("melanoma_TME_COX_MP_SKCM06") there

rm(list=ls())
outdir <- "~/Mypackage/MSofTimiGP/Fig6/Tirosh2016_melanoma_TME/"
data("melanoma_TME_COX_MP_SKCM06")
cox_res <- melanoma_TME_COX_MP_SKCM06

# a) Enrichment ----------------------------------------------------------------
data("CellType_Tirosh2016_melanoma_TME")
geneset <- CellType_Tirosh2016_melanoma_TME
geneset$Dataset <- "Tirosh2016"
cell_pair <- TimiCellPair(geneset = geneset,core = 20)
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
background <- TimiBG(marker.pair = row.names(cox_res))
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# b) Cell Chord Diagram---------------------------------------------------------
pdf(paste0(outdir,"circle.pdf"),width = 10,height =10,)
TimiCellChord(resdata = res,dataset = "Tirosh2016")
dev.off()

# c) Favorability score---------------------------------------------------------
# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p1 <- TimiFSBar(score,select = c(1:2,4:6))
p1

pdf(paste0(outdir,"score.pdf"),width = 6,height = 5)
print(p1)
dev.off()


#[2] Zheng2021_Tcell(Related to Figure 6D,E)####################################
# The tutorial is example/example04_Zheng2021_Tcell.R
# You can learn how to generate data(Tcell_COX_MP_SKCM06) there

rm(list=ls())
outdir <- "~/Mypackage/MSofTimiGPt/Fig6/Zheng2021_Tcell/"
data(Tcell_COX_MP_SKCM06)
cox_res <- Tcell_COX_MP_SKCM06

# a) Enrichment ----------------------------------------------------------------
data("CellType_Zheng2021_Tcell")
geneset <- CellType_Zheng2021_Tcell
cell_pair <- TimiCellPair(geneset = geneset,core = 20)
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
background <- TimiBG(marker.pair = row.names(cox_res))

res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)


# b) Cell Chord Diagram---------------------------------------------------------

pdf(paste0(outdir,"circle.pdf"),width = 20,height = 21)
TimiCellChord(resdata = res,dataset = "Zheng2021")
dev.off()
 
# c) Favorability score---------------------------------------------------------

# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p1 <-  TimiFSBar(score,select = c(1:4,6,7:10,12,16,29,32,38:40))+
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, 5),expand = c(0,0))
p1
pdf(paste0(outdir,"score.pdf"),width = 14,height = 10)
print(p1)
dev.off()

