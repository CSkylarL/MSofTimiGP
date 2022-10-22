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

#[2]Charoentong2017_Bindea2013_Xu2018_Immune####################################
# The tutorial is example/example02_Charoentong2017_Bindea2013_Xu2018_Immune.R
# You can learn how to generate data(Immune3_COX_MP_SKCM06) there

rm(list=ls())
data(Immune3_COX_MP_SKCM06)
cox_res <- Immune3_COX_MP_SKCM06

# a) TimiGP --------------------------------------------------------------------
data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
geneset <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune

# There are 3 genesets, find cell interaction separately

t <- names(cox_res)
res <- list()
for (i in 1:3) {
  cat("\n",t[i])
  GP <- rownames(cox_res[[t[i]]])[which(cox_res[[t[i]]]$QV<0.05)]
  background <- TimiBG(marker.pair = row.names(cox_res[[t[i]]]))
  se <- which(geneset$Dataset == t[i])
  cell_pair <- TimiCellPair(geneset = geneset[se,],core = 20)
  res[[t[i]]] <- TimiEnrich(gene = GP, background = background, 
                            geneset = cell_pair, p.adj = "BH",core=20)
}
Immune3_enrich_SKCM06 <- res
save(Immune3_enrich_SKCM06, 
     file =  "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda")
# b) Venn(Related to Figure 3E) ------------------------------------------------

myoutf1 <-"~/Mypackage/MSofTimiGP/Fig3/Venn_Immune.pdf"
library(RColorBrewer)
library(VennDiagram)
library(dplyr)

table(geneset$Dataset)
listInput <- list("Bindea2013" =  
                    unique(geneset$Gene[which(geneset$Dataset=="Bindea2013")]),               
                  "Xu2018" =  
                    unique(geneset$Gene[which(geneset$Dataset=="Xu2018")]), 
                  "Charoentong2017" =  
                    unique(geneset$Gene
                           [which(geneset$Dataset=="Charoentong2017")]))

col<-brewer.pal(9,"Set1")[c(1:4)]
venn.plot <- venn.diagram(
  listInput,
  filename = NULL,
  lwd = 3,
  col="white",
  main = paste0( "venn_3_annotations"),
  main.cex = 2,
  scaled=T,
  fill = col[1:3],
  alpha = 0.5,
  label.col = "white",
  cex = 1.7,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = col[1:3],
  cat.cex = 1.7,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.06,0.03),
  cat.pos = c(-20,20,180))

grid.newpage()
grid.draw(venn.plot)

pdf(myoutf1)
grid.newpage()
gl <- grid.layout(nrow=3, ncol=3, widths = c(0.1, 1, 0.1), 
                  heights = c(0.1, 1, 0.1))
# grid.show.layout(gl)
vp <- viewport(layout.pos.col=2, layout.pos.row=2) 
pushViewport(viewport(layout=gl))
pushViewport(vp)
grid.draw(venn.plot)
#popViewport()
dev.off()

# c) Compare 3 results (related to Figure 3F, S4) ------------------------------
rm(list = ls)
library(TimiGP)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(gridExtra)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"
myoutf2 <- "~/Mypackage/MSofTimiGP/Fig3/score_immune.pdf"
load(myinf1)  
res <- Immune3_enrich_SKCM06
# Figure S4, gene chord diagrams
pdf("~/Mypackage/MSofTimiGP/Fig3/T_Neutrophil_3_sets_gene_interactions.pdf",
    width = 12,height = 8)
par(mfrow=c(2,4)) 
TimiGeneChord(resdata = res$Bindea2013,select = 2)
TimiGeneChord(resdata = res$Charoentong2017,select = 1)
TimiGeneChord(resdata = res$Charoentong2017,select = 28)
TimiGeneChord(resdata = res$Xu2018,select = 3)

TimiGeneChord(resdata = res$Bindea2013,select = 8)
TimiGeneChord(resdata = res$Charoentong2017,select = 41)
TimiGeneChord(resdata = res$Charoentong2017,select = 59)
TimiGeneChord(resdata = res$Xu2018,select = 1)
dev.off()

# Charoentong2017 enriched markers of MDSC and Treg, 
# the results might be caused by the mutually exclusive setting

# Xu2018 result is great

# Favorability score(Figure 3F)
t <- names(res)
resdata <- data.frame()
for (i in 1:3){
  score <- TimiFS( res[[t[i]]])
  score$Dataset <- t[i]
  p <- TimiFSBar(score)
  resdata <- rbind(resdata,score)
}

# favorable cells

resdata <- resdata %>% group_by(Dataset) %>% 
  mutate(Frank=rank(-Favorable.Score,ties.method = "min"),
         Urank = rank(-Unfavorable.Score,ties.method = "min") ) %>% 
  ungroup()%>%data.frame()
rownames(resdata) <- paste0(resdata$Dataset,"_",resdata$Cell.Type)

select <- c("Bindea2013_Cytotoxic" , 
            "Charoentong2017_aCD8 T" ,"Charoentong2017_CD8 Tem" ,
            "Xu2018_CD8 T"   , 
            "Bindea2013_Th1" ,       
            "Charoentong2017_aCD4 T", "Charoentong2017_Th1" ,    
            "Xu2018_Th1",
            "Bindea2013_CD56dim NK" ,"Charoentong2017_NK" ,                               "Xu2018_NK"   )
faCell <- resdata[select,] %>%
  mutate(Favorable.Cell.Type=factor(select,levels = select),
         Dataset = factor(Dataset,levels = unique(Dataset))) %>%
  arrange(Favorable.Cell.Type)
faCell$group <- c(rep("CD8 T Cell",4),rep("CD4 T Cell",4),rep("NK Cell",3))
faCell$label <- paste0("#",faCell$Frank)
faCell
max <- max(faCell$Favorable.Score)+10

options(repr.plot.width = 1, repr.plot.height = 0.75) 
p1 <- ggplot(data = faCell,
             mapping = aes(x = Favorable.Cell.Type,y=Favorable.Score)) +
  geom_bar(mapping = aes(fill=Dataset),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  geom_text(aes(label = label), vjust = -0.2)+
  facet_grid( ~ group, space="free", scales="free",drop = TRUE) +
  scale_x_discrete(breaks = levels(faCell$Favorable.Cell.Type), labels= faCell$Cell.Type) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Cell Type",y="Favorable Score",
       title=paste0("#Favorable Cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(#legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    plot.margin = unit(c(1, 1, 1, 3), "lines") ,
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",family = "serif", size=10),
    legend.text= element_text(face="bold", color="black",family = "serif", size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=13, 
                               angle = 60,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=17),
    
    #  axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both")),
    axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(0, max), breaks=seq(0, max, 5),expand = c(0,0)) 


# unfavorable cells


select <- c(   "Bindea2013_Tcm"  ,"Charoentong2017_CD8 Tcm",
               "Charoentong2017_CD4 Tcm",
               "Bindea2013_Th2"  , "Charoentong2017_Th2", 
               "Charoentong2017_Th17","Xu2018_Th17",
               "Bindea2013_iDC","Charoentong2017_pDC"   ,
               "Charoentong2017_iDC" ,
               "Bindea2013_Neutrophil"  ,
               "Charoentong2017_Neutrophil" ,"Xu2018_Neutrophil" )

unfaCell <- resdata[select,] %>%
  mutate(Unfavorable.Cell.Type=factor(select,levels = select),
         Dataset = factor(Dataset,levels = unique(Dataset))) %>%
  arrange(Unfavorable.Cell.Type)
unfaCell$group <-c(rep("Tcm",3),rep("Th",4),
                   rep("DC",3),rep("Neutrophil",3))
unfaCell$group <- factor(unfaCell$group,levels=c("Tcm","Th","DC","Neutrophil"))
unfaCell$label <- paste0("#",unfaCell$Urank)
unfaCell
max <- max(unfaCell$Unfavorable.Score)+10

options(repr.plot.width = 1, repr.plot.height = 0.75) 
p2 <- ggplot(data = unfaCell,
             mapping = aes(x = Unfavorable.Cell.Type,y=Unfavorable.Score)) +
  scale_x_discrete(breaks = levels(unfaCell$Unfavorable.Cell), labels= unfaCell$Cell.Type) +
  geom_bar(mapping = aes(fill=Dataset),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  geom_text(aes(label = label), vjust = -0.2)+
  facet_grid( ~ group, space="free", scales="free",drop = TRUE) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Cell Type",y="Unfavorable Score",
       title=paste0("#Unfavorable Cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(#legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    plot.margin = unit(c(1, 1, 1, 3), "lines") ,
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",family = "serif", size=10),
    legend.text= element_text(face="bold", color="black",family = "serif", size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=13, 
                               angle = 60,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=17),
    
    #  axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both")),
    axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(0, max), breaks=seq(0, max, 5),expand = c(0,0)) 
p2

pdf(myoutf2,width = 6,height = 10)
do.call("grid.arrange", c(plotlist = list(p1,p2), nrow = 2))
dev.off() 

#[3] Newman2015(Related to Figure S5) ==========================================
# The tutorial is example/example03_Newman2015_LM22.R
# You can learn how to generate data("Newman2015_COX_MP_SKCM06") there

rm(list=ls())
outdir <- "~/Mypackage/MSofTimiGP/Fig3/Newman2015/"
dir.create(outdir)
data("Newman2015_COX_MP_SKCM06")
cox_res <- Newman2015_COX_MP_SKCM06

# a) Enrichment ----------------------------------------------------------------
data("CellType_Newman2015_LM22")
geneset <- CellType_Newman2015_LM22

cell_pair <- TimiCellPair(geneset = geneset,core = 20)
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
background <- TimiBG(marker.pair = row.names(cox_res))
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# b) Cell Chord Diagram---------------------------------------------------------
pdf(paste0(outdir,"circle.pdf"),width = 9,height = 9,)
TimiCellChord(resdata = res,dataset = "Newman2015")
dev.off()

# c) Favorability score---------------------------------------------------------

# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p <- TimiFSBar(score)
p
pdf(paste0(outdir,"score.pdf"),width = 8,height = 5)
print(p)
dev.off()