#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig4
# Date: 02/02/2023
# Author: Chenyang Skylar Li
# Note: 
# 1. Public validation sets are downloaded from GEO and 
# prepossessed to select metastatic melanoma with qualified survival statistics.
# 2. There are five public validation sets 
#   in "~/Mypackage/MSofTimiGP/Fig4/validation.dataset":  
#  Cirenajwis_GSE65904.rda,  Jayawardana_GSE54467.rda 
# Jonsson_GSE22155.rda,  Mann_GSE53118.rd,a  Xu_GSE8401.rda.
# 3. "Liu_phs000452.v3.p1" and "VanAllen_phs000452.v2.p1" have controlled access,
#  Please request them from the owners if you needed them.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [0] validation set ###########################################################
# Due to the size limitation of github
# The validation.dataset.rda has been splitted as belows to upload to github
# Please run below codes to re-combine the seperate files
mydir <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset/"
myfile <- list.files(path = mydir)
validation.dataset <- list()
for (i in 1: length(myfile)){
  nam <- gsub(myfile[i],replacement = "",pattern = ".rda")
  load(paste0(mydir,myfile[i]))
  validation.dataset[[nam]]$rna <- rna
  validation.dataset[[nam]]$info <- info
  rm(rna,info)
}

save(validation.dataset, 
     file = "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda")


# [1] Run TimiGP in public data (Validation data) ##############################
# Use CellType_Bindea2013_cancer
library(TimiGP)
library(dplyr)
rm(list=ls())

myinf1 <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda"
load(myinf1)


for (DS in names(validation.dataset)) {
  cat("\n",DS,"\n")
  myoutd <- paste0("~/Mypackage/MSofTimiGP/Fig4/",DS,"/")
  dir.create(myoutd)
  myoutf1 <- paste0(myoutd,"cox_res.rda")
  myoutf2 <- paste0(myoutd,"enrich_res.rda")
  myoutf3 <- paste0(myoutd,"cirle.pdf")
  myoutf4 <- paste0(myoutd,"score.pdf")
  
  
  # Data =======================================================================
  rna <- validation.dataset[[DS]]$rna
  info <- validation.dataset[[DS]]$info %>% select(e.surv, t.surv) 
  colnames(info) <- c("event","time")
  
  
  # TimiGP =====================================================================
  data("CellType_Bindea2013_cancer")
  geneset <- CellType_Bindea2013_cancer
  marker <- unique(geneset$Gene)
  
  se <- which(rownames(rna) %in% marker)
  rna <- rna[se,]
  
  mps <- TimiGenePair(rna)
  
  res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
  cox_res <- res$cox_res
  
  cell_pair <- TimiCellPair(geneset = geneset,core = 20)
  
  GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
  if (length(GP) == 0) {
    GP <- rownames(cox_res)[which(cox_res$PV<0.05)]
    cat(DS," PV < 0.05 No.GP=",length(GP))
  }else if ( length(GP) <1000){
    GP <- rownames(cox_res)[which(cox_res$PV<0.01)]
    cat(DS," PV < 0.01 No.GP=",length(GP))
  }
  
  
  background <- TimiBG(marker.pair = row.names(cox_res))
  
  enrich_res <- TimiEnrich(gene = GP, background = background, 
                           geneset = cell_pair, p.adj = "BH",core=20)
  # Export =====================================================================
  
  save(file = myoutf1,cox_res)
  save(file = myoutf2,enrich_res)
  
  pdf(myoutf3,width = 9,height = 11,)
  TimiCellChord(resdata = enrich_res,dataset = "Bindea2013_Cancer")
  dev.off()
  
  score <- TimiFS(enrich_res)
  p <- TimiFSBar(score)
  pdf(myoutf4,width = 10,height = 5)
  print(p)
  dev.off()
}

# [2] Run TimiGP in "Liu_phs000452.v3.p1" ######################################
# This data has controlled access. 
# Please request it from the owners if you needed it.
# Use CellType_Bindea2013_cancer

# control works
library(TimiGP)
library(dplyr)
rm(list=ls())

myinf1 <-  "~/Mydata/Liu_phs000452.v3.p1.rda"
load(myinf1)
DS <-  gsub(x = basename(myinf1),pattern = ".rda",replacement = "")

cat("\n",DS,"\n")
myoutd <- paste0("~/Mypackage/MSofTimiGP/Fig4/",DS,"/")
dir.create(myoutd)
myoutf1 <- paste0(myoutd,"cox_res.rda")
myoutf2 <- paste0(myoutd,"enrich_res.rda")
myoutf3 <- paste0(myoutd,"cirle.pdf")
myoutf4 <- paste0(myoutd,"score.pdf")


# Data =========================================================================
info <- info %>% select(OS.event,OS.time)  %>% TimiCheckEvent()
colnames(info) <- c("event","time")

data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- unique(geneset$Gene)

range(rna)
rna <- TimiPrePropress(marker = marker,rna = rna,
                       log = T,GMNorm = T,
                       cohort = rownames(info))

all(rownames(info) == colnames(rna))

# TimiGP =======================================================================
mps <- TimiGenePair(rna)

res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
cox_res <- res$cox_res



cell_pair <- TimiCellPair(geneset = geneset,core = 20)


GP <- rownames(cox_res)[which(cox_res$PV<0.05)]
cat(DS," PV < 0.05 No.GP=",length(GP))



background <- TimiBG(marker.pair = row.names(cox_res))

enrich_res <- TimiEnrich(gene = GP, background = background, 
                         geneset = cell_pair, p.adj = "BH",core=20)
# Export =======================================================================

save(file = myoutf1,cox_res)
save(file = myoutf2,enrich_res)

pdf(myoutf3,width = 9,height = 11,)
TimiCellChord(resdata = enrich_res,dataset = "Bindea2013_Cancer")
dev.off()

score <- TimiFS(enrich_res)
p <- TimiFSBar(score)
pdf(myoutf4,width = 10,height = 5)
print(p)
dev.off()



# [3] "VanAllen_phs000452.v2.p1" ###############################################
# This data has controlled access. 
# Please request it from the owners if you needed it.
# Use CellType_Bindea2013_cancer

library(TimiGP)
library(dplyr)
rm(list=ls())

myinf1 <-  "~/Mydata/VanAllen_phs000452.v2.p1.rda"
load(myinf1)
DS <-  gsub(x = basename(myinf1),pattern = ".rda",replacement = "")

cat("\n",DS,"\n")
myoutd <- paste0("~/Mypackage/MSofTimiGP/Fig4/",DS,"/")
dir.create(myoutd)
myoutf1 <- paste0(myoutd,"cox_res.rda")
myoutf2 <- paste0(myoutd,"enrich_res.rda")
myoutf3 <- paste0(myoutd,"cirle.pdf")
myoutf4 <- paste0(myoutd,"score.pdf")


# Data =========================================================================
info <- info %>% select(dead,overall_survival)  %>% TimiCheckEvent()
colnames(info) <- c("event","time")

data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- unique(geneset$Gene)

range(rna)
rna <- TimiPrePropress(marker = marker,rna = rna,
                       log = T,GMNorm = T,
                       cohort = rownames(info))

all(rownames(info) == colnames(rna))

# TimiGP =======================================================================
mps <- TimiGenePair(rna)

res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
cox_res <- res$cox_res



cell_pair <- TimiCellPair(geneset = geneset,core = 20)


GP <- rownames(cox_res)[which(cox_res$PV<0.05)]
cat(DS," PV < 0.05 No.GP=",length(GP))



background <- TimiBG(marker.pair = row.names(cox_res))

enrich_res <- TimiEnrich(gene = GP, background = background, 
                         geneset = cell_pair, p.adj = "BH",core=20)
# Export =======================================================================

save(file = myoutf1,cox_res)
save(file = myoutf2,enrich_res)

pdf(myoutf3,width = 9,height = 11,)
TimiCellChord(resdata = enrich_res,dataset = "Bindea2013_Cancer")
dev.off()

score <- TimiFS(enrich_res)
p <- TimiFSBar(score)
pdf(myoutf4,width = 10,height = 5)
print(p)
dev.off()


# [4] Figure 4A similarity #####################################################
rm(list=ls())
library(TimiGP)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
myoutf1 <-  "~/Mypackage/MSofTimiGP/Fig4/JC.pdf"
myinf <- list.files("~/Mypackage/MSofTimiGP/Fig4/",
                    recursive = T, full.names = T,
                    pattern = "enrich_res.rda") %>%.[c(1,3, 7, 2, 5, 6, 4)]
# load data --------------------------------------------------------
DS <- 
  myinf %>% gsub(pattern = "/home/cli15/Mypackage/MSofTimiGP/Fig4/",
                 replacement = "") %>%
  gsub(pattern = "/enrich_res.rda",replacement = "") %>% strsplit("_") %>%
  lapply("[[", 2) %>% unlist()

resdata <- data.frame()
for ( i in 1:length(DS)){
  load(myinf[i])
  score <- enrich_res %>% filter(Adjust.P.Value < 0.05) 
  score$Dataset <- DS[i]
  resdata <- rbind(resdata,score)
}

data("Bindea2013c_enrich")
enrich_res <- Bindea2013c_enrich [-11]
score <- enrich_res %>% filter(Adjust.P.Value < 0.05) 
score$Dataset <- "TCGA_SKCM06"

resdata <- rbind(score,resdata)%>%  
  mutate( Dataset = factor(Dataset,levels = unique(Dataset))) 


# Tversky index  --------------------------------------------------------
DS <- unique(resdata$Dataset)
p.data <- matrix(nrow = length(DS),ncol = length(DS))
colnames(p.data) <- rownames(p.data) <- DS
p.value <- p.data
# Tversky index 

DS <- unique(resdata$Dataset)
p.data <- matrix(nrow = length(DS),ncol = length(DS))
colnames(p.data) <- rownames(p.data) <- DS
for ( i in 1:nrow(p.data)){
  tmpi <- resdata %>% 
    filter(Dataset == rownames(p.data)[i]) %>% 
    pull(Cell.Interaction)
  for (j in 1:ncol(p.data)){
    tmpj <- resdata %>% 
      filter(Dataset == rownames(p.data)[j]) %>% 
      pull(Cell.Interaction)
    
    shared <- intersect(tmpi, tmpj) %>% length()
    union <- union(tmpi, tmpj) %>% length()
    p.data[i,j] <-  shared / min(length(tmpi),length(tmpj)) 
    
  }
}

round(p.data,2)

# p.value  --------------------------------------------------------
data("CellType_Bindea2013_cancer")
geneset  <- CellType_Bindea2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)
nn <- unique(cell_pair$Cell.Pair) %>% length()# total number of interaction
for ( i in 1:nrow(p.value)){
  tmpi <- resdata %>% 
    filter(Dataset == rownames(p.value)[i]) %>% 
    pull(Cell.Interaction)
  for (j in 1:ncol(p.value)){
    tmpj <- resdata %>% 
      filter(Dataset == rownames(p.value)[j]) %>% 
      pull(Cell.Interaction)
    
    shared <- intersect(tmpi, tmpj) %>% length()
    union <- union(tmpi, tmpj) %>% length()
    aa <- length(tmpi)
    bb <- length(tmpj)
    p.value[i,j] <-  sum(dhyper(shared:min(aa,bb), aa, nn - aa, bb)) # p-value
  }
}
round(p.value) 
max(p.value) # 0.0001150132. # all p-value <= 0.0001

# heatmap  --------------------------------------------------------
col_fun <- colorRamp2(c( 0, 1), c("white", "red"))
p1 <- Heatmap(p.data, name = "JC", cluster_rows = F, show_column_names = F,
              cluster_columns = F,
              show_row_names = T, column_title = "Tversky index ", 
              column_title_gp = gpar(fontsize = 15), col = col_fun,
              width = ncol(p.data)*unit(10, "mm"), 
              height = nrow(p.data)*unit(10, "mm"))

pdf(paste0(myoutf1),
    width = ncol(p.data)*1,  
    height = nrow(p.data)*0.7)
print(p1)
dev.off()

# [5] Figure B,C,D examples ####################################################

rm(list=ls())
library(TimiGP)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
myoutf1 <-  "~/Mypackage/MSofTimiGP/Fig4/FigBCD.pdf"
# load data --------------------------------------------------------
myinf <- list.files("~/Mypackage/MSofTimiGP/Fig4/",
                    recursive = T, full.names = T,
                    pattern = "enrich_res.rda") %>%.[c(1,3, 7, 2, 5, 6, 4)]

DS <- 
  myinf %>% gsub(pattern = "/home/cli15/Mypackage/MSofTimiGP/Fig4/",
                 replacement = "") %>%
  gsub(pattern = "/enrich_res.rda",replacement = "") %>% strsplit("_") %>%
  lapply("[[", 2) %>% unlist()
# score  --------------------------------------------------------
mycell <- c(          "Cytotoxic" ,
                      "Tumor")
resdata <- data.frame()
for ( i in 1:length(DS)){
  load(myinf[i])
  score <- TimiFS(enrich_res) %>% filter(Cell.Type %in% mycell) 
  score$Dataset <- DS[i]
  resdata <- rbind(resdata,score)
}

data("Bindea2013c_enrich")
enrich_res <- Bindea2013c_enrich [-11]
score <- TimiFS(enrich_res) %>% filter(Cell.Type %in% mycell) 
score$Dataset <- "TCGA_SKCM06"
resdata <- rbind(score,resdata) %>%  
  filter(Cell.Type %in% mycell)  %>%
  mutate( Dataset = factor(Dataset,levels = unique(Dataset))) 

max <- max(resdata$Favorable.Score)+10
p1 <-   ggplot() +
  geom_bar(data = resdata[which(resdata$Cell.Type %in% mycell[1]), ],
           mapping = aes(x = Dataset,y=Favorable.Score, fill=Dataset),
           stat="identity", alpha=0.5) +
  scale_fill_manual(values=c(brewer.pal(9,"Greys")[8],
                             rep(brewer.pal(9,"Set1")[1],3),
                             rep(brewer.pal(9,"Set1")[2],2),
                             brewer.pal(9,"Set1")[4],
                             brewer.pal(9,"Set1")[5])) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Dataset",y="Favorable Score",
       title=paste0("Cytotoxic Cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=13,
                                   angle = 60,  vjust = 1, hjust=1),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=17),
        axis.title.y = element_text(face="bold",color="black", size=17))+
  
  scale_y_continuous(limits=c(0, max), breaks=seq(0, max, 5),expand = c(0,0))  +
  guides(fill = "none")



max <- max(resdata$Unfavorable.Score)+10
p2 <-   ggplot() +
  geom_bar(data = resdata[which(resdata$Cell.Type %in% mycell[2]), ],
           mapping = aes(x = Dataset,y=Unfavorable.Score, fill=Dataset),
           stat="identity", alpha=0.5) +
  scale_fill_manual(values=c(brewer.pal(9,"Greys")[8],
                             rep(brewer.pal(9,"Set1")[1],3),
                             rep(brewer.pal(9,"Set1")[2],2),
                             brewer.pal(9,"Set1")[4],
                             brewer.pal(9,"Set1")[5])) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Dataset",y="Unfavorable Score",
       title=paste0("Tumor Cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=13,
                                   angle = 60,  vjust = 1, hjust=1),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=17),
        axis.title.y = element_text(face="bold",color="black", size=17))+
  
  scale_y_continuous(limits=c(0, max), breaks=seq(0, max, 5),expand = c(0,0)) +
  guides(fill = "none")

# dot plot Cytotoxic-->Tumor
p.data <- data.frame()
for ( i in 1:length(DS)){
  load(myinf[i])
  score <- enrich_res %>% filter(Cell.Interaction == "Cytotoxic_Tumor") 
  score$Dataset <- DS[i]
  p.data <- rbind(p.data,score)
}

data("Bindea2013c_enrich")
enrich_res <- Bindea2013c_enrich
score <- enrich_res %>% filter(Cell.Interaction == "Cytotoxic_Tumor") 
score$Dataset <- "TCGA_SKCM06"

p.data <- rbind(score,p.data)%>%  
  mutate( Dataset = factor(Dataset,levels = unique(Dataset))) 

max <- max(ceiling(p.data$Enrichment.Ratio)) + 1

p3 <- ggplot(p.data) + 
  geom_point(mapping = aes(x = Dataset, 
                           y = Enrichment.Ratio, 
                           color = Adjust.P.Value, 
                           size = No.Shared.IMGP)) + 
  scale_color_continuous(low = "#9970ab", high = "#5aae61", 
                         limits = c(0, 0.05), name = "Adjust P-Value", 
                         guide = guide_colorbar(reverse = TRUE)) + 
  scale_size(range = c(3, 8), name = "No. Shared\nMarker Pair") + 
  geom_hline(yintercept = 0, lty = 4, col = "black", lwd = 1) + 
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  labs(x = "Dataset", y = "Enrichment Ratio", 
       title = "Cytotoxic_Tumor", 
       base_size = 12, base_family = "serif", face = "bold") + 
  theme(legend.background = element_rect(linetype = 1, 
                                         size = 0.5, colour = 1), 
        legend.position = "right", 
        legend.box = "vertical", legend.direction = "vertical", 
        plot.margin = unit(c(1, 1, 1, 5), "lines"), 
        panel.grid = element_blank(), 
        legend.key.width = unit(0.5, "cm"), 
        legend.title = element_text(face = "bold", 
                                    color = "black", 
                                    family = "serif", size = 10), 
        legend.text = element_text(face = "bold", color = "black", 
                                   family = "serif", size = 10), 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 10, angle = 60, 
                                   vjust = 1, hjust = 1), 
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 13), 
        axis.title.x = element_text(face = "bold", 
                                    color = "black", size = 15), 
        axis.title.y = element_text(face = "bold", 
                                    color = "black", size = 15)) + 
  scale_y_continuous(limits = c(0,  max), 
                     breaks = seq(0, max, 2), expand = c(0, 0)) + 
  guides(shape = "none") + 
  annotate("text", x = -1, y = max, label = "Good\nPrognosis") + 
  coord_cartesian(xlim = c(1, nrow(p.data)), clip = "off")

library("cowplot")
pdf(myoutf1,width = 10,height = 10)
ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = .5, height = .45) +
  draw_plot(p2, x = .5, y = 0, width = .5, height = .45) +
  draw_plot(p3, x = 0.3, y = 0.5, width = 0.6, height = 0.45) 
dev.off()




# [6] Charoentong2017_Bindea2013_Xu2018_Immune##################################
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
     file =  "~/Mypackage/MSofTimiGP/Fig4/Immune3_enrich_SKCM06.rda")
# b) Venn(Related to Figure 4E) ------------------------------------------------

myoutf1 <-"~/Mypackage/MSofTimiGP/Fig4/Venn_Immune.pdf"
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

# c) Compare 3 results (related to Figure 4F, S4) ------------------------------
rm(list = ls)
library(TimiGP)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(gridExtra)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig4/Immune3_enrich_SKCM06.rda"
myoutf2 <- "~/Mypackage/MSofTimiGP/Fig4/score_immune.pdf"
load(myinf1)  
res <- Immune3_enrich_SKCM06
# Figure S4, gene chord diagrams
pdf("~/Mypackage/MSofTimiGP/Fig4/T_Neutrophil_3_sets_gene_interactions.pdf",
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

# Favorability score(Figure 4F)
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
            "Bindea2013_CD56dim NK" ,"Charoentong2017_NK" , "Xu2018_NK"   )
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
  scale_x_discrete(breaks = levels(faCell$Favorable.Cell.Type), 
                   labels= faCell$Cell.Type) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Cell Type",y="Favorable Score",
       title=paste0("#Favorable Cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(
    plot.margin = unit(c(1, 1, 1, 3), "lines") ,
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",
                                family = "serif", size=10),
    legend.text= element_text(face="bold", color="black",
                              family = "serif", size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=13, 
                               angle = 60,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=17),
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
  theme(
    plot.margin = unit(c(1, 1, 1, 3), "lines") ,
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",
                                family = "serif", size=10),
    legend.text= element_text(face="bold", color="black",
                              family = "serif", size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=13, 
                               angle = 60,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=17),
    
    axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(0, max), breaks=seq(0, max, 5),expand = c(0,0)) 
p2

pdf(myoutf2,width = 6,height = 10)
do.call("grid.arrange", c(plotlist = list(p1,p2), nrow = 2))
dev.off() 

#[7] Newman2015(Related to Figure S3) ==========================================
# The tutorial is example/example03_Newman2015_LM22.R
# You can learn how to generate data("Newman2015_COX_MP_SKCM06") there

rm(list=ls())
outdir <- "~/Mypackage/MSofTimiGP/Fig4/Newman2015/"
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