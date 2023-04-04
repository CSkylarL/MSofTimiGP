#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig7
# Date: 09/06/2022
# Author: Chenyang Skylar Li
# Note: 
# 1. Cell Type Annotation: CellType_Bindea2013_cancer.rda
# 2. "Bindea2013c_enrich_TCGA_pancancer.rda" 
#     and "Bindea2013c_COX_MP_TCGA_pancancer.rda"
#     are divided into csv by cancer types in the corresponding folder
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# [1]  Select Cancer type with high-quality survival statistics ################
library(TimiGP)
rm(list=ls())



myinf1 <- "~/geneset/TCGA/TCGA_Firehose_RNASeqV2_expr.rda"
# You can download the above files from https://gdac.broadinstitute.org/
myoutf1<- "~/Mypackage/MSofTimiGP/Fig7/clinical_survival_statistic.txt"
myoutf2<- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_COX_MP_TCGA_pancancer.rda"
## a) load  RNA-----------------------------------------------------------------
load(myinf1)
rna_all <- mydata
rm(mydata)
cancer_type <-
  colnames(rna_all) %>% 
  sub(pattern = "__.*", replacement = "") %>%
  unique()
dim(rna_all) # 20501 13048
length(cancer_type) #37

## b) clincal info -------------------------------------------------------------

info_all <- data.frame()
survival_event <- data.frame(type = character(),
                             alive_0 = numeric(),
                             dead_1 = numeric(),
                             total = numeric(),
                             event_percentage = numeric())

for (ct in 1:length(cancer_type)) {
  ## c) load clinical info -----------------------------------------------------
  type <- cancer_type[ct]
  myinf2 <- paste0("~/geneset/TCGA/Clinical/",type,"_Clincial_info.txt")
  # You can download the above files from https://gdac.broadinstitute.org/
  info <- read.table(myinf2, sep="\t", header=T, 
                     row.names=1, quote="",fill = TRUE)
  
  se <- c("vital_status", "days_to_death", "days_to_last_followup")
  
  if(  sum(colnames(info) %in% se) != 3 ) next
  info <- info[,se]
  
  xx <- info[, "vital_status"]
  e.surv <- ifelse(xx=="dead", 1, 0)
  e.surv[is.na(xx)] <- NA
  xx <- info[, "days_to_death"] 
  t.surv <- ifelse(!is.na(xx), xx, info[, "days_to_last_followup"])
  t.surv <- as.numeric(t.surv)
  info <- data.frame(event = e.surv,time = t.surv,row.names = rownames(info))
  
  info <- TimiCheckEvent(info)
  
  rownames(info) <- paste(type,rownames(info),sep="__")
  
  
  comSam <- intersect(row.names(info), colnames(rna_all))
  info <- info[comSam,]
  dim(info)
  
  info_all <- rbind(info_all,info)
  
  
  # d) statistical survial(Related to Table S7) --------------------------------
  ta <- table(info$event)
  alive_0<- ta[[1]]
  dead_1 <- ta[[2]]
  total <- alive_0+ dead_1
  event_percentage<- ta[[2]]/(ta[[1]]+ta[[2]]) *100
  tmp <- data.frame(type = type,
                    alive_0 = alive_0,
                    dead_1 =  dead_1,
                    total = total,
                    event_percentage = round(event_percentage,2))
  survival_event <- rbind(survival_event,tmp)
  
}

rm(info)

write.table(survival_event %>% arrange(-event_percentage,-total),
            myoutf1,sep = "\t",col.names = T,row.names = F,quote = F)


## e) select cancer type with high-quality survival data -----------------------
survival_event %>% arrange(event_percentage) %>% head()
# type alive_0 dead_1 total event_percentage
# 1 PRAD     489      8   497             1.61
# 2 THCA     482     14   496             2.82
# 3 PCPG     173      6   179             3.35
# 4 THYM     112      6   118             5.08
# cutoff
# 5 READ     116      9   125             7.20
# 6 UCEC     495     45   540             8.33
survival_event  %>% arrange(total) %>% head()
# type alive_0 dead_1 total event_percentage
# 1 CHOL      19     16    35            45.71
# 2 DLBC      40      5    45            11.11
# 3  UCS      24     32    56            57.14
# 4 KICH      56      8    64            12.50
# cutoff
# 5  ACC      54     25    79            31.65
# 6  UVM      67     13    80            16.25

# bad data due to event percentage: "PCPG","PRAD","THYM","THCA",

# bad data due to sample size: "CHOL","DLBC","KICH","UCS"

for( ct in cancer_type){
  se <- grep(pattern = paste0(ct,"__"),colnames(rna_all))
  if(length(se) == 0) {
    cat("\n",ct)}
}

# remove low-quality data & combination cohorts & liquid tumors
cancer_type <-survival_event$type
se <- which(cancer_type  %in%  
              c("PCPG","PRAD","THYM","THCA", # smallevent%
                "CHOL","DLBC","KICH","UCS", # small sample size
                "COADREAD","GBMLGG","KIPAN","STES", # combination cohorts
                "LAML","LCML","DLBC" # liquid tumors
                                 ))


cancer_type <-cancer_type[-se]

#[2] TimiCOX ###################################################################
data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- geneset$Gene %>% unique()

Bindea2013c_MPS_TCGA_pancancer <- list()
Bindea2013c_COX_MP_TCGA_pancancer <- list()
for (ct in 1:length(cancer_type)) {
  type <- cancer_type[ct]
  
  
  cat("\n",type)
  se <- grep(colnames(rna_all),pattern = paste0(type,"__"))
  if( length(se) == 0) next
  rna <- rna_all[,se]
  dim(rna)
  se <- grep(rownames(info_all),pattern = paste0(type,"__"))
  if( length(se) == 0) next
  info <- info_all[se,]
  dim(info)
  
  rna <- TimiPrePropress(marker = marker,rna = rna,
                         cohort = rownames(info),log = T,GMNorm = T)
  info <-info[colnames(rna),]
  all(colnames(rna) == rownames(info))
  
  mps_tmp <- TimiGenePair(rna)
  tmp_res <- TimiCOX(mps = mps_tmp,info = info,p.adj = "BH")
  Bindea2013c_COX_MP_TCGA_pancancer[[type]] <- tmp_res$cox_res
  Bindea2013c_MPS_TCGA_pancancer[[type]] <- tmp_res$mps
  
  xx <- tmp_res$cox_res
  cat("\n",type," PV<0.05#,%: ",sum(xx[,2]<0.05),sum(xx[,2]<0.05)/nrow(xx),
      "PV<0.01#,%: ",sum(xx[,2]<0.01),sum(xx[,2]<0.01)/nrow(xx),
      "QV<0.05#,%: ",sum(xx[,3]<0.05),sum(xx[,3]<0.05)/nrow(xx),
      "10000th,PV,QV: ",as.numeric(xx[10000,2:3]),
      "20000th,PV,QV: ",as.numeric(xx[20000,2:3]),
      "30000th,PV,QV: ",as.numeric(xx[30000,2:3]))
}


save(Bindea2013c_COX_MP_TCGA_pancancer,file=myoutf2)
# Due to the size limitation of github
# The .rda has been splited as belows to upload to github
mydir <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_COX_MP_TCGA_pancancer/"
cancer_type <- names(Bindea2013c_COX_MP_TCGA_pancancer)
for (i in cancer_type){
  write.csv(file = paste0(mydir,i,".csv"),
            Bindea2013c_COX_MP_TCGA_pancancer[[i]], 
            row.names = T,quote = F)
}

# You can choose above codes to generate "Bindea2013c_COX_MP_TCGA_pancancer.rda"
# Or you can use the re-combine the seperate files by following codes:
# mydir <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_COX_MP_TCGA_pancancer/"
# myfile <- list.files(path = mydir)
# Bindea2013c_COX_MP_TCGA_pancancer <- list()
# for (i in 1: length(myfile)){
#   nam <- gsub(myfile[i],replacement = "",pattern = ".csv")
#   Bindea2013c_COX_MP_TCGA_pancancer[[nam]] <- 
#     read.csv(file = paste0(mydir,myfile[i]),
#               row.names = 1)
# }
# save(Bindea2013c_COX_MP_TCGA_pancancer,file=myoutf2)

# [3] TimiEnrich ###############################################################

rm(list=ls())
library(dplyr)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_COX_MP_TCGA_pancancer.rda"
myoutd <- "~/Mypackage/MSofTimiGP/Fig7/"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_enrich_TCGA_pancancer.rda"
myoutf2 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_score_TCGA_pancancer.rda"

load(myinf1)
cancer_type <- names(Bindea2013c_COX_MP_TCGA_pancancer)
select_num_IMGP <- data.frame(type=cancer_type,
                              num_IMG0.05=numeric(length(cancer_type)),
                              per_IMG0.05=numeric(length(cancer_type)),
                              num_IMG0.01=numeric(length(cancer_type)),
                              per_IMG0.01=numeric(length(cancer_type)),
                              num_IMG0.001=numeric(length(cancer_type)),
                              per_IMG0.001=numeric(length(cancer_type)))


for (i in 1:nrow(select_num_IMGP)) {
  type <- select_num_IMGP$type[i]
  cat("\n",type)
 
  cox_res <-Bindea2013c_COX_MP_TCGA_pancancer[[type]]
  select_num_IMGP$num_IMG0.05[i] <- cox_res  %>% filter(PV <0.05) %>% nrow()
  select_num_IMGP$per_IMG0.05[i] <- select_num_IMGP$num_IMG0.05[i]/nrow(cox_res)
  select_num_IMGP$num_IMG0.01[i] <- cox_res  %>% filter(PV <0.01) %>% nrow()
  select_num_IMGP$per_IMG0.01[i] <- select_num_IMGP$num_IMG0.01[i]/nrow(cox_res)
  select_num_IMGP$num_IMG0.001[i] <- cox_res  %>% filter(PV <0.001) %>% nrow()
  select_num_IMGP$per_IMG0.001[i] <- select_num_IMGP$num_IMG0.001[i]/nrow(cox_res)
  
}
# [3.1] Choose the min value as the cutoff =====================================
IMGPcutoff <- min(select_num_IMGP$num_IMG0.05) 
IMGPcutoff #4352

# [3.2]Result ==================================================================
data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

res <- list()
score_all <- list()
for (ct in 1:length(cancer_type)) {
  
  type <- select_num_IMGP$type[ct]
  cat("\n",type)
  # a) Enrichment --------------------------------------------------------------
  cox_res <-Bindea2013c_COX_MP_TCGA_pancancer[[type]]
 
  GP <- cox_res %>% arrange(PV) %>% rownames() %>% .[1:IMGPcutoff]
  
  background <- TimiBG(marker.pair = row.names(cox_res))
  
  res[[type]] <- TimiEnrich(gene = GP, background = background, 
                            geneset = cell_pair, p.adj = "BH",core=20)
  
  
  # b) Calculate Favorability Score---------------------------------------------
  score <- TimiFS(res[[type]])
  head(score)
  score$Cancer.Type <- type
  score_all <- rbind(score_all,score)
  
  # c) Export PAAD resutls(related Figure S5)----------------------------------
  if(type == "PAAD"){
    # Cell Chord Diagram
    pdf(paste0(myoutd,"PAAD_circle.pdf"),width = 8,height = 9,)
    TimiCellChord(resdata = res[[type]],dataset = "Bindea2013_Cancer")
    dev.off()
    
    p <- TimiFSBar(score)
    pdf(paste0(myoutd,"PAAD_score.pdf"),width = 8,height = 5)
    print(p)
    dev.off()
  }
}

Bindea2013c_enrich_TCGA_pancancer <- res
save(Bindea2013c_enrich_TCGA_pancancer,file = myoutf1)

Bindea2013c_score_TCGA_pancancer <- score_all
save(Bindea2013c_score_TCGA_pancancer,file = myoutf2)

# Due to the size limitation of github
# The .rda has been splited as belows to upload to github
mydir <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_enrich_TCGA_pancancer/"
cancer_type <- names(Bindea2013c_enrich_TCGA_pancancer)
for (i in cancer_type){
  write.csv(file = paste0(mydir,i,".csv"),
            Bindea2013c_enrich_TCGA_pancancer[[i]], 
            row.names = F,quote = F)
}

# You can choose above codes to generate "Bindea2013c_enrich_TCGA_pancancer.rda"
# Or you can use the re-combine the seperate files by following codes:
# mydir <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_enrich_TCGA_pancancer/"
# myfile <- list.files(path = mydir)
# Bindea2013c_enrich_TCGA_pancancer <- list()
# for (i in 1: length(myfile)){
#   nam <- gsub(myfile[i],replacement = "",pattern = ".csv")
#   Bindea2013c_enrich_TCGA_pancancer[[nam]] <-
#     read.csv(file = paste0(mydir,myfile[i]))
# }
# save(Bindea2013c_enrich_TCGA_pancancer,file = myoutf1)

# [4] Figure 7 #################################################################
rm(list=ls())
library(dplyr)
library(reshape)
library(ggplot2)
library(scatterpie)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_score_TCGA_pancancer.rda"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_score_pie.pdf"
load(myinf1)
resdata <- Bindea2013c_score_TCGA_pancancer
rm(Bindea2013c_score_TCGA_pancancer)
resdata$total <- resdata$Favorable.Score+resdata$Unfavorable.Score
resdata$radius <- sqrt((resdata$total/pi/50))

# a) rank cell order:favorble-->unfavorable-------------------------------------

cel.order <- c("Cytotoxic", "T", "B","DC" , "NK" , 
               "CD8 T" ,"Tgd", "CD56dim NK",
               "Tem", "Tcm",
               "Th", "Tfh","Th1", "aDC" ,
                "CD56bright NK", "Eosinophil",
               "Mast","iDC", "Neutrophil", "Th2" ,"Macrophage","Tumor")
# b)rank: cancer order:Favorable.Score+Unfavorable.score high --> low-----------

not_immune <- c("Tumor")
can_order <-resdata[which(! resdata$Cell.Type %in% not_immune ),] %>% 
  mutate(n=Favorable.Score+Unfavorable.Score) %>% 
  group_by(Cancer.Type) %>% 
  summarise(infiltration = sum(n)) %>%
  arrange(-infiltration) %>%
  select(Cancer.Type)  

can_order <-can_order$Cancer.Type


# c) scatterpie plot -----------------------------------------------------------
resdata$x <- resdata$Cell.Type %>% 
  factor( levels=cel.order)%>% 
  as.numeric()
resdata$y <- resdata$Cancer.Type %>% 
  factor(levels =can_order  ) %>% 
  as.numeric()

colnames(resdata)


p1<- ggplot() + 
  geom_scatterpie(aes(x=x, y=y, r=radius), data=resdata,
                  cols=c("Favorable.Score","Unfavorable.Score"), 
                  legend_name = "Cell",
                  color=NA)  + 
  coord_equal() +
  scale_x_discrete(limits = factor(seq(1,length(cel.order))),
                   labels = cel.order) +
  scale_y_discrete(limits = factor(seq(1,length(can_order)+2)),
                   breaks = factor(seq(1,length(can_order))),
                   labels = c(can_order)) +
  scale_fill_manual(values=c("#f46d43","#3288bd")) +
  theme_bw(base_size = 15, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Cell Type",y="Cancer Type",
       title="TCGA 23 Solid Tumors",
       base_size = 30, base_family = "serif",face="bold") +
    legend.position = c(0.4,0.97),
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid.major = element_line(color = "#d9d9d9",
                                    size = 0.5,
                                    linetype = 2),
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",
                                family = "serif", size=30),
    legend.text= element_text(face="bold", color="black",
                              family = "serif", size=30),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=25, 
                               angle = 40,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=25),
    axis.title.x = element_text(face="bold", color="black", size=30),
    axis.title.y = element_text(face="bold",color="black", size=30))+
  geom_scatterpie_legend(resdata$radius, x=length(cel.order)-3, 
                         y=length(can_order)+1.5) 

pdf(myoutf1,
    width = 20, height = 20)
print(p1)


dev.off()

# [5]Figure S5 #################################################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(dplyr)
library(grid)
library(scales)


myinf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_enrich_TCGA_pancancer.rda"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig7/Bindea2013c_top_interaction_dotplot.pdf"
load(myinf1)

resdata <- NULL
for ( i in 1:length(Bindea2013c_enrich_TCGA_pancancer)) {
  tmp <- Bindea2013c_enrich_TCGA_pancancer[[i]]
  tmp$Cancer.Type  <- names(Bindea2013c_enrich_TCGA_pancancer)[i]
  resdata <- rbind(resdata,tmp)
}


# scale Adjust.P.Value for figure
resdata$padj.rescaled <- ifelse(resdata$Adjust.P.Value > 0.05, 
                                0.07+resdata$Adjust.P.Value/20, 
                                resdata$Adjust.P.Value)
sigCell <- resdata %>% filter(Adjust.P.Value < 0.05) %>% 
  select(Cell.Interaction) %>% table() %>% data.frame() %>% 
  arrange(-Freq)

# cell type
cel.order <- c("T",  "B",  "CD8 T", "Cytotoxic",  "DC",  "Tcm", 
               "CD56dim NK",  "Tgd", "Mast", 
               "Tem", "Tfh", "CD56bright NK", "Eosinophil",  
               "iDC", "NK",  "Th",  "aDC", "Th1",  
               "Neutrophil",  "Th2", "Tumor", "Macrophage")   

# cancer type order from tiny pie dotplot  

can_order <- c("ESCA", "SARC", "BRCA", "STAD", "LIHC", "HNSC",
               "BLCA", "UCEC", "COAD" ,
               "LUSC", "ACC",  "SKCM", "READ", "CESC", "OV",   
               "KIRP", "KIRC", "UVM"  ,
               "MESO", "GBM",  "LGG",  "LUAD", "PAAD")
# all cells related to 1 cell type.
# select top10 (#occurence in cancer type)
reference_selection <- NULL

for (cell in cel.order){
  cat("\n",cell)
  
  se <- which(resdata$Favorable.Cell.Type == cell | 
                resdata$Unfavorable.Cell.Type == cell)
  select <- resdata[se,] %>% filter(Adjust.P.Value < 0.05) %>% 
    select(Cell.Interaction) %>% table() %>% data.frame() %>% 
    arrange(-Freq)
  num <- 10
  reference_selection[[cell]] <- top_n(select,10,wt = Freq) 
}
reference_selection

p.list <- list()
ER_max <- ceiling(max(resdata$Enrichment.Ratio))

for (cell in c("Tumor",
               "Cytotoxic")){
  cat("\n",cell)
  
  
  se <- which(resdata$Favorable.Cell.Type == cell | 
                resdata$Unfavorable.Cell.Type == cell)
  select <- resdata[se,] %>% filter(Adjust.P.Value < 0.05) %>% 
    select(Cell.Interaction) %>% table() %>% data.frame() %>% 
    arrange(-Freq)
  num <- 10
  sel.int <- top_n(select,num,wt = Freq) 
  
  # select top 10
  se <-resdata$Cell.Interaction %in% sel.int$.
  
  p.list[[cell]]<- resdata[se,] %>% 
    mutate(Cell.Interaction=factor(Cell.Interaction, levels=sel.int$.),
           Cancer.Type = factor(Cancer.Type, levels=can_order)) %>%
    ggplot() +
    geom_point(
      mapping = aes(x=Cell.Interaction, y=Cancer.Type, color=padj.rescaled, 
                    size=Enrichment.Ratio)) +
    scale_colour_gradientn(colours=c("#9970ab","#5aae61","#525252","#969696"),
                           breaks = c(0,0.01,0.05,0.1),
                           labels = c(0,0.01,0.05,1),
                           limits = c( 0,0.1),
                           guide = guide_colorbar(barwidth = 1, barheight =8,
                                                  reverse=TRUE),
                           oob = scales::squish,
                           name = "Adjust\nP-Value") +
    
    scale_size(range = c(0,10), limits = c(0,ER_max),name="Enrichment\nRatio") +

    theme_bw(base_size = 15, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Cell Interaction",y="Cancer Type",
         title=paste0("TCGA-Pan-Cancer \n (",cell,") Top10"),
         base_size = 15, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, 
                                           size = 0.5, colour = 1),
          plot.margin = unit(c(1, 1, 1, 7), "lines"),
          legend.position="right",
          legend.box = "vertical",
          legend.direction= "vertical",
          panel.grid.major = element_line(color = "#d9d9d9",
                                          size = 0.5,
                                          linetype = 2),
          panel.grid=element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.title = element_text(face="bold", color="black",
                                      family = "serif", size=10),
          legend.text= element_text(face="bold", color="black",
                                    family = "serif", size=10),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black", size=10, 
                                     angle = 40,  vjust = 1, hjust=1),
          axis.text.y = element_text(face="bold", color="black", size=13),
          axis.title.x = element_text(face="bold", color="black", size=15),
          axis.title.y = element_text(face="bold",color="black", size=15))
  
  
}


pdf(myoutf1,
    width = 9, height = 8)
for (i in 1:length(p.list)){
  print(p.list[[i]])
}

dev.off()

