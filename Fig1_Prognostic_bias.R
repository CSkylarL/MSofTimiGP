#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig1
# Date: 08/09/2022
# Author: Chenyang Skylar Li

# Note: 
# 1. Immune marker gene
# 2. Immune cell infitration
# (1) Infiltration Estimation for all TCGA tumors
#     8 methods downloaded from 
#     http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz
#     immunedeconv with TIMER, CIBERSORT, quanTIseq, xCell, 
#     MCP-counter  and EPIC methods;
#     results of the cohort TCGA_SKCM06 were used in this analysis.
# 
# (2) CIBERSORTx_ABS|Ecotyper input: TPM normalized expression in SKCM06 cohort
#     ECotyper: default parameter(https://ecotyper.stanford.edu/carcinoma/)
#     CIBERSORTx: below(cibersortx.stanford.edu)
#       Date: 2022-08-09 16:06:37
#       Job type: Impute Cell Fractions
#       Signature matrix file: LM22.update-gene-symbols.txt 
#       Mixture file: SKCM06_TPM.txt
#       Batch correction: disabled
#       Disable quantile normalization: true
#       Run mode (relative or absolute): absolute
#       Permutations: 100
# (3) CIBERSORT and CIBERSORTx use absolute mode 
#     to ensure the result is comparable across samples
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[1] Immune Marker Gene results#################################################

#[1.1] Immune Marker Gene collected for this study==============================
rm(list=ls())
library(TimiGP) # The package in this study
# Install following the instruction:https://github.com/CSkylarL/TimiGP 

library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
data("Immune_Marker_n1293")
data("SKCM06info")
data("SKCM06rna")
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig1/IMG_COX.pdf"
myoutf1.1 <- "~/Mypackage/MSofTimiGP/Fig1/IMG_COX_Distribution.pdf"
myoutf1.2 <- "~/Mypackage/MSofTimiGP/Fig1/IMG_KM_eg.pdf"
myoutf1.4 <- "~/Mypackage/MSofTimiGP/Fig1/IMG_correlation_eg.pdf"
myoutf2 <- "~/Mypackage/MSofTimiGP/Fig1/IMG_Pearson.pdf"


info <- TimiCheckEvent(SKCM06info)
rna <- TimiPrePropress(marker = Immune_Marker_n1293,rna = SKCM06rna,
                       cohort = rownames(info),
                       log = T,GMNorm = T)

# a) survival(Related to Table S2) ---------------------------------------------
all((colnames(rna) == rownames(info)) == TRUE) 
mygrp <- rna

pval <- hr  <- rep(0, nrow(mygrp))
system.time({
  for(k in 1:nrow(mygrp)){
    cat("\r", k)
    mytag <- as.numeric(mygrp[k,])
    xx <- as.data.frame(cbind(mytag, info))
    mycox <- coxph(Surv(time, event)~mytag, xx) 
    mycox <- summary(mycox)
    pval[k] <- mycox$coefficients[5]
    tmp <- mycox$conf.int
    hr[k] <- tmp[1]
  }
})

QV <- p.adjust(pval, method="BH")

cox_img_res <- data.frame(HR=hr, PV=pval, QV=QV)
row.names(cox_img_res) <- row.names(mygrp)

# b) example(Related to Figure 1D and S2)----------------------------------------
good.gene <- rownames(resdata)[which(resdata$group == "Favorable")]
good.exp <- TimiPrePropress(marker = good.gene,rna = SKCM06rna,
                            cohort = rownames(info),
                            log = T,GMNorm = T)
mps <- TimiGenePair(rna = good.exp)
res <- TimiCOX(mps = mps,info = info,p.adj = "BH")

mps <- res$mps
cox_res <- res$cox_res %>% filter(PV<0.05)
cox_res$FG <- rownames(cox_res) %>% strsplit("_") %>% sapply("[[",1)
cox_res$UG <- rownames(cox_res) %>% strsplit("_") %>% sapply("[[",2)

select <- c( "IFNG_IDO2", "CXCL10_CD274","CXCL10_SIGLEC10",
             "CXCL10_IL10","CXCL11_IDO2","CXCL10_IL10","KLRK1_IL10","IFNG_IL10")

# gene KM curve
all(colnames(rna) == rownames(info))
p.sur <- list()
for( k in select%>%strsplit("_")%>%unlist()%>%unique())  {
  kmData <- info
  kmData$mytag <- as.numeric(rna[k,])
  cutoff <- median(kmData$mytag)
  # low (,median] high(median,)
  kmData$Exp <- as.factor(ifelse(test = kmData$mytag > cutoff, 
                                 "high", "low"))
  ta <- table(kmData$Exp)
  km_fit<- survfit(Surv(time,event) ~Exp, data = kmData)
  # log-rank test
  survdiff(Surv(time,event) ~Exp, data = kmData, rho = 0)
  
  if (k == "CD274") {
    k <- "PD-L1"
  }
  p.sur[[k]] <- ggsurvplot(
    km_fit,                     # survfit object with calculated statistics.
    pval = TRUE, 
    data = kmData,             # data used to fit survival curves.
    palette = c("red", "blue"),
    font.legend=25,
    legend = "top",
    legend.title="Expr.",
    xlab = "Years",   # customize X axis label.
    break.time.by = 365.25 *2,     # break X axis in time intervals by days.
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs = c(paste0("High (n=", ta[1], ")"),
                    paste0("Low (n=", ta[2], ")")),
    
    size=3,
    pval.method=T,
    pval.size=6,
    #  pval.coord = c(8000,0.6),
    Exp.table = TRUE,
    title = paste0(k,"(SKCM06)"),
    #  ggtheme = theme_bw()
    xscale ="d_y",
    ggtheme = theme_classic2(base_size=25, base_family = "serif"),
    font_family= "serif",
  )+
    guides(colour = guide_legend(nrow = 2))
}

p.list1 <- list()

for ( i in names(p.sur)){
  p.list1[[i]] <- p.sur[[i]]$plot
}

# Correlation between genes
scatterplot <- function(mydata,myx,myy){
  min <- min(mydata[myx])
  max <- max(mydata[myy]) 
  
  p<-ggplot(mydata,mapping = aes(x = get(myx),y=get(myy))) +
    geom_point(size=4,alpha=0.6) +
    geom_smooth(method=lm , color="blue", se=FALSE)+
    stat_cor(method = "pearson", label.x = min+0.1, label.y = max-0.1,
             p.accuracy = 0.001, r.accuracy = 0.01)+
    theme_bw(base_size = 12, base_family = "serif") + 
    labs(x=myx,y=myy,
         title=paste0("correlation_",myx,"_",myy)) +
    theme(legend.position="top",
          panel.grid=element_blank(),
          legend.title = element_blank(),
          # legend.title = element_text(face="bold", color="black",family = "serif", size=12),
          legend.text= element_text(face="bold", color="black",family = "serif", size=12),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black",family = "serif", size=10),
          axis.text.y = element_text(face="bold", color="black", family = "serif",size=10),
          axis.title.x = element_text(face="bold", color="black",family = "serif", size=12),
          axis.title.y = element_text(face="bold",color="black",family = "serif", size=12)) 
  return(p)
}



all(colnames(rna) == rownames(info))
p.cor <- list()
for( i in select)  {
  g1 <- strsplit(i,split = "_") %>% sapply("[[", 1)
  g2 <- strsplit(i,split = "_") %>% sapply("[[", 2)
  data <-  rna[c(g1,g2),] %>% t() %>% data.frame()
  
  
  if (g2 == "CD274") {
    g2 <- "PD-L1"
    colnames(data) <- c(g1,g2)
  }
  
  
  
  p.cor[[i]] <- scatterplot(data,g1,g2)
  
}

pdf(myoutf1.4,width = 20,height = 22)
do.call("grid.arrange", c(plotlist = p.cor, ncol=3))
dev.off() 

# Gene pair KM curve

p.sur <- list()
for( k in select)  {
  kmData <- info
  kmData$mytag <- mps[k,]
  ta <- table(kmData$mytag)
  km_fit<- survfit(Surv(time,event) ~mytag, data = kmData)
  # log-rank test
  survdiff(Surv(time,event) ~mytag, data = kmData, rho = 0)
  g1 <- strsplit(k,split = "_") %>% sapply("[[", 1)
  g2 <- strsplit(k,split = "_") %>% sapply("[[", 2)
  if (g1 == "CD274") {
    g1 <- "PD-L1"
  }
  if (g2 == "CD274") {
    g2 <- "PD-L1"
  }
  p.sur[[k]] <- ggsurvplot(
    km_fit,                     # survfit object with calculated statistics.
    pval = TRUE, 
    data = kmData,             # data used to fit survival curves.
    palette = c("red", "blue"),
    font.legend=25,
    legend = "top",
    legend.title="Diff.",
    xlab = "Years",   # customize X axis label.
    break.time.by = 365.25 *2,     # break X axis in time intervals by days.
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs = c(paste0(g1," < ",g2, "(n=", ta[1], ")"),
                    paste0(g1," > ",g2, "(n=", ta[2], ")")),
    
    size=3,
    pval.method=T,
    pval.size=6,
    Exp.table = TRUE,
    title = paste0(g1,"_",g2,"(SKCM06)"),
    ggtheme = theme_classic2(base_size=25, base_family = "serif"),
    font_family= "serif",
    xscale ="d_y")+
    guides(colour = guide_legend(nrow = 2))
}
p.list2 <- list()

for ( i in names(p.sur)){
  p.list2[[i]] <- p.sur[[i]]$plot
}

n <- length(p.list1)
nCol <- ceiling(sqrt(n))
pdf(myoutf1.2,width = 20,height = 22)
do.call("grid.arrange", c(plotlist = p.list1, ncol=nCol,top="KM_Gene"))
do.call("grid.arrange", c(plotlist = p.list2, ncol=nCol,top="KM_GP"))
dev.off() 

# c) Volcano plot(related to Figure 1B) -----------------------------------------
# Checkpoint
icg <- c("ADORA2A", "CD276", "VTCN1", "BTLA", "CTLA4", "IDO1", "KIR", "LAG3", "CYBB", "PDCD1", "CD274", "HAVCR2", "SIGLEC7", "VISTA", "VSIR", "C10orf54")

# inhibitor and stimulator
imm.inh <- c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2", "C10orf54")
imm.sti <- c("MICA", "MICB", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD70", "CD80", "CD86", "ICOS", "ICOSLG", "IL6", "IL6R", "PDCD1LG2", "TMEM173", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF13", "TNFSF13B", "TNFSF18", "TNFSF4", "TNFSF9", "TNFSF15", "TNFRSF25", "HHLA2", "TMIGD2", "BTNL2", "CD276", "CD48", "TNFSF14", "TNFRSF8", "PVR", "LTA",  "IL2RA", "ENTPD1", "NT5E", "CXCR4", "CXCL12", "KLRK1", "NKG2A", "RAET1E", "ULBP1")

# Volcano plot

resdata <- cox_img_res
resdata$label <- rownames(resdata)
resdata$group <- as.factor(
  ifelse(
    resdata$QV < 0.05,
    ifelse(resdata$HR < 1, "Favorable", "Unfavorable"),
    "NS"
  )
)

t <- table(resdata$group)
t
p1 <- ggplot(data = resdata, aes(label=label)) + 
  geom_point(aes(x=log10(HR), y=-log10(QV), color=group),alpha=0.3, size=3) +
  scale_color_manual(name="Prognosis",
                     limits=c( "Favorable","Unfavorable", "NS"),
                     values=c("#cb181d","#08519c", "#bdbdbd"),
                     labels=c(paste(names(t)[1],": ", t[1]),
                              paste(names(t)[3],": ", t[3]),
                              paste(names(t)[2],": ", t[2]))) +
  geom_vline(xintercept=log10(1), lty=4, col="#bdbdbd", lwd=2) + 
  geom_hline(yintercept=-log10(0.05), lty=4, col="#bdbdbd", lwd=2)+
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="log10 (Harzard Ratio)",y="-log10 (adj.P-Value)",
       title="Survival Analysis of IMG",
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 1, colour = 1),
        legend.position="top",
        legend.box = "horizontal",
        legend.direction= "horizontal",
        panel.grid=element_blank(),
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=15),
        legend.text= element_text(face="bold", color="black",family = "serif", size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=20),
        axis.text.y = element_text(face="bold", color="black", size=20),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))  +
  scale_y_continuous(limits=c(0, 3), breaks=seq(0,3, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(-1, 1), breaks=seq(-1, 1,0.2),expand = c(0,0))

p1.4 <- p1 +
  geom_label_repel(data = resdata[c("TIGIT","CD274","LAG3","IDO1") ,],  
                   aes(label = label,x=log10(HR), y=-log10(QV)),
                   color = "Blue",
                   segment.size  = 0.5,
                   segment.color = "grey50",
                   segment.linetype = 1,
                   min.segment.length = 0,
                   label.size = 0.5,
                   xlim = c(-0.2,NA),
                   ylim = c(1.3,NA),
                   direction = "both",
                   point.padding =0.5,
                   box.padding = 1,
                   max.overlaps = Inf,
                   seed = 1)   +
  
  geom_label_repel(data = resdata[c("KLRK1","CXCL10","CD40","CD80","CD86"), ],  
                   aes(label = label,x=log10(HR), y=-log10(QV)),
                   color = "Red",
                   segment.size  = 0.5,
                   segment.color = "grey50",
                   segment.linetype = 1,
                   min.segment.length = 0,
                   label.size = 0.5,
                   xlim = c(NA,-0.2),
                   ylim = c(1.3,NA),
                   direction = "both",
                   point.padding =0.5,
                   box.padding = 1,
                   max.overlaps = Inf,
                   seed = 1)  

pdf(myoutf1,width = 7.5,height = 8)
print(p1)
print(p1.2)
print(p1.3)
print(p1.4)
dev.off()

# Distribution of prognostic gene
distr <- data.frame(t)
distr <- distr[-2,] 
distr$Per <- distr$Freq /sum(distr$Freq) *100

p1.5 <- ggplot() +
  geom_bar(data = distr,
           mapping = aes(x = Var1,y=Per,fill=Var1),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity" ,alpha=0.7) +
  scale_fill_manual(values= c("#cb181d","#08519c"))+
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Prognostic Impact",y="Percentage of prognostic IMG(%)",
       title=paste0("Distribution of prognostic IMG"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        plot.margin = unit(c(1, 5, 1, 5), "lines"),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=13),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=17),
        axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(0, 100), breaks=seq(0,100, 10),expand = c(0,0)) +
  guides(fill ="none")

p1.5
pdf(myoutf1.1,width = 6,height = 6)
print(p1.5)
dev.off()



# d) correlation heatmap(Related to Figure 1C and S1)-----------------------------
# all gene
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
cor.all <- cor(t(rna), method ="pearson")
p2 <- Heatmap(cor.all, name = "PCC", show_column_names = F,
              show_row_names = F, 
              column_title = "Pearson Correlation of all IMG", 
              column_title_gp = gpar(fontsize = 20), col = col_fun)


# prognostic gene
se <- which(cox_img_res$QV < 0.05)
irg_sig <- rna[se, ]
cor.sig <- cor(t(irg_sig), method ="pearson")

immuGene <- unique(c(icg,imm.inh,imm.sti,"CXCL10"))
se <- which(rownames(irg_sig) %in% immuGene)
marker <- rowAnnotation("IRG" = anno_mark(at = se, 
                                          labels = rownames(irg_sig)[se]))
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
p2.1 <- Heatmap(cor.sig, name = "PCC", cluster_rows = T, 
                right_annotation = marker, show_column_names = F,
                show_row_names = F, 
                column_title = "Pearson Correlation of Prognostic IMG", 
                column_title_gp = gpar(fontsize = 20), col = col_fun)

pdf(myoutf2,width = 10,height = 9)
print(p2)
print(p2.1)
dev.off()
#[1.2] Immune Marker Gene in Current Methods====================================
rm(list=ls())
library(TimiGP)
library(ComplexHeatmap)
data("SKCM06info")
data("SKCM06rna")
info <- TimiCheckEvent(SKCM06info)
rna <- TimiPrePropress(marker = rownames(SKCM06rna),rna = SKCM06rna,
                       cohort = rownames(info),
                       log = T,GMNorm = T)

myinf <- list.files("~/Mypackage/MSofTimiGP/Fig1/MARKER/",
                    pattern = ".txt",full.names = T)

myoutf1 <- "~/Mypackage/MSofTimiGP/Fig1/MARKER/Methods_IMG_Pearson.pdf"

# a) prepare marker gene list---------------------------------------------------
myfile <- list()
for(i in 1:length(myinf)){
  nam <- basename(myinf[i]) %>% strsplit("_") %>% sapply("[[",1)
  myfile[[nam]] <- read.table(myinf[i],sep = "\t",header = T)
}


i=1
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $LM22_GENE %>% unique()
marker <- data.frame("IMG" = xx)
marker$Methods <- names(myfile[i])

i=2
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $Gene %>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])

marker <- rbind(marker,tmp)

i=3
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $Genesmarkers %>% strsplit(",",fixed = T) %>% 
  unlist() %>% na.omit()%>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])
marker <- rbind(marker,tmp)

i=4
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $Name %>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])

marker <- rbind(marker,tmp)

i=5
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $Symbol %>% strsplit(" /// ") %>% sapply("[[",1) %>% 
  na.omit()%>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])

marker <- rbind(marker,tmp)

i=6
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $GeneSymbol %>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])

marker <- rbind(marker,tmp)


i=7
names(myfile[i])
names(myfile[[i]])
xx <- myfile[[i]] $Gene %>% strsplit("|",fixed = T) %>% unlist() %>% 
  na.omit()%>% unique()
tmp <- data.frame("IMG" = xx)
tmp$Methods <- names(myfile[i])

marker <- rbind(marker,tmp)

# b) correlation & heatmap(Related to Figure S1) ----------------------------
t <- table(marker$Methods)
t <- names(t)
p <- list()
for ( i in t){
  cat("\r",i)
  se <- which(marker$Methods == i)
  mygen <- marker$IMG[se]
  
  com <- intersect(x = mygen, y=rownames(rna))
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  cor.all <- cor(t(rna[com,]), method ="pearson")
  p[[i]] <- Heatmap(cor.all, name = "PCC", show_column_names = F,
                    show_row_names = F, column_title = paste0("Pearson Correlation of " , i, " IMG"), 
                    column_title_gp = gpar(fontsize = 20), col = col_fun)
}

pdf(myoutf1)
for ( i in t){
  print(p[[i]])
}
dev.off()
#[2] Immune infiltration #######################################################
#[2.1] Integrate immune infiltration results ===================================
# a) TIMER2.0: 8 methods immune infiltration -----------------------------------
library(survival)

rm(list=ls())
myinf1 <- "~/Mypackage/MSofTimiGP/Fig1/infiltration_estimation_for_tcga.csv"

library(TimiGP)
library(ggrepel)


data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)
pid <- rownames(info)
# Select SKCM06 patients
data<- read.csv(myinf1, row.names = 1, header = T)
xx <- as.numeric(substr(rownames(data), 14, 15))
se <- which(xx==6)		## metastatic
length(se) # 394 samples
data <- data[se, ]
dim(data) # 394 119
rownames(data) <- substr(rownames(data), 1, 12)
imm.inf <- data[pid,]


# CIBERSORT.ABS, instead of CIBESORT, is comparable across sample,
# remove CIBERSORT results and rename "CIBERSORT.ABS" to CIBERSORT"
methods <- colnames(imm.inf) %>% strsplit("_") %>% sapply( "[[",2) 
mid <- unique(methods)
mid
se <- which(methods == "CIBERSORT")
imm.inf <- imm.inf[-se]
colnames(imm.inf) <- gsub(colnames(imm.inf),pattern = "CIBERSORT.ABS",replacement = "CIBERSORT")
# remove non-cell/ unkown cell type
xx <- c("immune.score_XCELL", "stroma.score_XCELL" ,
        "cytotoxicity.score_MCPCOUNTER"  ,
        "uncharacterized.cell_QUANTISEQ" ,
        "uncharacterized.cell_EPIC")
imm.inf <- imm.inf[-which(colnames(imm.inf) %in% xx)]
# b) CIBERSORTx Resut ----------------------------------------------------------
myinf2 <- "~/Mypackage/MSofTimiGP/Fig1/CIBERSORTx_ABS_LM22/CIBERSORTx_Job8_Results.txt"

data <- read.table(myinf2,header = T,row.names = 1,sep = "\t")

data <- data [1: (ncol(data)-4)]
colnames(data) <- paste0(colnames(data),"_CIBERSORTx" )
all(rownames(imm.inf) == rownames(data))
imm.inf <- cbind(imm.inf,data)
colnames(imm.inf)
# c) EcoTyper Resut ------------------------------------------------------------

xx <- list.dirs("~/Mypackage/MSofTimiGP/Fig1/ecotyper_output/") %>% 
  list.files(pattern = "Cell_State_Abundance.txt",full.names = T) 

for (i in 1:length(xx)){
  if (i==1){
    nam <- strsplit(basename(xx[i]),"_") %>% unlist()
    nam <- nam[1]
    file <- read.table(xx[i],header = T,row.names = 1)
    colnames(file) <- paste0(nam,".",colnames(file))
    data <- file
  } else {
    nam <- strsplit(basename(xx[i]),"_") %>% unlist()
    nam <- nam[1]
    file <- read.table(xx[i],header = T,row.names = 1)
    colnames(file) <- paste0(nam,".",colnames(file))
    if(all(rownames(file) == rownames(data))){
      data <- cbind(data,file)
    } else{
      warning(nam)
    }
    
    
  }
  
}
xx <- gsub(rownames(data),pattern = ".",replacement = "-",fixed = T)
all(xx == rownames(imm.inf))
rownames(data) <- xx

colnames(data) <- paste0(colnames(data),"_ECOTYPER" )

imm.inf <- cbind(imm.inf,data)
colnames(imm.inf)

# save--------------------------------------------------------------------------
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig1/Immune_Infiltration.rda"
save(imm.inf,file = myoutf1)

#[2.2] Calculate Correlation of each method (Related to Figure S1)==============
library(dplyr)
library(ComplexHeatmap)
library(circlize)
rm(list=ls())
myinf1 <- "~/Mypackage/MSofTimiGP/Fig1/Immune_Infiltration.rda"
load(myinf1)

methods <- colnames(imm.inf) %>% strsplit("_") %>% sapply( "[[",2) 
mid <- unique(methods)

for (i in 1:length(mid)){
  se <- which(methods == mid[i])
  data <- imm.inf[se]
  colnames(data) <- colnames(data) %>% strsplit("_") %>% sapply( "[[",1) 
  cor <- cor(data, method ="pearson")
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  p1 <- Heatmap(cor, name = "PCC", cluster_rows = T, show_column_names = F,
                show_row_names = T, column_title = paste0("Pearson Correlation of ",mid[i]), 
                column_title_gp = gpar(fontsize = 15), col = col_fun,
                width = ncol(cor)*unit(10, "mm"), 
                height = nrow(cor)*unit(10, "mm"))
  
  pdf(paste0("~/Mypackage/MSofTimiGP/Fig1/Pearson_Correlation_",mid[i],".pdf"),
      width = ncol(cor)*1,  
      height = nrow(cor)*0.7)
  print(p1)
  dev.off()
}


#[2.3] Survival Analysis========================================================
# a) cox regression(Related to Table S2, Figure 1A) ----------------------------
library(dplyr)
rm(list=ls())
myinf1 <- "~/Mypackage/MSofTimiGP/Fig1/Immune_Infiltration.rda"
load(myinf1)


library(TimiGP)


data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)

all((rownames(imm.inf) == rownames(info)) == TRUE)


methods <- colnames(imm.inf) %>% strsplit("_") %>% sapply( "[[",2) 
mid <- unique(methods)

res <- NULL
for (i in mid) {
  se <- which(methods == i)
  mygrp <- imm.inf[,se]
  colnames(mygrp) <- colnames(mygrp) %>% strsplit("_") %>% sapply( "[[",1) 
  ct <- pval <- hr  <- rep(0, ncol(mygrp))
  mt <- rep(i, ncol(mygrp))
  for (k in 1:ncol(mygrp)) {
    cat("\r", k)
    ct[k] <- names(mygrp)[k]
    mytag <- as.numeric(mygrp[, k])
    xx <- as.data.frame(cbind(mytag, info))
    mycox <- coxph(Surv(time, event)~mytag, xx) 
    mycox <- summary(mycox)
    pval[k] <- mycox$coefficients[5]
    tmp <- mycox$conf.int
    hr[k] <- tmp[1]
  }
  res <- rbind(res,data.frame(Method=mt, Cell=ct, HR=hr, PV=pval))
}

# b) save-----------------------------------------------------------------------
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig1/Immune_Infiltration_Cox.rda"
save(res,file = myoutf1)

# c) volcano plot --------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggrepel)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig1/Immune_Infiltration_Cox.rda"

myoutf1 <- "~/Mypackage/MSofTimiGP/Fig1/volcano_all_method.pdf"
myoutf2 <- "~/Mypackage/MSofTimiGP/Fig1/volcano_all_method_distribution.pdf"

load(myinf1)

# resdata <- res %>% filter(HR < 1e+6 & HR > 1e-6)
resdata <- res %>% mutate(HR = replace(HR, HR > 1e+6,  1e+6)) %>% 
  mutate(HR = replace(HR, HR < 1e-6, 1e-6))




resdata$label <- paste(resdata$Method, resdata$Cell, sep = "_")
resdata$group <- as.factor(
  ifelse(
    resdata$PV < 0.05,
    ifelse(resdata$HR < 1, "Favorable", "Unfavorable"),
    "NS"
  )
)
# Volcano plot
t <- table(resdata$group)

p <- ggplot(data = resdata, aes(label=Cell)) + 
  geom_point(aes(x=log10(HR), y=-log10(PV), color=group,shape=Method),alpha=0.7, size=5) +
  scale_shape_manual(name="Method",
                     limits=c("MCPCOUNTER","TIMER","ECOTYPER",
                              "CIBERSORTx","CIBERSORT","EPIC","XCELL","QUANTISEQ"),
                     values = c(15:18,0:3))+
  scale_color_manual(name="Prognosis",
                     limits=c( "Favorable","Unfavorable", "NS"),
                     values=c("#cb181d","#08519c", "#bdbdbd"),
                     labels=c(paste(names(t)[1],": ", t[1]),
                              paste(names(t)[3],": ", t[3]),
                              paste(names(t)[2],": ", t[2]))) +
  geom_vline(xintercept=log10(1), lty=4, col="#bdbdbd", lwd=2) + 
  geom_hline(yintercept=-log10(0.05), lty=4, col="#bdbdbd", lwd=2)+
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="log10 (Harzard Ratio)",y="-log10 (P-Value)",
       title="Cox Result of Immune Infiltration",
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_blank(),
        legend.position="right",
        legend.box = "vertical",
        legend.direction = "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=15),
        legend.text= element_text(face="bold", color="black",family = "serif", size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=20),
        axis.text.y = element_text(face="bold", color="black", size=20),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))  +
  scale_y_continuous(limits=c(0, 5), breaks=seq(0,5, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(-6.5, 6.5), breaks=seq(-6, 6,2),expand = c(0,0)) +
  guides(shape=guide_legend(order =1, nrow=8,byrow=TRUE),
         color=guide_legend(order =2, nrow=3,byrow=T)) +
  annotate("text",x=-5,y=3.5, color="#cb181d",family = "serif", size=10,
           label = paste0("Pronostic\nBias\n",
                          round(t[1]/(t[1]+t[3]) *100,1),"%"))


se1 <- which(resdata$label %in% c("TIMER_T.cell.CD8.","TIMER_Neutrophil",
                                  "XCELL_Plasmacytoid.dendritic.cell"))
                                

p2 <-p+
  geom_label_repel(data = resdata[se1,],  
                     aes(label = Cell,x=log10(HR), y=-log10(PV)),
                     color = "Red",
                     size=3,
                     segment.size  = 0.5,
                     segment.color = "grey50",
                     segment.linetype = 1,
                     min.segment.length = 0,
                     label.size = 0.5,
                     xlim = c(-6,-2),
                     ylim = c(1,3),
                     direction = "both",
                     # point.padding =0.5,
                     # box.padding = 1,
                     max.overlaps = 100,
                     force = 100,
                     seed = 1) 
  
p2

pdf(myoutf1,
    width = 10,height = 8)
print(p)
print(p1)
print(p2)
dev.off()

distr <- data.frame(t)
distr <- distr[-2,] 
distr$Per <- distr$Freq /sum(distr$Freq) *100

#d) Distribution barplot -------------------------------------------------------
p3 <- ggplot() +
  geom_bar(data = distr,
           mapping = aes(x = Var1,y=Per,fill=Var1),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity" ,alpha=0.7) +
  scale_fill_manual(values= c("#cb181d","#08519c"))+
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Prognostic Impact",y="Percentage of prognostic cell(%)",
       title=paste0("Distribution of prognostic cell"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        plot.margin = unit(c(1, 5, 1, 5), "lines"),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=13),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=17),
        axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(0, 100), breaks=seq(0,100, 10),expand = c(0,0)) +
  guides(fill ="none")

p3
pdf(myoutf2,width = 6,height = 6)
print(p3)
dev.off()


