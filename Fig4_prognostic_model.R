#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TimiGP manuscript Fig4
# Date: 08/31/2022
# Author: Chenyang Skylar Li
# Note: 
# 1. Validation sets are downloaded from GEO and 
#   prepossessed to select metastatic melanoma with qualified survival statistics.
# 2. There are five validation sets 
#   in "~/Mypackage/MSofTimiGP/Fig4/validation.dataset":  
#  Cirenajwis_GSE65904.rda,  Jayawardana_GSE54467.rda 
# Jonsson_GSE22155.rda,  Mann_GSE53118.rd,a  Xu_GSE8401.rda.
# 3. Training sets: TCGA_SKCM06.
#     Please follow the 
#     example/example02_Charoentong2017_Bindea2013_Xu2018_Immune.R
#     to generate Immune3_MPS_SKCM06.rda
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
# [1] LASSO for prognostic or enrich pairs (related to Figure 4B, Table S4)#####
rm(list = ls())
library(TimiGP)

library(glmnet)
library(doMC)
registerDoMC(cores = 50)

myinf1 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
# Please follow the example/example02_Charoentong2017_Bindea2013_Xu2018_Immune.R
# to generate Immune3_MPS_SKCM06

myinf2 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/lasso_res.rda"

# [1.1] LASSO ==================================================================
load(myinf1)
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)
all(rownames(info) == colnames(Immune3_MPS_SKCM06$Charoentong2017))

data("Immune3_COX_MP_SKCM06")
cox_res <- Immune3_COX_MP_SKCM06$Charoentong2017
prognosticGP <- rownames(cox_res)[which(cox_res$QV<0.05)]
length(prognosticGP)

mps_prognostic <- Immune3_MPS_SKCM06$Charoentong2017[prognosticGP,]

load(myinf2)
enrich_res <- Immune3_enrich_SKCM06$Charoentong2017
enrichGP <- enrich_res %>% filter(Adjust.P.Value < 0.05) %>% 
  pull(Shared.IMGP) %>%
  strsplit("/",fixed = T) %>% unlist() %>% unique()
length(enrichGP) 

mps_enrich <- Immune3_MPS_SKCM06$Charoentong2017[enrichGP,]


iter <- 1000
cv.lassoModel <- list()
lasso_res <- data.frame( "Seed" = numeric(length = iter),
                         "No.Predictor.Prognostic" = numeric(length = iter),
                         "Predictor.Prognostic" = character(length = iter),
                         "No.Predictor.Enrich" = numeric(length = iter),
                         "Predictor.Enrich" = character(length = iter),
                         "Cohort" = character(length = iter))


for (i in 1:iter) {
  cat("\n", i)

  # a) Random sampling 80% cohorts----------------------------------------------
  lasso_res$Seed[i] <- i
  set.seed(i)
  xx <- sample(x = 1:nrow(info), size = round(0.8*nrow(info)), replace = F)
  cohort <- rownames(info)[xx]

  lasso_res$Cohort[i] <- paste0(cohort,collapse = "/")

  # b) LASSO from prognostic pairs----------------------------------------------
  cv.lassoModel$prognostic[[i]] <- cv.glmnet(
    x=t(mps_prognostic[,xx]),
    y=Surv(info$time[xx], info$event[xx]),
    standardize=F,
    alpha=1.0,
    nfolds=10,
    family="cox",
    type.measure = "C",
    parallel=TRUE)

  lambda <- cv.lassoModel$prognostic[[i]] $lambda.1se 
  cv.lassoModel$prognostic[[i]] $co <- coef(cv.lassoModel$prognostic[[i]] , 
                                            s=lambda, exact=TRUE,
                                            x=t(mps_prognostic[,xx]),
                                            y=Surv(info$time[xx], 
                                                   info$event[xx]))
  predictor <- rownames(
    cv.lassoModel$prognostic[[i]] $co)[which(
      cv.lassoModel$prognostic[[i]] $co != 0)]
  lasso_res$No.Predictor.Prognostic[i] <- length(predictor)
  lasso_res$Predictor.Prognostic[i] <- paste0(predictor, collapse = "/")
  rm(lambda)
  rm(predictor)

  # c) LASSO from enriched pairs------------------------------------------------
  cv.lassoModel$enrich[[i]] <- cv.glmnet(
    x=t(mps_enrich[,xx]),
    y=Surv(info$time[xx], info$event[xx]),
    standardize=F,
    alpha=1.0,
    nfolds=10,
    family="cox",
    type.measure = "C",
    parallel=TRUE)

  lambda <- cv.lassoModel$enrich[[i]] $lambda.1se 
  cv.lassoModel$enrich[[i]] $co <- coef(object = cv.lassoModel$enrich[[i]] , s=lambda, exact=TRUE,
                                        x=t(mps_enrich[,xx]),
                                        y=Surv(info$time[xx], info$event[xx]))
  predictor <- rownames(cv.lassoModel$enrich[[i]] $co)[which(cv.lassoModel$enrich[[i]] $co != 0)]
  lasso_res$No.Predictor.Enrich[i] <- length(predictor)
  lasso_res$Predictor.Enrich[i] <- paste0(predictor, collapse = "/")
  rm(lambda)
  rm(predictor)
}

save(lasso_res, file = myoutf1)



# [1.2]  Frequency of predictors================================================
rm(list = ls())
myinf1 <- "~/Mypackage/MSofTimiGP/Fig4/lasso_res.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"

myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/LASSO_predictor.rda"

load(myinf1)

Predictor.Prognostic <- lasso_res$Predictor.Prognostic %>%
  strsplit(split = "/") %>%
  unlist() %>%
  table() %>%
  data.frame(row.names = 1)
names(Predictor.Prognostic) <- c( "Frequency")
Predictor.Prognostic$Percentage <- Predictor.Prognostic$Frequency/1000*100
Predictor.Prognostic <- Predictor.Prognostic[order(Predictor.Prognostic$Frequency, 
                                                   decreasing = T), ]
Predictor.Prognostic$Rank <- rank(-Predictor.Prognostic$Frequency,
                                  ties.method="min")


Predictor.Enrich <- lasso_res$Predictor.Enrich %>%
  strsplit(split = "/") %>%
  unlist() %>%
  table() %>%
  data.frame(row.names = 1)
names(Predictor.Enrich) <- c("Frequency")
Predictor.Enrich$Percentage <- Predictor.Enrich$Frequency/1000*100
Predictor.Enrich <- Predictor.Enrich[order(Predictor.Enrich$Frequency, 
                                           decreasing = T), ]
Predictor.Enrich$Rank <- rank(-Predictor.Enrich$Frequency,ties.method="min")

data("Immune3_COX_MP_SKCM06")
cox_res <- Immune3_COX_MP_SKCM06$Charoentong2017
prognosticGP <- rownames(cox_res)[which(cox_res$QV<0.05)]


load(myinf2)
enrich_res <- Immune3_enrich_SKCM06$Charoentong2017
enrichGP <- enrich_res %>% filter(Adjust.P.Value < 0.05) %>% 
  pull(Shared.IMGP) %>%
  strsplit("/",fixed = T) %>% unlist() %>% unique()


freq_res <- data.frame(IMGP = character(length(prognosticGP)),
                       Cell_interaction = character(length(prognosticGP)),
                       adjP_interaction = numeric(length(prognosticGP)),
                       adjP_survival = numeric(length(prognosticGP)),
                       Rank_enrich_lasso = numeric(length(prognosticGP)),
                       Freq_enrich_lasso = numeric(length(prognosticGP)),
                       Rank_prognostic_lasso = numeric(length(prognosticGP)),
                       Freq_prognostic_lasso = numeric(length(prognosticGP)))

for (i in 1:length(prognosticGP)) {
  cat("\r",i)
  IMGP <- prognosticGP[i]
  freq_res$IMGP[i] <- IMGP


  se1 <- grep(pattern = paste0("\\b",IMGP,"\\b"), x = enrich_res$Shared.IMGP)
  if(length(se1) == 0){
    freq_res$Cell_interaction[i] <- NA
    freq_res$adjP_interaction[i] <- NA
  } else{
    freq_res$Cell_interaction[i] <- paste0(unique(enrich_res$Cell.Interaction[se1]),
                                           collapse = ",")
    freq_res$adjP_interaction[i] <- paste0(unique(enrich_res$Adjust.P.Value[se1]),
                                           collapse = ",")
  }
  
  se2 <- which(rownames(cox_res) == IMGP)
  freq_res$adjP_survival[i] <- cox_res[se2,]$QV
  
  se3 <- which(rownames(Predictor.Prognostic) == IMGP)
  if(length(se3) == 0){
    freq_res$Freq_prognostic_lasso[i]  <- 0
    freq_res$Rank_prognostic_lasso[i]  <- NA

  } else{
    freq_res$Freq_prognostic_lasso[i]  <- Predictor.Prognostic[se3,]$Percentage
    freq_res$Rank_prognostic_lasso[i]  <- Predictor.Prognostic[se3,]$Rank
  }

  se4 <- which(rownames(Predictor.Enrich) == IMGP)
  if(length(se4) == 0){
    freq_res$Freq_enrich_lasso[i]  <- 0
    freq_res$Rank_enrich_lasso[i]  <- NA

  } else{
    freq_res$Freq_enrich_lasso[i]  <- Predictor.Enrich[se4,]$Percentage
    freq_res$Rank_enrich_lasso[i]  <- Predictor.Enrich[se4,]$Rank
  }
}


save(Predictor.Enrich,Predictor.Prognostic,freq_res, file = myoutf1)

# [1.3] Figure 4B ##############################################################
rm(list = ls())
myinf1  <- "~/Mypackage/MSofTimiGP/Fig4/LASSO_predictor.rda"
myoutf1  <- "~/Mypackage/MSofTimiGP/Fig4/Compare_LASSO.pdf"

load(myinf1)

library(dplyr)
library(ggplot2)
resdata <- freq_res
rownames(freq_res) <- freq_res$IMGP
colnames(resdata)

resdata$adjP_interaction <- resdata$adjP_interaction %>%
  strsplit(",") %>% lapply(as.numeric) %>% lapply(min) %>% as.numeric()
resdata$padj.rescaled <- ifelse(resdata$adjP_interaction > 0.05, 
                                0.05+resdata$adjP_interaction/20, 
                                resdata$adjP_interaction)

t.p <- nrow(resdata)
# prognostic pair； 14513
t.e <- resdata %>% filter(adjP_interaction < 0.05 ) %>% nrow()
# Enriched pair: 5766
both <- resdata %>% filter(Freq_enrich_lasso >0 & 
                             Freq_prognostic_lasso > 0) %>% nrow()
# both select: 544
s.p <- resdata %>% filter(Freq_prognostic_lasso > 0) %>% nrow()
# selected from prognostic IMGP: 2394
s.e <- resdata %>% filter(Freq_enrich_lasso > 0) %>% nrow()
# selected from enriched IMGP: 1124

h.e <- resdata %>% 
  filter(Freq_enrich_lasso >0 & Freq_prognostic_lasso > 0 ) %>%
  filter(Freq_enrich_lasso > Freq_prognostic_lasso)  %>% nrow()
# higher probability to be selected from enriched pair: 415

h.p <- resdata %>% 
  filter(Freq_enrich_lasso >0 & Freq_prognostic_lasso > 0 ) %>%
  filter(Freq_enrich_lasso < Freq_prognostic_lasso)  %>% nrow()

# higher probability to be selected from prognostic pair: 87

p <- ggplot() +
  geom_point(data = resdata,
             aes(x=Freq_prognostic_lasso, y=Freq_enrich_lasso,
                 color=adjP_interaction,size=adjP_survival),
             alpha=0.3) +
  scale_colour_gradientn(colours=c("#9970ab","#5aae61","#525252","#969696"),
                         breaks = c(0,0.01,0.05,0.1),
                         labels = c(0,0.01,0.05,1),
                         limits = c( 0,0.1),
                         guide = guide_colorbar(barwidth = 1, barheight =8,
                                                reverse=TRUE),
                         oob = scales::squish,
                         name = "Adjust\nP-Value\nof Cell\nInteraction") +
  scale_size_continuous(limits = c( 0,0.05),
                        name = "Adjust\nP-Value of\nPrognosis",
                        breaks = c(0,0.01,0.05),
                        labels = c(0,0.01,0.05),
                        range = c(7,2)) +
  geom_line(data = data.frame(x=c(0,100),
                              y=c(0,100)), aes(x=x,y=y), color="grey50", 
            lty=4, lwd=2) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x=paste0("Probability to be selected\n from Prognostic IMGP(%)\n","#",
                s.p,"/",t.p),
       y=paste0("Probability to be selected\n from Enriched IMGP(%)\n","#",
                s.e,"/",t.e),
       title=paste0(
         "Probability of selection by LASSO(#1000 iter) \n No. Selected by both: ",
                    both),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 1, colour = 1),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=15),
        legend.text= element_text(face="bold", color="black"
                                  family = "serif", size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=20),
        axis.text.y = element_text(face="bold", color="black", size=20),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))  +
  scale_y_continuous(limits=c(-3, 103), 
                     breaks=seq(0,100, 25),expand = c(0,0)) +
  scale_x_continuous(limits=c(-3, 103), 
                     breaks=seq(0,100, 25),expand = c(0,0)) +
  annotate("text", x = 30, y = 90,
           size = 8,color="black",
           family = "serif",
           label = paste0(h.e,
                          " IMGPs(",
                          round(h.e/both*100,1),
                          "%)")) +
  annotate("text", x = 75, y = 10,
           size = 8,color="black",
           family = "serif",
           label = paste0(h.p,
                          " IMGPs(",
                          round(h.p/both*100,1),
                          "%)"))

pdf(myoutf1, width = 9,height = 8)
print(p)
dev.off()


#[2] Examine single interaction TimiRS model(Related to Figure 4C,D) ###########
rm(list = ls())
myinf1 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda"
myinf3 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/cell.interaction.prognostic.value.rda"

# a) "Training" sets------------------------------------------------------------
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)


load(myinf1)
mps <- Immune3_MPS_SKCM06$Charoentong2017
all(rownames(info) == colnames(mps))

# b) validation sets -----------------------------------------------------------
load(myinf2)
data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
marker <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune %>%
  filter(Dataset == "Charoentong2017") %>% pull(Gene) %>% unique()

Dataset <- c("TCGA_SKCM06", # training
             "Cirenajwis_GSE65904",
             "Jonsson_GSE22155",
             "Jayawardana_GSE54467",
             "Mann_GSE53118",
             "Xu_GSE8401")

for (vDS in Dataset[-1]) {
  rna <- validation.dataset[[vDS]]$rna
  se <- which(rownames(rna) %in% marker)
  rna <- rna[se,]
  tmp1 <-  TimiGenePair(rna)
  tmp2 <- !tmp1
  rownames(tmp2) <- rownames(tmp1) %>% strsplit("_") %>%
    lapply(function(x) paste0(x[2],"_",x[1])) %>% unlist()
  validation.dataset[[vDS]]$mps<- rbind(tmp1,tmp2)
}

# c) result table (Related to Table S5)-----------------------------------------
load(myinf3)
resdata <- Immune3_enrich_SKCM06$Charoentong2017 %>% 
  filter(Adjust.P.Value < 0.05)


xx <- c(rep("Hazard.Ratio",length(Dataset)),
        rep("P.Value",     length(Dataset)),
        rep("Cindex",      length(Dataset)),
        rep("Cindex.se",   length(Dataset)))

res_col <- paste(Dataset,sep = "_", xx)
res_new <- matrix(nrow = nrow(resdata),ncol = length(res_col))
colnames(res_new) <- res_col
resdata <- data.frame(cbind(resdata,res_new),stringsAsFactors = F)
resdata$Miss.Predictor <- character(nrow(resdata))

for ( i in 1:nrow(resdata)){
  cat("\r",i)
  predictor <- resdata$Shared.IMGP[i]  %>%
    strsplit("/",fixed = T) %>% unlist() %>% unique()

  # Training
  # TimiRS=percentage of predictor=False
  TimiRS <- 100- colSums(mps[predictor,]*1)/ resdata$No.Shared.IMGP[i]*100
  if (all(names(TimiRS) == rownames(info)) == TRUE) {
    mycox <- coxph(Surv(info$time, info$event)~ as.numeric(TimiRS))
    mycox <- summary(mycox)
    resdata$TCGA_SKCM06_P.Value[i] <- mycox$coefficients[5]
    resdata$TCGA_SKCM06_Hazard.Ratio[i] <- mycox$conf.int %>% .[1]
    resdata$TCGA_SKCM06_Cindex[i] <- mycox$concordance[1]
    resdata$TCGA_SKCM06_Cindex.se[i] <- mycox$concordance[2]
    rm(mycox)
  } else {
    stop("TCGA_SKCM06 has different cohort ID in mps and info")
  }

  # Validation
  miss <- data.frame(row.names =Dataset[-1],
                     no=numeric(5),
                     mp=character(5))
  for (vDS in Dataset[-1]) {
    cat("\t",vDS)
    val.mps <- validation.dataset[[vDS]]$mps
    val.info <- validation.dataset[[vDS]]$info

    Com <- intersect(predictor,rownames(val.mps))
    val.rs <-  100- colSums(val.mps[Com,]*1)/ resdata$No.Shared.IMGP[i]*100

    if (all(names(val.mps) == rownames(val.info)) == TRUE) {
      mycox <- coxph(Surv(val.info$t.surv, val.info$e.surv)~ as.numeric(val.rs))
      mycox <- summary(mycox)
      resdata[i,paste0(vDS,"_P.Value")] <- mycox$coefficients[5]
      resdata[i,paste0(vDS,"_Hazard.Ratio")] <- mycox$conf.int %>% .[1]
      resdata[i,paste0(vDS,"_Cindex")] <- mycox$concordance[1]
      resdata[i,paste0(vDS,"_Cindex.se")] <- mycox$concordance[2]
      rm(mycox)
    } else {
      stop(vDS," has different cohort ID in mps and info")
    }

    # summarize miss.predictor

    miss[vDS,]$no <- length(predictor) - length(Com)

    miss[vDS,]$mp <- paste(predictor[! predictor %in% Com],collapse = "/")

    rm(val.info,val.mps,val.rs)
  }
  resdata$Miss.Predictor[i] <- paste0(rownames(miss),"(#",miss$no,":",miss$mp,")") %>%
    paste(collapse = ";")
}

save(resdata,file = myoutf1)

# d) Figure 4C,D ---------------------------------------------------------------

rm(list = ls())
library(RColorBrewer)
library(VennDiagram)
library(grid)
library(ggplot2)
library(dplyr)
myinf1 <- "~/Mypackage/MSofTimiGP/Fig4/cell.interaction.prognostic.value.rda"

myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/cell.interaction.prognostic.value.pdf"

load(myinf1)

Dataset<- c("TCGA_SKCM06","Cirenajwis_GSE65904","Jonsson_GSE22155",
            "Jayawardana_GSE54467","Mann_GSE53118","Xu_GSE8401")
sum <- data.frame("DS" = Dataset,
                  "Freq" = numeric(length(Dataset)),
                  "Percentage" = numeric(length(Dataset)))

for (i in 1:nrow(sum)) {
  ds <- sum$DS[i]
  xx <- paste0(ds,"_P.Value")
  sum$Freq[i] <- sum(resdata[,xx] <0.05)
  cat("\n",ds,"\t",sum(paste0(ds,"_Hazard.Ratio") < 1))
}

sum$Percentage <- sum$Freq/nrow(resdata) *100
sum <- sum %>% arrange(-Percentage) %>% mutate(DS=factor(DS,levels = Dataset))

listInput <- list()
for (i in Dataset) {
  xx <- paste0(i,"_P.Value")
  listInput[[i]] <- resdata$Cell.Interaction[which(resdata[,xx] <0.05)]

}


length(unique(unlist(listInput))) #141
length(unique(unlist(listInput)))/nrow(resdata)*100 #100%

sum <- rbind(sum,data.frame(
  "DS" = "At_least_1_set",
  "Freq" = length(unique(unlist(listInput))),
  "Percentage" = length(unique(unlist(listInput)))/nrow(resdata)*100)
) %>% arrange(-Percentage) %>% 
  mutate(DS=factor(DS,
                   levels = c("TCGA_SKCM06","At_least_1_set",
                              "Cirenajwis_GSE65904","Jonsson_GSE22155",
                              "Jayawardana_GSE54467","Mann_GSE53118",
                              "Xu_GSE8401")))
# Figure 4C
p <- ggplot() +
  geom_bar(data = sum,
           mapping = aes(x = DS,y=Percentage,fill = DS),
           stat="identity", alpha=0.5) +
  scale_fill_manual(values=c(brewer.pal(9,"Greys")[8],
                             brewer.pal(9,"Set1")[c(9,1,5,3,2,4)]))+
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Dataset",y="Percentage(%)",
       title=paste0("Prognostic Value of\n",nrow(resdata)," Functional Interaction"),
       base_size = 22, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=13,
                                   angle = 60,  vjust = 1, hjust=1),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=17),
        axis.title.y = element_text(face="bold",color="black", size=17),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))) +
  scale_y_continuous(limits=c(0, 105), breaks=seq(0, 100, 20),expand = c(0,0)) +
  guides(fill = "none")


# Figure 4D
col<-brewer.pal(9,"Set1")[c(1,5,3,2,4)]
venn.plot <- venn.diagram(
  listInput[-1],
  #  filename = "venn_lasso_CJJMX>5.png", imagetype = "png",
  filename = NULL,
  lwd = 3,
  col="white",
  main = paste0( "Prognostic Cell Functional Interaction"),
  main.cex = 2,
  scaled=T,
  fill = col,
  alpha = 0.5,
  label.col = "white",
  cex = 1.7,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = col,
  cat.cex = 1.7,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.2,
  cat.dist = c(0.2, 0.3,0.3, 0.25,0.25),
  cat.pos = c(0,-40,-120,145,25))

# The 6 interaction:
xx <- listInput[[1]]
for(i in 2:length(listInput)){
  yy <- listInput[[i]]
  xx <- intersect(xx,yy)
}
xx
# [1] "CD8 Tem_Neutrophil" "CD8 Tem_Tgd"        "CD8 Tem_CD8 Tcm"   
# "CD8 Tem_aDC"        "Th1_Tgd"            "CD8 Tem_CD4 Tem"   
# [7] "aCD8 T_Neutrophil"  "aCD8 T_CD56dim NK"  "CD8 Tem_iDC"      
# "CD8 Tem_Macrophage" "CD8 Tem_NKT"        "CD8 Tem_Th17"      
# [13] "iB_Th2"             "aB_Tgd"             "aB_Neutrophil"      "Tfh_Tgd" 
# "CD8 Tem_Eosinophil" "aCD8 T_Macrophage" 
# [19] "MDSC_NK"

# export
pdf(myoutf1,width = 7,height = 7)
print(p)

grid.newpage()
gl <- grid.layout(nrow=3, ncol=3, widths = c(0.1, 1, 0.1), heights = c(0.1, 1, 0.1))
vp <- viewport(layout.pos.col=2, layout.pos.row=2)
pushViewport(viewport(layout=gl))
pushViewport(vp)
grid.draw(venn.plot)
dev.off()

 
#[3] Develop TimiGP model based on TimiGP CD8 Tem ##############################
#[3.1]selected interactions(related to figure S6)===============================
rm(list = ls())

library(TimiGP)

myinf1 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"


myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/selected_interactions.pdf"

load(myinf1)
resdata <- Immune3_enrich_SKCM06$Charoentong2017 %>% 
  mutate(label = paste0("#",Rank))
select <-  resdata %>% 
  filter(Adjust.P.Value < 1e-4) %>%
  filter(Favorable.Cell.Type == "CD8 Tem" ) %>%
  pull(Index)

p <- TimiDotplot(resdata = resdata,select = select) +
  geom_text(data = resdata[select,],
            aes(x=Cell.Interaction, y=Enrichment.Ratio,label = label),
            vjust = -1.5)+
  scale_x_discrete(limits = resdata[select,]$Cell.Interaction,
                   labels = resdata[select,]$Unfavorable.Cell.Type,
                   name = " CD8 Tem -> Other Cells Interactions")
  
pdf(myoutf1, width = 8,height = 6)
print(p)
dev.off()
#[3.2]  LASSO for selected functional interactions =============================
rm(list = ls())

library(TimiGP)

library(glmnet)
library(doMC)
registerDoMC(cores = 50)
myinf1 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig3/Immune3_enrich_SKCM06.rda"


myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/lasso_res_selected.rda"
myoutf3 <- "~/Mypackage/MSofTimiGP/Fig4/selected_interactions.pdf"

load(myinf1)
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)
all(rownames(info) == colnames(Immune3_MPS_SKCM06$Charoentong2017))

# select interaction
load(myinf2)
selectGP <- Immune3_enrich_SKCM06$Charoentong2017 %>% 
  filter(Adjust.P.Value < 1e-4) %>%
  filter(Favorable.Cell.Type == "CD8 Tem" ) %>% 
  pull(Shared.IMGP) %>% 
  strsplit("/",fixed = T) %>% unlist() %>% unique()

mps <- Immune3_MPS_SKCM06$Charoentong2017
mps_select <- mps[selectGP,]


iter <- 1000
cv.lassoModel <- list()
lasso_res <- data.frame( "Seed" = numeric(length = iter),
                         "No.Predictor.select" = numeric(length = iter),
                         "Predictor.select" = character(length = iter),
                         "Cohort" = character(length = iter))


for (i in 1:iter) {
  cat("\n", i)
  
  # a) Random sampling 80% cohorts----------------------------------------------
  lasso_res$Seed[i] <- i
  set.seed(i)
  xx <- sample(x = 1:nrow(info), size = round(0.8*nrow(info)), replace = F)
  cohort <- rownames(info)[xx]
  
  lasso_res$Cohort[i] <- paste0(cohort,collapse = "/")
  
  # b) LASSO from selected pairs------------------------------------------------
  cv.lassoModel[[i]] <- cv.glmnet(
    x=t(mps_select[,xx]),
    y=Surv(info$time[xx], info$event[xx]),
    standardize=F,
    alpha=1.0,
    nfolds=10, 
    family="cox",
    type.measure = "C",
    parallel=TRUE)
  
  lambda <- cv.lassoModel[[i]] $lambda.1se
  cv.lassoModel[[i]] $co <- coef(object = cv.lassoModel[[i]] , s=lambda, 
                                 exact=TRUE,
                                 x=t(mps_select[,xx]),
                                 y=Surv(info$time[xx], info$event[xx]))
  predictor <- 
    rownames(cv.lassoModel[[i]] $co)[which(cv.lassoModel[[i]] $co != 0)] 
  lasso_res$No.Predictor.select[i] <- length(predictor)
  lasso_res$Predictor.select[i] <- paste0(predictor, collapse = "/")
  rm(lambda)
  rm(predictor)
}


lasso_res_selected <- lasso_res
cv.lassoModel_selected <- cv.lassoModel
save(lasso_res_selected, file = myoutf1)

#[3.3] select num(21) of predictor (related to Figure 4E-F,S7-8 Table S6)=======
rm(list = ls())
myinf1 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda"
myinf3 <-  "~/Mypackage/MSofTimiGP/Fig4/lasso_res_selected.rda"

myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/TimiGP_model.rda"

# a) "Training" set ------------------------------------------------------------
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)


load(myinf1)
mps <- Immune3_MPS_SKCM06$Charoentong2017
all(rownames(info) == colnames(mps))

# b) validation set-------------------------------------------------------------
load(myinf2)
data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
marker <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune %>%
  filter(Dataset == "Charoentong2017") %>% pull(Gene) %>% unique()

Dataset <- c("TCGA_SKCM06", # training
             "Cirenajwis_GSE65904",
             "Jonsson_GSE22155",
             "Jayawardana_GSE54467",
             "Mann_GSE53118",
             "Xu_GSE8401")

for (vDS in Dataset[-1]) {
  rna <- validation.dataset[[vDS]]$rna
  se <- which(rownames(rna) %in% marker)
  rna <- rna[se,]
  tmp1 <-  TimiGenePair(rna)
  tmp2 <- !tmp1
  rownames(tmp2) <- rownames(tmp1) %>% strsplit("_") %>%
    lapply(function(x) paste0(x[2],"_",x[1])) %>% unlist()
  validation.dataset[[vDS]]$mps<- rbind(tmp1,tmp2)
}



# c) predictor candidate--------------------------------------------------------
load(myinf3)
Predictor.select <- lasso_res_selected$Predictor.select %>%
  strsplit(split = "/") %>%
  unlist() %>%
  table() %>%
  data.frame(row.names = 1)
names(Predictor.select) <- c( "Frequency")
Predictor.select$Percentage <- Predictor.select$Frequency/1000*100
Predictor.select <- Predictor.select[order(Predictor.select$Frequency, 
                                           decreasing = T), ]
Predictor.select$Rank <- rank(-Predictor.select$Frequency,ties.method="min")

reference <- Predictor.select %>% filter(Predictor.select$Percentage > 20)


# d) result table for CD8 Tem model --------------------------------------------
nmin <- 10
nmax <- nrow(reference)
xx <- c(rep("Hazard.Ratio",length(Dataset)),
        rep("P.Value",     length(Dataset)),
        rep("Cindex",      length(Dataset)),
        rep("Cindex.se",   length(Dataset)))

res_col <- paste(Dataset,sep = "_", xx)
resdata <- matrix(nrow = nmax-nmin+1,ncol = length(res_col))
colnames(resdata) <- res_col
rownames(resdata) <- nmin:nmax
resdata <- data.frame(resdata)
resdata$Predictor <- character(nrow(resdata))
resdata$Miss.Predictor <- character(nrow(resdata))

for ( i in 1:nrow(resdata)){
  cat("\r",i)
  predictor <- reference %>% arrange(-Percentage) %>%
    row.names() %>% .[1:rownames(resdata)[i]]

  resdata$Predictor[i] <- paste0(predictor,collapse = "/")
  # Training
  # TimiRS=percentage of predictor=False
  
  se <- which(rownames(mps) %in% predictor)
  TimiRS <- 100- colSums(mps[se,]*1)/ length(predictor)*100
  
  if (all(names(TimiRS) == rownames(info)) == TRUE) {
    mycox <- coxph(Surv(info$time, info$event)~ as.numeric(TimiRS))
    mycox <- summary(mycox)
    resdata$TCGA_SKCM06_P.Value[i] <- mycox$coefficients[5]
    resdata$TCGA_SKCM06_Hazard.Ratio[i] <- mycox$conf.int %>% .[1]
    resdata$TCGA_SKCM06_Cindex[i] <- mycox$concordance[1]
    resdata$TCGA_SKCM06_Cindex.se[i] <- mycox$concordance[2]
    rm(mycox)
  } else {
    stop("TCGA_SKCM06 has different cohort ID in mps and info")
  }

  # Validation
  miss <- data.frame(row.names =Dataset[-1],
                     no=numeric(5),
                     mp=character(5))
  for (vDS in Dataset[-1]) {
    cat("\t",vDS)
    val.mps <- validation.dataset[[vDS]]$mps
    val.info <- validation.dataset[[vDS]]$info
    
    Com <- intersect(predictor,rownames(val.mps))
    
    val.se <- which(rownames(val.mps) %in% predictor)
    val.rs <-  100- colSums(val.mps[val.se,]*1)/ length(predictor)*100
    
    if (all(names(val.mps) == rownames(val.info)) == TRUE) {
      mycox <- coxph(Surv(val.info$t.surv, val.info$e.surv)~ as.numeric(val.rs))
      mycox <- summary(mycox)
      resdata[i,paste0(vDS,"_P.Value")] <- mycox$coefficients[5]
      resdata[i,paste0(vDS,"_Hazard.Ratio")] <- mycox$conf.int %>% .[1]
      resdata[i,paste0(vDS,"_Cindex")] <- mycox$concordance[1]
      resdata[i,paste0(vDS,"_Cindex.se")] <- mycox$concordance[2]
      rm(mycox)
    } else {
      stop(vDS," has different cohort ID in mps and info")
    }
    
    # summarize miss.predictor
    
    miss[vDS,]$no <- length(predictor) - length(Com)
    
    miss[vDS,]$mp <- paste(predictor[! predictor %in% Com],collapse = "/")
    
    rm(val.info,val.mps,val.rs)
  }
  resdata$Miss.Predictor[i] <- paste0(rownames(miss),"(#",miss$no,":",miss$mp,")") %>%
    paste(collapse = ";")

}

resdata.select <- resdata

save(resdata.select,file = myoutf1)


#[3.4] Develop control model with features selected from prognostic IMGP========
rm(list = ls())
myinf1 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda"
myinf3 <- "~/Mypackage/MSofTimiGP/Fig4/LASSO_predictor.rda"
myoutf1 <- "~/Mypackage/MSofTimiGP/Fig4/control.rda"

# a) "Training" ----------------------------------------------------------------
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)


load(myinf1)
mps <- Immune3_MPS_SKCM06$Charoentong2017
all(rownames(info) == colnames(mps))

# b) validation-----------------------------------------------------------------
load(myinf2)
data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
marker <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune %>%filter(Dataset == "Charoentong2017") %>% pull(Gene) %>% unique()

Dataset <- c("TCGA_SKCM06", # training
             "Cirenajwis_GSE65904",
             "Jonsson_GSE22155",
             "Jayawardana_GSE54467",
             "Mann_GSE53118",
             "Xu_GSE8401")

for (vDS in Dataset[-1]) {
  rna <- validation.dataset[[vDS]]$rna
  se <- which(rownames(rna) %in% marker)
  rna <- rna[se,]
  tmp1 <-  TimiGenePair(rna)
  tmp2 <- !tmp1
  rownames(tmp2) <- rownames(tmp1) %>% strsplit("_") %>%
    lapply(function(x) paste0(x[2],"_",x[1])) %>% unlist()
  validation.dataset[[vDS]]$mps<- rbind(tmp1,tmp2)
}



# c) predictor candidate--------------------------------------------------------
load(myinf3)

reference <- Predictor.Prognostic %>% 
  filter(Predictor.Prognostic$Percentage > 20)

# d) result table for control model --------------------------------------------
nmin <- 10
nmax <- min(nrow(reference),25)
xx <- c(rep("Hazard.Ratio",length(Dataset)),
        rep("P.Value",     length(Dataset)),
        rep("Cindex",      length(Dataset)),
        rep("Cindex.se",   length(Dataset)))

res_col <- paste(Dataset,sep = "_", xx)
resdata <- matrix(nrow = nmax-nmin+1,ncol = length(res_col))
colnames(resdata) <- res_col
rownames(resdata) <- nmin:nmax
resdata <- data.frame(resdata)
resdata$Predictor <- character(nrow(resdata))
resdata$Miss.Predictor <- character(nrow(resdata))

for ( i in 1:nrow(resdata)){
  cat("\r",i)
  predictor <- reference %>% arrange(-Percentage) %>%
    row.names() %>% .[1:rownames(resdata)[i]]
  
  resdata$Predictor[i] <- paste0(predictor,collapse = "/")
  # Training
  # TimiRS=percentage of predictor=False
  
  se <- which(rownames(mps) %in% predictor)
  TimiRS <- 100- colSums(mps[se,]*1)/ length(predictor)*100
  
  if (all(names(TimiRS) == rownames(info)) == TRUE) {
    mycox <- coxph(Surv(info$time, info$event)~ as.numeric(TimiRS))
    mycox <- summary(mycox)
    resdata$TCGA_SKCM06_P.Value[i] <- mycox$coefficients[5]
    resdata$TCGA_SKCM06_Hazard.Ratio[i] <- mycox$conf.int %>% .[1]
    resdata$TCGA_SKCM06_Cindex[i] <- mycox$concordance[1]
    resdata$TCGA_SKCM06_Cindex.se[i] <- mycox$concordance[2]
    rm(mycox)
  } else {
    stop("TCGA_SKCM06 has different cohort ID in mps and info")
  }
  
  # Validation
  miss <- data.frame(row.names =Dataset[-1],
                     no=numeric(5),
                     mp=character(5))
  for (vDS in Dataset[-1]) {
    cat("\t",vDS)
    val.mps <- validation.dataset[[vDS]]$mps
    val.info <- validation.dataset[[vDS]]$info
    
    Com <- intersect(predictor,rownames(val.mps))
    
    val.se <- which(rownames(val.mps) %in% predictor)
    val.rs <-  100- colSums(val.mps[val.se,]*1)/ length(predictor)*100
    
    if (all(names(val.mps) == rownames(val.info)) == TRUE) {
      mycox <- coxph(Surv(val.info$t.surv, val.info$e.surv)~ as.numeric(val.rs))
      mycox <- summary(mycox)
      resdata[i,paste0(vDS,"_P.Value")] <- mycox$coefficients[5]
      resdata[i,paste0(vDS,"_Hazard.Ratio")] <- mycox$conf.int %>% .[1]
      resdata[i,paste0(vDS,"_Cindex")] <- mycox$concordance[1]
      resdata[i,paste0(vDS,"_Cindex.se")] <- mycox$concordance[2]
      rm(mycox)
    } else {
      stop(vDS," has different cohort ID in mps and info")
    }
    
    # summarize miss.predictor
    
    miss[vDS,]$no <- length(predictor) - length(Com)
    
    miss[vDS,]$mp <- paste(predictor[! predictor %in% Com],collapse = "/")
    
    rm(val.info,val.mps,val.rs)
  }
  resdata$Miss.Predictor[i] <- paste0(rownames(miss),"(#",miss$no,":",miss$mp,")") %>%
    paste(collapse = ";")
  
}

resdata.prognostic <- resdata

save(resdata.prognostic,file = myoutf1)

#[3.4] Compare TimiGP model with control model =================================
# Both use top 21 IMGPs

rm(list = ls())
library(grid)
library(survival)
library(data.table) 
library(ggplot2)
library(survminer)
library(RColorBrewer)
library(survivalROC)
library(dplyr)
library(tidyr)
library(gridExtra)

myinf1 <- "~/Mypackage/MSofTimiGP/Fig4/TimiGP_model.rda"
myinf2 <- "~/Mypackage/MSofTimiGP/Fig4/control.rda"
myinf3 <- "~/Mypackage/MSofTimiGP/Fig4/LASSO_predictor.rda"
myinf4 <-"~/Mypackage/MSofTimiGP/Fig4/Immune3_MPS_SKCM06.rda"
myinf5 <- "~/Mypackage/MSofTimiGP/Fig4/validation.dataset.rda"

myoutf1 <-  "~/Mypackage/MSofTimiGP/Fig4/KM.pdf"
myoutf2 <-  "~/Mypackage/MSofTimiGP/Fig4/ROC.pdf"
myoutf3 <-  "~/Mypackage/MSofTimiGP/Fig4/Cindex.pdf"


mydir <- "~/Mypackage/MSofTimiGP/Fig4/network"
dir.create(mydir)
load(myinf1)
load(myinf2)
load(myinf3)
res <- rbind(resdata.prognostic["21",],resdata.select["21",])
rownames(res) <- res$Model <- c("Control","TimiGP")

# a) Network (related to Figure 4E) --------------------------------------------
xx <- res["TimiGP","Predictor"] %>% strsplit("/") %>% unlist()
se <- which(freq_res$IMGP %in% xx)
pre.net <- freq_res[se,]
rownames(pre.net) <- pre.net$IMGP
data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
geneset <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune %>% 
  filter(Dataset == "Charoentong2017")
NET <- TimiGeneNetwork(resdata = pre.net,select = rownames(pre.net),
                dataset = "Other",geneset = geneset,export = T,path = mydir )

# b) KM plot & ROC (related to figure 4F, S7,S8)----



# "Training" 
data("SKCM06info")
info <- TimiCheckEvent(SKCM06info)


load(myinf4)
mps <- Immune3_MPS_SKCM06$Charoentong2017
all(rownames(info) == colnames(mps))

# validation
load(myinf5)

data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
marker <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune %>%
  filter(Dataset == "Charoentong2017") %>% pull(Gene) %>% unique()

Dataset <- c("TCGA_SKCM06", # training
             "Cirenajwis_GSE65904",
             "Jonsson_GSE22155",
             "Jayawardana_GSE54467",
             "Mann_GSE53118",
             "Xu_GSE8401")

for (vDS in Dataset[-1]) {
  rna <- validation.dataset[[vDS]]$rna
  se <- which(rownames(rna) %in% marker)
  rna <- rna[se,]
  tmp1 <-  TimiGenePair(rna)
  tmp2 <- !tmp1
  rownames(tmp2) <- rownames(tmp1) %>% strsplit("_") %>%
    lapply(function(x) paste0(x[2],"_",x[1])) %>% unlist()
  validation.dataset[[vDS]]$mps<- rbind(tmp1,tmp2)
}

p.sur <- list()

roc_ext <- data.frame()

for ( i in 1:nrow(res)){
  cat("\r",i)
  predictor <- res[i,"Predictor"]  %>% strsplit("/") %>% unlist()
  # Training
  # TimiRS=percentage of predictor=False
  
  se <- which(rownames(mps) %in% predictor)
  TimiRS <- 100- colSums(mps[se,]*1)/ length(predictor)*100
  
  if (all(names(TimiRS) == rownames(info)) == TRUE) {
    # KM plot
    kmData <- info
    cutoff <- median(TimiRS)
    kmData$Risk <- as.factor(ifelse(test =  TimiRS > cutoff, 
                                    "high", "low"))
    ta <- table(kmData$Risk)
    km_fit<- survfit(Surv(time, event) ~Risk, data = kmData)
    
    p.sur[[rownames(res)[i]]][[Dataset[1]]] <- ggsurvplot(
      km_fit,                     # survfit object with calculated statistics.
      pval = TRUE, 
      data = kmData,             # data used to fit survival curves.
      palette = c("red", "blue"),
      font.legend=25,
      legend = "top",
      legend.title="TimiRS",
      xlab = "Years",   # customize X axis label.
      break.time.by = 365.25 *2,     # break X axis in time intervals by days.
      surv.median.line = "hv",  # add the median survival pointer.
      legend.labs = c(paste0("High (n=", ta[1], ")"),
                      paste0("Low (n=", ta[2], ")")),
      
      size=3,
      pval.method=T,
      pval.size=6,
      Exp.table = TRUE,
      title = paste0(rownames(res)[i],"_",Dataset[1]),
      ggtheme = theme_classic2(base_size=25, base_family = "serif"),
      font_family= "serif",
      xscale ="d_y")+
      guides(colour = guide_legend(nrow = 2))
    
    
    # Time-dependent ROC (3year)
    roc<- survivalROC(Stime        = info$time,
                      status       = info$event,
                      marker       = TimiRS,
                      predict.time = round(365.25*3,0),
                      method       = "NNE",
                      span = 0.25 * length(predictor)^(-0.20))
    roc_tmp <- as.data.frame(roc[1:3])
    roc_tmp$Dataset <- Dataset[1]
    roc_tmp$AUC <- roc$AUC
    roc_tmp$Model <- rownames(res)[i]
    
    roc_ext <- rbind(roc_ext,roc_tmp)
    
    rm(kmData,km_fit,TimiRS,roc,roc_tmp)
  } else {
    stop("TCGA_SKCM06 has different cohort ID in mps and info")
  }
  
  # Validation
  
  for (vDS in Dataset[-1]) {
    cat("\t",vDS)
    val.mps <- validation.dataset[[vDS]]$mps
    val.info <- validation.dataset[[vDS]]$info
    
    Com <- intersect(predictor,rownames(val.mps))
    
    val.se <- which(rownames(val.mps) %in% predictor)
    val.rs <-  100- colSums(val.mps[val.se,]*1)/ length(predictor)*100
    
    if (all(names(val.mps) == rownames(val.info)) == TRUE) {
      # KM plot
      kmData <- val.info
      cutoff <- median(val.rs)
      kmData$Risk <- as.factor(ifelse(test =  val.rs > cutoff, 
                                      "high", "low"))
      ta <- table(kmData$Risk)
      km_fit<- survfit(Surv(t.surv, e.surv) ~Risk, data = kmData)
      
      p.sur[[rownames(res)[i]]][[vDS]] <- ggsurvplot(
        km_fit,                     # survfit object with calculated statistics.
        pval = TRUE, 
        data = kmData,             # data used to fit survival curves.
        palette = c("red", "blue"),
        font.legend=25,
        legend = "top",
        legend.title="TimiRS",
        xlab = "Years",   # customize X axis label.
        break.time.by = 365.25 *2,     # break X axis in time intervals by days.
        surv.median.line = "hv",  # add the median survival pointer.
        legend.labs = c(paste0("High (n=", ta[1], ")"),
                        paste0("Low (n=", ta[2], ")")),
        
        size=3,
        pval.method=T,
        pval.size=6,
        Exp.table = TRUE,
        title = paste0(rownames(res)[i],"_",vDS),
        ggtheme = theme_classic2(base_size=25, base_family = "serif"),
        font_family= "serif",
        xscale ="d_y")+
        guides(colour = guide_legend(nrow = 2))
      
      # Time-dependent ROC (3year)
      roc<- survivalROC(Stime        = val.info$t.surv,
                        status       = val.info$e.surv,
                        marker       = val.rs,
                        predict.time = round(365.25*3,0),
                        method       = "NNE",
                        span = 0.25 * length(predictor)^(-0.20))
      roc_tmp <- as.data.frame(roc[1:3])
      roc_tmp$Dataset <- vDS
      roc_tmp$AUC <- roc$AUC
      roc_tmp$Model <- rownames(res)[i]
      roc_ext <- rbind(roc_ext,roc_tmp)
      
      rm(kmData,km_fit,val.rs,roc,roc_tmp)
    } else {
      stop(vDS," has different cohort ID in mps and info")
    }
    
    
  }
}

#plot ROC
p.roc <- list()
for ( i in 1:nrow(res)){
  line1 <- data.frame(x = c(0,1), y=c(0,1))
  auc <- roc_ext %>% filter(Model ==  rownames(res)[i]) %>% select(Dataset,AUC) %>% unique()
  p.roc[[ rownames(res)[i]]] <- ggplot() +
    geom_line(data = roc_ext[which(roc_ext$Model ==  rownames(res)[i]),], 
              mapping = aes(x = FP, y = TP, group = Dataset,color = Dataset), 
              size = 2, alpha=0.8) +
    scale_color_manual(name = "Dataset", 
                       values = c(brewer.pal(9,"Greys")[8],brewer.pal(9,"Set1")[c(1,5,3,2,4)]),
                       limits = Dataset,
                       labels =paste0(auc$Dataset,"(AUC = ", round(auc$AUC,3),")") )+
    geom_line(data = line1,
              mapping = aes(x = x, y = y), color = "grey50", size = 2,lty=4) +
    theme_bw(base_size = 20, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="False Positive Rate",y="True Positive Rate",
         title=paste0("Time-dependent ROC(3-year,", rownames(res)[i],")"),
         base_size = 22, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, size = 1, colour = 1),
          legend.position=c(0.68,0.2),
          legend.box = "vertical",
          legend.direction= "vertical",
          panel.grid=element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black", size=20),
          axis.text.y = element_text(face="bold", color="black", size=20),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))  +
    scale_y_continuous(limits=c(-0.01, 1.01), breaks=seq(0,1, 0.2),expand = c(0,0)) +
    scale_x_continuous(limits=c(-0.01, 1.01), breaks=seq(0, 1,0.2),expand = c(0,0))
}


# plot 

p.list1 <- list()
p.list2 <- list()

for ( i in 1:6){
  p.list1[[i]] <- p.sur[[1]][[i]]$plot
  p.list2[[i]] <- p.sur[[2]][[i]]$plot
}


pdf(myoutf1,width = 20,height = 15)
do.call("grid.arrange", c(plotlist = p.list1, ncol=3,top="Control"))
do.call("grid.arrange", c(plotlist = p.list2, ncol=3,top="TimiGP"))
dev.off() 

pdf(myoutf2,width = 20,height = 10)
do.call("grid.arrange", c(plotlist = p.roc, ncol=2))
dev.off() 
# c) C-index (related to figure 4F, S7)-----------------------------------------
colnames(res[13:18])

sum <- res[c(13:18,27)]
colnames(sum) <- colnames(sum) %>% strsplit("_C",fixed = T) %>% 
  lapply("[[", 1) %>% unlist()

sum <- pivot_longer(sum, cols=1:6, names_to = "DS", values_to = "Cindex") %>%
  data.frame() %>% 
  mutate(DS = factor(DS,levels = Dataset))
p.cindex <- list()
for (i in unique(sum$Model)) {
  p.cindex[[i]] <- sum %>% filter(Model == i) %>%
    ggplot() +
    geom_bar(
      mapping = aes(x = DS,y=Cindex,fill = DS),
      position="dodge", stat="identity", alpha=0.5) +
    scale_fill_manual(values=c(brewer.pal(9,"Greys")[8],
                               brewer.pal(9,"Set1")[c(1,5,3,2,4)]))+
    geom_hline(yintercept = 0.5,lty = 4, color = "grey",size =1.5) +
    theme_bw(base_size = 20, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Dataset",y="C-index",
         title=paste0(i," Model"),
         base_size = 22, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, 
                                           size = 0.5, colour = 1),
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
          axis.title.y = element_text(face="bold",color="black", size=17),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                         ends = "both"))) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2),expand = c(0,0)) +
    guides(fill = "none") 
  
}


pdf(myoutf3,width = 15,height = 7)
do.call("grid.arrange", c(plotlist = p.cindex, ncol=2))
dev.off() 
