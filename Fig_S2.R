# Create a directory to save figures and tables
mainDir <- "../Figures"
subDir <- "Figs_S2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_S2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")

library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(stringr)
library(ggplot2)
library(gtools)
library(ggrepel)
library(RColorBrewer)
library(circlize)



########################******************Fig S_2A  ***************################################
# HUMANN Pathways

# Import lme results from Metaphlan
lme_species <-  read.csv("../Tables/Tabs_1/Tab_LME_Species_1F.csv")
sig_species <-  unique(lme_species$Bac[lme_species$p_adj < 0.05])

search_str <-  paste0(paste0("s__",sig_species),collapse = "|")

# Import pathways:
path_dat <-  read.csv("../data/humann/humann3_pathabundance.txt",sep="\t")
names(path_dat)[1] <- "Pathway"

# Now remove the Unmapped and Unintegrated:
path_dat <- path_dat[!grepl("UNMAPPED|UNINTEGRATED",path_dat$Pathway),]
path_sp_dat <-  path_dat[grep(search_str,path_dat$Pathway),]
sel_pathways <- unique(gsub("\\|.*","",path_sp_dat$Pathway))
# Select only at the community level
path_sel_dat <- path_dat[!grepl("\\|",path_dat$Pathway),]

# Select only at the community level
path_sel_dat <- path_sel_dat[path_sel_dat$Pathway %in% sel_pathways,]
path_sel_dat$Pathway <- gsub("'","",path_sel_dat$Pathway)


# Import Pathways:
phy_path <-  readRDS("../Processed_data/phy_raw_path_mdr.rds")

phy_path <- prune_taxa(path_sel_dat$Pathway,phy_path)
sm_dt <- data.frame(sample_data(phy_path))
# Change the level of last time point as TC
levels(sm_dt$Visit)[6] <- "TC"
sample_data(phy_path) <-  sample_data(sm_dt)
phy_lme <- phy_path
# LME analysis:
# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_lme),
               MARGIN = ifelse(taxa_are_rows(phy_lme), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_lme))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.1
prev_cutoff <- prev_frac*nsamples(phy_lme) # Cut-off
abun_cutoff <- 0.00001
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]


phy_fil <-  prune_taxa(as.character(prevdf_fil$OTU),phy_lme)

melt_dt <-  psmelt(phy_fil)
taxa <- taxa_names(phy_fil)[1]
result_lme <- list()
for (taxa in taxa_names(phy_fil)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$PID <- factor(phy_m$PID)
  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$t_Abundance  <- log(phy_m$Abundance + 1)
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Sex + Age + Visit  ,random = ~1 |PID, phy_m)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Bac <- taxa
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only Visit
result_lme_dt <- result_lme_dt[grep("Visit",result_lme_dt$Var),]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]


p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]

# Save the LME results
write.csv(result_sig,paste(tab_folder,"Tab_Sig_LME_MDR.csv",sep = "/"))



unique(result_sig$Bac)
library(ggplot2)
ggplot(data = result_sig, aes(x = result_sig$Value)) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


result_sig <-  result_sig[result_sig$Value < -0.5,]
unique(result_sig$Bac)

# freq_dt <- stack(table(result_sig$Bac))
# freq_dt <-  freq_dt[freq_dt$values > 1,]
# result_sig <- result_sig[result_sig$Bac %in% freq_dt$ind,]

## Annotation of metacyc pathways:
ann_path <-  read.csv("../data/humann/map_metacyc-pwy_lineage.tsv",sep = "\t",header = F)
ann_path <-  unique(ann_path)
library(stringr)
sep_dt <-  data.frame(str_split_fixed(ann_path$V2,"\\|",n = 4))
ann_path <-  data.frame(ann_path,sep_dt)

sig_path <-  data.frame(Pathway = unique(result_sig$Bac),
                        ID = gsub(":.*","",unique(result_sig$Bac)))

table(sig_path$ID  %in% unique(ann_path$V1))


sig_path$ID[!sig_path$ID  %in% unique(ann_path$V1)]


ann_path <-  ann_path[ann_path$V1 %in% sig_path$ID,]
#ann_path <-  ann_path[!grepl("Super-Pathways",ann_path$V2),]

# Also remove the duplicates by choosing size of description:
ann_path$size <- nchar(ann_path$V2)
library(dplyr)
ann_path <- ann_path %>%
  group_by(V1) %>%
  slice_max(size) %>%
  ungroup %>% 
  group_by(V1) %>%
  sample_n(1) %>% data.frame()


sig_path <-  sig_path %>%
  left_join(ann_path,by  =  c("ID"="V1"))


# ann_path[grep("GLYOXYLATE",ann_path$V1)]
# tmp1 <- stack(table(ann_path$V1))
# tmp1 <- tmp1[tmp1$values > 1,]
# tmp2 <-  ann_path[ann_path$V1 %in% as.character(tmp1$ind),]


# Clean up the names 
sig_path <-  sig_path[,c("Pathway","ID","V2","X1","X2","X3")]
names(sig_path)[3:6] <-  c("Description","Class_1","Class_2","Class_3")

sig_path$Description[is.na(sig_path$Description)] <- "Other" 

sig_path$Class_1[is.na(sig_path$Class_1)] <- "Other"
table(sig_path$Class_1)

sig_path$Class_2[is.na(sig_path$Class_2)] <- "Other"
#tolower(unique(sig_path$Class_2))
library(stringr)
sig_path$Class_2 <- str_to_sentence(sig_path$Class_2)


# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_fil)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))
mixedsort(met$Visit)

library(ComplexHeatmap)

ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit


library(tidyverse)
#Make row annotation:
coeff_dt <-  result_sig[,c("Bac","Var","Value")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = Var,values_from = Value)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$Bac),]
#coeff_dt_w$VisitTC <-  NA

library(gtools)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(ps_sig)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

#col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

library(RColorBrewer)
library(circlize)
col_coef <- colorRamp2(c(-5,-2.5, 0, 2.5,5),col_grad) 

ha_left =  rowAnnotation(TwoWeeks = coeff_dt_w$Visit2.Wk,
                         OneMonth = coeff_dt_w$Visit1.Mo,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         TC =  coeff_dt_w$VisitTC,
                         col = list(TwoWeeks = col_coef,
                                    OneMonth = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef,
                                    TC =  col_coef),
                         show_legend = c(TRUE,F,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("2 Wk",  "1 Mo",  "2 Mo",  "6 Mo",  "TC"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        at =c(-5,-2.5, 0, 2.5,5),
                                                        legend_height = unit(10, "cm")))

colnames(mat) <- gsub("19300","",met$PID)

sig_path <- sig_path[match(rownames(mat),sig_path$Pathway),]
ht  =  Heatmap(log(mat+ 1),name = "Log(RPKM+1)",
               column_split = splitcols,
               row_split =  factor(sig_path$Class_2),
               #         top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T, heatmap_legend_param = list(title = "Log(RPKM + 1)", 
                                                       legend_width = unit(10, "cm"),
                                                       legend_direction = "horizontal",
                                                       title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                       labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 10)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 10)))
pdf(paste0(fig_folder,"/Fig_S2_A.pdf"),height = 8, width = 25,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()

###########################Fig S2 B ##############################################################

# Enriched :

p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]
unique(result_sig$Bac)
library(ggplot2)
ggplot(data = result_sig, aes(x = result_sig$Value)) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


result_sig <-  result_sig[result_sig$Value > 3,]
unique(result_sig$Bac)
unique(result_sig$Var)



## Annotation of metacyc pathways:
ann_path <-  read.csv("../data/humann/map_metacyc-pwy_lineage.tsv",sep = "\t",header = F)
ann_path <-  unique(ann_path)
library(stringr)
sep_dt <-  data.frame(str_split_fixed(ann_path$V2,"\\|",n = 4))
ann_path <-  data.frame(ann_path,sep_dt)

sig_path <-  data.frame(Pathway = unique(result_sig$Bac),
                        ID = gsub(":.*","",unique(result_sig$Bac)))

table(sig_path$ID  %in% unique(ann_path$V1))


sig_path$ID[!sig_path$ID  %in% unique(ann_path$V1)]


ann_path <-  ann_path[ann_path$V1 %in% sig_path$ID,]
#ann_path <-  ann_path[!grepl("Super-Pathways",ann_path$V2),]

# Also remove the duplicates by choosing size of description:
ann_path$size <- nchar(ann_path$V2)
library(dplyr)
ann_path <- ann_path %>%
  group_by(V1) %>%
  slice_max(size) %>%
  ungroup %>% 
  group_by(V1) %>%
  sample_n(1) %>% data.frame()


sig_path <-  sig_path %>%
  left_join(ann_path,by  =  c("ID"="V1"))


ann_path[grep("GLYOXYLATE",ann_path$V1)]
tmp1 <- stack(table(ann_path$V1))
tmp1 <- tmp1[tmp1$values > 1,]
tmp2 <-  ann_path[ann_path$V1 %in% as.character(tmp1$ind),]


# Clean up the names 
sig_path <-  sig_path[,c("Pathway","ID","V2","X1","X2","X3")]
names(sig_path)[3:6] <-  c("Description","Class_1","Class_2","Class_3")

sig_path$Description[is.na(sig_path$Description)] <- "Other" 

sig_path$Class_1[is.na(sig_path$Class_1)] <- "Other"
table(sig_path$Class_1)

sig_path$Class_2[is.na(sig_path$Class_2)] <- "Other"
library(stringr)
sig_path$Class_2 <- str_to_sentence(sig_path$Class_2)
tmp <- stack(table(sig_path$Class_2))


# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_fil)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))
mixedsort(met$Visit)

library(ComplexHeatmap)

ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit


library(tidyverse)
#Make row annotation:
coeff_dt <-  result_sig[,c("Bac","Var","Value")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = Var,values_from = Value)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$Bac),]
coeff_dt_w$VisitTC <-  NA

library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(ps_sig)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

#col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-5,-2.5, 0, 2.5,5),col_grad) 

ha_left =  rowAnnotation(TwoWeeks = coeff_dt_w$Visit2.Wk,
                         OneMonth = coeff_dt_w$Visit1.Mo,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         TC =  coeff_dt_w$VisitTC,
                         col = list(TwoWeeks = col_coef,
                                    OneMonth = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef,
                                    TC =  col_coef),
                         show_legend = c(TRUE,F,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("2 Wk",  "1 Mo",  "2 Mo",  "6 Mo",  "TC"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        at =c(-5,-2.5, 0, 2.5,5),
                                                        legend_height = unit(10, "cm")))

colnames(mat) <- gsub("19300","",met$PID)

sig_path <- sig_path[match(rownames(mat),sig_path$Pathway),]
ht  =  Heatmap(log(mat+ 1),name = "Log(RPKM+1)",
               column_split = splitcols,
               row_split =  factor(sig_path$Class_2),
               #         top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T, heatmap_legend_param = list(title = "Log(RPKM + 1)", 
                                                       legend_width = unit(10, "cm"),
                                                       legend_direction = "horizontal",
                                                       title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                       labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 10)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 10)))
pdf(paste0(fig_folder,"/Fig_S2_B.pdf"),height = 18, width = 25,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()



############################Fig S2 C #######################################################################

## Differential microbiome analysis:
phy_mic <- readRDS("../Processed_data/phy_mtph_hrze.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_mic))
levels(sm_dt$Visit)[c(4,5,6)] <- c("1 Mo","2 Mo","6 Mo")
sample_data(phy_mic) <-  sample_data(sm_dt)

############## Fig 2A  #########################
phy_lme <- phy_mic
phy_lme <- prune_taxa(taxa_sums(phy_lme)>0,phy_lme)

# LME analysis:

# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_lme),
               MARGIN = ifelse(taxa_are_rows(phy_lme), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_lme),
                    tax_table(phy_lme))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.05
prev_cutoff <- prev_frac*nsamples(phy_lme) # Cut-off
abun_cutoff <- 0.0
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]


phy_fil <-  prune_taxa(as.character(prevdf_fil$OTU),phy_lme)

melt_dt <-  psmelt(phy_fil)
taxa <- taxa_names(phy_fil)[1]
result_lme <- list()
for (taxa in taxa_names(phy_fil)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$PID <- factor(phy_m$PID)
  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Sex + Age + Visit  ,random = ~1 |PID, phy_m)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Bac <- taxa
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only Visit
result_lme_dt <- result_lme_dt[grep("Visit",result_lme_dt$Var),]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

# Save the LME results
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Species_HRZE.csv",sep = "/"))


# Filter results with pval
p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]
unique(result_sig$Bac)
# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_fil)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))

mixedsort(met$Visit)
ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit

library(tidyverse)
#Make row annotation:
coeff_dt <-  result_sig[,c("Bac","Var","Value")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = Var,values_from = Value)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$Bac),]


library(gtools)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(ps_sig)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_trt <- c("0 Day" = "grey50","1 Wk"= "#cb2314","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-0.5,-0.25, 0,0.25, 0.5), col_grad) 

ha_left =  rowAnnotation(OneWeek = coeff_dt_w$Visit1.Wk,
                         TwoWeeks = coeff_dt_w$Visit2.Wk,
                         OneMonth = coeff_dt_w$Visit1.Mo,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         col = list(OneWeek = col_coef,
                                    TwoWeeks = col_coef,
                                    OneMonth = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef),
                         show_legend = c(TRUE,F,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("1 Wk","2 Wk",  "1 Mo",  "2 Mo",  "6 Mo"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        at =c(-0.5,-0.25, 0,0.25, 0.5),
                                                        legend_height = unit(10, "cm")))
#lgd = Legend(col_fun = col_coef, title = "foo", direction = "horizontal")

tax_dt <- tax_table(ps_sig)@.Data %>% data.frame()
tax_dt <-  tax_dt[match(rownames(mat),rownames(tax_dt)),]
split_rows <- factor(tax_dt$Class)

colnames(mat) <- gsub("18300","",met$PID)
ht  =  Heatmap(log10(mat+ 0.00001),name = "Log10(RA)",
               column_split = splitcols,
               row_split =  split_rows,
               #         top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 14,fontface = "bold"),
               column_names_gp = gpar(fontsize = 10,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T,
               heatmap_legend_param = list(title = "RA", at = c(0,-1,-2,-3,-4,-5),
                                           legend_width = unit(10, "cm"),
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                           labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                           labels = c(expression("10"^0), expression("10"^-1), 
                                                      expression("10"^-2),expression("10"^-3),expression("10"^-4),expression("10"^-5))))

# row_names_max_width = max_text_width(rownames(mat),
#                                      gp = gpar(fontsize = 10)),
# column_names_max_height = max_text_width(colnames(mat),
#                                          gp = gpar(fontsize = 10)))
pdf(paste(fig_folder,"HRZE_LME.pdf",sep = "/"), height = 12, width = 20,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# Import pathways:
sig_species <- unique(result_sig$Bac)

search_str <-  paste0(paste0("s__",sig_species),collapse = "|")


path_dat <-  read.csv("../data/humann/humann3_pathabundance.txt",sep="\t")
names(path_dat)[1] <- "Pathway"

# Now remove the Unmapped and Unintegrated:
path_dat <- path_dat[!grepl("UNMAPPED|UNINTEGRATED",path_dat$Pathway),]
path_sp_dat <-  path_dat[grep(search_str,path_dat$Pathway),]
sel_pathways <- unique(gsub("\\|.*","",path_sp_dat$Pathway))
# Select only at the community level
path_sel_dat <- path_dat[!grepl("\\|",path_dat$Pathway),]

# Select only at the community level
path_sel_dat <- path_sel_dat[path_sel_dat$Pathway %in% sel_pathways,]
path_sel_dat$Pathway <- gsub("'","",path_sel_dat$Pathway)


# Import Pathways:
phy_path <-  readRDS("../Processed_data/phy_raw_path_hrze.rds")
#phy_path <- prune_taxa(path_sel_dat$Pathway,phy_path)
sm_dt <- data.frame(sample_data(phy_path))
# Change the level of last time point as TC
levels(sm_dt$Visit)[c(4,5,6)] <- c("1 Mo","2 Mo","6 Mo")
#sample_data(phy_tax) <-  sample_data(sm_dt)
sample_data(phy_path) <-  sample_data(sm_dt)
phy_lme <- phy_path
# LME analysis:
# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_lme),
               MARGIN = ifelse(taxa_are_rows(phy_lme), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_lme))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.1
prev_cutoff <- prev_frac*nsamples(phy_lme) # Cut-off
abun_cutoff <- 0.00001
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]


phy_fil <-  prune_taxa(as.character(prevdf_fil$OTU),phy_lme)

melt_dt <-  psmelt(phy_fil)
taxa <- taxa_names(phy_fil)[1]
result_lme <- list()
for (taxa in taxa_names(phy_fil)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$PID <- factor(phy_m$PID)
  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$t_Abundance  <- log(phy_m$Abundance + 1)
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Sex + Age + Visit  ,random = ~1 |PID, phy_m)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Bac <- taxa
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only Visit
result_lme_dt <- result_lme_dt[grep("Visit",result_lme_dt$Var),]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]


p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]
unique(result_sig$Bac)

# Save the LME results
write.csv(result_sig,paste(tab_folder,"Tab_LME_HRZE_Pathways.csv",sep = "/"))

# result_sig <-  result_sig[result_sig$Value < 0,]
# unique(result_sig$Bac)

# freq_dt <- stack(table(result_sig$Bac))
# freq_dt <-  freq_dt[freq_dt$values > 1,]
# result_sig <- result_sig[result_sig$Bac %in% freq_dt$ind,]

## Annotation of metacyc pathways:
ann_path <-  read.csv("../data/humann/map_metacyc-pwy_lineage.tsv",sep = "\t",header = F)
ann_path <-  unique(ann_path)
library(stringr)
sep_dt <-  data.frame(str_split_fixed(ann_path$V2,"\\|",n = 4))
ann_path <-  data.frame(ann_path,sep_dt)

sig_path <-  data.frame(Pathway = unique(result_sig$Bac),
                        ID = gsub(":.*","",unique(result_sig$Bac)))

table(sig_path$ID  %in% unique(ann_path$V1))


sig_path$ID[!sig_path$ID  %in% unique(ann_path$V1)]


ann_path <-  ann_path[ann_path$V1 %in% sig_path$ID,]
#ann_path <-  ann_path[!grepl("Super-Pathways",ann_path$V2),]

# Also remove the duplicates by choosing size of description:
ann_path$size <- nchar(ann_path$V2)
library(dplyr)
ann_path <- ann_path %>%
  group_by(V1) %>%
  slice_max(size) %>%
  ungroup %>% 
  group_by(V1) %>%
  sample_n(1) %>% data.frame()


sig_path <-  sig_path %>%
  left_join(ann_path,by  =  c("ID"="V1"))


# ann_path[grep("GLYOXYLATE",ann_path$V1)]
# tmp1 <- stack(table(ann_path$V1))
# tmp1 <- tmp1[tmp1$values > 1,]
# tmp2 <-  ann_path[ann_path$V1 %in% as.character(tmp1$ind),]


# Clean up the names 
sig_path <-  sig_path[,c("Pathway","ID","V2","X1","X2","X3")]
names(sig_path)[3:6] <-  c("Description","Class_1","Class_2","Class_3")

sig_path$Description[is.na(sig_path$Description)] <- "Other" 

sig_path$Class_1[is.na(sig_path$Class_1)] <- "Other"
table(sig_path$Class_1)

sig_path$Class_2[is.na(sig_path$Class_2)] <- "Other"
#tolower(unique(sig_path$Class_2))
library(stringr)
sig_path$Class_2 <- str_to_sentence(sig_path$Class_2)


tmp <- stack(table(sig_path$Class_2))


# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_fil)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))
mixedsort(met$Visit)

library(ComplexHeatmap)

ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit


library(tidyverse)
#Make row annotation:
coeff_dt <-  result_sig[,c("Bac","Var","Value")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = Var,values_from = Value)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$Bac),]
#coeff_dt_w$Visit6.Mo <-  NA

library(gtools)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(ps_sig)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_trt <- c("0 Day" = "grey50","1 Wk"= "#cb2314","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

#col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-2,-1, 0, 1,2),col_grad) 

ha_left =  rowAnnotation(OneWeek = coeff_dt_w$Visit1.Wk,
                         TwoWeeks = coeff_dt_w$Visit2.Wk,
                         OneMonth = coeff_dt_w$Visit1.Mo,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         col = list(OneWeek = col_coef,
                                    TwoWeeks = col_coef,
                                    OneMonth = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef),
                         show_legend = c(TRUE,F,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("1 Wk","2 Wk",  "1 Mo",  "2 Mo",  "6 Mo"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        at =c(-2,-1, 0, 1,2),
                                                        legend_height = unit(5, "cm")))

colnames(mat) <- gsub("18300","",met$PID)

sig_path <- sig_path[match(rownames(mat),sig_path$Pathway),]
ht  =  Heatmap(log(mat+ 1),name = "Log(RPKM+1)",
               column_split = splitcols,
               row_split =  factor(sig_path$Class_2),
               #         top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_names_gp = gpar(fontsize = 8,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T, heatmap_legend_param = list(title = "Log(RPKM + 1)", 
                                                       legend_width = unit(10, "cm"),
                                                       legend_direction = "horizontal",
                                                       title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                       labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 10)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 10)))
pdf(paste0(fig_folder,"/Fig_S2_C.pdf"),height = 5, width = 26,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


