#Load libraries
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(stringr)
library(ggplot2)
library(ggthemes)
library(msigdbr)
library(tidyr)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)
library(gtools)

# Create a directory to save figures and tables
mainDir <- "../Figures"
subDir <- "Fig_S3"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,  recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_S3"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


phy_gsva <- readRDS("../Processed_data/phy_rna_vst_mdr.rds")
#I think this should be TPM--or some other normalized count 
counts <- otu_table(phy_gsva) %>% as.data.frame %>% as.matrix() #or as.matrix() #TPM data here with gene rownames and sample colnames
#metadata from RNAseq phyloseq object
metadata <- data.frame(sample_data(phy_gsva))
#%>% column_to_rownames(var = "sample") #data frame with sample metadata with rownames = colnames of counts df


#gene lists
msigdbr_species()
m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == "H" )
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_t2g <- m_t2g[,c("gene_symbol","gs_name")]
#m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
hall_set <- unstack(m_t2g)

rownames(counts) <- gsub(".*\\^","",rownames(counts))

set.seed(1057)
gsva_mdr <- gsva(counts, hall_set,
                 min.sz=5, max.sz=500,
                 method = "gsva",
                 ssgsea.norm = F,
                 kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1
)

write.csv(gsva_mdr,paste0(fig_folder,"/GSVA_all_Hallmark.csv"))

phy_gsva <- phyloseq(otu_table(gsva_mdr,taxa_are_rows = T),sample_data(metadata))
saveRDS(phy_gsva,paste0(fig_folder,"/phy_GSVA_all_hallmark.rds"))


phy_gsva
met <- data.frame(sample_data(phy_gsva))
mat <- data.frame(otu_table(phy_gsva))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(make.names(rownames(met)),colnames(mat))]

colnames(mat) == make.names(rownames(met))

status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))

mixedsort(met$Visit)
ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit

#colfunc <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))

library(circlize)
colfunc <- colorRampPalette(c("#313695", "white", "#A50026"))
col_coef <- colorRamp2(c(-1,-0.5, 0,0.5, 1), colfunc(5)) 

min(mat)
max(mat)
rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))

library(tibble)
Hallmark_pathway_info <- read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv")
rownames(Hallmark_pathway_info) <- Hallmark_pathway_info$Hallmark.Name

setdiff(rownames(mat) ,rownames(Hallmark_pathway_info))

Hallmark_pathway_info <-  Hallmark_pathway_info[match(rownames(mat),
                                                      rownames(Hallmark_pathway_info)),]
split_rows <- Hallmark_pathway_info$Process.Category

rownames(mat) <- gsub("_"," ",rownames(mat))
rownames(mat) <- str_to_sentence(rownames(mat))


colnames(mat) <-  gsub("19300","",met$PID)

ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              col = col_coef,
              cluster_column_slices = F,
              row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 14,fontface = "bold"),
              row_title_rot = 0,
              column_title_gp = gpar(fontsize = 14,fontface = "bold"),
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = F, 
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
              show_column_names = T, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              border = T,
              heatmap_legend_param = list(title = "NES", 
                                          legend_width = unit(10, "cm"),
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                          labels_gp = gpar(fontsize = 12,fontface = "bold" )),

              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
pdf(paste0(fig_folder,"/GSVA_deseq_norm_hallmark.pdf"),width = 16,height = 8)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()

phy_fil <- phy_gsva
# Perform  LME on the gsva data to 
melt_dt <-  psmelt(phy_fil)
taxa <- taxa_names(phy_fil)[1]
result_lme <- list()
for (taxa in taxa_names(phy_fil)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  #phy_m$PID <- factor(phy_m$PID)
#  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$t_Abundance  <- phy_m$Abundance 
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
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Species_Hallmark.csv",sep = "/"))


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
mat <-  mat[,match(make.names(rownames(met)),colnames(mat))]
colnames(mat) == make.names(rownames(met))

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

library(scico)
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-0.5,-0.25, 0,0.25, 0.5), col_grad) 

ha_left =  rowAnnotation(TwoWeeks = coeff_dt_w$Visit2.Wk,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         TC =  coeff_dt_w$VisitTC,
                         col = list(TwoWeeks = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef,
                                    TC =  col_coef),
                         show_legend = c(TRUE,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("2 Wk",  "2 Mo",  "6 Mo",  "TC"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        legend_height = unit(10, "cm"),
                                                        at = c(-0.5,-0.25, 0,0.25, 0.5)))


library(circlize)
colfunc <- colorRampPalette(c("#313695", "white", "#A50026"))
col_mat <- colorRamp2(c(-.8,-0.4, 0,0.4, .8), colfunc(5)) 


rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))

library(tibble)
Hallmark_pathway_info <- read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv")
rownames(Hallmark_pathway_info) <- Hallmark_pathway_info$Hallmark.Name

setdiff(rownames(mat) ,rownames(Hallmark_pathway_info))

Hallmark_pathway_info <-  Hallmark_pathway_info[match(rownames(mat),
                                                      rownames(Hallmark_pathway_info)),]
split_rows <- Hallmark_pathway_info$Process.Category

rownames(mat) <- gsub("_"," ",rownames(mat))
rownames(mat) <- str_to_sentence(rownames(mat))


colnames(mat) <-  gsub("19300","",met$PID)

max(mat)
min(mat)

ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              left_annotation = ha_left,
              col = col_mat,
              cluster_column_slices = F,
              row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 14,fontface = "bold"),
              row_title_rot = 0,
              column_title_gp = gpar(fontsize = 14,fontface = "bold"),
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = F, 
             border = TRUE,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
              show_column_names = T, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              heatmap_legend_param = list(title = "NES", 
                                          legend_width = unit(10, "cm"),
                                          at = c(-.8,-0.4, 0,0.4, .8),
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                          labels_gp = gpar(fontsize = 12,fontface = "bold" )),
              
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
pdf(paste0(fig_folder,"/Fig_S3_A_hallmark_lme.pdf"),width = 18,height = 8)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# Do similar for TB Signature:
library(TBSignatureProfiler)
## List all signatures in the profiler
names(TBsignatures)
# Add Dupnik_20 tb sig
Dupnik_20 <-  c("ANXA3","CACNA1E","CR1","CREB5","DYSF","ENTPD1","GK","GYG1","LIMK2","MAPK14",
                "NAIP","OSM","PFKFB3","PGLYRP1","PHTF1","RAB20","SIPA1L2","SOCS3","SRPK1","WDFY3")

TBsignatures$Dupnik_20 <- Dupnik_20
## We can use all of these signatures for further analysis
siglist_mdr <- names(TBsignatures)


phy_gsva <-  readRDS("../Processed_data/phy_rna_vst_mdr.rds")

#I think this should be TPM--or some other normalized count 
counts <- otu_table(phy_gsva) %>% as.data.frame %>% as.matrix() #or as.matrix() #TPM data here with gene rownames and sample colnames

#metadata from RNAseq phyloseq object
metadata <- data.frame(sample_data(phy_gsva))
#%>% column_to_rownames(var = "sample") #data frame with sample metadata with rownames = colnames of counts df

#gene lists
rownames(counts) <- gsub(".*\\^","",rownames(counts))

tbsig_list <- TBsignatures
# gsvaRes_ssgsea <- GSVA::gsva(counts, custom.gene.sets, 
#       method = "ssgsea", parallel.sz = 2) #can alter method here:("gsva", "ssgsea", "zscore", "plage")
# gsvaRes_ssgsea
set.seed(1057)
gsva_mdr <- gsva(counts, tbsig_list,
                 min.sz=2, max.sz=500,
                 method = "gsva",
                 ssgsea.norm = F,
                 kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1
)

write.csv(gsva_mdr,paste0(fig_folder,"/GSVA_all_TBSig.csv"))
#write.csv(gsvaRes_ssgsea,paste0(fig_folder,"/GSVA_all.csv"))


phy_gsva <- phyloseq(otu_table(gsva_mdr,taxa_are_rows = T),sample_data(metadata))
saveRDS(phy_gsva,paste0(fig_folder,"/phy_GSVA_all_TBSig.rds"))


phy_gsva
met <- data.frame(sample_data(phy_gsva))
mat <- data.frame(otu_table(phy_gsva))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(make.names(rownames(met)),colnames(mat))]

colnames(mat) == make.names(rownames(met))

status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))

mixedsort(met$Visit)
ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit

brewer.pal(n = 10, name = "RdYlBu")
library(circlize)
colfunc <- colorRampPalette(c("#313695", "white", "#A50026"))
col_coef <- colorRamp2(c(-1,-0.5, 0,0.5, 1), colfunc(5)) 


library(tibble)
colnames(mat) <-  gsub("19300","",met$PID)

ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              col = col_coef,
              cluster_column_slices = F,
            #  row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 14,fontface = "bold"),
              row_title_rot = 0,
              column_title_gp = gpar(fontsize = 14,fontface = "bold"),
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = F, 
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
              show_column_names = T, show_row_names = T,
              border = TRUE,
              column_names_side = "bottom",na_col="white",
              heatmap_legend_param = list(title = "NES", 
                                          legend_width = unit(10, "cm"),
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                          labels_gp = gpar(fontsize = 12,fontface = "bold" )),
              
              show_row_dend = F)


pdf(paste0(fig_folder,"/GSVA_deseq_norm_TBSig.pdf"),width = 15,height = 12)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()
# GSVA on TB signature profiles


phy_fil <- phy_gsva
# Perform  LME on the gsva data to 
melt_dt <-  psmelt(phy_fil)
taxa <- taxa_names(phy_fil)[1]
result_lme <- list()
for (taxa in taxa_names(phy_fil)){
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  #phy_m$PID <- factor(phy_m$PID)
  #  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$t_Abundance  <- phy_m$Abundance 
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
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Species_TBSig.csv",sep = "/"))


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
mat <-  mat[,match(make.names(rownames(met)),colnames(mat))]
colnames(mat) == make.names(rownames(met))

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

library(scico)
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-1.5,-0.75, 0,0.75, 1.5), col_grad) 

ha_left =  rowAnnotation(TwoWeeks = coeff_dt_w$Visit2.Wk,
                         TwoMonths =  coeff_dt_w$Visit2.Mo,
                         SixMonths = coeff_dt_w$Visit6.Mo,
                         TC =  coeff_dt_w$VisitTC,
                         col = list(TwoWeeks = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef,
                                    TC =  col_coef),
                         show_legend = c(TRUE,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("2 Wk",  "2 Mo",  "6 Mo",  "TC"),
                         
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        legend_height = unit(10, "cm"),
                                                        at =c(-1.5,-0.75, 0,0.75, 1.5)))

min(mat)
max(mat)
library(circlize)
colfunc <- colorRampPalette(c("#313695", "white", "#A50026"))
col_mat <- colorRamp2(c(-1,-0.5, 0,0.5, 1), colfunc(5)) 


colnames(mat) <-  gsub("19300","",met$PID)

max(mat)
min(mat)

ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              left_annotation = ha_left,
              col = col_mat,
              cluster_column_slices = F,
            #  row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 14,fontface = "bold"),
              row_title_rot = 0,
              column_title_gp = gpar(fontsize = 14,fontface = "bold"),
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = F, 
              border = TRUE,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
              show_column_names = T, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              heatmap_legend_param = list(title = "NES", 
                                          legend_width = unit(10, "cm"),
                                          at = c(-1,-0.5, 0,0.5, 1),
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                          labels_gp = gpar(fontsize = 12,fontface = "bold" )),
              
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
pdf(paste0(fig_folder,"/Fig_S3_B_TBSig_lme.pdf"),width = 16,height = 12)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()




