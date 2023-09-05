# Create a directory to save figures and tables
mainDir <- "../Figures"
subDir <- "Fig_S1"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_S1"
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


########## Read Microbiome data #################################################
phy_mic_mdr <- readRDS("../Processed_data/phy_mtph_mdr.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_mic_mdr))
# Change the level of last time point as TC
levels(sm_dt$Visit)[6] <- "TC"


# Now select only modeled variables:
sm_dt <-  sm_dt[,c("Sample","PID","Visit","Sex","Age")]
sm_dt$Treatment <-  "MDR"
sample_data(phy_mic_mdr) <-  sample_data(sm_dt)

# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_mic_mdr),
               MARGIN = ifelse(taxa_are_rows(phy_mic_mdr), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_mic_mdr),
                    tax_table(phy_mic_mdr))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.05
prev_cutoff <- prev_frac*nsamples(phy_mic_mdr) # Cut-off
abun_cutoff <- 0.0
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]


phy_fil_mdr <-  prune_taxa(as.character(prevdf_fil$OTU),phy_mic_mdr)

##########*********HRZE###################### 
phy_mic_hrze <- readRDS("../Processed_data/phy_mtph_hrze.rds")
# Sample data
sm_dt <- data.frame(sample_data(phy_mic_hrze))
levels(sm_dt$Visit)[c(4,5,6)] <- c("1 Mo","2 Mo","6 Mo")

# Now select only modeled variables:
sm_dt <-  sm_dt[,c("Sample","PID","Visit","Sex","Age")]
sm_dt$Treatment <-  "HRZE"
sample_data(phy_mic_hrze) <-  sample_data(sm_dt)

# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_mic_hrze),
               MARGIN = ifelse(taxa_are_rows(phy_mic_hrze), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_mic_hrze),
                    tax_table(phy_mic_hrze))
prevdf$OTU <-  rownames(prevdf)

prevdf <- prevdf[order(prevdf$TotalAbundance,decreasing = T), ]
prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)

prev_frac<-0.05
prev_cutoff <- prev_frac*nsamples(phy_mic_hrze) # Cut-off
abun_cutoff <- 0.0
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff & prevdf$Prop >= abun_cutoff, ]
phy_fil_hrze <-  prune_taxa(as.character(prevdf_fil$OTU),phy_mic_hrze)


phy_mer <-  merge_phyloseq(phy_fil_mdr, phy_fil_hrze)
sm_dt <-  data.frame(sample_data(phy_mer))

# Remove Eukaryota 
phy_mer <- subset_taxa(phy_mer,Kingdom %in% "Bacteria")

table(sm_dt$Visit,sm_dt$Treatment)

phy_mer <-  subset_samples(phy_mer, !Visit %in% c("TC","1 Wk"))
phy_mer <- prune_taxa(taxa_sums(phy_mer)>0,phy_mer)

melt_dt <-  psmelt(phy_mer)
taxa <- taxa_names(phy_mer)[1]
result_lme <- list()
for (taxa in taxa_names(phy_mer)){
  print(taxa)
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$PID <- factor(phy_m$PID)
  phy_m$Sex <- factor(phy_m$Sex, levels = c("Male","Female"))
  phy_m$Treatment <- factor(phy_m$Treatment, levels = c("HRZE","MDR"))
  phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Sex + Age + Visit*Treatment  ,random = ~1 |PID, phy_m)
  sum_mod <-  summary(mod_bac)
  library(emmeans)
  comp_trt <-  emmeans(mod_bac, pairwise ~ Treatment | Visit)
  sum_mod_dt <- data.frame(comp_trt$contrasts)
  #sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Bac <- taxa
  #sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[7] <- "pval"
# Keep only Visit
#result_lme_dt <- result_lme_dt[grep("Visit",result_lme_dt$Var),]
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

# Save the LME results
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Species_S1.csv",sep = "/"))


# Filter results with pval
p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]
unique(result_sig$Bac)

# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_mer)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
trt_col <-  c("HRZE" = "#ee6c4d","MDR" = "#293241")
mixedsort(met$Visit)
ha_column = HeatmapAnnotation(Treatment = met$Treatment,
                              Sex =  met$Sex,
                              col=list(Treatment = trt_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit

library(tidyverse)
#Make row annotation:
coeff_dt <-  result_sig[,c("Bac","Visit","estimate")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = Visit,values_from = estimate)%>%
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
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$estimate)
min(coeff_dt$estimate)

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-0.8,-0.4, 0,0.4, 0.8), col_grad) 

ha_left =  rowAnnotation(Day_0 = coeff_dt_w$X0.Day,
                         TwoWeeks = coeff_dt_w$X2.Wk,
                         OneMonth =  coeff_dt_w$X1.Mo,
                         TwoMonths =  coeff_dt_w$X2.Mo,
                         SixMonths = coeff_dt_w$X6.Mo,
                         col = list(Day_0 = col_coef,
                                    TwoWeeks = col_coef,
                                    OneMonth = col_coef,
                                    TwoMonths = col_coef,
                                    SixMonths = col_coef),
                         show_legend = c(TRUE,F,F,F,F),
                         show_annotation_name = c(TRUE,TRUE,TRUE,TRUE,TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("Day 0",  "2 Wks","1 Mo",  "2 Mo",  "6 Mo"),
                         
                         annotation_legend_param = list(title = "Coeff",
                           title_gp = gpar(fontsize = 14,fontface = "bold" ),
                           labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                         legend_height = unit(10, "cm"),
                           at = c(-0.8,-0.4, 0,0.4, 0.8)))
#lgd = Legend(col_fun = col_coef, title = "foo", direction = "horizontal")

tax_dt <- tax_table(ps_sig)@.Data %>% data.frame()
tax_dt <-  tax_dt[match(rownames(mat),rownames(tax_dt)),]
split_rows <- factor(tax_dt$Class)

colnames(mat) <- gsub("19300|18300","",met$PID)
ht  =  Heatmap(log10(mat+ 0.00001),name = "Log10(RA)",
               column_split = splitcols,
               row_split =  split_rows,
                top_annotation = ha_column,
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
pdf(paste(fig_folder,"Fig_S1.pdf",sep = "/"), height = 18, width = 30,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()

