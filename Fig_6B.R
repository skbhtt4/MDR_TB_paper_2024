# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_6B"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_6"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


# Import phyloseq for metaphlan:
library(phyloseq)
phy_tax <-  readRDS("../Processed_data/phy_mtph_fmt_rel.rds")
sm_dt <-  data.frame(sample_data(phy_tax))
tax_dt <-  tax_table(phy_tax)@.Data
# Visualize the data:

#phy_hc <-  subset_samples(phy_tax, FMT %in% c("HC"))
phy_hc <-  subset_samples(phy_tax, !key %in% c("Input","FMTGD1","PFGD1"))
phy_hc <-  prune_taxa(taxa_sums(phy_hc)> 0,phy_hc)


library(vegan)
library(ggplot2)
set.seed(1057)
#phy.log <-  transform_sample_counts(phy_mic,function(x) log(x + 1))
out.pcoa<- ordinate(phy_hc,  method = "PCoA", distance = "bray")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_hc, out.pcoa,axes = c(1,2)) 

data_pca <- p$data

col_shades <-  c("HC_BDQ" = "#293241", "HC_VEH" = "#c0c7d1","TC_BDQ" = "#ee6c4d", "TC_VEH" = "#f9d2c9")
data_pca[1,]

data_pca$key <- factor(data_pca$key,levels = c("BDQD1",  "BDQD3" , "BDQD5"  ,"BDQD7" , "PBD2"))

data_pca$Time <- data_pca$key
levels(data_pca$Time) <- c("D0","D2","D4","D6","D8")

data_pca$Shades <-  paste0(data_pca$FMT,"_",data_pca$BDQorVehicle)

data_pca$label <- paste0(data_pca$MouseID)

# Add geom line to the D0 and D6 
library(tidyverse)
line_dt <-  data_pca %>%
  filter(Time %in% c("D0","D6"))


library(ggprism)
library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,fill = Shades,shape = Time))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .9,size=3) + 
  scale_color_manual(name = "TRT",values=col_shades) +
  scale_shape_manual(name = "Time", values = c("D0"= 21,"D2" = 22  ,"D4"= 23,"D6"=24,"D8"= 25))+
#  facet_wrap(~FMT)+
  #, linetype = 2
# geom_text_repel(aes(label = label))+
  scale_fill_manual(name = "FMT",values=col_shades) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
 # facet_wrap(~key,ncol = 8)+
  guides(fill = guide_legend(override.aes = list( shape = 21),title="TRT"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y)

print(p_pca)



# library(patchwork)
# pdf(paste0(results_folder,"/PCA_BDQ_plots.pdf"),height = 8, width = 8,useDingbats = F)
# print(p_pca)
# dev.off()

# Perform permanova test:
#Not  a pairwise comparisons 
source("perm_rep_test.R")

sm.dist <- vegdist(t(data.frame(otu_table(phy_hc))), "euclidean")
class(sm.dist)
rownames(sm.dist)

# Repeated measures permanova
rnames <-  rownames(sm.dist %>% as.matrix())
sm_dt <-  data.frame(sample_data(phy_hc))


sm_dt$key <- factor(sm_dt$key,levels = c("BDQD1",  "BDQD3" , "BDQD5"  ,"BDQD7" , "PBD2"))
sm_dt$Time <- sm_dt$key
levels(sm_dt$Time) <- c("D0","D2","D4","D6","D8")

sm_dt$Group <-  paste0(sm_dt$FMT,"_",data_pca$BDQorVehicle)

sm_dt <-  sm_dt[order(sm_dt$MouseID),]
block_sub_dt <-  sm_dt %>% filter(Time == "D0") %>% select(MouseID,WeightGrams)%>% unique()
rownames(block_sub_dt) <- block_sub_dt$MouseID
block_sub_dt$MouseID <-  NULL

meta_sub_dt <-  sm_dt %>% select(MouseID, Group)
rownames(meta_sub_dt) == rnames
meta_sub_dt <-meta_sub_dt[rnames,]
rownames(meta_sub_dt) == rnames
meta_sub_dt$MouseID <- NULL


# From OMNIBUS script
# https://bitbucket.org/biobakery/hmp2_analysis/src/master/overview/src/
set.seed(2057)
permanova <- PERMANOVA_repeat_measures(
  D = sm.dist, permutations=999,
  permute_within=meta_sub_dt,
    blocks=factor(sm_dt$MouseID,levels = unique(sm_dt$MouseID) ), block_data=block_sub_dt)

test_res <-  permanova$aov.tab
test_res

# Do a pairwise comparisons
phy <-  phy_hc
meta_var <-  "Group"

# Separate timepoints into four categories
# Pre , Early treatment, 6 Months , TC
sample_data(phy)$Group <-  paste0(sample_data(phy_hc)$FMT,"_",sample_data(phy_hc)$BDQorVehicle)

# Custom fuctions to do group comparisons in a loop
pairwise_rep_permanova <- function(phy,meta_var, dist = "bray", adj = "BH", perm = 999) {
  
  require(vegan)
  sample_dt <-  data.frame(sample_data(phy))
  group_var <- sample_dt[,c(meta_var)]
  
  ## list contrasts
  group_var <- as.character(group_var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, p_value = NA
  )
  
  for (i in seq(nrow(contrasts))) {
    grp1 <-  contrasts$group1[i]
    grp2 <-  contrasts$group2[i]
    
    # Subset grp1 and grp2 only:
    
    keep_idx = as.character(get_variable(phy, meta_var)) %in% c(grp1,grp2)
    phy_sub_clr <- prune_samples(keep_idx,phy)
    
    sm.dist <- vegdist(t(data.frame(otu_table(phy_sub_clr))), "bray")
    #    class(sm.dist)
    #    rownames(sm.dist)
    
    # Repeated measures permanova
    rnames <-  rownames(sm.dist %>% as.matrix())
    sm_dt <-  data.frame(sample_data(phy_sub_clr))
    
    sm_dt$key <- factor(sm_dt$key,levels = c("BDQD1",  "BDQD3" , "BDQD5"  ,"BDQD7" , "PBD2"))
    sm_dt$Time <- sm_dt$key
    levels(sm_dt$Time) <- c("D0","D2","D4","D6","D8")
    
    
    sm_dt <-  sm_dt[order(sm_dt$MouseID),]
    block_sub_dt <-  sm_dt %>% filter(Time == "D0") %>% select(MouseID,WeightGrams)%>% unique()
    rownames(block_sub_dt) <- block_sub_dt$MouseID
    block_sub_dt$MouseID <-  NULL
    
    meta_sub_dt <-  sm_dt %>% select(MouseID, Group)
    rownames(meta_sub_dt) == rnames
    meta_sub_dt <-meta_sub_dt[rnames,]
    rownames(meta_sub_dt) == rnames
    meta_sub_dt$MouseID <- NULL
    
    
    # From OMNIBUS script
    # https://bitbucket.org/biobakery/hmp2_analysis/src/master/overview/src/
    set.seed(2057)
    permanova <- PERMANOVA_repeat_measures(
      D = sm.dist, permutations=1999,
      permute_within=meta_sub_dt,
      blocks=factor(sm_dt$MouseID,levels = unique(sm_dt$MouseID) ), block_data=block_sub_dt)
    
    aov.tab <-  as.data.frame(permanova$aov.tab)
    
    contrasts$R2[i] <- round(aov.tab$R2[1], digits = 4)
    contrasts$F_value[i] <- round(aov.tab$F.Model[1], digits = 4)
    contrasts$p_value[i] <- aov.tab$`Pr(>F)`[1]
  }
  
  ## adjust p-values for multiple comparisons
  contrasts$p_adj <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)
  
  return(list(
    contrasts = contrasts, 
    "p-value adjustment" = adj, 
    permutations = perm
  ))
}


perm_rep <- pairwise_rep_permanova(phy, "Group")
ctr_mat <-  perm_rep$contrasts
test_res

# Save the LME results
write.csv(ctr_mat,paste(tab_folder,"Tab_Permanova_Group_Comparison_6B.csv",sep = "/"))
write.csv(test_res,paste(tab_folder,"Tab_Permanova_6B.csv",sep = "/"))


# Now add the variation in PCoA using ggside:

data_pca$Time_grp <-  paste0(data_pca$Shades,"_",data_pca$Time)


library(ggside)
p_plot <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .9,size=3,aes(fill = Shades,shape = Time)) + 
  scale_color_manual(name = "TRT",values=col_shades) +
  scale_shape_manual(name = "Time", values = c("D0"= 21,"D2" = 22  ,"D4"= 23,"D6"=24,"D8"= 25))+
  #  facet_wrap(~FMT)+
  #, linetype = 2
  # geom_text_repel(aes(label = label))+
  scale_fill_manual(name = "FMT",values=col_shades) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  # facet_wrap(~key,ncol = 8)+
  guides(fill = guide_legend(override.aes = list( shape = 21),title="TRT"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y) + 
  coord_fixed(ratio=1)+
  ggside::geom_xsideboxplot(aes( y = Shades,fill = Shades), orientation = "y",outlier.colour = NA, alpha = 0.3) +
  
  ggside::geom_ysideboxplot(aes( x = Time_grp,fill = Shades), orientation = "x",outlier.colour = NA, alpha = 0.3) +
  # 
  # ggside::geom_xsideboxplot(aes( y = Shades,fill = Shades), orientation = "y",outlier.colour = NA, alpha = 0.3) +
  # ggside::geom_ysideboxplot(aes( x = Shades, fill = Shades), orientation = "x",outlier.colour = NA, alpha = 0.3) +
  theme(ggside.panel.scale.x = .1,
        ggside.panel.scale.y = .5,
        ggside.axis.text.x = element_text(size=rel(0.5) ),
        ggside.axis.text.y = element_text(size=rel(0.5))) +
  ggside::scale_xsidey_discrete() +
  ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) 


pdf(paste(results_folder,"Fig_6B.pdf",sep = "/"), height = 10,width = 10,useDingbats = F)
print(p_plot)
dev.off()


# Now explain variations in each PCs using LME:

# Draw boxplot separately:
library(ggrepel)
library(lemon)

lme_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,MouseID,Time,Shades)
lme_pca_dt$Time <- as.numeric(gsub("D","",lme_pca_dt$Time))
names(lme_pca_dt)[5] <- "Treatment"


# Model PCoA1 and PCoA2
library(nlme)
lme_pca_dt$PID <- factor(lme_pca_dt$MouseID)
names(lme_pca_dt)
sum_pc1  <-lme( Axis.1 ~  Time + Treatment,
                random = ~ 1|PID,
                data = lme_pca_dt)

sum_pc2  <-lme( Axis.2  ~  Time + Treatment,
                random = ~ 1|PID,
                data = lme_pca_dt)

# DO pairwise comparison using emmeans
library(emmeans)
emmodel_1 <- emmeans(sum_pc1, pairwise ~  Treatment)
con_dt1 <-  data.frame(emmodel_1$contrasts)
con_dt1$Axis <-  "Axis.1"
emmodel_2 <- emmeans(sum_pc2, pairwise ~ Treatment)
con_dt2 <-  data.frame(emmodel_2$contrasts)
con_dt2$Axis <-  "Axis.2"
con_dt <-  rbind(con_dt1,con_dt2)

# Contrasts across both PCs 
# Save the LME results
write.csv(con_dt,paste(tab_folder,"Tab_LME_PCs_Treatment_6B.csv",sep = "/"))


# Within treatments compare time points:
lme_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,MouseID,Time,Shades)
trt <- "HC_BDQ"

lme_list <-  list()
for(trt in unique(lme_pca_dt$Shades)){
  
  library(nlme)
  sel_lme_dt <- lme_pca_dt %>%
                filter(lme_pca_dt$Shades %in% trt)
  sel_lme_dt$PID <- factor(sel_lme_dt$MouseID)
  
  mod_pc1  <-lme( Axis.1 ~  Time,
                  random = ~ 1|PID,
                  data = sel_lme_dt)
  lme_dt_pc1 <- summary(mod_pc1)
  lme_dt_pc1 <-  data.frame(lme_dt_pc1$tTable )
  lme_dt_pc1$Treatment <- trt
  lme_dt_pc1$Variable <- rownames(lme_dt_pc1)
  lme_dt_pc1$PC <- "PC1"
  
  
  mod_pc2  <-lme( Axis.2 ~  Time,
                  random = ~ 1|PID,
                  data = sel_lme_dt)
  lme_dt_pc2 <- summary(mod_pc2)
  lme_dt_pc2 <-  data.frame(lme_dt_pc2$tTable )
  lme_dt_pc2$Treatment <- trt
  lme_dt_pc2$Variable <- rownames(lme_dt_pc2)
  lme_dt_pc2$PC <- "PC2"
  
  comb_lme_dt <-  rbind(lme_dt_pc1,lme_dt_pc2)
  lme_list[[trt]] <- comb_lme_dt
  
}

lme_final_dt <-  do.call("rbind",lme_list)
lme_final_dt <-  lme_final_dt[grepl("Time",lme_final_dt$Variable),]

write.csv(lme_final_dt,paste(tab_folder,"Tab_LME_PCs_Time_6B.csv",sep = "/"))
