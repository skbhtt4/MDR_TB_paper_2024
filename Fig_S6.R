# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Figs_S6"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_S6"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


# Import phyloseq for metaphlan:
library(phyloseq)
phy_tax <-  readRDS("../Processed_data/phy_mtph_fmt_rel.rds")
sm_dt <-  data.frame(sample_data(phy_tax))
tax_dt <-  tax_table(phy_tax)@.Data
# Visualize the data:

phy_pre_fmt <-  subset_samples(phy_tax, key %in% c("FMTGD1"))
phy_pre_fmt <-  prune_taxa(taxa_sums(phy_pre_fmt)> 0,phy_pre_fmt)


library(vegan)
library(ggplot2)
set.seed(1057)
#phy.log <-  transform_sample_counts(phy_mic,function(x) log(x + 1))
out.pcoa<- ordinate(phy_pre_fmt,  method = "PCoA", distance = "bray")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_pre_fmt, out.pcoa,axes = c(1,2)) 

data_pca <- p$data

col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")
data_pca[1,]

data_pca$group <-  paste0(data_pca$FMT,"_",data_pca$BDQorVehicle)
library(ggprism)
library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,fill = FMT))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .9,size=3, shape = 21) + 
  scale_color_manual(name = "TRT",values=col_trt) +
  #facet_wrap(~phase,scales = "free")+
  #, linetype = 2
#  geom_text_repel(aes(label = BDQorVehicle))+
  scale_fill_manual(name = "FMT",values=col_trt) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  
  guides(fill = guide_legend(override.aes = list( shape = 21),title="TRT"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y)

print(p_pca)

# Perform PERMANOVA analysis 
sm.dist <- vegdist(t(data.frame(otu_table(phy_pre_fmt))), "bray")
class(sm.dist)
rownames(sm.dist)

meta <- data.frame(sample_data(phy_pre_fmt))

library(vegan)
permanova <- adonis(sm.dist ~ FMT + BDQorVehicle,
                    data = meta, permutations=999, method = "bray")
permanova$aov.tab

perm_tab <- data.frame(permanova$aov.tab)

write.csv(perm_tab,paste(tab_folder,"Tab_Permanova_S6A.csv",sep = "/"))

# Now draw PCA plot with variance explained by variables
library(ggside)
p_plot <- p_pca + 
  coord_fixed(ratio=2)+
  
  ggside::geom_xsideboxplot(aes( y = FMT), orientation = "y",outlier.colour = NA, alpha = 0.3) +
  ggside::geom_ysideboxplot(aes( x = FMT), orientation = "x",outlier.colour = NA, alpha = 0.3) +
  theme(ggside.panel.scale = .3) +
  ggside::scale_xsidey_discrete() +
  ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) 


# pdf(paste(fig_folder,"Fig_S6A.pdf",sep = "/"), height =6,width = 6,useDingbats = F)
# print(p_plot)
# dev.off()

library(tidyverse)
library(rstatix)
test_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,FMT)
test_pca_dt <- test_pca_dt %>%
               pivot_longer(cols = Axis.1:Axis.2,names_to = "PC",values_to = "val") %>%
               group_by(PC) %>%
               wilcox_test(val ~ FMT)%>%
               adjust_pvalue(method = "BH") %>%
               add_significance("p.adj")
# Save the Wilcox Test
write.csv(test_pca_dt,paste(tab_folder,"Tab_Wilcox_test_PCs_S6A.csv",sep = "/"))





plot_pca_pre <- p_plot
# Run MAaslin on the abundance data across treatment
library(Maaslin2)
otu_data = data.frame(otu_table(phy_pre_fmt))
metadata = data.frame(sample_data(phy_pre_fmt))

metadata$FMT <- factor(metadata$FMT,levels = c("HC","TC"))

fit_data = Maaslin2(
  input_data = otu_data, 
  input_metadata = metadata, 
  output = paste0(fig_folder,"/Masslin_pre_fmt_output"),
  min_prevalence = 0,
  normalization = "NONE",
  transform = "AST",
  fixed_effects = c("FMT"),
  reference = "HC",
  standardize = "FALSE")


# Volcano plots:
library(ggthemes)
res_dt <-  fit_data$results

write.csv(res_dt,paste(tab_folder,"Tab_Maaslin2_S6B.csv",sep = "/"))

res_dt$L2FC <-  log2(exp(res_dt$coef))
res_dt$dir <-  ifelse(res_dt$L2FC > 0, "TC","HC")


q_cut_off <- 0.05  
data_pvalue_less <- res_dt[which(res_dt$qval< q_cut_off),]

# p_vol  <-ggplot() + 
#   geom_point(data=res_dt[which(res_dt$qval>= q_cut_off),],
#              aes(y=-log10(qval),x=L2FC),
#              size=2,color="black",alpha=0.1)+
#   geom_text_repel(data=data_pvalue_less,
#                   aes(y=-log10(qval),x=L2FC,label = feature))+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(qval),x=L2FC),
#              color="blue", size=2,alpha=1)+
#   #scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(q_cut_off),size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change (TC/HC)")+
#   ylab("-Log10 qval")+
#   xlim(-1,1)+
#   ylim(0,8)
# 
# print(p_vol)

library(tidyverse)
p_vol  <-ggplot() + 
  geom_point(data=res_dt[which(res_dt$qval>= q_cut_off),],
             aes(y=-log10(qval),x=L2FC),
             size=3,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less %>% filter(abs(L2FC)< 0.1),
             aes(y=-log10(qval),x=L2FC,fill = dir),color = "black",shape = 21,size=3,alpha=0.1)+
  geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)< 0.1),
                  aes(y=-log10(qval),x=L2FC,label = feature), alpha = 0.1)+
  geom_point(data= data_pvalue_less %>% filter(abs(L2FC)>= 0.1),
             aes(y=-log10(qval),x=L2FC,fill = dir),color = "black",shape = 21,size=3,alpha=1)+
  geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)>= 0.1),
                  aes(y=-log10(qval),x=L2FC,label = feature))+
  scale_fill_manual(name="FMT",values = col_trt) +
  geom_hline(yintercept = -log10(q_cut_off),linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_prism(border = TRUE)+
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 7)))+
  theme(legend.position =  "node")+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change (TC/HC)")+
  ylab("-Log10 qval")+
  xlim(-1,1)+
  ylim(0,11)

print(p_vol)

p_vol_pre <-  p_vol

library(patchwork)

plot_pca_pre + p_vol_pre


# Differential analysis post fmt gavaage before BDQ
phy_post_fmt <-  subset_samples(phy_tax, key %in% c("BDQD1"))
phy_post_fmt <-  prune_taxa(taxa_sums(phy_post_fmt)> 0,phy_post_fmt)

# Looking at samples postFMT comparing HC and TC before BDQ
# PCoA plots:

library(vegan)
set.seed(1057)
#phy.log <-  transform_sample_counts(phy_mic,function(x) log(x + 1))
out.pcoa<- ordinate(phy_post_fmt,  method = "PCoA", distance = "bray")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_post_fmt, out.pcoa,axes = c(1,2)) 

data_pca <- p$data

col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")

data_pca[1,]

library(ggprism)
library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,fill = FMT))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .9,size=3, shape = 21) + 
  scale_color_manual(name = "TRT",values=col_trt) +
  #facet_wrap(~phase,scales = "free")+
  #, linetype = 2
#  geom_text_repel(aes(label = BDQorVehicle))+
  scale_fill_manual(name = "TRT",values=col_trt) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  
  guides(fill = guide_legend(override.aes = list( shape = 21),title="TRT"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y)

print(p_pca)


# Perform PERMANOVA analysis 
sm.dist <- vegdist(t(data.frame(otu_table(phy_post_fmt))), "bray")
class(sm.dist)
rownames(sm.dist)

meta <- data.frame(sample_data(phy_post_fmt))

library(vegan)
permanova <- adonis(sm.dist ~ FMT + BDQorVehicle,
                    data = meta, permutations=999, method = "bray")
permanova$aov.tab

perm_tab <- data.frame(permanova$aov.tab)

write.csv(perm_tab,paste(tab_folder,"Tab_Permanova_S6C.csv",sep = "/"))


# Now draw PCA plot with variance explained by variables
library(ggside)
p_plot <- p_pca + 
  coord_fixed(ratio=2)+
  
  ggside::geom_xsideboxplot(aes( y = FMT), orientation = "y",outlier.colour = NA, alpha = 0.3) +
  ggside::geom_ysideboxplot(aes( x = FMT), orientation = "x",outlier.colour = NA, alpha = 0.3) +
  theme(ggside.panel.scale = .3) +
  ggside::scale_xsidey_discrete() +
  ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) 

# 
# pdf(paste(fig_folder,"Fig_S6A.pdf",sep = "/"), height =6,width = 6,useDingbats = F)
# print(p_plot)
# dev.off()

library(tidyverse)
library(rstatix)
test_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,FMT)
test_pca_dt <- test_pca_dt %>%
  pivot_longer(cols = Axis.1:Axis.2,names_to = "PC",values_to = "val") %>%
  group_by(PC) %>%
  wilcox_test(val ~ FMT)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Save the Wilcox Test
write.csv(test_pca_dt,paste(tab_folder,"Tab_Wilcox_test_PCs_S6C.csv",sep = "/"))


plot_pca_post <-  p_plot
# Run MAaslin on the abundance data across treatment
library(Maaslin2)
otu_data = data.frame(otu_table(phy_post_fmt))
metadata = data.frame(sample_data(phy_post_fmt))

metadata$FMT <- factor(metadata$FMT,levels = c("HC","TC"))

fit_data = Maaslin2(
  input_data = otu_data, 
  input_metadata = metadata, 
  output = paste0(fig_folder,"/Masslin_post_fmt_output"),
  min_prevalence = 0,
  normalization = "NONE",
  transform = "AST",
  fixed_effects = c("FMT"),
  reference = "HC",
  standardize = "FALSE")


res_dt <-  fit_data$results
write.csv(res_dt,paste(tab_folder,"Tab_Maaslin2_S6D.csv",sep = "/"))

# Volcano plots:
library(ggthemes)
res_dt <-  fit_data$results

res_dt$L2FC <-  log2(exp(res_dt$coef))

res_dt$dir <-  ifelse(res_dt$L2FC > 0, "TC","HC")

q_cut_off <- 0.05  
data_pvalue_less <- res_dt[which(res_dt$qval< q_cut_off),]

col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")


library(tidyverse)
p_vol  <-ggplot() + 
  geom_point(data=res_dt[which(res_dt$qval>= q_cut_off),],
             aes(y=-log10(qval),x=L2FC),
             size=3,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less %>% filter(abs(L2FC)< 0.1),
             aes(y=-log10(qval),x=L2FC,fill = dir),color = "black",shape = 21,size=3,alpha=0.1)+
  geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)< 0.1),
                  aes(y=-log10(qval),x=L2FC,label = feature), alpha = 0.1)+
  geom_point(data= data_pvalue_less %>% filter(abs(L2FC)>= 0.1),
             aes(y=-log10(qval),x=L2FC,fill = dir),color = "black",shape = 21,size=3,alpha=1)+
  geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)>= 0.1),
                  aes(y=-log10(qval),x=L2FC,label = feature))+
  scale_fill_manual(name="FMT",values = col_trt) +
  geom_hline(yintercept = -log10(q_cut_off),linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_prism(border = TRUE)+
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 7)))+
  theme(legend.position =  "node")+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change (TC/HC)")+
  ylab("-Log10 qval")+
  xlim(-1,1)+
  ylim(0,11)

print(p_vol)

plot_vol_post <-  p_vol


library(patchwork)
pdf(paste0(fig_folder,"/Fig_S6.pdf"),height = 15, width = 15,useDingbats = F)
(plot_pca_pre + p_vol_pre)/(plot_pca_post + plot_vol_post)
dev.off()

