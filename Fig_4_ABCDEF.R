# List files 
library(stringr)
library(tidyverse)
library(phyloseq)
library(ape)

# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_4"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")


mainDir <- "../Tables"
subDir <- "Tabs_4"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")

library(phyloseq)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(gtools)
library(ggrepel)

##########*********MDR###################### 


phy_mdr <-  readRDS("../Processed_data/phy_mtph_mdr.rds")  
sm_dt <-  data.frame(sample_data(phy_mdr))  

# Loop through species list
sp_list <- c("101338","102279")
sp <-  sp_list[2]

list_div_test <- list()
list_dist_pc <-  list()
list_permanova <-  list()

for(sp in sp_list){
  
  
  # Import consensus tree
  tree <-  read.tree(paste0("../data/Midas/MDR/consensus_",sp,".tree"))  
  
  # Change the tip label to PID and time point:
  meta_dt <- data.frame(sample_data(phy_mdr))
  meta_dt <-  meta_dt[meta_dt$Sample %in% tree$tip.label,]
  
  meta_dt <-  meta_dt[,c("Sample","PID","Visit")]
  meta_st <-  data.frame(Sample = meta_dt$Sample, ID = paste0(meta_dt$PID,"_",meta_dt$Visit))
  meta_st <- meta_st[match(tree$tip.label,meta_st$Sample),]
  tree$tip.label <-  meta_st$ID
  
  
  # Import diversity:
  div_dt <-  read.csv(paste0("../data/Midas/MDR/snp_diversity_",sp,".txt"),sep = "\t")
  div_dt <- div_dt %>%
    left_join(meta_st,by = c("sample_id" = "Sample"))
  div_dt$Time <-  factor(gsub(".*_","",div_dt$ID))
  levels(div_dt$Time)[2] <- "TC"
  div_dt$PID <-  gsub("_.*","",div_dt$ID)
  library(rstatix)
  stat.test <- div_dt %>% wilcox_test(snps_kb ~ Time)
  wil_test <-  data.frame(stat.test)
  wil_test$Species <- sp
  
  list_div_test[[sp]] <- wil_test  
  
  library(ggplot2)
  library(ggprism)
  set.seed(100)
  p_div <- ggplot(div_dt, aes(x =Time, y = snps_kb,color = Time, fill = Time)) +
    #geom_line()+
    geom_violin(alpha = 0.1)+
    geom_jitter(height = 0, width = 0.1,size = 3, alpha = 0.9)+
    # geom_point(size = 2, alpha =0.6, color = "blue")+
    scale_color_manual(values = c("0 Day" = "#bf7c2a","TC" = "#663695"))+
    scale_fill_manual(values = c("0 Day" = "#bf7c2a","TC" = "#663695"))+
    ggtitle(paste0("SNPs diversity (p_val = ",stat.test$p,")"))+
    ylab("SNPs/kb")+
    #  ggrepel::geom_text_repel( show.legend = FALSE , alpha = 0.7 )+
    theme_prism(border = TRUE)+
    theme(legend.title = element_text(size=12, face = "bold"),
          legend.text = element_text(size=10,face = "bold"))+
    theme(legend.position="none")
  # p_div
  pdf(paste0(results_folder,"/SNP_diversity_",gsub(" ","_",sp),".pdf"),width = 5,height = 5)
  print(p_div)
  dev.off()
  
  # Import SNPs frequency tab
  # Allele frequency:
  freq_dt <-  read.csv(paste0("../data/Midas/MDR/snps_freq_",sp,".txt"),sep ="\t")
  rownames(freq_dt) <- freq_dt$site_id
  freq_dt$site_id <-  NULL
  freq_dt <- freq_dt[,match(make.names(meta_st$Sample),names(freq_dt))]
  library(ComplexHeatmap)
  names(freq_dt) <- meta_st$ID
  library(circlize)
  
  # Filter freq_dt
  freq_dt <-  freq_dt[rowSums(freq_dt)>0,]
  
  prev_dt <- stack(apply(freq_dt, 1, function(c)sum(c>0)))
  head(prev_dt)
  # PCA plot
  set.seed(200)
  pca_allele <-  prcomp(t(freq_dt),scale. = F)
  #screeplot(pca_allele, type = "line", main = "Scree plot")
  # Proportion of variance
  PoV <- pca_allele$sdev^2/sum(pca_allele$sdev^2)
  PoV <-  signif(100*PoV[1:2], digits = 3)
  
  df_pcoa <-  data.frame(pca_allele$x[,c(1,2,3,4)])
  df_pcoa$sample <-  rownames(df_pcoa)
  df_pcoa$Time <- gsub(".*_","",df_pcoa$sample)
  df_pcoa$Time <-  factor(df_pcoa$Time)
  levels(df_pcoa$Time)[2] <- "TC"
  library(ggplot2)
  
  centroid <- df_pcoa %>%
    group_by(Time)%>%
    summarize(PC1_mean = mean(PC1),PC2_mean = mean(PC2))
  
  df_pcoa <-  df_pcoa %>%
    left_join(centroid)
  library(ggside)
  
  p_pca_all <- ggplot(df_pcoa, aes(x =PC1, y = PC2, color = Time, fill = Time )) +
    geom_point(size = 3,shape = 21,stroke = 1,color = "black",)+
    stat_ellipse(data=df_pcoa,aes(color=Time),level = 0.95) +
    scale_color_manual(values = c("0 Day" = "#bf7c2a","TC" = "#663695"))+
    scale_fill_manual(values = c("0 Day" = "#bf7c2a","TC" = "#663695"))+
    ggtitle("PCA allele frequency")+
    # Add the median 
    geom_point(data = centroid,aes(x =PC1_mean, y = PC2_mean,color = Time ),
               size = 3,shape = 24,alpha = 1,stroke = 1,fill = "white")+
    geom_segment(data =  df_pcoa,aes(x=PC1_mean, y=PC2_mean, xend=PC1, yend=PC2,color = Time),alpha = 0.2)+
    theme_prism(border = TRUE)+
    theme(legend.title = element_text(size=12, face = "bold"),
          legend.text = element_text(size=10,face = "bold"))+
    theme(legend.position="none")+
    coord_fixed(ratio=1.5)+
    ggside::geom_xsideboxplot(aes( y = Time), orientation = "y",outlier.colour = NA, alpha = 0.3) +
    ggside::geom_ysideboxplot(aes( x = Time), orientation = "x",outlier.colour = NA, alpha = 0.3) +
    theme(ggside.panel.scale.x = .2,ggside.panel.scale.y = .2) +
    ggside::scale_xsidey_discrete() +
    ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) +
    xlab(paste0("Component 1 ","[",PoV[1],"%]"))+
    ylab(paste0("Component 2 ","[",PoV[2],"%]"))
  print(p_pca_all)
  
  pdf(paste0(results_folder,"/PCA_",gsub(" ","_",sp),".pdf"),width = 6,height = 6)
  print(p_pca_all)
  dev.off()
  
  library(vegan)
  # Perform PERMANOVA analysis 
  sm.dist <- vegdist(t(freq_dt), "bray")
  class(sm.dist)
  
  
  meta <- meta_dt
  meta$Sample <-  paste0(meta$PID,"_",meta$Visit)
  library(vegan)
  permanova <- adonis(sm.dist ~ Visit,
                      data = meta, permutations=999, method = "bray")
  permanova$aov.tab
  
  perm_tab <- data.frame(permanova$aov.tab)
  perm_tab$Species <- sp
  
  list_permanova[[sp]] <- perm_tab
  
  
  library(tidyverse)
  library(rstatix)
  test_pca_dt <-  df_pcoa %>% select(PC1,PC2,Time)
  test_pca_dt <- test_pca_dt %>%
    pivot_longer(cols = PC1:PC2,names_to = "PC",values_to = "val") %>%
    group_by(PC) %>%
    wilcox_test(val ~ Time)%>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  test_pca_dt$Species <- sp
  
  list_dist_pc[[sp]] <- test_pca_dt
  
  
  
  # PIE Chart
  
  # Import variable importance
  imp_dt <-  read.csv(paste0("../data/Midas/MDR/VIMP_",sp,".csv"))
  
  snps_dt <-  read.csv(paste0("../data/Midas/MDR/SNPs_",sp,".csv"))
  snps_dt$X <- NULL
  snps_dt <- snps_dt %>%
    left_join(imp_dt[,c(6:ncol(imp_dt))])
  # Change the tip label to PID and time point:
  meta_dt <- data.frame(sample_data(phy_mdr))
  meta_dt <-  meta_dt[meta_dt$Sample %in% unique(snps_dt$Sample),]
  
  meta_dt <-  meta_dt[,c("Sample","PID","Visit")]
  
  snps_dt <- snps_dt %>%
    left_join(meta_dt)
  
  
  library(scatterpie)
  counts <- snps_dt
  counts$Visit <-  as.character(counts$Visit)
  counts$site_id <- factor(counts$site_id,levels = imp_dt$site_id)
  
  
  counts_sum <- counts %>%
    group_by(Visit,site_id)%>%
    summarise(A = sum(count_a), C = sum(count_c),G = sum(count_g), T = sum(count_t))%>%
    as.data.frame()
  
  counts_sum <- counts_sum %>% 
    mutate(Visit_num = as.numeric(as.factor(Visit)),
           Site_num = as.numeric(site_id))
  
  
  library(ggprism)
  pdf(paste0(results_folder,"/",sp,"_SNPs_profile.pdf"),width = 15,height = 10)
  p_profile <-    ggplot()+
    geom_scatterpie(aes(x=Visit_num, y=Site_num, r=0.4),
                    cols= c("A","C","G","T"), 
                    color= "black", data= counts_sum)+
    ggtitle(paste0(" (",sp,")"))+
    scale_x_continuous(breaks=c(1,2), labels=c("Day 0", "TC")) + 
    scale_y_continuous(breaks=1:nrow(imp_dt), labels = imp_dt$Function) + 
    labs(x="Visit", y="RF imp site") + 
    coord_fixed()+
    theme_prism()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(p_profile)     
  dev.off() 
  
}

perm_comb_dt <-  do.call("rbind",list_permanova)
write.csv(perm_comb_dt,paste(tab_folder,"/Tab_Permanova_Fig_4_B_D.csv",sep = "/"))

dist_pc_comb_dt <-  do.call("rbind",list_dist_pc)
write.csv(dist_pc_comb_dt,paste(tab_folder,"/Tab_PC_dist_test_Fig_4_B_D.csv",sep = "/"))

div_test_comb_dt <-  do.call("rbind",list_div_test)
write.csv(div_test_comb_dt,paste(tab_folder,"/Tab_Diversity_test_Fig_4_A_C.csv",sep = "/"))
