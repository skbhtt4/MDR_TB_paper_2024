# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_6CDE"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_6"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")

# Import phyloseq for metaphlan:
library(phyloseq)
library(tidyverse)
phy_tax <-  readRDS("../Processed_data/phy_mtph_fmt_rel.rds")
sm_dt <-  data.frame(sample_data(phy_tax))
tax_dt <-  tax_table(phy_tax)@.Data
# Visualize the data:


# Subset samples:
phy_bdq <-  subset_samples(phy_tax, key %in% c("PFGD1","BDQD1","BDQD3","BDQD5","BDQD7","PBD2"))
phy_bdq <-  subset_samples(phy_bdq, FMT %in% c("HC"))
phy_bdq <-  prune_taxa(taxa_sums(phy_bdq)> 0,phy_bdq)

# Now Compute distance
phy_pca <-  phy_bdq
out.pcoa<- ordinate(phy_pca,  method = "PCoA", distance = "bray")

bd <- phyloseq::distance(phy_bdq,"bray")

bd <- as.matrix(bd)
sm_dt <-  data.frame(sample_data(phy_bdq))
# Do it for every time point:

# Within group distances
df_group_sim <- list()
time <-  unique(sm_dt$key)[2]
for(time in unique(sm_dt$key)){
  
  # Within groups:
  sub_dist <- list()
  sub_sm_dt <-  sm_dt %>% filter(key %in% time)
  groups_all <- unique(sm_dt$BDQorVehicle)
  group <- "VEH"
  for (group in groups_all) { 
    sample_group <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group]
    sub_dist[[group]] <- bd[ sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
  }
  braygroups<- reshape::melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  df.bray$Time <- time
  df_group_sim[[time]] <-  df.bray
}  

df_sim_dt <-  do.call("rbind",df_group_sim)



# Across group distance
# Within group distances
df_acc_sim <- list()
time <-  unique(sm_dt$key)[1]
for(time in unique(sm_dt$key)){
  
  # Across groups:
  sub_dist <-  list()
  sub_sm_dt <-  sm_dt %>% filter(key %in% time)
  groups_all <- unique(sm_dt$BDQorVehicle)
  group_1 <- "VEH"
  group_2 <- "BDQ"
  sample_group_1 <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group_1]
  sample_group_2 <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group_2]
  sub_dist[[paste0(group_1,"-",group_2)]] <- bd[ sample_group_1, sample_group_2]
  
  braygroups<- reshape::melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  df.bray$Time <- time
  df_acc_sim[[time]] <-  df.bray
}  

df_acc_dt <-  do.call("rbind",df_acc_sim)
comb_dt_hc <-  rbind(df_sim_dt,df_acc_dt)


# TC 
# Subset samples:
phy_bdq <-  subset_samples(phy_tax, key %in% c("PFGD1","BDQD1","BDQD3","BDQD5","BDQD7","PBD2"))
phy_bdq <-  subset_samples(phy_bdq, FMT %in% c("TC"))
#phy_bdq <-  subset_samples(phy_bdq, key %in% c("BDQD3"))
phy_bdq <-  prune_taxa(taxa_sums(phy_bdq)> 0,phy_bdq)

# Now Compute distance
phy_pca <-  phy_bdq
out.pcoa<- ordinate(phy_pca,  method = "PCoA", distance = "bray")

bd <- phyloseq::distance(phy_bdq,"bray")

bd <- as.matrix(bd)
sm_dt <-  data.frame(sample_data(phy_bdq))
# Do it for every time point:

# Within group distances
df_group_sim <- list()
time <-  unique(sm_dt$key)[1]
for(time in unique(sm_dt$key)){
  
  # Within groups:
  sub_dist <- list()
  sub_sm_dt <-  sm_dt %>% filter(key %in% time)
  groups_all <- unique(sm_dt$BDQorVehicle)
  group <- "VEH"
  for (group in groups_all) { 
    sample_group <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group]
    sub_dist[[group]] <- bd[ sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
  }
  braygroups<- reshape::melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  df.bray$Time <- time
  df_group_sim[[time]] <-  df.bray
}  

df_sim_dt <-  do.call("rbind",df_group_sim)



# Across group distance

# Within group distances
df_acc_sim <- list()
time <-  unique(sm_dt$key)[1]
for(time in unique(sm_dt$key)){
  
  # Across groups:
  sub_dist <-  list()
  sub_sm_dt <-  sm_dt %>% filter(key %in% time)
  groups_all <- unique(sm_dt$BDQorVehicle)
  group_1 <- "VEH"
  group_2 <- "BDQ"
  sample_group_1 <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group_1]
  sample_group_2 <- sub_sm_dt$SeqID[sub_sm_dt$BDQorVehicle ==  group_2]
  sub_dist[[paste0(group_1,"-",group_2)]] <- bd[ sample_group_1, sample_group_2]
  
  braygroups<- reshape::melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  df.bray$Time <- time
  df_acc_sim[[time]] <-  df.bray
}  

df_acc_dt <-  do.call("rbind",df_acc_sim)
comb_dt_tc <-  rbind(df_sim_dt,df_acc_dt)

comb_dt_tc$FMT <-  "TC"
comb_dt_hc$FMT <-  "HC"

final_comb_dt <-  rbind(comb_dt_hc,comb_dt_tc)

head(final_comb_dt)

final_comb_dt$Time <-  factor(final_comb_dt$Time,levels = c("PFGD1","BDQD1","BDQD3","BDQD5","BDQD7","PBD2"))
levels(final_comb_dt$Time) <-  c("Post FMT","D1","D3","D5","D7","D9")

final_comb_dt$L1 <-  factor(final_comb_dt$L1,levels = c("VEH","BDQ","VEH-BDQ"))
#col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")
# Fit a regression line:
reg_dt  <- final_comb_dt
levels(reg_dt$Time) <- c("-10","0","2","4","6","8")
reg_dt$Time <-  as.numeric(as.character(reg_dt$Time))
col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")
head(reg_dt)

library(ggprism)
library(lemon)

model_dt <-  reg_dt %>% filter(Time != c(-10))

sum_dt <- reg_dt %>%
  dplyr::group_by(L1,FMT,Time) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  

library(ggprism)
library(lemon)
p_dist <- ggplot() +
  #  geom_boxplot(position=position_dodge(.9),alpha = 0.4) +
  annotate("rect", xmin = 1, xmax = 7, ymin = 0, ymax = Inf,
           alpha = .2)+
  geom_smooth(data = model_dt,aes(x=Time,y =value,color = FMT),method = "lm",se = F, alpha = 0.2)+
  geom_errorbar(data = sum_dt,aes(x=Time, ymin=mean-se, ymax=mean+se,color = FMT),width = .1, 
                position=position_dodge(width=0),
                alpha=0.9) +
  geom_point(data = sum_dt,aes(x=Time, y=mean,fill = FMT),position=position_dodge(width=0),
             size = 3,shape = 21,color = "black", alpha =1,stroke = 0.8)+
  geom_point(data = reg_dt,aes(x=Time,y =  value,fill = FMT),
             position=position_dodge(0),size=2, shape = 21,color = "black", alpha =0.1,stroke = 0.8) +
  
  scale_fill_manual(name = "FMT",values = col_trt)+
  scale_color_manual(name = "FMT",values = col_trt)+
  theme_prism()+
  facet_rep_wrap(~L1, scales = "free_x",repeat.tick.labels = T) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ggtitle(paste0("Distance Metric = bray"))+
  ylab("Bray-Curtis distance")+
  scale_x_continuous(breaks = c(-10,0,2,4,6,8))+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 60,hjust=1),
        #axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10,face="bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),
        axis.text.y = element_text(size=10,face="bold"))+
  ggtitle("")
print(p_dist)

# Compare the distance between each post fmt samples to
pdf(paste0(results_folder,"/Fig_6CDE.pdf"),height = 5, width = 10,useDingbats = F)
print(p_dist)
dev.off()

# Linear mixed effect model:
lme_dt <-  model_dt
lme_dt$PID <-  factor(paste0(lme_dt$X1,"_",lme_dt$X2))
grp <- "BDQ"
result_lme <-  list()
for(grp in unique(lme_dt$L1)){
  library("nlme")
  sub_lme_dt <-  lme_dt %>% filter(L1 %in% grp)
  
  mod_bac <- lme(value ~ Time*FMT  ,random = ~1 |PID, sub_lme_dt)
  sum_mod_dt <- data.frame(summary(mod_bac)$tTable)
  sum_mod_dt$L1 <- grp
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[grp]] <-  sum_mod_dt
}

result_lme_dt <-  do.call("rbind",result_lme)
# Only select interaction
result_lme_dt <-  result_lme_dt[grep(":",result_lme_dt$Var),]

# Save the LME results
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Distance_6CDE.csv",sep = "/"))

