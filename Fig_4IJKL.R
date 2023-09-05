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

##########*********HRZE###################### 
phy_mic_hrze <- readRDS("../Processed_data/phy_mtph_hrze.rds")
# Sample data
sm_dt <- data.frame(sample_data(phy_mic_hrze))
levels(sm_dt$Visit)[c(4,5,6)] <- c("1 Mo","2 Mo","6 Mo")

# Now select only modeled variables:
sm_dt <-  sm_dt[,c("Sample","PID","Visit","Sex","Age")]
sm_dt$Treatment <-  "HRZE"
sample_data(phy_mic_hrze) <-  sample_data(sm_dt)


# Select 
# Gemmiger formicilis and Phascolarcto

grep("Gemmiger_formicilis",taxa_names(phy_mic_hrze),value = T)
grep("Phascolarctobacterium_succinatutens",taxa_names(phy_mic_hrze),value = T)

sel_taxa <-  c("Gemmiger_formicilis","Phascolarctobacterium_succinatutens")

phy_rif <-  prune_taxa(sel_taxa,phy_mic_hrze)

melt_dt <-  psmelt(phy_rif)

# Now lets get the samples for Gemmiger
sp <- "100084"
freq_dt <-  read.csv(paste0("../data/Midas/HRZE/snps_freq_Gemmiger_formicilis.txt"),sep ="\t")
rownames(freq_dt) <- freq_dt$site_id
freq_dt$site_id <-  NULL

# Subset only samples 
sm_sel_dt <-sm_dt[make.names(sm_dt$Sample) %in% names(freq_dt), ]
sm_sel_dt <-sm_sel_dt[sm_sel_dt$Visit %in% c("6 Mo"), ]

phy_gem <- melt_dt[melt_dt$PID %in% unique(sm_sel_dt$PID),]
phy_gem <- phy_gem[phy_gem$OTU %in% c("Gemmiger_formicilis"),]

phy_gem$PID <-  gsub("18300","",phy_gem$PID)

lme_dt <- phy_gem
# Perform LME compared to TC
lme_dt$PID <- factor(lme_dt$PID)
lme_dt$t_Abundance  <- asin(sqrt(lme_dt$Abundance ))
lme_dt$Visit <- relevel(lme_dt$Visit,ref = "6 Mo")
library("nlme")
mod_bac <- lme(t_Abundance ~ Visit  ,random = ~1 |PID, lme_dt)
sum_mod <-  summary(mod_bac)
sum_mod_dt <- data.frame(sum_mod$tTable)
sum_mod_dt$Var <-  rownames(sum_mod_dt)
sum_mod_dt$Species <- "Gemmiger_Formicilis"

lme_gemmiger <- sum_mod_dt

write.csv(lme_gemmiger,paste0(tab_folder,"/Fig_4J_LME_abundance_Gemmiger.csv"))


# Including all the samples
library(tidyverse)
sum_dt_gem <-  phy_gem %>%
  group_by(Visit,OTU)%>%
  summarise(mean_det = mean(Abundance),
            se = sd(Abundance)/sqrt(n()),
            Total_samples = n()) 


library(ggprism)
p_gem <-  ggplot() +
  geom_pointrange(data = sum_dt_gem, aes(x=Visit, y=mean_det, ymin=mean_det-se, ymax=mean_det+se),
                  colour="black", alpha=1,fill = "white",shape = 21,stroke = 1 ,size=1)+
  facet_wrap(~OTU,scales = "free")+
  ylab("Relative Abundance")+
  theme_prism()+
  scale_y_continuous(guide = "prism_offset_minor")


pdf(paste0(results_folder,"/Fig_4_J_Gemmiger_lineplot.pdf"),width = 5,height = 5)
print(p_gem)
dev.off()


# Print variable importance
imp_dt <- read.csv(paste0("../data/Midas/HRZE","/VIMP_",sp,".csv"))


# Plot importance:s
imp_dt$Function <- factor(imp_dt$Function,levels = imp_dt$Function)
p_imp <- ggplot(imp_dt, aes(x= Function, y=MeanDecreaseAccuracy))+
  geom_bar(aes(y = MeanDecreaseAccuracy, x = Function, alpha = Freq ),
           stat="identity",fill = "darkblue",color = "black",
           data = imp_dt)+
  scale_alpha_continuous(limits = c(0,1),range = c(0,1))+
  #geom_point(size = 3, shape = 21, stroke = 1,aes( alpha = Freq), fill =  "darkblue")+
  coord_flip()+
  ylab("Variable Importance")+
  xlab("")+theme_bw()+
  ggtitle("Permutation variable importance")+
  theme_prism()

pdf(paste0(results_folder,"/Fig_4_L_Gemmiger_Importance.pdf"),width = 10,height = 4)
print(p_imp)
dev.off()




# Now lets get the samples for Phasco
sp <- "101365"
freq_dt <-  read.csv(paste0("../data/Midas/HRZE/snps_freq_Phascolarcto_succinatutens.txt"),sep ="\t")
rownames(freq_dt) <- freq_dt$site_id
freq_dt$site_id <-  NULL

# Subset only samples 
sm_sel_dt <-sm_dt[make.names(sm_dt$Sample) %in% names(freq_dt), ]
sm_sel_dt <-sm_sel_dt[sm_sel_dt$Visit %in% c("6 Mo"), ]

phy_gem <- melt_dt[melt_dt$PID %in% unique(sm_sel_dt$PID),]
phy_gem <- phy_gem[phy_gem$OTU %in% c("Phascolarctobacterium_succinatutens"),]

phy_gem$PID <-  gsub("18300","",phy_gem$PID)

lme_dt <- phy_gem
# Perform LME compared to TC
lme_dt$PID <- factor(lme_dt$PID)
lme_dt$t_Abundance  <- asin(sqrt(lme_dt$Abundance ))
lme_dt$Visit <- relevel(lme_dt$Visit,ref = "6 Mo")
library("nlme")
mod_bac <- lme(t_Abundance ~ Visit  ,random = ~1 |PID, lme_dt)
sum_mod <-  summary(mod_bac)
sum_mod_dt <- data.frame(sum_mod$tTable)
sum_mod_dt$Var <-  rownames(sum_mod_dt)
sum_mod_dt$Species <- "Phascolarctobacterium_succinatutens"

lme_phasco <- sum_mod_dt
write.csv(lme_phasco,paste0(tab_folder,"/Fig_4_I_LME_abundance_Phasco.csv"))
# Including all the samples
library(tidyverse)
sum_dt_phasco <-  phy_gem %>%
  group_by(Visit,OTU)%>%
  summarise(mean_det = mean(Abundance),
            se = sd(Abundance)/sqrt(n()),
            Total_samples = n()) 

p_phasco <-  ggplot() +
  geom_pointrange(data = sum_dt_phasco, aes(x=Visit, y=mean_det, ymin=mean_det-se, ymax=mean_det+se),
                  colour="black", alpha=1,fill = "white",shape = 21,stroke = 1 ,size=1)+
  facet_wrap(~OTU,scales = "free")+
  ylab("Relative Abundance")+
  theme_prism()+
  scale_y_continuous(guide = "prism_offset_minor")


pdf(paste0(results_folder,"/Fig_4_I_Phascolarcto_lineplot.pdf"),width = 5,height = 5)
print(p_phasco)
dev.off()



# Print variable importance
imp_dt <- read.csv(paste0("../data/Midas/HRZE","/VIMP_",sp,".csv"))


# Plot importance:s
imp_dt$Function <- factor(imp_dt$Function,levels = imp_dt$Function)
p_imp <- ggplot(imp_dt, aes(x= Function, y=MeanDecreaseAccuracy))+
  geom_bar(aes(y = MeanDecreaseAccuracy, x = Function, alpha = Freq ),
           stat="identity",fill = "darkblue",color = "black",
           data = imp_dt)+
  scale_alpha_continuous(limits = c(0,1),range = c(0,1))+
  #geom_point(size = 3, shape = 21, stroke = 1,aes( alpha = Freq), fill =  "darkblue")+
  coord_flip()+
  ylab("Variable Importance")+
  xlab("")+theme_bw()+
  ggtitle("Permutation variable importance")+
  theme_prism()

pdf(paste0(results_folder,"/Fig_4_K_Phasco_Importance.pdf"),width = 10,height = 4)
print(p_imp)
dev.off()
