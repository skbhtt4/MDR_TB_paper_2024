library(tidyverse)
# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_S7"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_S7"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")

# Import phyloseq for metaphlan:
library(phyloseq)
phy_tax <-  readRDS("../Processed_data/phy_mtph_fmt_rel.rds")
sm_dt <-  data.frame(sample_data(phy_tax))
tax_dt <-  tax_table(phy_tax)@.Data
# Visualize the data:

# For all
phy_bdq <-  subset_samples(phy_tax, key %in% c("BDQD1","BDQD7"))

# Filter taxa that are present in both HC and TC at BDQD1
phy_bdq_0 <- subset_samples(phy_tax,key %in% "BDQD1")
phy_bdq_0 <- prune_taxa(taxa_sums(phy_bdq)> 0,phy_bdq_0)
melt_dt_w <- psmelt(phy_bdq_0) %>%
  select(Abundance, OTU, Trxt)%>%
  group_by(OTU, Trxt)%>%
  summarise(mean = mean(Abundance))%>%
  pivot_wider(id_cols = OTU,names_from = Trxt,values_from = mean)

mat_presence <- melt_dt_w[,2:5]
mat_presence[mat_presence > 0] <- 1
mat_presence$Sum <- rowSums(mat_presence)
mat_presence$Species <-  melt_dt_w$OTU

#Select the ones that is present in all four groups
mat_sel <- mat_presence[mat_presence$Sum == 4,]
phy_bdq <- prune_taxa(mat_sel$Species,phy_bdq)

melt_dt <-  psmelt(phy_bdq)
melt_dt$MouseID <- paste0("M",melt_dt$MouseID)

grep("Flavoni",taxa_names(phy_bdq))
taxa <- taxa_names(phy_bdq)[1]

result_lme <- list()
for (taxa in taxa_names(phy_bdq)){
  print(taxa)
  phy_m <-  melt_dt[melt_dt$OTU %in% taxa,]
  phy_m$PID <- factor(phy_m$MouseID)
  phy_m$Treatment <- factor(phy_m$BDQorVehicle, levels = c("VEH","BDQ"))
  phy_m$t_Abundance  <- asin(sqrt(phy_m$Abundance ))
  #  phy_m$FMT <- factor(phy_m$FMT, levels = c("TC","HC"))
  phy_m$FMT <- factor(phy_m$FMT, levels = c("HC","TC"))
  
  phy_m$Time <-  factor(phy_m$key, levels = c("BDQD1","BDQD7"))
  levels(phy_m$Time) <- c(1,7)  
  phy_m$Time <-  as.numeric(as.character(phy_m$Time)) -1
  library("nlme")
  mod_bac <- lme(t_Abundance ~ Time*Treatment*FMT  ,random = ~1 |PID, phy_m)
  library(emmeans)
  em_mod <-  emmeans(mod_bac, pairwise ~ Time*Treatment|FMT)
  contr_dt <- contrast(em_mod, interaction = TRUE, "pairwise", adjust="mvt")
  contr_dt <-  data.frame(contr_dt)
  sum_mod_dt <- contr_dt
  sum_mod_dt$Bac <- taxa
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[8] <- "pval"

result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")

result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

# Save the LME results
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Volcano_plot_S7.csv",sep = "/"))



# First create a volcano plots:
res_dt <-  result_lme_dt
head(res_dt)

res_dt$L2FC <-  log2(exp(res_dt$estimate))

#res_dt$dir <-  ifelse(res_dt$L2FC > 0, "D6","D0")
res_dt$FMT <-  factor(res_dt$FMT,levels = c("HC","TC"))

# Input
phy_ip <-  subset_samples(phy_tax, key %in% c("Input"))
phy_ip <-  prune_taxa(taxa_sums(phy_ip)> 0,phy_ip)

res_dt$Input <-  ifelse(res_dt$Bac %in% taxa_names(phy_ip),"Input","Mouse")
#res_dt$Bac <-  res_dt$Var

q_cut_off <- 0.1  
data_pvalue_less <- res_dt[which(res_dt$p_adj< q_cut_off),]

col_trt <- c("TC" = "#ee6c4d","HC" = "#293241")

library(ggrepel)
library(tidyverse)
library(ggprism)
library(lemon)
p_vol  <-ggplot() + 
  geom_point(data=res_dt[which(res_dt$pval>= q_cut_off),],
             aes(y=-log10(p_adj),x=L2FC),
             size=3,color="black",alpha=0.1)+
  # geom_point(data= data_pvalue_less %>% filter(abs(L2FC)< 0.1),
  #            aes(y=-log10(pval),x=L2FC,fill = dir),color = "black",shape = 21,size=3,alpha=0.1)+
  # geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)< 0.1),
  #                 aes(y=-log10(pval),x=L2FC,label = Bac), alpha = 0.1)+
  geom_point(data= data_pvalue_less %>% filter(abs(L2FC)>= 0),
             aes(y=-log10(p_adj),x=L2FC,fill = FMT,shape = Input),color = "black",size=3,alpha=0.8)+
  geom_text_repel(data=data_pvalue_less %>% filter(abs(L2FC)>= 0),
                  aes(y=-log10(p_adj),x=L2FC,label = Bac),
                  size = 2,
                  #fill = "white",xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                  box.padding = .5,max.overlaps = Inf)+
  scale_fill_manual(name="FMT",values = col_trt) +
  scale_shape_manual(name = "Input",values = c("Input" = 23, "Mouse" = 21))+
  geom_hline(yintercept = -log10(q_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_prism(border = TRUE)+
  facet_rep_wrap(~FMT,repeat.tick.labels = T, nrow=2)+
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 7)))+
  theme(legend.position =  "node")+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change (D6/D0)")+
  ylab("-Log10 qval")+
  xlim(-1,1)+
  ylim(0,3)

print(p_vol)

pdf(paste0(results_folder,"/Fig_S7_A.pdf"),width = 5 ,height = 5)
print(p_vol)
dev.off()

# Now  line plots for the input species:

sel_sp <-  data_pvalue_less$Bac[data_pvalue_less$Input == "Input"]
library(tidyverse)
#phy_bdq <-  subset_samples(phy_tax, key %in% c("BDQD1","BDQD3","BDQD5","BDQD7"))
phy_bdq <-  prune_taxa(taxa_sums(phy_bdq)> 0,phy_bdq)
phy_sp_sel <-  prune_taxa(sel_sp,phy_bdq)

# Make a line plots:
melt_dt <-  psmelt(phy_sp_sel)
melt_dt$MouseID <- paste0("M",melt_dt$MouseID)
#melt_dt$Time <-  as.character(gsub("BDQD","",melt_dt$key))

melt_dt$Time <-  factor(melt_dt$key,levels = c("BDQD1","BDQD7"))
library(ggplot2)
library(lemon)

col_shades <-  c("HC_BDQ" = "#293241", "HC_VEH" = "#c0c7d1","TC_BDQ" = "#ee6c4d", "TC_VEH" = "#f9d2c9")
melt_dt$Shades <-  paste0(melt_dt$FMT,"_",melt_dt$BDQorVehicle)

levels(melt_dt$Time) <- c("0","6")
melt_dt$Time <- as.numeric(as.character(melt_dt$Time))

#melt_dt <- melt_dt %>% filter(OTU %in% c("Flavonifractor_plautii"))


line_p <- ggplot(melt_dt,aes(y = Abundance +10^-4,x = Time )) +
  geom_point(data = melt_dt, aes(y = Abundance+10^-4,x = Time,fill = Shades,group = MouseID),colour="black", alpha=1, size=3,
             shape = 21,stroke = 1) +
  geom_line(data =melt_dt, aes(y = Abundance +10^-4,x = Time,color = Shades,group = MouseID,linetype = BDQorVehicle),
            alpha = 1)+
  facet_rep_wrap(OTU~FMT,ncol = 2,repeat.tick.labels = T,scales = "free_x")+
  #scale_x_continuous(breaks = c(-15,1,7, 14, 28, 56, 168), labels = c("Screening","1","7", "14", "28", "56", "168"))+
  #facet_wrap(~name~TRT.norm,scales = "free")+
  scale_fill_manual(name = "TRT",values =col_shades)+
  scale_color_manual(name = "TRT",values =col_shades)+
  ylab("Log10(Abundance+0.0001)")+
  scale_y_log10()+
  # xlab("Log10(Cytokines + 0.001)")+
  theme_prism()+
  theme(axis.text.x=element_text(size=10,face="bold",angle = 60,hjust=1),
        #axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10,face="bold"),
        #strip.background = element_rect(colour="blue", fill="white"), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),
        axis.text.y = element_text(size=10,face="bold"))+
  ggtitle("")

pdf(paste0(results_folder,"/Fig_S7_B.pdf"),width = 8,height = 6)
print(line_p)
dev.off()

