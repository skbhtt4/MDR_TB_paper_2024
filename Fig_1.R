# Create a directory to save figures and tables
mainDir <- "../Figures"
subDir <- "Figs_1"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_1"
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
phy_mic <- readRDS("../Processed_data/phy_mtph_counts.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_mic))

############## Fig 1 B  #########################

ttp_dt <-  sm_dt %>%
  select(PID,Visit,Sex, Age,BMI, Av_TTP,ind,Time)
# Remove data rows where TTP is NA
ttp_dt <-  ttp_dt[!is.na(ttp_dt$Av_TTP),]

# Change 24 months to Treatment Completion
levels(ttp_dt$Visit)[6] <- "TC"
# Remove missing TTP values
#ttp_dt <-  ttp_dt[ttp_dt$ind == 0,]

# Model TTP 
library(nlme)
lme_dt <-  ttp_dt
lme_dt$Sex <- factor(lme_dt$Sex)
lme_dt$PID <- factor(lme_dt$PID)
names(lme_dt)
fit_ttp <-lme( log(Av_TTP) ~ Sex + Age + Visit ,
               random = ~ 1|PID,
               data = lme_dt)

sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)

# Save the LME results
write.csv(dt_sum,paste(tab_folder,"Tab_LME_TTP_1B.csv",sep = "/"))

# Plot TTP
library(ggprism)
library(ggplot2)
p <- ggplot()+
  geom_point(data=ttp_dt, aes(x = Visit, y= Av_TTP,
                              group = PID), 
             # position=position_dodge(.5),
             size=3, shape = 21, alpha =1,fill = "#ff2f17",stroke = 1)+
  geom_line(data=ttp_dt, aes(x = Visit, y=Av_TTP, group = PID), 
            #position=position_dodge(.5), 
            alpha = 0.4)+
  theme_prism()+
  xlab("Day")+
  ylab("Time to positivity (TTP)")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=18,face = "bold"),
        axis.title.y = element_text(size=18,face = "bold"),
        axis.text.x = element_text(size=15,angle=45, hjust=1,vjust = 1,face = "bold"),
        axis.text.y = element_text(size=15,face = "bold"))

pdf(paste(fig_folder,"Fig_1B.pdf",sep = "/"), height = 5,width = 5,useDingbats = F)
print(p)
dev.off()


############## Fig 1 C  #########################
# PCoA plots:

library(microViz)
clr_trans <-  phy_mic %>%
  tax_transform(trans = "clr",rank = "Species")
phy_clr <- clr_trans

#phy_clr <- transform_sample_counts(phy_mic,function(x) log(x + 1))
library(vegan)
set.seed(1057)
#phy.log <-  transform_sample_counts(phy_mic,function(x) log(x + 1))
out.pcoa<- ordinate(phy_clr,  method = "PCoA", distance = "euclidean")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_clr, out.pcoa,axes = c(1,2)) 

data_pca <- p$data

#col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo" = "#51BBFE","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "lightgreen")

col_trt <-  c("0 Day" = "#8783c2","2 Wk"= "#a7aed8","1 Mo" = "#96e0e1","2 Mo" = "#97dab8","6 Mo" = "#bcdf92","TC" = "#ece97c")

data_pca[1,]

levels(data_pca$Visit)[6] <- "TC"

library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,fill = Visit))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .7,size=3, shape = 21) + 
  #stat_ellipse(data=data_pca,aes(color=treatment),  level = 0.8) +
  scale_color_manual(name = "Visit",values=col_trt) +
  #facet_wrap(~phase,scales = "free")+
  #, linetype = 2
  # geom_text_repel(aes(label = mouse))+
  scale_fill_manual(name = "Visit",values=col_trt) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  
  guides(fill = guide_legend(override.aes = list( shape = 21),title="Visit"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y)
# print(p_pca)
# 
# pdf(paste(fig_folder,"Fig_1C.pdf",sep = "/"), height = 4,width = 6,useDingbats = F)
# print(p_pca)
# dev.off()


# Draw boxplot separately:
library(ggrepel)
library(lemon)

lme_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,PID,Visit)
# Model PCoA1 and PCoA2
library(nlme)
lme_pca_dt$PID <- factor(lme_pca_dt$PID)
names(lme_pca_dt)
sum_pc1  <-lme( Axis.1 ~ Visit ,
               random = ~ 1|PID,
               data = lme_pca_dt)
sum_pc2  <-lme( Axis.2 ~ Visit ,
                random = ~ 1|PID,
                data = lme_pca_dt)

# DO pairwise comparison using emmeans
library(emmeans)
emmodel_1 <- emmeans(sum_pc1, pairwise ~ Visit)
con_dt1 <-  data.frame(emmodel_1$contrasts)
con_dt1$Axis <-  "Axis.1"
emmodel_2 <- emmeans(sum_pc2, pairwise ~ Visit)
con_dt2 <-  data.frame(emmodel_2$contrasts)
con_dt2$Axis <-  "Axis.2"
con_dt <-  rbind(con_dt1,con_dt2)

# Contrasts across both PCs 
# Save the LME results
write.csv(con_dt,paste(tab_folder,"Tab_LME_PCs_1C.csv",sep = "/"))


library(ggside)
p_plot <- p_pca + 
  coord_fixed(ratio=1)+

ggside::geom_xsideboxplot(aes( y = Visit), orientation = "y",outlier.colour = NA, alpha = 0.3) +
ggside::geom_ysideboxplot(aes( x = Visit), orientation = "x",outlier.colour = NA, alpha = 0.3) +
theme(ggside.panel.scale = .4) +
ggside::scale_xsidey_discrete() +
ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) 
  

pdf(paste(fig_folder,"Fig_1C.pdf",sep = "/"), height =9,width = 9,useDingbats = F)
print(p_plot)
dev.off()

# Now do the repeated measure permanova
# Permanova
#sm.dist<-vegdist(, method='euclidean')

#Not  a pairwise comparisons 
source("perm_rep_test.R")

sm.dist <- vegdist(t(data.frame(otu_table(phy_clr))), "euclidean")
class(sm.dist)
rownames(sm.dist)

# Repeated measures permanova
rnames <-  rownames(sm.dist %>% as.matrix())
sm_dt <-  data.frame(sample_data(phy_clr))
levels(sm_dt$Visit)[6] <- "TC"

sm_dt <-  sm_dt[order(sm_dt$PID),]
block_sub_dt <-  sm_dt %>% select(PID,Sex, Age,Height) %>% unique()
# Replace height with mean height:
block_sub_dt$Height[is.na(block_sub_dt$Height)] <- mean(block_sub_dt$Height,na.rm = T)
rownames(block_sub_dt) <- block_sub_dt$PID
block_sub_dt$PID <-  NULL

meta_sub_dt <-  sm_dt %>% select(PID, Visit)
rownames(meta_sub_dt) == rnames
meta_sub_dt <-meta_sub_dt[rnames,]
rownames(meta_sub_dt) == rnames
meta_sub_dt$PID <- NULL


# From OMNIBUS script
# https://bitbucket.org/biobakery/hmp2_analysis/src/master/overview/src/
set.seed(2057)
permanova <- PERMANOVA_repeat_measures(
  D = sm.dist, permutations=999,
  permute_within=meta_sub_dt,
  blocks=factor(sm_dt$PID,levels = unique(sm_dt$PID) ), block_data=block_sub_dt)

test_res <-  permanova$aov.tab
test_res

# Do a pairwise comparisons
phy <-  phy_clr
meta_var <-  "TimeID"

# Separate timepoints into four categories
# Pre , Early treatment, 6 Months , TC
sample_data(phy)$TimeID <- ifelse(sample_data(phy)$Visit == "0 Day","Pre",
                                  ifelse(sample_data(phy)$Visit == "6 Mo","6 Mo",
                                         ifelse(sample_data(phy)$Visit == "24 Mo","TC","Early_treat")))

# Custom fuctions to do group comparisons in a loop
pairwise_rep_permanova <- function(phy,meta_var, dist = "euclidean", adj = "BH", perm = 999) {
  
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
    
    sm.dist <- vegdist(t(data.frame(otu_table(phy_sub_clr))), "euclidean")
#    class(sm.dist)
#    rownames(sm.dist)
    
    # Repeated measures permanova
    rnames <-  rownames(sm.dist %>% as.matrix())
    sm_dt <-  data.frame(sample_data(phy_sub_clr))
    
    sm_dt <-  sm_dt[order(sm_dt$PID),]
    block_sub_dt <-  sm_dt %>% select(PID,Sex, Age,Height) %>% unique()
    # Replace height with mean height:
    block_sub_dt$Height[is.na(block_sub_dt$Height)] <- mean(block_sub_dt$Height,na.rm = T)
    rownames(block_sub_dt) <- block_sub_dt$PID
    block_sub_dt$PID <-  NULL
    
    meta_sub_dt <-  sm_dt %>% select(PID, {{meta_var}})
    rownames(meta_sub_dt) == rnames
    meta_sub_dt <-meta_sub_dt[rnames,]
    rownames(meta_sub_dt) == rnames
    meta_sub_dt$PID <- NULL
    
    
    # From OMNIBUS script
    # https://bitbucket.org/biobakery/hmp2_analysis/src/master/overview/src/
    set.seed(2057)
    permanova <- PERMANOVA_repeat_measures(
      D = sm.dist, permutations=perm,
      permute_within=meta_sub_dt,
      blocks=factor(sm_dt$PID,levels = unique(sm_dt$PID) ), block_data=block_sub_dt)
    
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


perm_rep <- pairwise_rep_permanova(phy, "TimeID")
ctr_mat <-  perm_rep$contrasts
test_res

# Save the LME results
write.csv(ctr_mat,paste(tab_folder,"Tab_Permanova_Group_Comparison_1C.csv",sep = "/"))
write.csv(test_res,paste(tab_folder,"Tab_Permanova_1C.csv",sep = "/"))


############## Fig 1 D  #########################

# PCA plots for metabolic pathways:
phy_path <- readRDS("../Processed_data/phy_raw_path_mdr.rds")
tax_dt <- data.frame(Path1 = taxa_names(phy_path), Path2 = taxa_names(phy_path) )
rownames(tax_dt) <- taxa_names(phy_path)
tax_table(phy_path) <-  tax_table(as.matrix(tax_dt))

phy_path <- transform_sample_counts(phy_path,function(x) ceiling(x))
#phy_path <- transform_sample_counts(phy_path,function(x) log(x + 1))

clr_trans <-  phy_path %>%
  tax_transform(trans = "clr",rank = NA)
phy_clr <- clr_trans

library(vegan)
set.seed(1057)
out.pcoa<- ordinate(phy_clr,  method = "PCoA", distance = "euclidean")
evals <- out.pcoa$values[,1]
phyloseq::plot_scree(out.pcoa) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

p <- plot_ordination(phy_clr, out.pcoa,axes = c(1,2)) 

data_pca <- p$data

#col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo" = "#51BBFE","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "lightgreen")
col_trt <-  c("0 Day" = "#8783c2","2 Wk"= "#a7aed8","1 Mo" = "#96e0e1","2 Mo" = "#97dab8","6 Mo" = "#bcdf92","TC" = "#ece97c")

data_pca[1,]

levels(data_pca$Visit)[6] <- "TC"

library(ggrepel)
library(lemon)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,fill = Visit))+ 
  geom_point(color = "black",stroke = 1,
             alpha = .7,size=3, shape = 21) + 
  #stat_ellipse(data=data_pca,aes(color=treatment),  level = 0.8) +
  scale_color_manual(name = "Visit",values=col_trt) +
  #facet_wrap(~phase,scales = "free")+
  #, linetype = 2
  # geom_text_repel(aes(label = mouse))+
  scale_fill_manual(name = "Visit",values=col_trt) +
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  
  guides(fill = guide_legend(override.aes = list( shape = 21),title="Visit"),
         shape = guide_legend(override.aes = list(fill = "black") ))+
  xlab(p$labels$x)+
  ylab (p$labels$y)

# print(p_pca)

# pdf(paste(fig_folder,"Fig_1D.pdf",sep = "/"), height = 9,width = 9,useDingbats = F)
# print(p_pca)
# dev.off()

library(ggrepel)
library(lemon)

lme_pca_dt <-  data_pca %>% select(Axis.1,Axis.2,PID,Visit)
# Model PCoA1 and PCoA2
library(nlme)
lme_pca_dt$PID <- factor(lme_pca_dt$PID)
names(lme_pca_dt)
sum_pc1  <-lme( Axis.1 ~ Visit ,
                random = ~ 1|PID,
                data = lme_pca_dt)
sum_pc2  <-lme( Axis.2 ~ Visit ,
                random = ~ 1|PID,
                data = lme_pca_dt)

# DO pairwise comparison using emmeans
library(emmeans)
emmodel_1 <- emmeans(sum_pc1, pairwise ~ Visit)
con_dt1 <-  data.frame(emmodel_1$contrasts)
con_dt1$Axis <-  "Axis.1"
emmodel_2 <- emmeans(sum_pc2, pairwise ~ Visit)
con_dt2 <-  data.frame(emmodel_2$contrasts)
con_dt2$Axis <-  "Axis.2"
con_dt <-  rbind(con_dt1,con_dt2)

# Contrasts across both PCs 
# Save the LME results
write.csv(con_dt,paste(tab_folder,"Tab_LME_PCs_1D.csv",sep = "/"))


library(ggside)
p_plot <- p_pca + 
  coord_fixed(ratio=1.5)+
  
  ggside::geom_xsideboxplot(aes( y = Visit), orientation = "y",outlier.colour = NA, alpha = 0.3) +
  ggside::geom_ysideboxplot(aes( x = Visit), orientation = "x",outlier.colour = NA, alpha = 0.3) +
  theme(ggside.panel.scale = .3) +
  ggside::scale_xsidey_discrete() +
  ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) 


pdf(paste(fig_folder,"Fig_1D.pdf",sep = "/"), height =9,width = 9,useDingbats = F)
print(p_plot)
dev.off()



# Now do the repeated measure permanova
# Permanova
#sm.dist<-vegdist(, method='euclidean')

#Not  a pairwise comparisons 
source("perm_rep_test.R")

sm.dist <- vegdist(t(data.frame(otu_table(phy_clr))), "euclidean")
class(sm.dist)
rownames(sm.dist)

# Repeated measures permanova
rnames <-  rownames(sm.dist %>% as.matrix())
sm_dt <-  data.frame(sample_data(phy_clr))
levels(sm_dt$Visit)[6] <- "TC"

sm_dt <-  sm_dt[order(sm_dt$PID),]
block_sub_dt <-  sm_dt %>% select(PID,Sex, Age,Height) %>% unique()
# Replace height with mean height:
block_sub_dt$Height[is.na(block_sub_dt$Height)] <- mean(block_sub_dt$Height,na.rm = T)
rownames(block_sub_dt) <- block_sub_dt$PID
block_sub_dt$PID <-  NULL

meta_sub_dt <-  sm_dt %>% select(PID, Visit)
rownames(meta_sub_dt) == rnames
meta_sub_dt <-meta_sub_dt[rnames,]
rownames(meta_sub_dt) == rnames
meta_sub_dt$PID <- NULL


# From OMNIBUS script
# https://bitbucket.org/biobakery/hmp2_analysis/src/master/overview/src/
set.seed(2057)
permanova <- PERMANOVA_repeat_measures(
  D = sm.dist, permutations=999,
  permute_within=meta_sub_dt,
  blocks=factor(sm_dt$PID,levels = unique(sm_dt$PID) ), block_data=block_sub_dt)

test_res <-  permanova$aov.tab
test_res

# Do a pairwise comparisons
phy <-  phy_clr
meta_var <-  "TimeID"

# Separate timepoints into four categories
# Pre , Early treatment, 6 Months , TC
sample_data(phy)$TimeID <- ifelse(sample_data(phy)$Visit == "0 Day","Pre",
                                  ifelse(sample_data(phy)$Visit == "6 Mo","6 Mo",
                                         ifelse(sample_data(phy)$Visit == "24 Mo","TC","Early_treat")))

# Custom fuctions to do group comparisons in a loop
perm_rep <- pairwise_rep_permanova(phy, "TimeID")
ctr_mat <-  perm_rep$contrasts
test_res

# Save the LME results
write.csv(ctr_mat,paste(tab_folder,"Tab_Permanova_Group_Comparison_1D.csv",sep = "/"))
write.csv(test_res,paste(tab_folder,"Tab_Permanova_1D.csv",sep = "/"))



# Diversity plots:
library(yingtools2)
library(purrr)
#alpha <- estimate_richness(phy_mic) #%>% mutate(sample=row.names(.))
div_dt <- get.samp(phy_mic,stats = T) %>% as.data.frame()

lme_div_dt <-  div_dt %>%
  select(PID,Visit,Sex, Age,BMI, Av_TTP,ind,Time,InvSimpson)
# Remove data rows where TTP is NA
lme_div_dt <-  lme_div_dt[!is.na(lme_div_dt$InvSimpson),]

# Change 24 months to Treatment Completion
levels(lme_div_dt$Visit)[6] <- "TC"

lme_div_dt$Sex <-  factor(lme_div_dt$Sex)
lme_div_dt$PID <- factor(lme_div_dt$PID)

# Inv Simpson 
fit_div <-lme( InvSimpson ~ Sex + Age + Visit ,
               random = ~ 1| PID,
               data = lme_div_dt)
sum_fit <- summary(fit_div)
dt_sum <-  data.frame(sum_fit$tTable)
dt_sum
# Save the lme results
write.csv(dt_sum,paste(tab_folder,"LME_InvSimpson_div_1_E.csv",sep = "/"))

# Do pairwise comparisons too
emmodel_div <- emmeans(fit_div, pairwise ~ Visit)
con_div_dt <-  data.frame(emmodel_div$contrasts)
# Save the lme results
write.csv(con_div_dt,paste(tab_folder,"LME_InvSimpson_div_pairwise_comp_1_E.csv",sep = "/"))




library(ggprism)
p <- ggplot()+
  geom_point(data=lme_div_dt, aes(x = Visit, y= InvSimpson,
                                  group = PID), 
             # position=position_dodge(.5),
             size=2, shape = 21, alpha =1,fill = "#ff2f17",stroke = 1)+
  geom_line(data=lme_div_dt, aes(x = Visit, y=InvSimpson, group = PID), 
            #position=position_dodge(.5), 
            alpha = 0.4)+
  stat_summary(data=lme_div_dt, aes(x = Visit, y=InvSimpson, group = Visit), fun = median,
               geom = "crossbar", width = 0.2)+
  theme_prism()+
  xlab("Day")+
  ylab("Inv Simpson")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=18,face = "bold"),
        axis.title.y = element_text(size=18,face = "bold"),
        axis.text.x = element_text(size=15,angle=45, hjust=1,vjust = 1,face = "bold"),
        axis.text.y = element_text(size=15,face = "bold"))

pdf(paste(fig_folder,"Fig_1E.pdf",sep = "/"), height = 5,width = 5,useDingbats = F)
print(p)
dev.off()



############ Fig 1F ###################################

# Read Microbiome data
phy_mic <- readRDS("../Processed_data/phy_mtph_mdr.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_mic))
# Change the level of last time point as TC
levels(sm_dt$Visit)[6] <- "TC"
sample_data(phy_mic) <-  sample_data(sm_dt)

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
  mod_bac <- lme(t_Abundance ~ Sex  + Age +  Visit  ,random = ~1 |PID, data = phy_m)
  # Optimize if it doesnt work
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
write.csv(result_lme_dt,paste(tab_folder,"Tab_LME_Species_1F.csv",sep = "/"))


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
                                                        legend_height = unit(10, "cm"),
                                                        at = c(-0.5,-0.25, 0,0.25, 0.5)))
#lgd = Legend(col_fun = col_coef, title = "foo", direction = "horizontal")

tax_dt <- tax_table(ps_sig)@.Data %>% data.frame()
tax_dt <-  tax_dt[match(rownames(mat),rownames(tax_dt)),]
split_rows <- factor(tax_dt$Class)

colnames(mat) <- gsub("19300","",met$PID)
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
pdf(paste(fig_folder,"Fig_1F.pdf",sep = "/"), height = 18, width = 20,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


## Correlation plot 1 E

sm_dt <-  read.csv("../data/metadata/Sample_metadata_MDR_04aug2023.csv")
sm_dt <-  sm_dt %>%
  select(PID,Days_after_Last_Dose_Date_that_stool_was_collected) %>%
  unique()

sm_dt$Days_after_Last_Dose_Date_that_stool_was_collected <-  as.numeric(sm_dt$Days_after_Last_Dose_Date_that_stool_was_collected)
sm_dt <-  sm_dt[complete.cases(sm_dt),]

# Extract diversity measures
phy_mic <- readRDS("../Processed_data/phy_mtph_counts.rds")

# Diversity plots:
library(yingtools2)
library(purrr)
library(tidyverse)
div_dt <- get.samp(phy_mic,stats = T) %>% as.data.frame()
div_dt <-  div_dt %>%
  select(PID, Visit, InvSimpson)%>%
  filter(PID %in% c(sm_dt$PID)) %>%
  filter(Visit %in% c("6 Mo","24 Mo"))%>%
  pivot_wider(id_cols = PID,names_from = Visit,values_from = InvSimpson) %>%
  data.frame()%>%
  mutate(L2FC = log2(X24.Mo/X6.Mo))

cor_dt <- sm_dt %>%
  left_join(div_dt)

names(cor_dt)[2] <- "Days_after_last_dose"



cor_dt$PID <- gsub("19300","",cor_dt$PID)

lm_model <- lm(data = cor_dt,formula = L2FC ~ Days_after_last_dose)
summary(lm_model)


cor_res <- cor.test(cor_dt$L2FC, cor_dt$Days_after_last_dose,method = "spearman")
pval <-  cor_res$p.value
r <-  as.numeric(cor_res$estimate)
label_test <- paste0( paste0("R = ", round(r,2)) , " pval = ",round(pval,4)  )




library(ggprism)
library(ggpmisc)
p <- ggplot(data =  cor_dt,aes(y = L2FC, x= Days_after_last_dose))+
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "adj.R2","p", "n")),label.x = "right",label.y = .75) +
  geom_point(size=2, shape = 21, alpha =1,fill = "blue",stroke = 1)+
  ggrepel::geom_text_repel(aes(label = PID))+
  theme_prism()+
  xlab("Days after last dose Abx")+
  ylab("L2FC(TC/6 Mo)")+
  #  annotate("label", y=3, x=50, label=label_test, hjust=0,size = 5) +
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=18,face = "bold"),
        axis.title.y = element_text(size=18,face = "bold"),
        axis.text.x = element_text(size=15,angle=45, hjust=1,vjust = 1,face = "bold"),
        axis.text.y = element_text(size=15,face = "bold"))


pdf(paste(fig_folder,"Fig_Cor_1_E.pdf",sep = "/"), height = 5,width = 5,useDingbats = F)
print(p)
dev.off()



