#Load libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(gtools)
library(reshape2)
library(ggplot2)
library(edgeR)
library(limma)

# Create a directory to save figures and tables
mainDir <- "../Figures"
subDir <- "Fig_2A"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,  recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")
# Import microbiome data 
# Clinical Trial
# RNASeq
phy_mtph <-  readRDS("../Processed_data/phy_mtph_mdr.rds")
phy_mtph <- prune_taxa(taxa_sums(phy_mtph)>0,phy_mtph)

# Now import RNASeq data for MDR:
phy_rna_mdr <-  readRDS("../Processed_data/phy_rna_mdr.rds")
sample_names(phy_rna_mdr) <- make.names(sample_names(phy_rna_mdr))

# Microbiome
phy_mic <- prune_taxa(taxa_sums(phy_mtph)>0,phy_mtph)
# Only keep samples that pairs with RNASeq data
phy_mic <-  prune_samples(sample_names(phy_rna_mdr), phy_mic)
phy_mic <-  prune_taxa(taxa_sums(phy_mic)>0, phy_mic)



# Only keep the needed info:
# Time, TTP, SampleID, PatientID
sm_dt <-  data.frame(sample_data(phy_rna_mdr))
sm_dt <- sm_dt[,c("PID","Visit","Av_TTP")]
sm_dt$SampleID <- rownames(sm_dt)

# Remove samples that don't have TTP value
sm_dt <- sm_dt[!is.na(sm_dt$Av_TTP),]

sm_dt$Time <- ifelse(sm_dt$Visit == "0 Day",0,
                     ifelse(sm_dt$Visit == "2 Wk",14,
                            ifelse(sm_dt$Visit == "2 Mo",60,
                                   ifelse(sm_dt$Visit == "6 Mo",180,720))))

# Assign sample data back to phyloseq obj
sample_data(phy_mic) <- sample_data(sm_dt) 
sample_data(phy_rna_mdr) <- sample_data(sm_dt)

phy_mic
phy_rna_mdr

sample_names(phy_mic) <-  make.names(paste(sample_data(phy_mic)$PID,"_D",
                                           sample_data(phy_mic)$Time,sep=""))
sample_names(phy_rna_mdr) <-  make.names(paste(sample_data(phy_rna_mdr)$PID,"_D",
                                               sample_data(phy_rna_mdr)$Time,sep=""))

# GSVA matrix
phy_gsva <-  readRDS("../Figures/Fig_S3/phy_GSVA_all_hallmark.rds")
sample_names(phy_gsva) <- make.names(sample_names(phy_gsva))
#phy_gsva <-  readRDS("../results_MDR/GSVA/phy_GSVA_all_hallmark.rds")
phy_gsva <- prune_samples(sm_dt$SampleID,phy_gsva)
sample_data(phy_gsva) <-  sample_data(sm_dt)
sample_names(phy_gsva) <-  make.names(paste(sample_data(phy_gsva)$PID,"_D",
                                            sample_data(phy_gsva)$Time,sep=""))

phy_mic
phy_rna_mdr
phy_gsva

# Lets look a the timeline for these patients:
time_dt <-  data.frame(sample_data(phy_mic))
length(unique(time_dt$PID))
p_time <- ggplot(time_dt, aes(x = Visit,y = PID,
                              group = PID, label = Av_TTP))+
  geom_point(size = 4,alpha = 0.7,fill = "blue", shape = 21, stroke = 1)+
  geom_text(hjust=0.5, vjust=-0.5)+
  geom_line()+
  theme_bw()+
  ylab("Subjects")+
  theme(axis.text.x=element_text(size=14,angle = 45, hjust=1,face="bold"),
        axis.text.y=element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14,face="bold"))

print(p_time)

pdf(paste(fig_folder,paste0("Timeline_MDR",".pdf"),sep="/"),height = 7, width = 6)
print(p_time)
dev.off()



# Filter species that are prevalent 
perc_samples <-  0.15
phy_fil = filter_taxa(phy_mic, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)

saveRDS(phy_fil,paste(fig_folder,paste0("phy_mic",".rds"),sep="/"))

library(tidyverse)
####### DATA IMPORT AND PREP #####
phy_genes <- phy_gsva
phy_micr <- phy_fil

library(R6)
library(tidyverse)
####### MERF #####
# define an R6 class to store results
ml_res <- R6Class("ml_results",
                  public = list(
                    X_train = NULL,
                    Y_train = NULL,
                    cluster_train = NULL,
                    Z_train = NULL,
                    X_test = NULL,
                    Y_test = NULL,
                    Z_test = NULL,
                    cluster_test = NULL, 
                    importance = NULL, 
                    Yhat = NULL,
                    Yhat_val = NULL,
                    Yhat_rf_base = NULL,
                    Yhat_val_rf_base = NULL,
                    merf_model = NULL
                  )
)

ml_res_holder <-  R6Class("results_container",
                          public = list(
                            pathway = NULL,
                            Xdata_otu = NULL,
                            Ydata_otu = NULL,
                            clusterdata_otu = NULL,
                            clusterdata_otu_numeric = NULL,
                            Zdata_otu = NULL,
                            Tdata_otu = NULL,
                            rf_results = NULL
                          )
)

#adding some basic pred vars
pathways <- taxa_names(phy_genes) # gene ids
u_patient <- unique(sample_data(phy_micr)$PID)
# for every Pathway build a longitudinal rf model 
ipath <- 1


ml_res_holder_list <- list()
for (ipath in seq(1,length(pathways))){
  #for (iotu in seq(1,2)){  
  ml_res_holder_iotu <- NULL
  ml_res_holder_iotu <- ml_res_holder$new()
  
  # set the matrices for the modeling
  Ydata <- otu_table(phy_genes)[pathways[ipath],] ## THIS IS THE GENE
  Ydata <- as.data.frame(Ydata)
  Ydata <- t(Ydata)
  Xdata <- data.frame(sample_data(phy_micr))
  Xdata <- rownames_to_column(Xdata)
  
  # add microbes to Xdata
  xtmp <- data.frame(t(otu_table(phy_micr)))
  xtmp <- rownames_to_column(xtmp)
  Xdata <- merge(Xdata, xtmp, by = "rowname")
  rownames(Xdata) <- Xdata$rowname 
  Xdata$rowname <- NULL
  Xdata <-  Xdata[,6:ncol(Xdata)]
  
  # Order Xdata as Ydata
  Xdata <-  Xdata[match(rownames(Ydata),rownames(Xdata)),]
  
  clusterdata <- data.frame(sample_data(phy_micr)[,c("PID","Av_TTP","Time")])
  clusterdata <-  clusterdata[match(rownames(Xdata),rownames(clusterdata)),]
  Xdata <- cbind(Xdata,TTP = log(clusterdata$Av_TTP))
  
  
  
  Zdata <- rep(1,length(clusterdata$PID))
  Tdata <- clusterdata$Time 
  
  clusterdata_numeric <- as.character(clusterdata$PID)
  
  # assign to object
  ml_res_holder_iotu$Xdata_otu<-Xdata
  ml_res_holder_iotu$Ydata_otu<-Ydata
  ml_res_holder_iotu$clusterdata_otu<-clusterdata
  ml_res_holder_iotu$clusterdata_otu_numeric<-clusterdata_numeric
  ml_res_holder_iotu$Zdata_otu<-Zdata
  ml_res_holder_iotu$Tdata_otu<-Tdata
  ml_res_holder_iotu$pathway <- pathways[ipath]
  
  mlo_k <- NULL
  mlo_k <- ml_res$new()
  
  # train data
  Y_train <- Ydata
  X_train <- as.matrix(Xdata) 
  X_train <- as.matrix(as.data.frame(X_train))
  Z_train <- as.matrix(Zdata)
  cluster_train <-clusterdata_numeric 
  Time_train <- Tdata
  
  library(vita)
  library(LongituRF)
  
  if (all(!is.na(Y_train)) ){
    # mixed effect random forest
    set.seed(1057)
    model_merf<-MERF(X = X_train, Y = Y_train, id = cluster_train, 
                     Z =  Z_train, time = Time_train, sto = "none", 
                     mtry = round((1/3)*ncol(X_train)), ntree = 3000, iter = 100)
    mlo_k$merf_model <- model_merf
    
    # permutated importance plot
    set.seed(300)
    pimp <- PIMP(X_train, Y_train, model_merf$forest, 
                 parallel=TRUE, ncores = 20, seed = 300)
    pimp_test <- PimpTest(pimp)
    pimp_all <- data.frame(orig =  pimp_test$VarImp, pimp_test$pvalue)
    imp_df<- pimp_all     
    imp_df$Predictors <- rownames(imp_df)
    rownames(imp_df) <- NULL
    
    # predict training data
    # Yhat <- predict.merf(merf_sto = model_merf, X = X_train, 
    #                      Z = Z_train, id = cluster_train, time = Time_train)
    mlo_k$X_train <- X_train
    mlo_k$Y_train <- Y_train
    mlo_k$Z_train <- Z_train
    mlo_k$cluster_train <- cluster_train
    mlo_k$importance <- imp_df
    #mlo_k$Yhat <- Yhat
    
  }else{
    mlo_k <- NA
  }
  ml_res_holder_iotu$rf_results<-mlo_k
  ml_res_holder_list[[ipath]]<-ml_res_holder_iotu
}
#ml_res_holder_list_bup <- ml_res_holder_list

saveRDS(ml_res_holder_list, paste0(fig_folder,"/","ml_res_holder_list.rds"))
#save.image(paste0(fig_folder,"/","ml_res.rdata"))



#######################Summary of MERF model ####################################


df_imp<-c()
i<-1
for (i in seq(1:length(ml_res_holder_list))){
  instance <- ml_res_holder_list[[i]]
  #if (is.na(instance$rf_results)==F){
  temp<-instance$rf_results$importance
  temp$iter<-i  
  temp$pathway<- instance$pathway
  df_imp<-rbind(df_imp,temp)
  #}
}
df_imp<-na.omit(df_imp) 
df_imp_filtered <- df_imp[df_imp$VarImp>0,] 
unique(df_imp_filtered$pathway) 

# Widen the data frame:
head(df_imp_filtered)
library(tidyverse)

imp_dt_w <- df_imp_filtered %>%
  select(Predictors,pathway,VarImp)%>%
  pivot_wider(names_from = pathway,values_from = VarImp)


# Directionality
# Use LME to get the slope or the relation between dependent variables:
imp_bugs <- imp_dt_w$Predictors
#phy_mic <- readRDS("../Figures/Merf_Hallmark/phy_mic.rds")
phy_fil_mic <-  prune_taxa(imp_bugs,phy_mic)
sample_names(phy_fil_mic)

met_dt <-  data.frame(sample_data(phy_fil_mic))
mic_dt <- psmelt(phy_fil_mic)
mic_dt$Sample_ID <-  mic_dt$Sample
library(tidyverse)
mic_dt_w <-  mic_dt %>%
  dplyr::select(OTU,Sample_ID,Abundance)%>%
  pivot_wider(names_from = OTU, values_from = Abundance)

phy_gene_mod <- readRDS("../Figures/Fig_S3/phy_GSVA_all_hallmark.rds")
sample_names(phy_gene_mod) <-  make.names(paste(sample_data(phy_gene_mod)$PID,"_D",
                                                sample_data(phy_gene_mod)$Time,sep=""))
phy_gene_mod <-  prune_samples(sample_names(phy_fil_mic),phy_gene_mod)

met_dt <-  data.frame(sample_data(phy_gene_mod))
gene_dt <- psmelt(phy_gene_mod)
gene_dt$Sample_ID <-  gene_dt$Sample


merge_mic_dt <-  gene_dt %>%  inner_join(mic_dt_w)
length(unique(merge_mic_dt$Sample_ID))
total_taxa <-  taxa_names(phy_fil_mic)
# For Microbiome
result_lme <-  list()
gene_modules <- taxa_names(phy_gene_mod)
gene_m <- gene_modules[1]
for(gene_m in gene_modules){
  
  print(gene_m)
  lme_dt <-  merge_mic_dt %>% filter(OTU == gene_m)
  lme_dt$PID <- factor(lme_dt$PID)
  lme_dt$Time 
  # Remove data with NAs
  #lme_dt <-  lme_dt[!is.na(lme_dt$Average.TTP),]
  lme_dt$Time <-  scale(lme_dt$Time,center = F)
  library("nlme")
  
  result_lme_sp <-  list()
  sp <- total_taxa[1]
  for(sp in total_taxa){
    
    
    sum_mod_dt <- tryCatch({
      print(sp)
      sp_dt <- data.frame(lme_dt)
      sp_dt$t_Abundance <-  asin(sqrt(sp_dt[,c(make.names(sp))]))
      
      library("nlme")
      mod_bac <- lme(Abundance ~ t_Abundance + Time, random = ~1 |PID, sp_dt)
      sum_mod <-  summary(mod_bac)
      sum_mod_dt <- data.frame(sum_mod$tTable)
      sum_mod_dt$Module <- gene_m
      sum_mod_dt$Var <-  rownames(sum_mod_dt)
      sum_mod_dt$Species <-  sp
      result_lme_sp[[sp]] <-  sum_mod_dt
      
    }, error = function(err){
      sum_mod_dt <- NA
      sum_mod_dt
    }
    )
    result_lme_sp[[sp]] <-  sum_mod_dt
    
  }
  result_lme_dt_sp <-  do.call("rbind",result_lme_sp)
  # result_lme_dt_sp$adj_pval <-  p.adjust(result_lme_dt_sp$p.value,method = "BH")
  # Remove intercept
  #  result_lme_dt_sp <-  result_lme_dt_sp[result_lme_dt_sp$Var %in% "t_Abundance",]
  
  result_lme[[gene_m]] <- result_lme_dt_sp
  
}

res_dt <- do.call("rbind",result_lme)

#res_dt <-  res_dt[res_dt$adj_pval < 0.1,]
res_dt$adj_pval <-  p.adjust(res_dt$p.value,method = "BH")
# Remove intercept
res_dt <-  res_dt[res_dt$Var %in% "t_Abundance",]
# res_dt <-  res_dt[res_dt$adj_pval < 0.05,]
# table(res_dt$Module)
#res_dt$Relation <-  ifelse(res_dt$Value > 0, 
res_dt_w <-  res_dt %>%
  select(Species,Value,Module)%>%
  pivot_wider(names_from = Module, values_from = Value)


# Include Sex and TTP
# For TTP
result_lme <-  list()
#merge_dt
gene_m
for(gene_m in gene_modules){
  
  lme_dt <-  gene_dt %>% filter(OTU == gene_m)
  # Remove data with NAs
  lme_dt <-  lme_dt[!is.na(lme_dt$Av_TTP),]
  lme_dt$Time <-  scale(lme_dt$Time,center = F)
  library("nlme")
  mod_bac <- lme(Abundance ~ log( Av_TTP) +  Time ,random = ~1 |PID, lme_dt)
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Module <- gene_m
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[gene_m]] <-  sum_mod_dt  
}

res_dt_TTP <- do.call("rbind",result_lme)
res_dt_TTP$adj_pval <-  p.adjust(res_dt_TTP$p.value,method = "BH")
# Remove intercept
#res_dt_TTP <-  res_dt_TTP[!res_dt$Var %in% "(Intercept)",]
res_dt_TTP <-  res_dt_TTP[res_dt_TTP$Var %in% c("log(Av_TTP)"),]

#res_dt[res_dt$adj_pval < 0.05,]
res_dt_w_TTP <-  res_dt_TTP %>%
  select(Var,Value,Module)%>%
  pivot_wider(names_from = Module, values_from = Value)



names(res_dt_w)[1] <- "Var"
comb_res_dt <- rbind(res_dt_w,res_dt_w_TTP)
comb_res_dt$Var[comb_res_dt$Var == "log(Av_TTP)"] <- "TTP"


dir_dt <- comb_res_dt %>%
  pivot_longer(!Var,names_to = "Modules", values_to = "dir")


comp_imp <- df_imp_filtered
# Top predictor 
var_sel <- comp_imp %>%
  group_by(pathway) %>%
  filter(VarImp == max(VarImp))

# Select only top 10
sel_imp <- comp_imp %>%
  dplyr::group_by(pathway) %>%
  dplyr::arrange(desc(VarImp)) %>%
  dplyr::slice_max(n= 10,VarImp)%>% data.frame()

# sel_imp <- sel_imp %>%
#            left_join(imp_slope, by = c("Predictors" = "pred","pathway" = "Genes"))

sel_imp <- sel_imp %>%
  left_join(dir_dt, by = c("Predictors" = "Var","pathway" = "Modules"))


# sel_imp$value <- 1
# sel_imp$value[sel_imp$Slope < 0 ] <- -1 
# sel_imp <- sel_imp[,c("Predictors","pathway","value")]
sel_imp$value <- 1
sel_imp$value[sel_imp$dir < 0 ] <- -1 
sel_imp <- sel_imp[,c("Predictors","pathway","value")]


library(tidyr)

comp_imp_w <- sel_imp %>%
  pivot_wider(names_from = pathway, values_from = value,values_fill = list(value = 0))%>%
  data.frame()

rownames(comp_imp_w) <-  comp_imp_w$Predictors
comp_imp_w$Predictors <- NULL


library(RColorBrewer)
library(circlize)
# Lets only draw the immune pathways from lme analysis
mat <- as.matrix(comp_imp_w) 
#col_binary <-  colorRampPalette(c('#7fc97f','white','#d95f02'))(3)

col_binary <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdYlBu")))

col_binary <-  c("#91BFDB","white","#FC8D59")



#sel_cols = colorRampPalette(c("white", "darkblue"))(2)
col_bar = colorRamp2(c(0,1), c("white", "darkblue"))


# Display All 
# Loop across the pathways
p_ann <-  read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv")
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
immune_path <- p_ann$Hallmark.Name[p_ann$Process.Category %in% "immune"]

#lme_pathway 
sel_imm_path <- p_ann

mat_imm <- mat[,colnames(mat) %in% sel_imm_path$Hallmark.Name]


sel_p_ann <- p_ann[p_ann$Hallmark.Name %in% colnames(mat_imm),]
sel_p_ann <- sel_p_ann[match(colnames(mat_imm),sel_p_ann$Hallmark.Name),]
# # Match it with the colnames of the matrix
# rsq_dt <-  rsq_dt[match(colnames(mat),rsq_dt$Comp),]
# 
sp_cols <- factor(sel_p_ann$Process.Category)

mat_imm <- mat_imm[apply(mat_imm, 1, function(x) !all(x==0)),]


# Find the index of  top predictor for each pathway
head(var_sel)

imm_var_sel <- var_sel[var_sel$pathway %in% colnames(mat_imm),]


var <- imm_var_sel$pathway[1]

idx_c <- c()
idx_r <- c()
for (var in imm_var_sel$pathway){
  
  spec <- imm_var_sel$Predictors[imm_var_sel$pathway %in% var]
  idx_c <- c(idx_c,which(colnames(mat_imm) == var))
  idx_r <- c(idx_r,which(rownames(mat_imm) == spec))
}



# Add tax dt info
library(phyloseq)
phy_ht <- readRDS("../Processed_data/phy_mtph_mdr.rds")
phy_ht <-  prune_taxa(rownames(mat_imm),phy_ht)

library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])
topN <- length(unique(tax_dt$Order))
match_tax_dt <- data.frame(tax_table(phy_ht))
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]
rownames(match_tax_dt)


phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
                "Spirochaetes" = "grey50",
                "Candidatus Melainabacteria" = "#999999")

dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
shades_num <- match_tax_dt %>%
  group_by(Phylum) %>%
  mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col )
tax_sel <-  unique(match_tax_dt[,c("Phylum","Order")])

#shades_num$phy_col <- as.vector(phylum_col[order(names(phylum_col))])
order_col <- mapply(function(x, y){ shades(color = x, ncolor = y, variation = 1)}, x= shades_num$phylum_col ,y = shades_num$unique_types,
                    SIMPLIFY = F)

order_col <- as.vector(unlist(order_col))
names(order_col) <- tax_sel$Order

mycol <- order_col

# Add a row for TTP
match_tax_dt[nrow(match_tax_dt)+1,] <- NA
rownames(match_tax_dt)[nrow(match_tax_dt)] <- "TTP"
match_tax_dt <- match_tax_dt[match(rownames(mat_imm),rownames(match_tax_dt)),]

library(ComplexHeatmap)

ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")

split_rows <-  match_tax_dt$Phylum
split_rows[is.na(split_rows)] <- "TTP"
split_rows <- factor(split_rows, levels= c("TTP","Actinobacteria","Firmicutes","Bacteroidetes","Proteobacteria"))


# Replace dots with space
rownames(mat_imm)  <- gsub("\\."," ",rownames(mat_imm))
rownames(mat_imm)  <- gsub("X","",rownames(mat_imm))

colnames(mat_imm) <- gsub("HALLMARK_","",colnames(mat_imm))
colnames(mat_imm) <- gsub("_"," ",colnames(mat_imm))

# Find the index where t


ht1  = Heatmap(mat_imm,name = "Relation",
               column_split = sp_cols,
               row_split =  split_rows, row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 90,
               cluster_row_slices = F,
               left_annotation = ha1,
               show_parent_dend_line = F,
               border = T,
               rect_gp = gpar(col = "gainsboro"),
               row_names_side = "left",
               show_row_dend = F,
               show_column_dend = F,
               col = col_binary,
               row_names_max_width = max_text_width(rownames(mat_imm),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_imm),
                                                        gp = gpar(fontsize = 12)),
               heatmap_legend_param = list(
                 title = "Relation", at = c(-1, 1), 
                 labels = c("Neg", "Pos")
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 #if(mat_imm[i, j] > 0)
                 for (k in 1:length(idx_r)){
                   if(i == idx_r[k] & j == idx_c[k]){
                     print(paste0(i," ",j))
                     grid.points(x, y, pch = 16, size = unit(4, "mm"))
                   }
                 }
               }
               # layer_fun = function(j, i, x, y, w, h, fill) {
               #   ind_mat = restore_matrix(j, i, x, y)
               #   
               #   
               #   ind = ind_mat[1, ]
               #   grid.points(x[ind], y[ind], pch = 16, size = unit(4, "mm"))
               # }
               
)


pdf(paste0(fig_folder,"/Fig_2A.pdf"),height = 15, width = 15,useDingbats = F)
draw(ht1)
dev.off()













