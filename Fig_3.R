# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_3"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_3"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


# Import butyrate encoding gene's profile across MDR samples
CARD_genes <-  read.csv("../data/CARD/CARD_MDR.csv")
colnames(CARD_genes) <- gsub("_S.*","",colnames(CARD_genes))
rownames(CARD_genes) <- CARD_genes$X
CARD_genes$X <- NULL
CARD_genes

library(stringr)
gene_id_dt <- data.frame(str_split_fixed(rownames(CARD_genes),"\\|",4))
gene_id_dt <-  data.frame(gene_id = rownames(CARD_genes),gene_id_dt)
gene_id_dt <-  gene_id_dt[,c(1,4,5)]
names(gene_id_dt)[c(2,3)] <- c("ARO_id","Gene") 
gene_id_dt$ARO_id <- gsub("_",":",gene_id_dt$ARO_id)
## Metadata for CARD:

# Read the metadata for the genes
# genes metadata
genes_met_cat_file <- "../data/metadata_CARD/aro_categories.tsv"
genes_met_cat <- read.csv(genes_met_cat_file, sep ="\t")

genes_met_cat_index_file <-"../data/metadata_CARD/aro_categories_index.tsv"
genes_met_cat_index <- read.csv(genes_met_cat_index_file, sep ="\t")

genes_met_index_file <- "../data/metadata_CARD/aro_index.tsv"
genes_met_index <- read.csv(genes_met_index_file, sep="\t")


length(unique(genes_met_index$ARO.Accession))
length(unique(gene_id_dt$ARO_id))
table(unique(gene_id_dt$ARO_id) %in% unique(genes_met_index$ARO.Accession))

head(genes_met_index)
head(genes_met_cat_index)

library(tidyverse)
meta_genes <-  gene_id_dt %>%
  inner_join(genes_met_index, by = c("ARO_id"= "ARO.Accession"))
rownames(meta_genes) <-  meta_genes$gene_id

# Read homolog fasta files:
# for all the ARG look where they came from and add
library(ShortRead)
fasta_files <-  list.files("../data/metadata_CARD",pattern = ".fasta",full.names = T)
fasta_list <-  lapply(fasta_files,readAAStringSet)

comb_fasta <-  do.call(c,fasta_list)
seq_name = names(comb_fasta)
# sequence = paste(fastaFile)
fasta_df <- data.frame(seq_name, data.frame(str_split_fixed(seq_name,"\\|",4)))
fasta_df$microbe <- sub(".*\\[(.*)\\].*", "\\1", fasta_df$seq_name, perl=TRUE)
fasta_df <- fasta_df[,c("seq_name","X3","microbe")]
names(fasta_df)[2] <- "ARO_id"

table(meta_genes$ARO_id %in% fasta_df$ARO_id)
meta_genes <-  meta_genes %>%
                left_join(fasta_df)


library(phyloseq)
# Import phyloseq for metaphlan:
phy_tax <-  readRDS("../Processed_data/phy_mtph_mdr.rds")
sm_dt <- data.frame(sample_data(phy_tax))
# Change the level of last time point as TC
levels(sm_dt$Visit)[6] <- "TC"
sample_data(phy_tax) <-  sample_data(sm_dt)
tax_table(phy_tax) <- NULL
otu_table(phy_tax) <-  otu_table(as.matrix(CARD_genes),taxa_are_rows = T)
phy_tax <- prune_taxa(taxa_sums(phy_tax) > 0,phy_tax )

# Visualize the data:
# Composition over time:
# Only Post :
phy_ht <- phy_tax
# Remove less prevalent species across all the samples:
# Overall abundance of ASVs
sm_dt <-  data.frame(sample_data(phy_ht))
phy_ht <-  prune_taxa(taxa_sums(phy_ht)>0,phy_ht)


# Relative abundance:
mat <- data.frame(otu_table(phy_ht))
met <- data.frame(sample_data(phy_ht))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))

library(ComplexHeatmap)
library(tidyverse)
library(gtools)
library(RColorBrewer)
library(circlize)

mixedsort(met$Visit)
ha_column = HeatmapAnnotation(Status =  met$Visit,
                              Sex =  met$Sex,
                              col=list(Status = status_col,
                                       Sex = c("Male" = "white","Female" = "black")))

split_cols<-  met$Visit


col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(phy_ht)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_cfu <- scico(5, palette = 'lajolla')
col_time <- scico(2,palette = "devon")
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- scico(10, palette = 'tofino')
col_mat <- scico(10, palette = 'lajolla')
library(ComplexHeatmap)
col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
colnames(mat) <- gsub("19300","",met$PID)




ht  =  Heatmap(log(mat+ 1),name = "Log(RPKM+1)",
               column_split = splitcols,
               # row_split =  split_rows,
               #         top_annotation = ha_column,
               # left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T, heatmap_legend_param = list(title = "Log(RPKM + 1)", 
                                                       legend_width = unit(10, "cm"),
                                                       legend_direction = "horizontal",
                                                       title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                       labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 10)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 10)))
# pdf(paste0(results_folder,"/Heatmap_CARD.pdf"),height = 30, width = 25,useDingbats = F)
# draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
# dev.off()




# Differential analysis over time compared to day 0:
phy_lme <- phy_tax
# LME analysis:
# Keep only prevalence species:
### filter by prevalence
prevdf = apply(X = otu_table(phy_lme),
               MARGIN = ifelse(taxa_are_rows(phy_lme), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy_lme))
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
  phy_m$t_Abundance  <- log(phy_m$Abundance + 1)
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

# Remove intercept
result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
result_lme_dt <- result_lme_dt[grep("Visit",result_lme_dt$Var),]
#result_lme_dt <- result_lme_dt[result_lme_dt$Var != "(Intercept)",]

result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

result_lme_dt_sub <- result_lme_dt
#result_lme_dt_sub <- result_lme_dt[grep('6 Mo',rownames(result_lme_dt)),]

p_cutoff <-  0.05
result_sig <- result_lme_dt_sub[result_lme_dt_sub$p_adj < p_cutoff ,]
unique(result_sig$Bac)

# Save the LME results
write.csv(result_sig,paste(tab_folder,"Tab_LME_CARD.csv",sep = "/"))

# Out of 113 ABX gene:
# neg_res <-  result_sig[result_sig$Value < 0,]
# unique(neg_res$Bac)

# Filter by negative coefficient (Lower in Post compared to baseline)
ps_sig <- prune_taxa(unique(result_sig$Bac),phy_fil)

# Relative abundance:
mat <- data.frame(otu_table(ps_sig))
met <- data.frame(sample_data(ps_sig))

met <-  met[order(met$Visit,met$PID),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)

library(RColorBrewer)
status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
names(status_col) <- unique(unique(met$Visit))
mixedsort(met$Visit)

library(ComplexHeatmap)

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



two_week_num <-  table(coeff_dt_w$Visit2.Wk > 0)
one_month_num <-  table(coeff_dt_w$Visit1.Mo > 0)
two_month_num <-  table(coeff_dt_w$Visit2.Mo > 0)
six_month_num <-  table(coeff_dt_w$Visit6.Mo > 0)
TC_num <-  table(coeff_dt_w$VisitTC > 0)

gene_num_dt <-  data.frame(rbind(two_week_num,one_month_num,two_month_num,six_month_num,TC_num))
names(gene_num_dt) <-  c("Depressed","Enriched")
gene_num_dt$time <-  gsub("_num","",rownames(gene_num_dt))

gene_num_dt <- pivot_longer(gene_num_dt, cols = c(Enriched, Depressed))
gene_num_dt$time <-  factor(gene_num_dt$time,levels = c("two_week","one_month","two_month","six_month","TC"))
levels(gene_num_dt$time) <- c("2 Wk","1 Mo","2 Mo","6 Mo","TC")

# Save the LME results
write.csv(gene_num_dt,paste(tab_folder,"Num_ARGS_Time.csv",sep = "/"))


library(ggrepel)
library(ggprism)
# Line plot :
p_line <- ggplot() +
  geom_point(data  = gene_num_dt,aes(x = time,y = value,fill = name),
             alpha = 0.7,shape = 21, stroke = 1, color = "black", size = 3) +
  geom_line(data  = gene_num_dt,aes(x = time,y = value,group = name,color = name),
            alpha = 0.7)+
  geom_text_repel(data  = gene_num_dt,aes(x = time,y = value,label = value))+
  ylab("Number of ARGs")+
  xlab("Time")+
  scale_fill_manual(name = "Pattern",values = c("Enriched" = "violetred","Depressed" = "darkblue"))+
  scale_color_manual(name = "Pattern",values = c("Enriched" = "violetred","Depressed" = "darkblue"))+
  theme_prism()+
  theme(axis.text.x=element_text(size=10,angle = 45, vjust = 1, hjust=1,face="bold"),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14))
pdf(paste0(results_folder,"/Fig_3B.pdf"), width = 6, height = 6)
print(p_line)
dev.off()



library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(ps_sig)

library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
dist_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6",
                 "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3","#808000","#ffd8b1",
                 "#000080")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



library(scico)
col_cfu <- scico(5, palette = 'lajolla')
col_time <- scico(2,palette = "devon")
col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","1 Mo"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple","TC" = "purple")
splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))

col_mat <- scico(10, palette = 'tofino')
col_mat <- scico(10, palette = 'lajolla')
library(ComplexHeatmap)
col_mat <- viridis::viridis(10)
#ha_column = HeatmapAnnotation(Treatment = met$treatment)
max(coeff_dt$Value)
min(coeff_dt$Value)


col_grad <- colorRampPalette(c("blue", "white", "red"))(5)

col_coef <- colorRamp2(c(-6,-3, 0, 3,6),col_grad) 

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
                                                        at =c(-6,-3, 0, 3,6),
                                                        legend_height = unit(10, "cm")))

colnames(mat) <- gsub("19300","",met$PID)

sel_meta_dt <-  meta_genes[meta_genes$gene_id %in% rownames(mat), ]
# Order sel_meta_dt by mat
sel_meta_dt <- sel_meta_dt[match(rownames(mat),sel_meta_dt$gene_id),]
# drug_class_dt <- data.frame(str_split_fixed(sel_meta_dt$Drug.Class,";",n = 15))
drug_class_dt <- str_split_fixed(sel_meta_dt$Drug.Class,";",n = 15)

vec_class <- unique(as.vector(drug_class_dt))
# Remove empty elements
vec_class <-  vec_class[!vec_class %in% ""]
# Make a data frame of unique antibiotics class:

drug_ann_dt <-  data.frame(matrix(data = 0,nrow = nrow(drug_class_dt), ncol = length(vec_class) ) )
names(drug_ann_dt) <- vec_class
row_num <-  40
for(row_num in 1:nrow(drug_class_dt)){
    # Vector in each row
  vec_comp <-  as.character(drug_class_dt[row_num,])
  drug_ann_dt[row_num,which(vec_class %in% vec_comp)] <- 1
}


# Combine columns in dataframe into one
t_drug_ann_dt <-  data.frame(t(drug_ann_dt))
t_drug_ann_dt$sub_class <- rownames(t_drug_ann_dt)
# Meta ann:
meta_ann <-  data.frame(
              Drug_Class =  c("Tetracycline",
                              "Tetracycline",
                              "Aminoglycoside",
                              "Nitroimidazole",
                              "Fluoroquinolones",
                              "Aminocoumarin",
                              "Diaminopyrimidine",
                              "Fosfomycin",
                              "Rifamycin",
                              "Macrolide",
                              "Nucleoside",
                              "Sulfonamides",
                              "Sulfonamides",
                              "Beta-Lactams",
                              "Beta-Lactams",
                              "Beta-Lactams",
                              "Beta-Lactams",
                              "Beta-Lactams",
                              "Beta-Lactams",
                              "Oxazolidinones",
                              "LSPP",
                              "LSPP",
                              "LSPP",
                              "LSPP",
                              "Nitrofuran",
                              "Other",
                              "Other",
                              "Other",
                              "Other"),
              sub_class = c("tetracycline antibiotic","glycylcycline",
                            "aminoglycoside antibiotic",
                            "nitroimidazole antibiotic",
                            "fluoroquinolone antibiotic",
                            "aminocoumarin antibiotic",
                            "diaminopyrimidine antibiotic",
                            "fosfomycin",
                            "rifamycin antibiotic",
                            "macrolide antibiotic",
                            "nucleoside antibiotic",
                            "sulfone antibiotic","sulfonamide antibiotic",
                            "carbapenem","cephalosporin","cephamycin","monobactam","penem","penam",
                            "oxazolidinone antibiotic",
                            "lincosamide antibiotic","streptogramin antibiotic","phenicol antibiotic","pleuromutilin antibiotic",
                            "nitrofuran antibiotic",
                            "acridine dye","peptide antibiotic","disinfecting agents and intercalating dyes","triclosan")
)


merge_ann_dt <-  t_drug_ann_dt %>%
                  left_join(meta_ann)
merge_ann_dt$sub_class <- NULL 
# Now aggregate by drug class
agg_ann_dt  <-  aggregate(.~Drug_Class ,merge_ann_dt,sum)
rownames(agg_ann_dt) <- agg_ann_dt$Drug_Class
agg_ann_dt$Drug_Class <-  NULL

agg_ann_dt <- data.frame(t(agg_ann_dt))
agg_ann_dt[agg_ann_dt > 0] <-  1
col_drug_class <-  c("0" = "white","1" = "#2048a6")

# Right annotations
ha_right <- rowAnnotation(Aminocoumarin = agg_ann_dt$Aminocoumarin,
                          Aminoglycoside = agg_ann_dt$Aminoglycoside,
                          Beta_Lactams = agg_ann_dt$Beta.Lactams,
                          Diaminopyrimidine = agg_ann_dt$Diaminopyrimidine,
                          Fluoroquinolones = agg_ann_dt$Fluoroquinolones,
                          Fosfomycin = agg_ann_dt$Fosfomycin,
                          Oxazolidinones = agg_ann_dt$Oxazolidinones,
                          LSPP = agg_ann_dt$LSPP,
                          Macrolide =  agg_ann_dt$Macrolide,
                          Nitrofuran =  agg_ann_dt$Nitrofuran,
                          Nitroimidazole = agg_ann_dt$Nitroimidazole,
                          Nucleoside =  agg_ann_dt$Nucleoside,
                          Rifamycin = agg_ann_dt$Rifamycin,
                          Sulfonamides = agg_ann_dt$Sulfonamides,
                          Tetracycline = agg_ann_dt$Tetracycline,
                          Other = agg_ann_dt$Other,
                          col = list(Aminocoumarin = col_drug_class,
                                     Aminoglycoside = col_drug_class,
                                     Beta_Lactams = col_drug_class,
                                     Diaminopyrimidine = col_drug_class,
                                     Fluoroquinolones = col_drug_class,
                                     Fosfomycin = col_drug_class,
                                     Oxazolidinones = col_drug_class,
                                     LSPP =col_drug_class,
                                     Macrolide = col_drug_class,
                                     Nitrofuran =  col_drug_class,
                                     Nitroimidazole = col_drug_class,
                                     Nucleoside =  col_drug_class,
                                     Rifamycin = col_drug_class,
                                     Sulfonamides = col_drug_class,
                                     Tetracycline = col_drug_class,
                                     Other = col_drug_class),
                          show_legend = c(TRUE,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
                          annotation_name_gp = gpar(fontface = "bold" ),
                          annotation_name_rot = 75,
                          annotation_legend_param = list(title = "Presence",
                                                         title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                         labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                         legend_height = unit(8, "cm")))


micr_genes <-  sel_meta_dt$microbe
micr_genes <-  gsub("Escherichia coli.*","Escherichia coli",micr_genes)
micr_genes <-  gsub("Klebsiella pneumoniae.*","Klebsiella pneumoniae",micr_genes)
micr_genes <-  gsub("Enterococcus faecalis.*","Enterococcus faecalis",micr_genes)
micr_genes <-  gsub("Vibrio cholerae.*","Vibrio cholerae",micr_genes)
micr_genes <-  gsub("Salmonella enterica.*","Salmonella enterica",micr_genes)
micr_genes <-  gsub(" .*","",micr_genes)
unique(micr_genes)
split_rows <- factor(micr_genes)

mat_aro <- mat 
rownames(mat_aro) <-  sel_meta_dt$ARO_id
ht  =  Heatmap(log(mat_aro+ 1),name = "Log(RPKM + 1)",
               column_split = splitcols,
               #row_split =  split_rows,
               row_split =  split_rows,
               right_annotation = ha_right,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T, heatmap_legend_param = list(title = "Log(RPKM + 1)", 
                                                       legend_width = unit(10, "cm"),
                                                       legend_direction = "horizontal",
                                                       title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                       labels_gp = gpar(fontsize = 12,fontface = "bold" )))
pdf(paste0(results_folder,"/Fig_3A.pdf"),height = 25, width = 24,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# Now lets make a piechart of the increasing ARGs
mat_meta <-  agg_ann_dt
mat_meta <-  data.frame(ARG_name =  rownames(mat), mat_meta)

# For each time point find all the ARGs that go up and sum their groups
timepoints <-  names(coeff_dt_w)[2:6]
time_comp <-  timepoints[1]

row_list <-  list()
for(time_comp in timepoints){

sel_coeff_dt <-  coeff_dt_w %>%
                 filter((!!as.symbol(time_comp))  > 0)
sel_args <-  sel_coeff_dt$Bac
sel_mat_meta <-  mat_meta[mat_meta$ARG_name  %in% sel_args,]
counts_sum <- sel_mat_meta %>%
  dplyr::summarise_at(vars(-ARG_name), sum)  
row_list[[time_comp]] <-  counts_sum
}

pie_dt <-  do.call("rbind",row_list)
pie_dt <-  data.frame(time = rownames(pie_dt),pie_dt)
pie_dt <-  pie_dt[,c(names(pie_dt)[!names(pie_dt) %in% c("Other")],"Other")]


pie_dt_l <-  pie_dt %>%
             pivot_longer(!time,names_to = "Cat",values_to = "count")
pie_dt_l$time <-  factor(pie_dt_l$time, levels = c("Visit2.Wk","Visit1.Mo",  "Visit2.Mo","Visit6.Mo", "VisitTC"))
levels(pie_dt_l$time) <- c("2_Wks","1 Mo","2 Mo","6 Mo", "TC")

# Use plain ggplot2:
# Load ggplot2
library(ggplot2)
library(dplyr)

# Compute the position of labels
pie_dt_l <- pie_dt_l %>%
   group_by(time)%>%
  mutate(prop = count / sum(count) *100) %>%
  mutate(ypos = 0.5*prop + lead(cumsum(prop), 1),
         ypos=  if_else(is.na(ypos), prop/2, ypos),#ypos = cumsum(prop)- 0.5*prop,
         labels = paste0(round(prop, 1), "%"))

# Save the PieChart
write.csv(pie_dt_l,paste(tab_folder,"Tab_PieChart.csv",sep = "/"))


# Only display top 8 categories in each group
library(tidyverse)
disp_dt <-  pie_dt_l %>% 
              group_by(time) %>% 
              arrange(desc(prop)) %>%
              dplyr::slice(1:8)
disp_dt$labels_mod <-  disp_dt$labels

pie_dt_merge <-  pie_dt_l %>%
              left_join(disp_dt)
pie_dt_merge$labels_mod[is.na(pie_dt_merge$labels_mod)] <- ""
# Basic piechart

library(hues)
set.seed(1057)
cat_abx <- names(pie_dt)[2:(ncol(pie_dt) )]
col_cat <-  iwanthue(length(cat_abx), plot=TRUE)
names(col_cat) <- c(cat_abx[!cat_abx %in% "Other"],"Other")
col_cat[length(col_cat)] <- "grey80"






pie_chart <-  ggplot(pie_dt_merge, aes(x= "", y=prop, fill=Cat)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  geom_text(data = pie_dt_merge,  aes(x = 1.75,label = labels_mod),position = position_stack(vjust=0.5), color = "black", size= 3) +
  theme(legend.position="none") +
  scale_fill_manual(name = "Categories",values = col_cat)+
 # scale_x_discrete(limits = c("2_Wks","1 Mo","2 Mo","6 Mo", "TC")) 
  facet_grid(.~ time) +theme_void()+
  theme(strip.text = element_text(size = 14,face = "bold"),
        legend.text=element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 14,face = "bold"))
  
pdf(paste0(results_folder,"/Fig_3C.pdf"),height = 15, width = 15,useDingbats = F)
print(pie_chart)
dev.off()


