# Import genome ids to species:
library(ggplot2)
library(tidyverse)

# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_5"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")

# Read instrain results:

library(readxl)
sheet_names <- excel_sheets("../data/InStrain_mutations/InStrain_mutations.xlsx")           # Get sheet names
sheet_names   

list_all <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel("../data/InStrain_mutations/InStrain_mutations.xlsx", sheet = x)) } )
names(list_all) <- sheet_names     


# atpC mutations:

atpC_mut <-  list_all$InStrain_atpC_mutations
atpC_mut_mtb_map <- atpC_mut[!is.na(atpC_mut$dorea_ref_position),]
# Only display the ones that have matching seq WT Dorea
atpC_mut_mtb_map <- atpC_mut_mtb_map[atpC_mut_mtb_map$dorea_align == atpC_mut_mtb_map$aa_align,]


# Make a heatmap of mutation species and mutation site:
plot_ht_dt <-  atpC_mut_mtb_map %>%
  select(SNP_ID,mutation,gene,Lineage,dorea_align,dorea_ref_position)
# Extract the last residue:
plot_ht_dt$residue <-  str_sub(plot_ht_dt$mutation, start= -1)

mat <-  plot_ht_dt %>%
  ungroup() %>%
  select(gene,dorea_ref_position,residue)%>%
  group_by(gene,dorea_ref_position)%>%
  summarise(residue = paste(residue, collapse=","))%>%
  pivot_wider(id_cols = gene,names_from = dorea_ref_position,values_from = residue,values_fill = "NO") %>%
  data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL

library(gtools)
mat <- mat[,mixedsort(colnames(mat))]

char_vec <-unique(unlist(mat))

library(hues)
set.seed(200)
cols_mat <- iwanthue(length(char_vec), plot=T)
cols_mat[1] <- "white"
names(cols_mat) <- char_vec
#namescolors = structure(1:length(char_vec), names = char_vec)

# Dorea position
colnames(mat) <- gsub("X","",colnames(mat))

pos_dt <- atpC_mut_mtb_map %>% select(dorea_ref_position,dorea_align) %>% unique() %>% arrange(dorea_ref_position)
colnames(mat) <- paste0(colnames(mat),":",pos_dt$dorea_align )

# Row names:
rownames(mat) <- str_extract(rownames(mat), "[^_]*_[^_]*")

sub_genome_dt <-  atpC_mut_mtb_map %>% select(gene,Lineage) %>% unique()
sub_genome_dt$Genome <-  str_extract(sub_genome_dt$gene, "[^_]*_[^_]*")
table(sub_genome_dt$Genome ==  rownames(mat))

tax_dt <-  str_split_fixed(sub_genome_dt$Lineage,";",n = 7)
tax_dt <- data.frame( gsub(".*__","",tax_dt))
names(tax_dt) <-  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
idx_sp <-  which(tax_dt$Species == "")
tax_dt$Species[idx_sp] <-  paste0(tax_dt$Genus[idx_sp]," uncl_Sp")


rownames(mat) <- paste0(tax_dt$Species," ",gsub("GUT_","",sub_genome_dt$Genome))

colors <-   structure(1:length(char_vec), names = char_vec)
library(ComplexHeatmap)
ht <- Heatmap(mat, name = "mat", col = colors,
              column_title = paste0(unique(atpC_mut_mtb_map$gene_name)," Mutations across genomes"),
              #  cluster_rows = F,
              #  cluster_columns = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mat[i, j] != "NO"){
                  
                  grid.rect(x, y, width, height, gp = gpar(fill = "lightblue", col = "lightblue"))
                  grid.text(sprintf(mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }else{
                  grid.rect(x, y, width, height, gp = gpar(fill = "grey50", col = "grey50"))
                  grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                  
                }
              },
              row_names_side = "left",
              border = T,
              show_heatmap_legend = F,
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_max_width = max_text_width(rownames(mat),
                                                   gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                                       gp = gpar(fontsize = 10)))


# Make a heatmap 
pdf(paste0(results_folder,"/","Fig_5_A_atpC","_MDR.pdf"),width = 6,height = 2)
draw(ht,heatmap_legend_side = c("right"),annotation_legend_side="right",legend_grouping = "original")
dev.off()



# glpK mutations:

glpK_mut <-  list_all$InStrain_glpK_mutations

# Only display the ones that align with MTB gene and have same amino acid in the WT
glpK_mut_mtb_map <- glpK_mut[!is.na(glpK_mut$mtb_ref_position),]
# Only display the ones that have matching seq WT MTB
glpK_mut_mtb_map <- glpK_mut_mtb_map[glpK_mut_mtb_map$mtb_align == glpK_mut_mtb_map$aa_align,]


# Make a heatmap of mutation species and mutation site:
plot_ht_dt <-  glpK_mut_mtb_map %>%
  select(SNP_ID,mutation,gene,Lineage,mtb_align,mtb_ref_position)
# Extract the last residue:
plot_ht_dt$residue <-  str_sub(plot_ht_dt$mutation, start= -1)

mat <-  plot_ht_dt %>%
  ungroup() %>%
  select(gene,mtb_ref_position,residue)%>%
  group_by(gene,mtb_ref_position)%>%
  summarise(residue = paste(residue, collapse=","))%>%
  pivot_wider(id_cols = gene,names_from = mtb_ref_position,values_from = residue,values_fill = "NO") %>%
  data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL

library(gtools)
mat <- mat[,mixedsort(colnames(mat))]

char_vec <-unique(unlist(mat))

library(hues)
set.seed(200)
cols_mat <- iwanthue(length(char_vec), plot=T)
cols_mat[1] <- "white"
names(cols_mat) <- char_vec

# MTB position
colnames(mat) <- gsub("X","",colnames(mat))

pos_dt <- glpK_mut_mtb_map %>% select(mtb_ref_position,mtb_align) %>% unique() %>% arrange(mtb_ref_position)
colnames(mat) <- paste0(colnames(mat),":",pos_dt$mtb_align )

# Row names:
rownames(mat) <- str_extract(rownames(mat), "[^_]*_[^_]*")

sub_genome_dt <-  glpK_mut_mtb_map %>% select(gene,Lineage) %>% unique()
sub_genome_dt$Genome <-  str_extract(sub_genome_dt$gene, "[^_]*_[^_]*")
table(sub_genome_dt$Genome ==  rownames(mat))


tax_dt <-  str_split_fixed(sub_genome_dt$Lineage,";",n = 7)
tax_dt <- data.frame( gsub(".*__","",tax_dt))
names(tax_dt) <-  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
idx_sp <-  which(tax_dt$Species == "")
tax_dt$Species[idx_sp] <-  paste0(tax_dt$Genus[idx_sp]," uncl_Sp")


rownames(mat) <- paste0(tax_dt$Species," ",gsub("GUT_","",sub_genome_dt$Genome))

colors <-   structure(1:length(char_vec), names = char_vec)
library(ComplexHeatmap)
ht <- Heatmap(mat, name = "mat", col = colors,
              column_title = "glpK Mutations across genomes",
              #  cluster_rows = F,
              #  cluster_columns = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mat[i, j] != "NO"){
                  
                  grid.rect(x, y, width, height, gp = gpar(fill = "lightblue", col = "lightblue"))
                  grid.text(sprintf(mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }else{
                  grid.rect(x, y, width, height, gp = gpar(fill = "grey50", col = "grey50"))
                  grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                  
                }
              },
              row_names_side = "left",
              border = T,
              show_heatmap_legend = F,
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_max_width = max_text_width(rownames(mat),
                                                   gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                                       gp = gpar(fontsize = 10)))


# Make a heatmap 
pdf(paste0(results_folder,"/Fig_5_B_glpK_Mutations_MDR.pdf"),width = 6,height = 2)
draw(ht,heatmap_legend_side = c("right"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# atpE mutations:

atpE_mut <-  list_all$InStrain_atpE_mutations
atpE_mut_mtb_map <- atpE_mut[!is.na(atpE_mut$mtb_ref_position),]

# Make a heatmap of mutation species and mutation site:
plot_ht_dt <-  atpE_mut_mtb_map %>%
  select(SNP_ID,mutation,gene,Lineage,aa_align,mtb_align,mtb_ref_position)
# Extract the last residue:
plot_ht_dt$residue <-  str_sub(plot_ht_dt$mutation, start= -1)
plot_ht_dt$residue <-  paste0(plot_ht_dt$aa_align,"->",plot_ht_dt$residue)


mat <-  plot_ht_dt %>%
  ungroup() %>%
  select(gene,mtb_ref_position,residue)%>%
  group_by(gene,mtb_ref_position)%>%
  summarise(residue = paste(residue, collapse=","))%>%
  pivot_wider(id_cols = gene,names_from = mtb_ref_position,values_from = residue,values_fill = "NO") %>%
  data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL

library(gtools)
mat <- mat[,mixedsort(colnames(mat))]

char_vec <-unique(unlist(mat))

library(hues)
set.seed(200)
cols_mat <- iwanthue(length(char_vec), plot=T)
cols_mat[1] <- "white"
names(cols_mat) <- char_vec

# MTB position
colnames(mat) <- gsub("X","",colnames(mat))


pos_dt <-  atpE_mut_mtb_map %>% select(mtb_ref_position,mtb_align) %>% unique() %>% arrange(mtb_ref_position)
colnames(mat) <- paste0(colnames(mat),":",pos_dt$mtb_align )

# Row names:
rownames(mat) <- str_extract(rownames(mat), "[^_]*_[^_]*")

sub_genome_dt <- atpE_mut_mtb_map %>% select(gene,Lineage) %>% unique()
sub_genome_dt$Genome <- str_extract(sub_genome_dt$gene, "[^_]*_[^_]*")

table(sub_genome_dt$Genome ==  rownames(mat))



tax_dt <-  str_split_fixed(sub_genome_dt$Lineage,";",n = 7)
tax_dt <- data.frame( gsub(".*__","",tax_dt))
names(tax_dt) <-  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
idx_sp <-  which(tax_dt$Species == "")
tax_dt$Species[idx_sp] <-  paste0(tax_dt$Genus[idx_sp]," uncl_Sp")


rownames(mat) <- paste0(tax_dt$Species," ",gsub("GUT_","",sub_genome_dt$Genome))

colors <-   structure(1:length(char_vec), names = char_vec)
library(ComplexHeatmap)
ht <- Heatmap(mat, name = "mat", col = colors,
              column_title = "atpE mutations",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mat[i, j] != "NO"){
                  
                  grid.rect(x, y, width, height, gp = gpar(fill = "lightblue", col = "lightblue"))
                  grid.text(sprintf(mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }else{
                  grid.rect(x, y, width, height, gp = gpar(fill = "grey50", col = "grey50"))
                  grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                  
                }
              },
              row_names_side = "left",
              border = T,
              show_heatmap_legend = F,
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_max_width = max_text_width(rownames(mat),
                                                   gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                                       gp = gpar(fontsize = 10)))


# Make a heatmap 
pdf(paste0(results_folder,"/Fig_5_C_atpE_MDR.pdf"),width = 6,height = 2)
draw(ht,heatmap_legend_side = c("right"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# GyrA mutations:

gyrA_mut <- list_all$InStrain_gyrA_mutations
# Only display the ones that align with MTB  and have same amino acid in the WT
gyr_mut_mtb_map <- gyrA_mut[!is.na(gyrA_mut$mtb_ref_position),]
# Only display the ones that mutations in multiple species:
freq_mut <-  stack(table(gyr_mut_mtb_map$mtb_ref_position))

sel_pos_mtb <-  as.character(freq_mut$ind[freq_mut$values > 1])
gyr_mut_mtb_map <- gyr_mut_mtb_map[gyr_mut_mtb_map$mtb_ref_position %in% sel_pos_mtb,]

# Make a heatmap of mutation species and mutation site:
plot_ht_dt <-  gyr_mut_mtb_map %>%
  select(SNP_ID,mutation,gene,Lineage,aa_align,mtb_align,mtb_ref_position)
# Extract the last residue:
plot_ht_dt$residue <-  str_sub(plot_ht_dt$mutation, start= -1)
plot_ht_dt$residue <-  paste0(plot_ht_dt$aa_align,"->",plot_ht_dt$residue)

mat <-  plot_ht_dt %>%
  ungroup() %>%
  select(gene,mtb_ref_position,residue)%>%
  group_by(gene,mtb_ref_position)%>%
  summarise(residue = paste(residue, collapse=","))%>%
  pivot_wider(id_cols = gene,names_from = mtb_ref_position,values_from = residue,values_fill = "NO") %>%
  data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL

library(gtools)
mat <- mat[,mixedsort(colnames(mat))]

char_vec <-unique(unlist(mat))

library(hues)
set.seed(200)
cols_mat <- iwanthue(length(char_vec), plot=T)
cols_mat[1] <- "white"
names(cols_mat) <- char_vec

# MTB position
colnames(mat) <- gsub("X","",colnames(mat))

pos_dt <- gyr_mut_mtb_map %>% select(mtb_ref_position,mtb_align) %>% unique() %>% arrange(mtb_ref_position)
colnames(mat) <- paste0(colnames(mat),":",pos_dt$mtb_align )


sub_genome_dt <- gyr_mut_mtb_map %>% select(gene,Lineage) %>% unique()
sub_genome_dt$Genome <- str_extract(sub_genome_dt$gene, "[^_]*_[^_]*")


table(sub_genome_dt$gene ==  rownames(mat))

tax_dt <-  str_split_fixed(sub_genome_dt$Lineage,";",n = 7)
tax_dt <- data.frame( gsub(".*__","",tax_dt))
names(tax_dt) <-  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
idx_sp <-  which(tax_dt$Species == "")
tax_dt$Species[idx_sp] <-  paste0(tax_dt$Genus[idx_sp]," uncl_Sp")


rownames(mat) <- make.unique(paste0(tax_dt$Species," ",gsub("GUT_","",sub_genome_dt$Genome)))

colors <-   structure(1:length(char_vec), names = char_vec)
library(ComplexHeatmap)
ht <- Heatmap(mat, name = "mat", col = colors,
              column_title = "gyrA Mutations across genomes",
              #  cluster_rows = F,
              #  cluster_columns = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mat[i, j] != "NO"){
                  
                  grid.rect(x, y, width, height, gp = gpar(fill = "lightblue", col = "lightblue"))
                  grid.text(sprintf(mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }else{
                  grid.rect(x, y, width, height, gp = gpar(fill = "grey50", col = "grey50"))
                  grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                  
                }
              },
              row_names_side = "left",
              border = T,
              show_heatmap_legend = F,
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_max_width = max_text_width(rownames(mat),
                                                   gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                                       gp = gpar(fontsize = 10)))


# Make a heatmap 
pdf(paste0(results_folder,"/Fig_5_E_gyrA_Mutations_MDR.pdf"),width = 6,height = 4)
draw(ht,heatmap_legend_side = c("right"),annotation_legend_side="right",legend_grouping = "original")
dev.off()

# rpoB mutations:


rif_mut <-  list_all$InStrain_rpoB_mutations
# Only display the ones that align with MTB rpoB gene and have same amino acid in the WT
rif_mut_mtb_map <- rif_mut[!is.na(rif_mut$mtb_ref_position),]
# Only display the ones that have matching seq WT MTB
rif_mut_mtb_map <- rif_mut_mtb_map[rif_mut_mtb_map$mtb_align == rif_mut_mtb_map$aa_align,]

# Make a heatmap of mutation species and mutation site:
plot_ht_dt <-  rif_mut_mtb_map %>%
  select(SNP_ID,mutation,gene,Lineage,mtb_align,mtb_ref_position)
# Extract the last residue:
plot_ht_dt$residue <-  str_sub(plot_ht_dt$mutation, start= -1)

mat <-  plot_ht_dt %>%
  ungroup() %>%
  select(gene,mtb_ref_position,residue)%>%
  group_by(gene,mtb_ref_position)%>%
  summarise(residue = paste(residue, collapse=","))%>%
  pivot_wider(id_cols = gene,names_from = mtb_ref_position,values_from = residue,values_fill = "NO") %>%
  data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL

mat <- mat[,sort(colnames(mat))]

char_vec <-unique(unlist(mat))

library(hues)
set.seed(200)
cols_mat <- iwanthue(length(char_vec), plot=T)
cols_mat[1] <- "white"
names(cols_mat) <- char_vec


colors = structure(1:length(char_vec), names = char_vec)

# MTB position
colnames(mat) <- gsub("X","",colnames(mat))

pos_dt <- rif_mut_mtb_map %>% select(mtb_ref_position,mtb_align) %>% unique() %>% arrange(mtb_ref_position)
colnames(mat) <- paste0(colnames(mat),":",pos_dt$mtb_align )

sub_genome_dt <- rif_mut_mtb_map %>% select(gene,Lineage) %>% unique()
sub_genome_dt$Genome <- str_extract(sub_genome_dt$gene, "[^_]*_[^_]*")

# Row names:
rownames(mat) <- str_extract(rownames(mat), "[^_]*_[^_]*")

table(sub_genome_dt$Genome ==  rownames(mat))
tax_dt <-  str_split_fixed(sub_genome_dt$Lineage,";",n = 7)
tax_dt <- data.frame( gsub(".*__","",tax_dt))
names(tax_dt) <-  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
idx_sp <-  which(tax_dt$Species == "")
tax_dt$Species[idx_sp] <-  paste0(tax_dt$Genus[idx_sp]," uncl_Sp")


rownames(mat) <- paste0(tax_dt$Species," ",gsub("GUT_","",sub_genome_dt$Genome))

library(ComplexHeatmap)
ht <- Heatmap(mat, name = "mat", col = colors,
              column_title = "rpoB Mutations across genomes",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(mat[i, j] != "NO"){
                  
                  grid.rect(x, y, width, height, gp = gpar(fill = "lightblue", col = "lightblue"))
                  grid.text(sprintf(mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }else{
                  grid.rect(x, y, width, height, gp = gpar(fill = "grey50", col = "grey50"))
                  grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                  
                }
              },
              row_names_side = "left",
              border = T,
              show_heatmap_legend = F,
              column_names_gp = gpar(fontsize = 10,fontface = "bold"),
              row_names_gp = gpar(fontsize = 8,fontface = "bold"),
              row_names_max_width = max_text_width(rownames(mat),
                                                   gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                                       gp = gpar(fontsize = 10)))


# Make a heatmap 
pdf(paste0(results_folder,"/Fig_5_G_rpoB_Mutations.pdf"),width = 7,height = 10)
draw(ht,heatmap_legend_side = c("right"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


