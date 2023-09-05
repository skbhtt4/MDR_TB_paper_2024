library("phyloseq")
library("gtools")
library("stringr")

# Create a directories to store results from the analysis
mainDir <- "../Processed_data"
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

###############################*********************************Import Metaphlan*****************************######################
# Read metaphlan and skip the first row
mtph_dat <- read.csv("../data/humann/20220825_merged_mph.tsv",sep="\t",skip = 1)
# Only select the taxa at the species level
mtph_dat <-  mtph_dat[grep("s__",mtph_dat$clade_name),]

# Split into taxonomy table
tax_dt <-  str_split_fixed(mtph_dat$clade_name,pattern = "\\|", n = 7) %>% data.frame()
colnames(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove the *__  from all taxa levels
tax_dt[] <- lapply(tax_dt, function(x) gsub("^.{0,3}", "", x)) 
rownames(tax_dt) <-  tax_dt$Species

# Relative abundance table
rel_dt <- mtph_dat[,3:ncol(mtph_dat)] 
rownames(rel_dt) <-  tax_dt$Species
colnames(rel_dt) <- gsub("_mph","",colnames(rel_dt))

# Import metadata for MDR
meta_mdr <-  readxl::read_xlsx("../data/metadata/corrected_MDR_metadata_11feb2022.xlsx") %>% data.frame()
head(meta_mdr)
# Associate sample names with Processing number
# Attach metadata 
fname <- "../data/metadata/mdr_stool_delivered_02dec2021 .xlsx"
sheets <- readxl::excel_sheets(fname)
meta_1 <- readxl::read_excel(fname, sheet = sheets[1])
meta_2 <- readxl::read_excel(fname, sheet = sheets[2])
meta_3 <-  readxl::read_excel(fname, sheet =  sheets[3])

library(tidyverse)
meta_dt <- meta_2 %>%
  inner_join(meta_3)%>% data.frame()
rownames(meta_dt) <- gsub(" ","_",meta_dt$submission.sample.number)
meta_dt <-  meta_dt[,c("Global.Specimen.ID","submission.sample.number")]

head(meta_dt)
meta_mdr_dt <-  meta_dt %>%
  inner_join(meta_mdr) %>% data.frame()

# Correct sample labels
meta_mdr_dt$Global.Specimen.ID[meta_mdr_dt$Global.Specimen.ID == "0487-0093XC00-002"] <-"0487-0093XC00-001"
rownames(meta_mdr_dt) <- make.names(meta_mdr_dt$Global.Specimen.ID)
all(rownames(meta_mdr_dt) %in% colnames(rel_dt))

meta_mdr_dt$submission.sample.number <- NULL
meta_mdr_dt$sample.number.for.processing <- NULL
meta_mdr_dt$Visit.. <- NULL

head(meta_mdr_dt)


names(meta_mdr_dt) <-  gsub("Calculated.","",names(meta_mdr_dt))
names(meta_mdr_dt)[1] <- "Sample"

# 1930019T didn't participate in the study:
# Remove it
meta_mdr_dt <-  meta_mdr_dt[meta_mdr_dt$PID != "1930019T",]

meta_mdr_dt$Visit <- factor(meta_mdr_dt$Visit,c("0 Day","2 Wk" ,"1 Mo" , "2 Mo" , "6 Mo","24 Mo"))
# Repeat age and Sex per patient ID.
meta_mdr_dt <- meta_mdr_dt[order(meta_mdr_dt$PID, meta_mdr_dt$Visit),]

library(tidyverse)
meta_mdr_dt <-  meta_mdr_dt %>%
  group_by(PID)%>%
  fill(Sex, Age,Height,Prior.TB.Treatment.Check,Prior.Treatment.Date,Prior.Treatment.Regimen,Prior.Treatment.Outcome,HIV.Status)%>%
  data.frame()

rownames(meta_mdr_dt) <- make.names(meta_mdr_dt$Sample)

names(rel_dt)
names(rel_dt)[!names(rel_dt) %in% rownames(meta_mdr_dt)]


# Fix the missing TTP and fill all the 24 months and six months samples with max TTP
# Make an indicator column to identify the imputed samples:

names(meta_mdr_dt)[11] <- "Av_TTP"

# If Any TTP measurement is negative , replace that for the average
table(meta_mdr_dt$TTP1)
meta_mdr_dt$Av_TTP[meta_mdr_dt$TTP1 %in% c("Negative")] <- "Negative"
meta_mdr_dt$Av_TTP[meta_mdr_dt$TTP2 %in% c("Negative")] <- "Negative"

# Set 6 Mo and 24 Mo as Max TTP
meta_mdr_dt$Av_TTP[meta_mdr_dt$Av_TTP %in% c("Negative")] <- "1000"
meta_mdr_dt$Av_TTP <-  as.numeric(meta_mdr_dt$Av_TTP)
#meta_mdr_dt$Av_TTP[meta_mdr_dt$Visit %in% c("2 Mo") & meta_mdr_dt$PID %in% c("1930004T")] <- 1000
# Indicator for substitution:
meta_mdr_dt$ind <- 0
meta_mdr_dt$ind[is.na(meta_mdr_dt$Av_TTP) & meta_mdr_dt$Visit %in% c("6 Mo","24 Mo")] <-  1

meta_mdr_dt$Av_TTP[meta_mdr_dt$Visit %in% c("6 Mo","24 Mo")] <- 1000

# Remove samples that don't have TTP value
meta_mdr_dt$Time <- ifelse(meta_mdr_dt$Visit == "0 Day",0,
                     ifelse(meta_mdr_dt$Visit == "2 Wk",14,
                            ifelse(meta_mdr_dt$Visit == "1 Mo",30,
                            ifelse(meta_mdr_dt$Visit == "2 Mo",60,
                                   ifelse(meta_mdr_dt$Visit == "6 Mo",180,720)))))

phy_mtph_mdr <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                          tax_table(as.matrix(tax_dt)),
                          sample_data(meta_mdr_dt))

sample_sums(phy_mtph_mdr)
phy_mtph_mdr <-  prune_taxa(taxa_sums(phy_mtph_mdr)>0, phy_mtph_mdr)
phy_mtph_mdr <- transform_sample_counts(phy_mtph_mdr, function(x) x/sum(x))

saveRDS(phy_mtph_mdr,paste0(results_folder,"/phy_mtph_mdr.rds"))

# Add counts for the alpha diversity measures 
# Skip the first row
mtph_counts <- read.csv("../data/humann/mtph_counts.txt",sep="\t",skip = 1)
# Only select the taxa at the species level
mtph_counts <-  mtph_counts[grep("s__",mtph_counts$clade_name),]
# Relative abundance table
counts_dt <- mtph_counts[,3:ncol(mtph_counts)] 
rownames(counts_dt) <-  gsub(".*s__","",mtph_counts$clade_name)
colnames(counts_dt) <- gsub("_mph","",colnames(counts_dt))

phy_mtph_counts <- phyloseq(otu_table(counts_dt,taxa_are_rows = T),
                            tax_table(as.matrix(tax_dt)),
                            sample_data(meta_mdr_dt))
phy_mtph_counts <-  prune_taxa(taxa_sums(phy_mtph_counts)>0,phy_mtph_counts)
saveRDS(phy_mtph_counts,paste0(results_folder,"/phy_mtph_counts.rds"))
  

##################################
# Import pathways:
path_dat <-  read.csv("../data/humann/humann3_pathabundance.txt",sep="\t")
names(path_dat)[1] <- "Pathway"

# Select only at the community level
path_dat <- path_dat[!grepl("\\|",path_dat$Pathway),]

# Now remove the Unmapped and Unintegrated:
path_dat <- path_dat[!grepl("UNMAPPED|UNINTEGRATED",path_dat$Pathway),]
rownames(path_dat) <-  path_dat$Pathway
path_dat$Pathway <-  NULL

#colnames(path_dat) <- gsub("d10s","_",colnames(path_dat))
colnames(path_dat) <- gsub("_Abundance","",colnames(path_dat))

phy_path <-  phyloseq(otu_table(path_dat,taxa_are_rows = T),
                      sample_data(meta_mdr_dt))
phy_path <-  prune_taxa(taxa_sums(phy_path)>0, phy_path)


saveRDS(phy_path,paste0(results_folder,"/phy_raw_path_mdr.rds"))

phy_path <- transform_sample_counts(phy_path, function(x) x/sum(x))
sample_sums(phy_path)
saveRDS(phy_path,paste0(results_folder,"/phy_rel_path_mdr.rds"))

# phy_path <- prune_taxa(taxa_sums(phy_path)>0,phy_path)


############## Import CARD data ###################
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


library(phyloseq)
# Import phyloseq for metaphlan:
phy_card <- phy_mtph_mdr
sm_dt <-  data.frame(sample_data(phy_card))
tax_table(phy_card) <- NULL
otu_table(phy_card) <-  otu_table(as.matrix(CARD_genes),taxa_are_rows = T)
tax_table(phy_card) <- tax_table(as.matrix(meta_genes))

saveRDS(phy_card,paste0(results_folder,"/phy_card.rds"))




############################# Import RNASeq data #################################
# Now import RNASeq data for MDR:
rna_counts_1 <-  read.csv("../data/RNASeq/raw.counts.csv")
head(rna_counts_1)

snames <-  names(rna_counts_1)

rna_counts_2 <-  read.csv("../data/RNASeq/raw_counts_30MDR_06jun2022.csv")

all(rna_counts_1$Gene == rna_counts_2$Gene)
dim(rna_counts_1)
dim(rna_counts_2)


# Match with the one 
rna_counts_2 <- rna_counts_2[match(rna_counts_1$Gene,rna_counts_2$Gene),]
rna_counts <- cbind(rna_counts_1,rna_counts_2[,2:ncol(rna_counts_2)])
snames <-  names(rna_counts)

library(stringr)
dt_split <- str_split_fixed(snames,pattern = "\\.",n = 2) %>% data.frame()
# Remove the first row
dt_split <- dt_split[-1,]
names(dt_split) <- c("PID","Visit")
dt_split$PID <-  gsub("X","",dt_split$PID)
dt_split$Visit <- factor(dt_split$Visit)
levels(dt_split$Visit) <-  c("0 Day", "2 Mo","24 Mo","6 Mo","2 Wk")
dt_split$Visit <-  as.character(dt_split$Visit)

sm_dt <-  data.frame(sample_data(phy_mtph_mdr))
library(dplyr)
rna_meta <- dt_split %>%
  left_join(sm_dt,by = c("PID","Visit"))
rownames(rna_meta) <-  rna_meta$Sample
rownames(rna_counts) <-  rna_counts$Gene
rna_counts$Gene <-  NULL
rownames(rna_counts)
colnames(rna_counts) <- rna_meta$Sample

rna_meta$Visit <- factor(rna_meta$Visit, levels = c("0 Day" ,"2 Wk", "2 Mo"  ,"6 Mo","24 Mo"))
# Change the level of last time point as TC
levels(rna_meta$Visit)[5] <- "TC"

phy_rna <-  phyloseq(otu_table(rna_counts,taxa_are_rows = T),
                     sample_data(rna_meta))
phy_rna <- prune_taxa(taxa_sums(phy_rna)>0,phy_rna)
taxa_names(phy_rna)

saveRDS(phy_rna,paste0(results_folder,"/phy_rna_mdr.rds"))

# Combine metadata:
meta_mdr <-  data.frame(sample_data(phy_rna))
meta_mdr$PID <- factor(meta_mdr$PID )
#meta_mdr$Visit <- factor(meta_mdr$Visit , levels = c("0 Day",   "2 Wk", "2 Mo",  "6 Mo", "24 Mo"  ))
meta_mdr$Sex <-  factor(meta_mdr$Sex,levels = c("Male","Female") )
#meta_mdr <-  meta_mdr[,c("PID","Visit","Sex","Age")]
# meta_mdr$Treatment <-  "MDR"
sample_data(phy_rna) <-  sample_data(meta_mdr)

library(DESeq2)
phy_gene_sel <-  phy_rna
# Prune out low abundant transcripts
phy_gene_sel <- prune_taxa(taxa_sums(phy_gene_sel)>10,phy_gene_sel)

dds <- phyloseq_to_deseq2(phy_gene_sel , ~ PID + Visit) #replace this with any sample variable(s)
#dds <- estimateSizeFactors(dds,"poscounts")
dds <- DESeq(dds,parallel = T)

# VST phyloseq for Heatmap later
phy_vst_rna <- phy_gene_sel
vst_dt <- getVarianceStabilizedData(dds)
otu_table(phy_vst_rna) <-  otu_table(vst_dt, taxa_are_rows = T)
saveRDS(phy_vst_rna,paste0(results_folder,"/phy_rna_vst_mdr.rds"))




################# Import HRZE data #################################

# Import metaphlan
# Read metaphlan and skip the first row
mtph_dat <- read.csv("../data/humann/20220825_merged_mph.tsv",sep="\t",skip = 1)
# Only select the taxa at the species level
mtph_dat <-  mtph_dat[grep("s__",mtph_dat$clade_name),]

# Split into taxonomy table
tax_dt <-  str_split_fixed(mtph_dat$clade_name,pattern = "\\|", n = 7) %>% data.frame()
colnames(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove the *__  from all taxa levels
tax_dt[] <- lapply(tax_dt, function(x) gsub("^.{0,3}", "", x)) 
rownames(tax_dt) <-  tax_dt$Species

# Relative abundance table
rel_dt <- mtph_dat[,3:ncol(mtph_dat)] 
rownames(rel_dt) <-  tax_dt$Species
colnames(rel_dt) <- gsub("_mph","",colnames(rel_dt))

meta_hrze <-  readxl::read_xlsx("../data/metadata/HRZE_cohort_metadata_12aug2022.xlsx") %>% data.frame()
head(meta_hrze)

meta_hrze <-  meta_hrze[,c("Global.Spec.ID","PID","Visit","Age","Sex","Weight","Height","BMI","AFB.Smear",
                           "AFB.Date","GeneXpert","GeneXpert.Date","HIV","TTP1","TTP2","Average.TTP")]
# Change 14 Day to 2 Wk
meta_hrze$Visit[meta_hrze$Visit == "14 Day"] <- "2 Wk"
table(meta_hrze$Visit)
meta_hrze$Visit <- factor(meta_hrze$Visit,c("0 Day","1 Wk" ,"2 Wk" ,"4 Wk", "8 Wk" ,"24 Wk"))
# Repeat age and Sex per patient ID.
meta_hrze <- meta_hrze[order(meta_hrze$PID, meta_hrze$Visit),]
head(meta_hrze)
library(tidyverse)
meta_hrze <-  meta_hrze %>%
  group_by(PID)%>%
  fill(Sex, Age,Height,HIV)%>%
  data.frame()

names(meta_hrze)[1] <- "Sample"
rownames(meta_hrze) <- make.names(meta_hrze$Sample)

names(rel_dt)
names(rel_dt)[!names(rel_dt) %in% rownames(meta_hrze)]

phy_mtph_hrze <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                           tax_table(as.matrix(tax_dt)),
                           sample_data(meta_hrze))

sample_sums(phy_mtph_hrze)
phy_mtph_hrze <- transform_sample_counts(phy_mtph_hrze, function(x) x/sum(x))
phy_mtph_hrze <-  prune_taxa(taxa_sums(phy_mtph_hrze)>0, phy_mtph_hrze)

saveRDS(phy_mtph_hrze,paste0(results_folder,"/phy_mtph_hrze.rds"))


# Add counts for the alpha diversity measures 
# Skip the first row
hrze_counts <- read.csv("../data/humann/mtph_counts.txt",sep="\t",skip = 1)
# Only select the taxa at the species level
hrze_counts <-  hrze_counts[grep("s__",hrze_counts$clade_name),]
# Relative abundance table
counts_dt <- hrze_counts[,3:ncol(hrze_counts)] 
rownames(counts_dt) <-  gsub(".*s__","",hrze_counts$clade_name)
colnames(counts_dt) <- gsub("_mph","",colnames(counts_dt))

phy_hrze_counts <- phyloseq(otu_table(counts_dt,taxa_are_rows = T),
                            tax_table(as.matrix(tax_dt)),
                            sample_data(meta_hrze))
phy_hrze_counts <-  prune_taxa(taxa_sums(phy_hrze_counts)>0,phy_hrze_counts)
saveRDS(phy_hrze_counts,paste0(results_folder,"/phy_hrze_counts.rds"))




################Import CARD for HRZE #######################

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


library(phyloseq)
# Import phyloseq for metaphlan:
phy_card <- phy_mtph_hrze
sm_dt <-  data.frame(sample_data(phy_card))
tax_table(phy_card) <- NULL
otu_table(phy_card) <-  otu_table(as.matrix(CARD_genes),taxa_are_rows = T)
tax_table(phy_card) <- tax_table(as.matrix(meta_genes))

saveRDS(phy_card,paste0(results_folder,"/phy_card_hrze.rds"))

#######################Import Pathways for HRZE#################################
# Import pathways:
path_dat <-  read.csv("../data/humann/humann3_pathabundance.txt",sep="\t")
names(path_dat)[1] <- "Pathway"

# Select only at the community level
path_dat <- path_dat[!grepl("\\|",path_dat$Pathway),]

# Now remove the Unmapped and Unintegrated:
path_dat <- path_dat[!grepl("UNMAPPED|UNINTEGRATED",path_dat$Pathway),]
rownames(path_dat) <-  path_dat$Pathway
path_dat$Pathway <-  NULL

#colnames(path_dat) <- gsub("d10s","_",colnames(path_dat))
colnames(path_dat) <- gsub("_Abundance","",colnames(path_dat))

phy_path <-  phyloseq(otu_table(path_dat,taxa_are_rows = T),
                      sample_data(meta_hrze))
phy_path <-  prune_taxa(taxa_sums(phy_path)>0, phy_path)


saveRDS(phy_path,paste0(results_folder,"/phy_raw_path_hrze.rds"))

phy_path <- transform_sample_counts(phy_path, function(x) x/sum(x))
sample_sums(phy_path)
saveRDS(phy_path,paste0(results_folder,"/phy_rel_path_hrze.rds"))