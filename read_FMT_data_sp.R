library("phyloseq")
library("gtools")
library("stringr")

# Create a directories to store results from the analysis
mainDir <- "../Processed_data"
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

#*********************************Import Metaphlan*****************************
# Read metaphlan and skip the first row
mtph_dat <- read.csv("../data/metaphlan_fmt/merged_mph_rel_ab.tsv",sep="\t",skip = 1)
# Only select the taxa at the species level
mtph_dat <-  mtph_dat[grep("s__",mtph_dat$clade_name),]

# Split into taxonomy table
tax_dt <-  str_split_fixed(mtph_dat$clade_name,pattern = "\\|", n = 7) %>% data.frame()
colnames(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove the *__  from all taxa levels
tax_dt[] <- lapply(tax_dt, function(x) gsub("^.{0,3}", "", x)) 
rownames(tax_dt) <-  tax_dt$Species

# Relative abundance table
rel_dt <- mtph_dat[,2:ncol(mtph_dat)] 

#colSums(rel_dt)
rownames(rel_dt) <-  tax_dt$Species
colnames(rel_dt) <- gsub("_mph","",colnames(rel_dt))


meta_dt <-  read.csv("../data/metadata_fmt/expm24_metadata_switch.csv") %>% data.frame()
head(meta_dt)
rownames(meta_dt) <-  meta_dt$SeqID
meta_dt$DaysSinceAbx <-  as.character(meta_dt$DaysSinceAbx)

day_dt <- data.frame(key = c("FMTGD1","PFGD1","BDQD1","BDQD3","BDQD5","BDQD7","PBD2"),
                     DaysSinceAbx = as.character(c(0,3,13,15,17,19,20)))
library(tidyverse)                              
meta_dt <- meta_dt %>%
  left_join(day_dt)
meta_dt$key[meta_dt$Trxt == "FMT"] <- "Input" 
rownames(meta_dt) <- meta_dt$SeqID
meta_dt[meta_dt == ""] <- NA

# Remove samples:
# E24MBD13
# E24MBD15
# Remove the blank ones:
meta_dt <- meta_dt[!is.na(meta_dt$FMT),]
phy_mtph_fmt <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                          tax_table(as.matrix(tax_dt)),
                          sample_data(meta_dt))

sample_sums(phy_mtph_fmt) 
table(sample_sums(phy_mtph_fmt) > 0)

# Remove E24MBD15 , sample has no reads
phy_mtph_fmt <-  prune_samples(sample_sums(phy_mtph_fmt)>0, phy_mtph_fmt)
phy_mtph_fmt <-  prune_taxa(taxa_sums(phy_mtph_fmt)>0, phy_mtph_fmt)

sample_names(phy_mtph_fmt) <-  paste0("E24M",sample_data(phy_mtph_fmt)$MouseID,"D",sample_data(phy_mtph_fmt)$DaysSinceAbx)

sm_dt_final <-  data.frame(sample_data(phy_mtph_fmt))
sm_dt_final$SeqID <-  rownames(sm_dt_final)
sample_data(phy_mtph_fmt) <- sample_data(sm_dt_final)


phy_mtph_fmt <- transform_sample_counts(phy_mtph_fmt, function(x) x/sum(x))
phy_mtph_fmt <-  prune_taxa(taxa_sums(phy_mtph_fmt)>0, phy_mtph_fmt)

saveRDS(phy_mtph_fmt,paste0(results_folder,"/phy_mtph_fmt_rel.rds"))


#*********************************Import Metaphlan Counts*****************************
# Read metaphlan and skip the first row
mtph_dat <- read.csv("../data/metaphlan_fmt/merge_mph_counts.tsv",sep="\t",skip = 1)
# Only select the taxa at the species level
mtph_dat <-  mtph_dat[grep("s__",mtph_dat$clade_name),]

# Split into taxonomy table
tax_dt <-  str_split_fixed(mtph_dat$clade_name,pattern = "\\|", n = 7) %>% data.frame()
colnames(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove the *__  from all taxa levels
tax_dt[] <- lapply(tax_dt, function(x) gsub("^.{0,3}", "", x)) 
rownames(tax_dt) <-  tax_dt$Species

#mtph_sp <-  read.csv("../data/metaphlan_FMT/merged_abundance_table_species.txt",sep="\t")

# Relative abundance table
rel_dt <- mtph_dat[,3:ncol(mtph_dat)] 

#colSums(rel_dt)
rownames(rel_dt) <-  tax_dt$Species
colnames(rel_dt) <- gsub("_mph","",colnames(rel_dt))

phy_mtph_fmt <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                          tax_table(as.matrix(tax_dt)),
                          sample_data(meta_dt))

sample_sums(phy_mtph_fmt) 
table(sample_sums(phy_mtph_fmt) > 0)

# Remove E24MBD15 , sample has no reads
phy_mtph_fmt <-  prune_samples(sample_sums(phy_mtph_fmt)>0, phy_mtph_fmt)
phy_mtph_fmt <-  prune_taxa(taxa_sums(phy_mtph_fmt)>0, phy_mtph_fmt)

sample_names(phy_mtph_fmt) <-  paste0("E24M",sample_data(phy_mtph_fmt)$MouseID,"D",sample_data(phy_mtph_fmt)$DaysSinceAbx)
sm_dt_final <-  data.frame(sample_data(phy_mtph_fmt))
sm_dt_final$SeqID <-  rownames(sm_dt_final)
sample_data(phy_mtph_fmt) <- sample_data(sm_dt_final)

saveRDS(phy_mtph_fmt,paste0(results_folder,"/phy_mtph_fmt_counts.rds"))


