# List files 
library(stringr)
library(tidyverse)
library(phyloseq)
library(ape)

# Create a directories to store results from the analysis
mainDir <- "../Figures"
subDir <- "Fig_S5"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,subDir,sep="/")


phy_mdr <-  readRDS("../Processed_data/phy_mtph_mdr.rds")  
sm_dt <-  data.frame(sample_data(phy_mdr))  

# List of files:
sp_list <- c("102506","102538")


sp <-  sp_list[1]
for(sp in sp_list){
  
  
  flist <- list.files(path = paste0("../data/Midas/MDR_EColi_Kleb/",sp),recursive = T,full.names = T)
  myfiles = lapply(flist, read.csv)
  visit_time <- gsub("VIMP_X","",basename(flist))
  visit_time <- gsub("_.*","",visit_time)
  names(myfiles) <- visit_time
  visit_t <- visit_time[1]
  imp_dt_list <- list()
  for(visit_t in visit_time){
    
    imp_dt <- myfiles[[visit_t]]   
    # Print variable importance
    # Plot importance:
    imp_dt <- imp_dt[order(imp_dt$MeanDecreaseAccuracy,decreasing = F),]
    imp_dt$Function <-  as.character(imp_dt$Function)
    imp_dt$Function <-  factor(imp_dt$Function , levels =  imp_dt$Function)
    imp_dt$Visit <-  visit_t
    imp_dt_list[[visit_t]] <-  imp_dt[4:ncol(imp_dt)]
  }
  
  imp_comb_dt <-  do.call("rbind",imp_dt_list)
  
  snps_ann_dt <-  read.csv(paste0("../data/Midas/MDR_EColi_Kleb/SNPs_",sp,".csv"))
  
  meta_dt <- data.frame(sample_data(phy_mdr))
  meta_dt <-  meta_dt[meta_dt$Sample %in% unique(snps_ann_dt$Sample),]
  meta_dt <-  meta_dt[,c("Sample","PID","Visit")]
  snps_ann_dt <- snps_ann_dt %>%
    left_join(meta_dt)
  
  
  
  
  # Only select upto 6 months
  table(snps_ann_dt$Visit)
  snps_ann_dt  <-  snps_ann_dt[snps_ann_dt$Visit %in% c("0 Day"  ,"2 Wk"  ,"1 Mo"  ,"2 Mo"  ,"6 Mo"),]
  
  # Now summarize it across time points:
  counts <- snps_ann_dt
  
  counts_sum <- counts %>%
    group_by(Visit,site_id)%>%
    summarise(A = sum(count_a), C = sum(count_c),G = sum(count_g), T = sum(count_t))%>%
    as.data.frame()
  
  
  # sids 
  sids <-  unique(counts_sum$site_id)
  
  time_dt_w <- imp_comb_dt %>%
    select(site_id,Visit)%>%
    unique()%>%
    mutate(val = 1)%>%
    pivot_wider(names_from = Visit,values_from = val)%>%
    data.frame()
  
  
  pval_dt_long <-  time_dt_w %>%
    pivot_longer(cols = starts_with("X"),names_to = "Visit",values_to = "pval")
  
  pval_dt_long$Visit <- factor(pval_dt_long$Visit)
  levels(pval_dt_long$Visit) <-  c("1 Mo", "2 Mo","2 Wk", "6 Mo")
  #levels(pval_dt_long$Visit) <-  c("2 Mo", "6 Mo")
  
  # p_mat <- p_mat %>%
  #   filter( p_mat$X1.Mo < 0.05 | p_mat$X2.Mo < 0.05 | p_mat$X6.Mo < 0.05 )
  counts_fil_sum <- counts_sum[counts_sum$site_id %in% time_dt_w$site_id,] %>%
    left_join(unique(counts[,c("site_id","ref_allele")]))%>%
    left_join(unique(imp_comb_dt[,c("site_id","Function")]))
  
  fil_dt <- counts_fil_sum %>%
    left_join(pval_dt_long)
  fil_dt$pval[fil_dt$Visit == "0 Day"] <- 1
  #fil_dt$pval[is.na(fil_dt$pval)] <- 10
  
  head(fil_dt)
  
  fil_dt_sel <-  fil_dt
  
  fil_dt_sel$Visit <- factor(as.character(fil_dt_sel$Visit),levels =c("0 Day"  ,"2 Wk"  ,"1 Mo"  ,"2 Mo"  ,"6 Mo"))
  
  fil_dt_sel <- fil_dt_sel %>% 
    mutate(Visit_num = as.numeric(Visit),
           Site_num = as.numeric(factor(site_id)))
  
  site_names_dt <-  fil_dt_sel %>% select(site_id,Site_num,Function)%>%unique() %>%
    arrange(Site_num)
  
  
  library(scatterpie)
  library(ggprism)
  
  fil_dt_sel <- fil_dt_sel %>% filter(!is.na(pval))
  
  p_profile <-    ggplot()+
    geom_scatterpie(aes(x=Visit_num, y=Site_num, r=0.4),
                    cols= c("A","C","G","T"), 
                    color= "black", data= fil_dt_sel)+
    scale_x_continuous(breaks=c(1,2,3,4,5), labels=c("0 Day","2 Wk","1 Mo","2 Mo","6 Mo")) + 
    
    scale_y_continuous(breaks=1:nrow(site_names_dt), labels = site_names_dt$Function) + 
    labs(x="Visit", y="Site IDs") + 
    coord_fixed()+
    theme_prism()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  pdf(paste0(results_folder,"/SNPs_profile_",sp,".pdf"),width = 15,height = 20)
  print(p_profile)     
  dev.off() 
  
}



