rm(list=ls())
setwd("C:\\Users\\Desktop\\179Lipids_Rcode")

library(yulab.utils)
library(ieugwasr)
library(gwasglue)
library(gwasvcf)
library(ieugwasr)
library(MRInstruments)
library(TwoSampleMR)
library(VariantAnnotation)
library("pacman")
p_load(VariantAnnotation)
p_load(data.table,ggplot2,purrr,MendelianRandomization,dplyr) 
library("TwoSampleMR","MRInstruments","gwasvcf","gwasglue")
library("pacman")
p_load(data.table,TwoSampleMR,purrr) 
library("ieugwasr")
library(MRInstruments)
library(plyr)
library(dplyr)

folder <- "Lipids_1e-5/Clumpexposure3"
dir.create(folder)
if(!dir.exists(folder)){
  dir.create(folder)
}
workdir <- getwd() 
workdir <- paste0(getwd(), "/Lipids_1e-5")  
workdir
files <- list.files(path = workdir,pattern=".csv")  
print(files)  
files
result <- data.frame()  

setwd("C:\\Users\\Desktop\\179Lipids_Rcode\\Lipids_1e-5")
for (i in files) {
  tryCatch({
    #data=fread(files[1])
    head(data)
    expo_data <- read_exposure_data(
      filename = i,  
      sep = ",", 
      phenotype_col = "Phenotype",
      snp_col = "SNP",  
      beta_col = "beta",  
      se_col = "se",  
      effect_allele_col = "effect_allele",  
      other_allele_col = "other_allele",  
      pval_col = "pval",  
      eaf_col = "eaf", 
      chr_col = "chr",  
      #pos_col = "pos",  
      samplesize_col = "samplesize" 
    )
    expo_data <- expo_data %>%
      distinct(SNP, .keep_all=TRUE)  
    expo_data <- expo_data[expo_data$SNP != "", ]  
    expo_data <- expo_data[expo_data$pval.exposure < 1e-5, ] 
    biof_iv <- expo_data[,c("SNP","pval.exposure")] 
    colnames(biof_iv) <- c("rsid","pval")
    clump_dat <- ld_clump_local(dat = biof_iv,  
                                clump_kb = 10000,  
                                clump_r2 = 0.001,  
                                clump_p = 1,  
                                bfile = "C:\\Users\\Desktop\\179Lipids_Rcode\\data_maf0.01_rs_ref\\data_maf0.01_rs_ref",  
                                plink_bin = "C:\\Users\\Desktop\\179Lipids_Rcode\\plink_win64_20231018\\plink"  
    )
    expo_data <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] 
    write.table(expo_data, file = paste0("C:/Users/Desktop/179Lipids_Rcode/Lipids_1e-5/Clumpexposure3/", basename(i)), row.names = FALSE, sep = ",", quote = FALSE)
    result <- rbind(result, expo_data) 
  }, error = function(e) {   
    cat("Error occurred for file:", i, "\n") 
    cat("Error message:", conditionMessage(e), "\n")
  })
}

library("pacman")
p_load(data.table,ggplot2,TwoSampleMR,purrr) 
library(MRInstruments)
library(dplyr)
folder <- "Clumpexposure3/FClumpexposure4"
if(!dir.exists(folder)){
  dir.create(folder)
}
orig_dir <- getwd()
workdir <- getwd() 
workdir <- paste0(getwd(),"/Clumpexposure3/")
workdir
files <- list.files(path = workdir,pattern=".csv", full.names = TRUE)
print(files)
files 
seq_along(files)
for (i in seq_along(files)) {
  file = files[i]
  expo_data <- fread(file)
  numerator <- 2 * expo_data$beta.exposure^2 * expo_data$eaf.exposure * (1 - expo_data$eaf.exposure)
  denominator <- numerator + 2 * expo_data$se.exposure^2 * expo_data$samplesize.exposure * expo_data$eaf.exposure * (1 - expo_data$eaf.exposure)
  expo_data$R2 <- numerator / denominator
  expo_data$F <- expo_data$R2 * (expo_data$samplesize.exposure - 2) / (1 - expo_data$R2)
  output=expo_data[as.numeric(expo_data$F)>10,]
  fwrite(output, paste0("./Clumpexposure3/FClumpexposure4/", basename(file)))
}

library("pacman")
p_load(data.table,MendelianRandomization,purrr,readr,TwoSampleMR)
library(MendelianRandomization)   
orig_dir <- getwd()
workdir <- getwd() 
workdir <- paste0(getwd(),"/Clumpexposure3/FClumpexposure4/")
workdir
files <- list.files(path = workdir,pattern=".csv", full.names = TRUE)
print(files)
print(files)
files
if(!dir.exists("PHENOdata5")){
  dir.create("PHENOdata5")
}

for (i in files) {
  expo_data <- read_exposure_data(
    filename = i, 
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure", 
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = "eaf.exposure",
    chr_col = "chr.exposure",
    pos_col = "pos.exposure",
    samplesize_col = "samplesize.exposure")
  grp_size <- 100 
  grps <- split(expo_data$SNP, ceiling(seq_along(expo_data$SNP)/grp_size))
  results <- list() 
  for(grp in grps) {
    res <- phenoscanner(snpquery = grp,  
                        catalogue = "GWAS",
                        pvalue = 1e-05,   
                        proxies = "None", 
                        r2 = 0.8,         
                        build = 38)$results  
    results[[length(results)+1]] <- res
  }
  results <- do.call(rbind, results)
  filename <- basename(sub("\\.csv$","",i))
  fwrite(results, file = paste0("./PHENOdata5/",filename))
}

library("pacman")
p_load(data.table,dplyr) 
df_list <- list()
files <- list.files(pattern="*noBMI.txt")
files
for(f in files){
  df <- fread(f, header=TRUE)
  df$id <- f
  df$id <- gsub(".txt", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
combined_df=combined_df[,-c("id")]
combined_df <- combined_df[!duplicated(combined_df),]
fwrite(combined_df,"allpheno6.csv")

library("pacman")
p_load(data.table,dplyr,TwoSampleMR,MRInstruments,purrr) 
workdir <- getwd()
workdir <- paste0(getwd(), "/FClumpexposure4/")
workdir
files <- list.files(path = workdir,pattern=".csv", full.names = TRUE)
files 
if (!dir.exists("PHONEresult7")) {
  dir.create("PHONEresult7")
}
phdata <- fread("allpheno6.csv")
confounders <- c("smoking")
pattern <- paste0("^(", paste(confounders, collapse = "|"), ")$")
remove_snps <- phdata %>%
  filter(grepl(pattern, trait))
remove_snps
get_print_filepath <- function(files) {
  filepath <- paste0("./PHONEresult7/", basename(files))
  cat("File path:", filepath, "\n")
  return(filepath)
}
get_print_filepath
for (i in seq_along(files)) {
  expo_data <- read_exposure_data(
    filename = files[i],
    sep = ",", 
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure", 
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = "eaf.exposure",
    chr_col = "chr.exposure",
    pos_col = "pos.exposure",
    samplesize_col = "samplesize.exposure"
  )
  if (nrow(remove_snps) > 0) {
    remove_snps$phenotype <- files[i]
    expo_data <- expo_data[!(expo_data$SNP %in% remove_snps$snp), ]
  }
  filepath <- get_print_filepath(files[i])
  fwrite(expo_data, filepath, na = "", quote = FALSE)
}

library("pacman")
p_load(data.table,ggplot2,TwoSampleMR,purrr) 
library(MRInstruments)
library(plyr)
library(dplyr)
orig_dir <- getwd()
workdir <- getwd() 
workdir <- paste0(getwd(),"/PHONEresult7/")
workdir
files <- list.files(path = workdir,pattern=".csv", full.names = TRUE)
print(files)
print(files)
files
outcomefile <- "gen_samplesize_phenocode-151.tsv"
if(!dir.exists("ORdata9")){
  dir.create("ORdata9")
}
if(!dir.exists("Pleiotropydata9")){
  dir.create("Pleiotropydata9")
}
if(!dir.exists("Pressodata9")){
  dir.create("Pressodata9")
}
if(!dir.exists("Heterogeneity9")){
  dir.create("Heterogeneity9")
}
#result=data.frame()
for (i in files) {
  expo_data <- read_exposure_data(
    filename = i, 
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure", 
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = "eaf.exposure",
    chr_col = "chr.exposure",
    pos_col = "pos.exposure",
    samplesize_col = "samplesize.exposure")
  outc_data <- read_outcome_data(
    snps = expo_data$SNP,
    filename = outcomefile,
    sep = ",", 
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    #eaf_col = "eaf",
    chr_col = "#chrom" ,
    pval_col = "pval",
    pos_col = "pos",
    samplesize_col = "samplesize"
  )
  harm_data <- harmonise_data(exposure_dat = expo_data, outcome_dat = outc_data,action = 2)
  mr_result <- mr(harm_data)
  result_or <- generate_odds_ratios(mr_result)
  filename <- basename(sub("\\.txt$","",i))
  dir.create(filename )
  write.table(harm_data, file = paste0(filename,"/harmonise.csv"),
              row.names = F, sep = "\t", quote = F)
  write.table(result_or[,5:ncol(result_or)], 
              file = paste0(filename,"/OR.csv"),
              row.names = F, sep = "\t", quote = F)
  write.table(result_or[, 5:ncol(result_or)], 
              file = paste0("./ORdata9/",filename, "_OR.csv"), 
              row.names = FALSE, sep = "\t", quote = FALSE)
  p1 <- mr_scatter_plot(mr_result, harm_data)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), 
         width=8, height=8)
  pleiotropy <- mr_pleiotropy_test(harm_data)
  write.table(pleiotropy, file = paste0(filename,"/pleiotropy.csv"),
              sep = "\t", quote = F)
  write.table(pleiotropy, 
              file = paste0("./Pleiotropydata9/",filename, "_pleiotropy.csv"), 
              sep = "\t", quote = F)
  heterogeneity <- mr_heterogeneity(harm_data)
  write.table(heterogeneity, file = paste0(filename,"/heterogeneity.csv"),
              sep = "\t", quote = F)
  write.table(heterogeneity, 
              paste0("./Heterogeneity9/",filename,"_heterogeneity.csv"),
              sep = "\t", quote = F)
# MR-PRESSO
  presso <- run_mr_presso(harm_data, NbDistribution = 1000)
  capture.output(presso, file = paste0(filename,"/presso.csv"))
  write.table(presso[[1]]$`Main MR results`, 
              paste0("./Pressodata9/",filename,"_pressodata.csv"),
              sep = "\t", quote = F)
  singlesnp_res <- mr_singlesnp(harm_data)
  singlesnpOR <- generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR, file=paste0(filename,"/singlesnpOR.csv"),
              row.names = F, sep = "\t", quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
# Leave-one-out
  sen_res <- mr_leaveoneout(harm_data)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), 
         width=8, height=8)
# Funnel plot
  res_single <- mr_singlesnp(harm_data)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
}

library("pacman")
p_load(data.table,dplyr) 
df_list <- list()
files <- list.files(pattern="*_heterogeneity.csv")
files
for(f in files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub("_heterogeneity.csv", "", df$id)
  df_list <- c(df_list, list(df))
}
df
#
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
combined_df$method
ivw_df <-combined_df %>% 
  filter(method =="Inverse variance weighted")
fwrite(combined_df,"allheterogeneity10.csv")
fwrite(ivw_df,"all_ivw_heterogeneity10.csv")
##2)OR
library("pacman")
library("pacman")
p_load(data.table,dplyr) 
df_list <- list()
files <- list.files(pattern="*_OR.csv")
files
for(f in files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub("_OR.csv", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allOR10.csv")
ivw_data <- combined_df %>% 
  filter(method=="Inverse variance weighted")
fwrite(ivw_data,"ivw_allOR10.csv")
##3)Pleiotropydata
library("pacman")
p_load(data.table,dplyr) 
df_list <- list()
files <- list.files(pattern="*_pleiotropy.csv")
files
for(f in files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub("_pleiotropy.csv", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allplei10.csv")
##4)Pressodata
library("pacman")
p_load(data.table,dplyr) 
df_list <- list()
files <- list.files(pattern="*_pressodata.csv")
files
for(f in files){
  df <- fread(f)[1,]
  df$id <- f
  df$id <- gsub("_pressodata.csv", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allpresso10.csv")