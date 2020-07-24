#########################################
#MR pQTL analysis script -- Jie Zheng #
#########################################

###install the Two sample MR package (just need do this once) 
##source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")

##to update the R package (once there is a )
#library(devtools)
#install_github("MRCIEU/TwoSampleMR")

#example of use the older version of the package
#devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")

##call necessary libraries
library(TwoSampleMR)
#library(MRInstruments)
library("readxl") #plesae install the package
rm(list=ls(all=TRUE)) 

##setup your working folder
setwd("/your/working/folder/")

##read in the exposure data from a file
Ins<-data <- read_excel("./data/Instruments.xlsx",1)

Ins<-format_data(Ins, type = "exposure", header = TRUE,
                 phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                 se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele", pval_col = "pval")


##read in the outcome data
ao<-available_outcomes(access_token=NULL)

ids<-as.character(unlist(read.table("./data/outcome.id.txt",header=F)))

outcome_dat<- extract_outcome_data(
  snps = Ins$SNP,
  outcomes = ids)

##harmonise the exposure and outcome data
dat <- NULL
dat <- harmonise_data(
  exposure_dat = Ins, 
  outcome_dat = outcome_dat
)

##run the MR and sensitivity analyses 
mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))  # main MR analysis
mr_hetero <- mr_heterogeneity(dat) # heterogeneity test across instruments
mr_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test  
try(mr_single <- mr_singlesnp(dat)) #single SNP MR using Wald ratio

##save the MR results 
exposure <- "pQTL-combined"
result_file0 <- paste0("./your/results/folder/",exposure,".harmonise.txt")
result_file <- paste0("./your/results/folder/",exposure,".mr.txt")
result_file2 <- paste0("./your/results/folder/",exposure,".mr_hetero.txt")
result_file3 <- paste0("./your/results/folder/",exposure,".mr_pleio.txt")
result_file4 <- paste0("./your/results/folder/",exposure,".mr_single.txt")
if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}


