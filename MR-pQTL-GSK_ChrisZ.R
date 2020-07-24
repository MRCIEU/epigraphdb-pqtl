#########################################
#MR pQTL analysis script -- Chris Zheng #
#########################################
# source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
#update the R package
library(devtools)
install_github("MRCIEU/TwoSampleMR")
#use the elder version of the package
#devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")


library(ggplot2)
#link to the test version
library(TwoSampleMR)
#library(MRInstruments)
library("readxl")
rm(list=ls(all=TRUE)) 
toggle_dev("elastic")
toggle_dev("test")


#revoke_mrbase_access_token()
#get_mrbase_access_token()
#ao<-available_outcomes()

setwd("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2/projects/pQTL-MR/4.GSK/pQTL-all-diseases/")

Ins<-data <- read_excel("data/Instruments-Suhre.xlsx",2)
Ins<-Ins[Ins$SNP!="rs35233997",]
Ins<-Ins[Ins$SNP!="rs11292716",]
Ins<-Ins[Ins$SNP!="rs5784945",]
Ins<-Ins[Ins$SNP!="rs35400166",]
Ins<-Ins[Ins$SNP!="rs35365539",]
Ins<-Ins[Ins$SNP!="rs11548618",] #
#Ins<-Ins[Ins$SNP!="rs6036478",] 
Ins<-Ins[Ins$SNP!="rs45602433",] #

#for non-MHC snps
#Ins<-Ins[Ins$`MHC/NK` == "NO",]

#test snps with issues
#Ins<-data <- read_excel("data/Instruments-Suhre.xlsx",2)
#Ins<-Ins[Ins$SNP=="rs11548618",] #
#Ins<-Ins[Ins$SNP=="rs6036478",] #
#Ins<-Ins[Ins$SNP=="rs45602433",] #
#Ins <- clump_data(Ins)


Ins<-format_data(Ins, type = "exposure", header = TRUE,
                 phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                 se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele", pval_col = "pval",
                 gene_col = "gene")
Ins <- clump_data(Ins)
save(Ins,file="data/Chris/mr_pqtl_instruments_Suhre2.R")
save(Ins,file="data/Chris/mr_pqtl_instruments_Suhre_noMHC.R")
# load("data/Chris/mr_pqtl_instruments_Suhre2.R")
# load("data/Chris/mr_pqtl_instruments_Suhre_noMHC.R")
# write.table(Ins,"data/Chris/instruments-clump.txt",sep="\t",col.names=T,row.names=F,quote=F)


ao<-available_outcomes()
id.dis<-ao$id
#ao.dis<-ao[ao$category=="Disease",]
#ao.dis<-ao.dis[!is.na(ao.dis$ncase),]
#ao.dis<-ao.dis[ao.dis$mr==1,]
#ao.dis<-ao.dis[order(ao.dis$trait,ao.dis$nsnp,decreasing=T),]
# write.table(ao.dis,"data/study.txt",sep="\t",col.names=T,row.names=F,quote=F)
#ao.dis[,c("trait","ncase","nsnp","id")]
#id.dis<-ao.dis$id

#metastroke and Primary sclerosing cholangitis analysis
#ao<-available_outcomes()
#ao.dis<-ao[ao$filename=="ipscsg2016.result.combined.full.txt.tab.all_pos",]
#ao.dis<-ao[ao$filename=="metastroke.all.chr.bp.tab",]
#id.dis<-c(1108,1109,1110,1111,1112)


outcome_dat<- extract_outcome_data(
  snps = Ins$SNP,
  outcomes = id.dis)

dat <- harmonise_data(
  exposure_dat = Ins, 
  outcome_dat = outcome_dat
)

Res <- mr(dat)
Res_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
Res_hetero <- mr_heterogeneity(dat)
Res_single <- mr_singlesnp(dat)

write.table(Res,"results/Chris/MR_results_Suhre.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_pleio,"results/Chris/MR_pleio_Suhre.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_hetero,"results/Chris/MR_hetero_Suhre.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_single,"results/Chris/MR_single_Suhre.txt",sep="\t",col.names=T,row.names=F,quote=F)
save(list = ls(all.names = TRUE),file="results/Chris/mr_pqtl_Suhre.R")
# load("results/Chris/mr_pqtl_Suhre.R")

# Scatter plots 95Q3tc	33
id.exp<-"95Q3tc"
id.out<-"33"
Res_select<-Res[Res$id.exposure==id.exp & Res$id.outcome==id.out,]
Dat_select<-dat[dat$id.exposure==id.exp & dat$id.outcome==id.out,]
p1 <- mr_scatter_plot(Res_select, Dat_select)
p1

# Forest plot 
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]


#########################################
#MR pQTL analysis script -- Chris Zheng #
#########################################
# source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
#update the R package
library(devtools)
install_github("MRCIEU/TwoSampleMR")

library(ggplot2)
#link to the test version
library(TwoSampleMR)
#library(MRInstruments)
library("readxl")
rm(list=ls(all=TRUE)) 

setwd("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2/projects/pQTL-MR/4.GSK/pQTL-all-diseases/")

Ins<-data <- read_excel("data/instruments-Sun.xlsx",2)
Ins<-Ins[Ins$SNP!="rs536720776",]
Ins<-Ins[Ins$SNP!="rs397780227",]
Ins<-Ins[Ins$SNP!="rs142201367",]
Ins<-Ins[Ins$SNP!="rs34983651",]
Ins<-Ins[Ins$SNP!="rs369160772",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs11434877",]
Ins<-Ins[Ins$SNP!="rs201022770",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs140142103",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="chr12:89825226",]
Ins<-Ins[Ins$SNP!="rs140142103",]
Ins<-Ins[Ins$SNP!="rs139130389",]
Ins<-Ins[Ins$SNP!="rs75077631",]
Ins<-Ins[Ins$SNP!="rs145827860",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="rs142201367",]
Ins<-Ins[Ins$SNP!="rs200997557",]
Ins<-Ins[Ins$SNP!="rs3833490",]
Ins<-Ins[Ins$SNP!="rs371314787",]
Ins<-Ins[Ins$SNP!="rs11390840",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="rs147678983",]
Ins<-Ins[Ins$SNP!="chr5:54327845",]
Ins<-Ins[Ins$SNP!="rs202135714",]
Ins<-Ins[Ins$SNP!="rs10624573",]
Ins<-Ins[Ins$SNP!="rs75077631",]
Ins<-Ins[Ins$SNP!="rs113216780",]
Ins<-Ins[Ins$SNP!="rs148410779",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs3032928",]
Ins<-Ins[Ins$SNP!="rs113358888",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs3216676",]
Ins<-Ins[Ins$SNP!="rs146682150",]
Ins<-Ins[Ins$SNP!="rs200357845",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs34983651",]
Ins<-Ins[Ins$SNP!="rs200436316",]
Ins<-Ins[Ins$SNP!="rs34983651",]
Ins<-Ins[Ins$SNP!="rs140142103",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="rs200357845",]
Ins<-Ins[Ins$SNP!="rs68055275",]
Ins<-Ins[Ins$SNP!="rs147859411",]
Ins<-Ins[Ins$SNP!="rs35157100",]
Ins<-Ins[Ins$SNP!="rs3917539",]
Ins<-Ins[Ins$SNP!="rs139130389",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="rs139121778",]
Ins<-Ins[Ins$SNP!="rs113022368",]
Ins<-Ins[Ins$SNP!="rs3917539",]
Ins<-Ins[Ins$SNP!="rs33944729",]
Ins<-Ins[Ins$SNP!="rs11447348",]
Ins<-Ins[Ins$SNP!="rs34983651",]
Ins<-Ins[Ins$SNP!="chr4:187154272",]
Ins<-Ins[Ins$SNP!="rs372065719",]
Ins<-Ins[Ins$SNP!="rs34813609",]
Ins<-Ins[Ins$SNP!="rs28480494",] #
Ins<-Ins[Ins$SNP!="rs11342894",] #
Ins<-Ins[Ins$SNP!="rs3079633",] #
Ins<-Ins[Ins$SNP!="rs35179000",] #
Ins<-Ins[Ins$SNP!="rs138942288",] #
Ins<-Ins[Ins$SNP!="rs10569394",] #
Ins<-Ins[Ins$SNP!="rs540781594",] #
Ins<-Ins[Ins$SNP!="rs145766673",] #
Ins<-Ins[Ins$SNP!="rs117103342",]
Ins<-Ins[Ins$SNP!="rs144979264",]
Ins<-Ins[Ins$SNP!="rs55942822",]
Ins<-Ins[Ins$SNP!="rs57736976",]
Ins<-Ins[Ins$SNP!="rs2229331",]
Ins<-Ins[Ins$SNP!="rs557663478",]
Ins<-Ins[Ins$SNP!="rs139531404",]
Ins<-Ins[Ins$SNP!="rs59765118",]
Ins<-Ins[Ins$SNP!="rs574174001",]
Ins<-Ins[Ins$SNP!="rs112433249",]
Ins<-Ins[Ins$SNP!="rs35814191",]
Ins<-Ins[Ins$SNP!="rs565102642",]
Ins<-Ins[Ins$SNP!="rs562022020",]
Ins<-Ins[Ins$SNP!="rs193280350",]
Ins<-Ins[Ins$SNP!="rs189821701",]
Ins<-Ins[Ins$SNP!="rs72743610",]
Ins<-Ins[Ins$SNP!="rs34963977",]
Ins<-Ins[Ins$SNP!="rs28461437",]
Ins<-Ins[Ins$SNP!="rs34963977",]
Ins<-Ins[Ins$SNP!="rs62358364",]
Ins<-Ins[Ins$SNP!="rs201263987",]
Ins<-Ins[Ins$SNP!="rs11291564",]
Ins<-Ins[Ins$SNP!="rs587729126",]
Ins<-Ins[Ins$SNP!="rs62037084",]
Ins<-Ins[Ins$SNP!="rs199605734",]
Ins<-Ins[Ins$SNP!="rs573431210",]
Ins<-Ins[Ins$SNP!="rs375375234",]
Ins<-Ins[Ins$SNP!="rs3917549",]
Ins<-Ins[Ins$SNP!="rs36082205",]
save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")

#Ins<-Ins[Ins$`MHC/NK` == "NO",]


#split into three jobs
Ins1<-Ins[1:301,]
Ins2<-Ins[302:603,]
Ins3<-Ins[604:910,]
Ins4<-Ins[911:1213,]
Ins5<-Ins[1214:1516,]
Ins6<-Ins[1517:1822,]

load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins1<-Ins1[Ins1$SNP!="rs41313926",]
Ins1<-Ins1[Ins1$SNP!="rs373214284",]
Ins1<-Ins1[Ins1$SNP!="rs41313926",]
Ins1<-Ins1[Ins1$SNP!="rs373214284",]
Ins1<-Ins1[Ins1$SNP!="rs191448232",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins1 <- clump_data(Ins1))


load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins2<-Ins2[Ins2$SNP!="rs144147445",]
Ins2<-Ins2[Ins2$SNP!="rs150359707",]
Ins2<-Ins2[Ins2$SNP!="rs62358361",]
Ins2<-Ins2[Ins2$SNP!="rs7789303",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins2 <- clump_data(Ins2))


load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins3<-Ins3[Ins3$SNP!="rs41313926",]
Ins3<-Ins3[Ins3$SNP!="rs35223184",]
Ins3<-Ins3[Ins3$SNP!="rs367804438",]
Ins3<-Ins3[Ins3$SNP!="rs151288400",]
Ins3<-Ins3[Ins3$SNP!="rs532540191",]
Ins3<-Ins3[Ins3$SNP!="rs11548618",]
Ins3<-Ins3[Ins3$SNP!="rs35457250",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins3 <- clump_data(Ins3))

load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins4<-Ins4[Ins4$SNP=="rs527942",]
Ins4 <- clump_data(Ins4)

load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins4<-Ins4[Ins4$SNP!="rs77320622",]
Ins4<-Ins4[Ins4$SNP!="rs2116308",]
Ins4<-Ins4[Ins4$SNP!="rs35505705",]
Ins4<-Ins4[Ins4$SNP!="rs8176643",]
Ins4<-Ins4[Ins4$SNP!="rs532540191",]
Ins4<-Ins4[Ins4$SNP!="rs12972028",]
Ins4<-Ins4[Ins4$SNP!="rs199696982",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins4 <- clump_data(Ins4))


load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins5<-Ins5[Ins5$SNP=="rs5757973",]
Ins5 <- clump_data(Ins5)

load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins5<-Ins5[Ins5$SNP!="rs35143187",]
Ins5<-Ins5[Ins5$SNP!="rs140647145",]
Ins5<-Ins5[Ins5$SNP!="rs5757973",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins5 <- clump_data(Ins5))


load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins6<-Ins6[Ins6$SNP=="rs6889164",]
Ins6 <- clump_data(Ins6)

load("data/Chris/mr_pqtl_instruments_Sun.R")
Ins6<-Ins6[Ins6$SNP!="rs398062996",]
Ins6<-Ins6[Ins6$SNP!="rs1799895",]
Ins6<-Ins6[Ins6$SNP!="rs35039172",]
Ins6<-Ins6[Ins6$SNP!="rs139574809",]
Ins6<-Ins6[Ins6$SNP!="rs6889164",]

save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
try(Ins6 <- clump_data(Ins6))


Ins<-format_data(Ins, type = "exposure", header = TRUE,
                 phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                 se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele", pval_col = "pval",
                 gene_col = "gene")
Ins <- clump_data(Ins)
save(list = ls(all.names = TRUE),file="data/Chris/mr_pqtl_instruments_Sun.R")
# load("data/Chris/mr_pqtl_instruments_Sun.R")
# write.table(Ins,"data/Chris/instruments-clump.txt",sep="\t",col.names=T,row.names=F,quote=F)


ao<-available_outcomes()
id.dis<-ao$id

id.dis1<-id.dis[1:200]
id.dis2<-id.dis[201:400]
id.dis3<-id.dis[401:600]
id.dis4<-id.dis[601:800]
id.dis5<-id.dis[801:1000]
id.dis6<-id.dis[1001:1102]

Ins_GP1BA <- Ins3[Ins3$exposure=="GP1BA",]

suffix <- "GP1BA" 

outcome_dat<- extract_outcome_data(
  snps = Ins_GP1BA$SNP,
  outcomes = "279")

dat <- harmonise_data(
  exposure_dat = Ins_GP1BA, 
  outcome_dat = outcome_dat
)

Res <- mr(dat)
Res_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
Res_hetero <- mr_heterogeneity(dat)
Res_single <- mr_singlesnp(dat)

write.table(Res,paste0("results/Chris/MR_results_Sun_",suffix),sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_pleio,paste0("results/Chris/MR_pleio_Sun_",suffix),sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_hetero,paste0("results/Chris/MR_hetero_Sun_",suffix),sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_single,paste0("results/Chris/MR_single_Sun_",suffix),sep="\t",col.names=T,row.names=F,quote=F)
save(list = ls(all.names = TRUE),file=paste0("results/Chris/mr_pqtl_Sun.R"))
# load(paste0("results/Chris/mr_pqtl_Sun.R"))


# Scatter plots 95Q3tc	33
id.exp<-"95Q3tc"
id.out<-"33"
Res_select<-Res[Res$id.exposure==id.exp & Res$id.outcome==id.out,]
Dat_select<-dat[dat$id.exposure==id.exp & dat$id.outcome==id.out,]
p1 <- mr_scatter_plot(Res_select, Dat_select)
p1

# Forest plot 
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

##################forest plot for GP1BA results#####################
library(ggplot2)
library(readxl)
rm(list=ls(all=TRUE)) 

setwd("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2/projects/pQTL-MR/manuscript/")
data <- read_excel("Tables.xlsx",3)
#data <- read_excel("prelimiary-results-v3.xlsx",4)

label <-data$Outcome
LogOR <-data$Beta 
lower<-data$Beta-1.96*data$SE
upper<-data$Beta+1.96*data$SE
Category <- data$Category



df <- data.frame(label, LogOR, lower, upper,Category)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label))

ggplot(data=df, aes(x=label, y=LogOR, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  geom_point(aes(size=LogOR, fill=Category), colour="black",shape=21) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  #xlab("") + ylab("SD change in BMD (95%CI) / log(OR) of fracture per unit change in Sclerostin") +
  xlab("") + ylab("log(OR) of diseases per unit change in GP1BA level") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=15)) 

#ggsave("~/Google Drive/working space/Postdoc_year1_year2/projects/SOST-George/manuscript/forest_plot_sclerostin_bone.png", width=12, height=8)
ggsave("~/Google Drive/working space/Postdoc_year1_year2/projects/pQTL-MR/manuscript/forest_plot_GP1BA_selected_diseases.png", width=12, height=8)









