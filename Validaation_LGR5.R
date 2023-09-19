setwd("~/Documents/Validation_LGR5_Ican")
library(data.table)
data=fread("ukb40011_histology.pcs.csv")
data=fread("ukb40006_cancer_icd10.pcs.csv")
data1=as.matrix(data.frame(data))
data2=data1[rowSums(is.na(data1))>0,]


sed 's/","/'115

setwd("/Users/anilkuma/Documents/Validation_LGR5_Ican/file/")
library(data.table)
data=fread("finngen_R9_finngen_R9_analysis_data_summary_stats_release_finngen_R9_C3_SMALL_INTESTINE_NEUROENDOCRINE_EXALLC")

library(qqman)







snpsOfInterest=c("cg19693031","cg01178710","cg04673737","cg11024682","cg00574958")
pdf("manhattan_plot_manhattanagexcellcompepc12genepc12.pdf",width=8)
manhattan(
 data1,
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  col = c("red","blue","pink","black","cyan","coral1","darkorchid","deepskyblue", "gray60","blue4","orange3","darkolivegreen1","aquamarine","bisque","blue2","darkgreen","brown3","darkorchid","darkred","green","hotpink","navy","mediumblue"),
  suggestiveline = -log10(1e-05),
  genomewideline = -log10(5e-08),
 highlight = snpsOfInterest,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE)
 dev.off()

 data1=fread("temporary_file_sharing_C3_SMALL_INTESTINE_NEUROENDOCRINE_EXALLC_meta_C3_SMALL_INTESTINE_NEUROENDOCRINE_EXALLC_meta_out_filtered.tsv")
 data1=data1[data1$all_inv_var_meta_p<0.01,]
 snpsOfInterest=c("rs193194921","rs1019866379","rs200138614","rs2660861","rs8009633","rs6043873","rs2328547")
 
 pdf("manhattan_plot_meta.pdf",width=25)
 manhattan(
   data1,
   chr = "#CHR",
   bp = "POS",
   p = "all_inv_var_meta_p",
   snp = "rsid",
   col = c("blue","red","pink","black","cyan","coral1","darkorchid","deepskyblue","gray10", "gray60","blue4","orange3","darkolivegreen1","aquamarine"),
   chrlabs = NULL,
   suggestiveline = -log10(1e-05),
   genomewideline = -log10(5e-08),
   highlight = snpsOfInterest,
   logp = TRUE,
   annotatePval = NULL,
   annotateTop = TRUE)
 dev.off()
 
data_s1=data[data$`20001-0.0`=="1027",]
data_s2=data[data$`20001-0.1`=="1027",]
data_s3=data[data$`20001-0.2`=="1027",]
data_s4=data[data$`20001-0.3`=="1027",]
data_s5=data[data$`20001-0.4`=="1027",]
data_s6=data[data$`20001-0.5`=="1027",]
data_s7=data[data$`20001-0.6`=="1027",]
data_s8=data[data$`20001-0.7`=="1027",]
data_s9=data[data$`20001-0.8`=="1027",]
data_s10=data[data$`20001-0.9`=="1027",]
data_s11=data[data$`20001-1.0`=="1027",]

data_s12=data[data$`20001-1.1`=="1027",]
data_s13=data[data$`20001-1.2`=="1027",]
data_s14=data[data$`20001-1.3`=="1027",]
data_s15= data[data$`20001-1.5`=="1027",]

data_s5=data[data$`20001-0.5`=="1027",]
data_s6=data[data$`20001-0.6`=="1027",]
data_s7=data[data$`20001-0.7`=="1027",]
data_s8=data[data$`20001-0.8`=="1027",]
data_s9=data[data$`20001-0.9`=="1027",]
data_s10=data[data$`20001-1.0`=="1027",]


data1 <- data[data$`40011-0.0`=="8936",]
data2 <- data[data$`40011-1.0`=="8936",]
data3 <- data[data$`40011-2.0`=="8936",]
data4 <- data[data$`40011-3.0`=="8936",]
data5 <- data[data$`40011-4.0`=="8936",]
data6 <- data[data$`40011-5.0`=="8936",]
data7 <- data[data$`40011-6.0`=="8936",]
data8 <- data[data$`40011-7.0`=="8936",]
data9 <- data[data$`40011-8.0`=="8936",]
data10 <- data[data$`40011-9.0`=="8936",]
data11 <- data[data$`40011-10.0`=="8936",]
data12 <- data[data$`40011-11.0`=="8936",]
data13 <- data[data$`40011-12.0`=="8936",]
data14 <- data[data$`40011-13.0`=="8936",]

data15 <- data[data$`40011-14.0`=="8936",]
data16 <- data[data$`40011-15.0`=="8936",]
data17 <- data[data$`40011-16.0`=="8936",]
data18 <- data[data$`40011-17.0`=="8936",]
data19 <- data[data$`40011-18.0`=="8936",]
data20 <- data[data$`40011-19.0`=="8936",]
data21 <- data[data$`40011-20.0`=="8936",]
data22 <- data[data$`40011-21.0`=="8936",]
data23 <- data[data$`40011-22.0`=="8936",]


data24 <- data[data$`40011-23.0`=="8936",]
data25 <- data[data$`40011-24.0`=="8936",]
data26 <- data[data$`40011-25.0`=="8936",]
data27 <- data[data$`40011-26.0`=="8936",]
data28 <- data[data$`40011-27.0`=="8936",]
data29 <- data[data$`40011-28.0`=="8936",]
data30 <- data[data$`40011-29.0`=="8936",]
data31 <- data[data$`40011-30.0`=="8936",]
data32 <- data[data$`40011-31.0`=="8936",]

data_all=data[data$eid%in%common,]

datanew=rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,data21,data22,data23,data24,data25,data26,data27,data28,data29,data30,data31,data32)
sample=unique(datanew$eid)
sample_8936=cbind (sample,rep("8936",length(sample)))
write.csv(sample_8936,"sample_8936.csv")

sample_1432=cbind(sample_1432,sample_1432[,2])
sampleall=rbind(sample_1432,sample_8240,sample_8241,sample_8244,sample_8246,sample_8249)

data_comb=c()
files=list.files(pattern = "UKBB_SI_pgen_")
for (i in 1:length(files )){
  data=fread(files[[i]])
  data1=data[data$TEST=="ADD",]
  data1=data1[data1$P<0.05,]
  data_comb=rbind(data_comb, data1)
  print (i)
}
write.csv( data_comb,"data_comb.csv")


pdf("manhattan_NET_UKBB.pdf")
manhattan(data1,col=c("red","blue","green","#EE1289","#9400D3","#1C86EE","#6A5ACD","#00FF7F","#00F5FF","#EE5C42","#6B8E23","#8B8682","#EE7942","#912CEE","#FFFF00","#FFBBFF","#FF34B3","#87CEEB","#54FF9F","#0000FF","#B8860B","#DEB887","#458B00","#CD3333"),suggestiveline=-log10(1e-5),genomewideline=-log10(1e-7))
dev.off()