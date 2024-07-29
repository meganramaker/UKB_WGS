################################################################################
#UKB WGS VCF Bed File Creation
#Megan E Ramaker, PhD
#Svati Shah Lab
#Duke University
#July 29, 2024
################################################################################
#R script for looping through individual files containing the start position for 
#each UKB DRAGEN WGS VCF file to annotate the chr, start and stop position for each file
#in one bed file for easier WGS variant filtering in the UKB
#The UKB_VCF_INDEX directory contains an individual bed file for each DRAGEN WGS VCF 
#file in the UKB RAP created using the following DX command in a bash script for each chromosome
#Example of file contents: chrX 156021119
################################################################################
#Bash script
################################################################################
#!/bin/bash

# #set this to the project
# project="HF_48785_Dec23"
# #set this to the wgs directory that you want 
# wgs_file_dir="/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, pVCF format [500k release]/chrX/"
# #set this to the wgs data field for your release
# data_field="ukb24310"
# #set this to the output directory
# data_file_dir="/Megan/MYH6_Enhancer_Analysis/WGS_VCF_INDEX/"
# run_vcf_index="bcftools query -f'%CHROM %POS\n' ${data_field}_cX_b${i}_v1.vcf.gz | head -n1 > ${data_field}_cX_b${i}_v1_pos.txt"
# dx run swiss-army-knife -iin="${wgs_file_dir}/${data_field}_cX_b${i}_v1.vcf.gz" \
# -icmd="${run_vcf_index}" --tag="wgs_vcf_index" --instance-type "mem1_ssd1_v2_x8"\
# --destination="${project}:/${data_file_dir}" --brief --yes
################################################################################
##R code for creating final bed file
library(data.table)
options(scipen=999)
system("./fetchChromSizes hg38 > ~/Documents/UKB/hg38.chrom.sizes")
chrom_lengths<-fread("~/Documents/UKB/hg38.chrom.sizes")
setwd("~/Documents/UKB/WGS_VCF_INDEX/")
files<-list.files(getwd())[-which(file.size(list.files(getwd()))==0)]
chrs<-unique(unlist(lapply(1:length(files), function(x) strsplit(files[x], "_")[[1]][2])))
#Need to put chrs in order for bed file here
chrs_ideal<-c(paste0("c",1:22),"cX")
chrs<-chrs_ideal[which(chrs_ideal%in%chrs)]
bed_list<-list()
for(x in 1:length(chrs)){
  chr<-chrs[x]
  print(chr)
  chr_files<-files[grep(paste0("_",chr,"_"),files)]
  file_list<-list()
  for(i in 1:length(chr_files)){
    file<-chr_files[i]
    file_list[[paste0(file)]]<-data.frame(fread(file))
  }
  file_db<-data.frame(do.call(rbind,file_list))
  file_db[,3]<-NA
  file_db<-file_db[order(file_db$V2),]
  for(j in 1:nrow(file_db)-1){
    file_db[j,3]<-file_db[j+1,2]-1
  }
  file_db[,4]<-rownames(file_db)
  file_db[nrow(file_db),3]<-chrom_lengths[which(chrom_lengths[,1]==file_db[nrow(file_db),1]),2]
  bed_list[[paste0(chr)]]<-file_db
}
tmp<-data.frame(do.call(rbind, bed_list))
rownames(tmp)<-NULL
write.table(tmp, file = "~/Documents/UKB/UKB_DRAGEN_WGS_VCF.bed", sep = "\t", row.names = F, col.names = F, quote = F)
