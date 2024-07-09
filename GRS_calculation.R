library(dplyr)

args = commandArgs(TRUE)
argsLen <- length(args)

if(argsLen != 3) stop('error: right arguments (data_dir[dosage & vcf] loci_list updated_or_not(Y/N)')

# arguments ####
data_dir = args[1] # directory where dosage and vcf files exit 
loci_list = args[2] # list of collected associated loci (updated or not)
updated = args[3] # updated list or not?

# output ####
all_ds = c()
# data_dir = "/kimlab_wd/rhdfyd/IBD_study/IBD_GRS/GRS/GRS_result/Control"
# loci_list = "/kimlab_wd/rhdfyd/IBD_study/IBD_GRS/GRS/loci/final_loci/IBD_grs_input_final.txt"
# updated = "N"

# check chromosome in loci_list [21.03.04 updated]
temp_list = read.table(loci_list, header=T,stringsAsFactors=F)
chromosomes = temp_list[,"chr"] %>% sort %>% unique

#

# making single dosage file with allele informaiton ####
#for (i in 1:22){ [21.03.04 updated]
for (i in chromosomes){

  # updated (*_2nd) or not (*_1st) ? ##
  if (updated == "Y"){
  input_prefix = paste(data_dir,"/chr",i,"_loci_2nd", sep="")
  }
  else if (updated == "N"){
  input_prefix = paste(data_dir,"/chr",i,"_loci_1st", sep="")
  }
  
  # dosage + allele ##
  # dosage file
  temp_ds = read.table(file=paste(input_prefix,".DS.FORMAT",sep=""), header=T, stringsAsFactors=F) 
  # vcf file for allele information
  temp_vcf = read.table(file=paste(input_prefix,".recode.vcf",sep=""), header=F, stringsAsFactors=F) 
  temp_vcf = temp_vcf[,1:5]; colnames(temp_vcf) = c("chr","pos","SNP","ref","alt")
  # attaching alleles to dosage
  temp_ds = cbind(temp_vcf, temp_ds[,3:ncol(temp_ds)])
  
  # merging into single dosage file
  all_ds = rbind(all_ds, temp_ds)
}

# loci_list = read.table(loci_list, header=T,stringsAsFactors=F)

# making dosage + effect file with matched alleles ####
all_ds_with_OR = inner_join(all_ds, temp_list, by = c("chr","pos")) %>% mutate(check = case_when(ref == A1 & alt == A2 ~ "OK", ref == A2 & alt == A1 ~ "OK")) %>% # matching alleles
                 filter(check == "OK") %>% # filtering unmatched pairs
                 mutate(effect_matched = case_when(ref == A2 & alt == A1 ~ log(OR), 
                                                   ref == A1 & alt == A2 ~ -log(OR))) # converting OR to beta

# calculating GRS ####
sample_size = ncol(temp_ds) - 5
GRS = all_ds_with_OR[,6:(sample_size+5)] * all_ds_with_OR$effect_matched
GRS = colSums(GRS) %>% data.frame(ID = names(.), GRS = .)

# saving results ####
setwd(data_dir)

write.table(all_ds_with_OR, file="Final_dos_in_loci_with_OR.txt",sep="\t", quote=F,row.names=F,col.names=T)
save(all_ds_with_OR, file="Final_dos_in_loci_with_OR.RData")

write.table(GRS, file="Final_GRS.txt",sep="\t", quote=F,row.names=F,col.names=T)
save(GRS, file="Final_GRS.RData")
