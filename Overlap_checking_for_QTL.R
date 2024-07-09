#! /usr/bin/Rscript-4.1

args = commandArgs(TRUE)
argsLen <- length(args)

if(argsLen != 2) stop('error: wrong number of arguments (5e-08_loci_file, OUTDIR)')

library(dplyr)
library(bigreadr)

# arguments #####
loci_file = args[1] #All.results.5e-8.txt
outdir=args[2]


mbQTL_loci = fread2(loci_file)

ax_related_loci = fread2("anxiety_GWAS.txt")

dp_related_loci_1 = fread2("depression_GWAS_1.txt")
dp_related_loci_2 = fread2("depression_GWAS_2.txt")

IBD_related_loci_meta = fread2("IBD_GWAS_meta.txt")
IBD_related_loci_EUR = fread2("IBD_GWAS_EUR.txt")

UC_related_loci_meta = fread2("UC_GWAS_meta.txt")
UC_related_loci_EUR = fread2("UC_GWAS_EUR.txt")

CD_related_loci_meta = fread2("CD_GWAS_meta.txt")
CD_related_loci_EUR = fread2("CD_GWAS_EUR.txt")

ax_overlap <- subset(ax_related_loci, CHR == mbQTL_loci$CHR & POS>= mbQTL_loci$START_pos & POS <= mbQTL_loci$END_pos)
dp_overlap_1 <- subset(dp_related_loci_1, CHR == mbQTL_loci$CHR & BP>= mbQTL_loci$START_pos & BP <= mbQTL_loci$END_pos)
dp_overlap_2 <- subset(dp_related_loci_2, CHR == mbQTL_loci$CHR & POS>= mbQTL_loci$START_pos & POS <= mbQTL_loci$END_pos)

IBD_overlap_meta <- subset(IBD_related_loci_meta, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)
IBD_overlap_EUR <- subset(IBD_related_loci_EUR, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)

UC_overlap_meta <- subset(UC_related_loci_meta, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)
UC_overlap_EUR <- subset(UC_related_loci_EUR, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)

CD_overlap_meta <- subset(CD_related_loci_meta, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)
CD_overlap_EUR <- subset(CD_related_loci_EUR, chromosome == mbQTL_loci$CHR & base_pair_location >= mbQTL_loci$START_pos & base_pair_location <= mbQTL_loci$END_pos)


file_save <- function(result, disease){
   if (nrow(result) > 0) {
      write.table(result, file= paste0(outdir,"/",disease,"_physical_overlapped.result.txt"),sep="\t",col.names = T,row.names = F,quote = F)
   }
}

file_save(ax_overlap,"Ax")
file_save(dp_overlap_1,"Dp_MDD_sig")
file_save(dp_overlap_2,"Dp_33859377")
file_save(IBD_overlap_meta,"IBD_meta")
file_save(IBD_overlap_EUR,"IBD_EUR")
file_save(UC_overlap_meta,"UC_meta")
file_save(UC_overlap_EUR,"UC_EUR")
file_save(CD_overlap_meta,"CD_meta")
file_save(CD_overlap_EUR,"CD_EUR")
