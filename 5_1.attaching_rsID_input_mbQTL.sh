#!/bin/bash 

if [ $# -ne 1 ]; then
	echo "USAGE: ./1_attaching_rsID_input_mbQTL WD "
exit 1
else

# arguments ##

WD=$1 # working directory 

# dbsnp information
dbsnp=hg19_avsnp150.txt 
echo -e "rsid\ta1\ta2\tN\tbeta\tSE" > $WD/header

mkdir ${WD}/LDSC
outdir=${WD}/LDSC

cat ${WD}/autosome_result_sorted | cut -f 1,2  > ${WD}/SNP_info_for_grepping_dbsnp # CHR POS for "grep"
fgrep -wf ${WD}/SNP_info_for_grepping_dbsnp ${dbsnp} | cut -f 1,2,4,5,6 | sort -rk 1,2 > ${WD}/SNP_info_rsID_for_down # rsID of variants[chr pos ref alt ID]

# R for redundant rsID and matching alleles ##
Rscript matching_markername_with_rsID.R ${WD}/SNP_info_rsID_for_down ${WD}/autosome_result_sorted ${outdir}

fi