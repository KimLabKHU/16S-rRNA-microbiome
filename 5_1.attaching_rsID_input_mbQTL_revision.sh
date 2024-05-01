#!/bin/bash 

if [ $# -ne 1 ]; then
	echo "USAGE: ./1_attaching_rsID_input_mbQTL WD "
exit 1
else

# arguments ##

WD=$1 # working directory 
# count=(`wc -l $WD/taxa_list`)

# count=(`wc -l $WD/Temp_taxa_list`)

# dbsnp information
dbsnp=/home1/eha/packages/annovar/humandb/hg19_avsnp150.txt 
echo -e "rsid\ta1\ta2\tN\tbeta\tSE" > $WD/header

# for t_count in $(seq 1 $count); do

# taxa_name=`head -n${t_count} $WD/taxa_list | tail -n1`
mkdir ${WD}/LDSC
outdir=${WD}/LDSC

cat ${WD}/autosome_result_sorted | cut -f 1,2  > ${WD}/SNP_info_for_grepping_dbsnp # CHR POS for "grep"
fgrep -wf ${WD}/SNP_info_for_grepping_dbsnp ${dbsnp} | cut -f 1,2,4,5,6 | sort -rk 1,2 > ${WD}/SNP_info_rsID_for_down # rsID of variants[chr pos ref alt ID]
# R for redundant rsID and matching alleles ##
Rscript /home1/rhdfyd/code/Popcorn/matching_markername_with_rsID.R ${WD}/SNP_info_rsID_for_down ${WD}/autosome_result_sorted ${outdir}

fi