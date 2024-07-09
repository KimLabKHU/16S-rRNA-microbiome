#!/bin/bash 

if [ $# -ne 3 ]; then
	echo "USAGE: ./0_mbQTL_input mbQTL_result_dir out_dir disease"
exit 1
else

# arguments ##
WD=$1
out_dir=$2
disease=$3

mkdir ${out_dir}

# ls $WD | grep k__ > $out_dir/taxa_list
ls -l $WD | grep "^d" | awk '{print $9}' > $out_dir/taxa_list
count=(`wc -l $out_dir/taxa_list`)

for t_count in $(seq 1 $count); do

taxa_name=`head -n${t_count} $out_dir/taxa_list | tail -n1`
# taxa_dir=`find $WD -name '${taxa_name}'`

mkdir $out_dir/${taxa_name}

# Rvtest
zcat ${WD}/"${taxa_name}"/IBD_Genus_rank_inv_DESEQ_chr*.result.meta_score.MetaScore.assoc.gz |  grep -v "#" | grep -v "POS" | grep -v "NA" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > autosome_$taxa_name.input
done

fi
