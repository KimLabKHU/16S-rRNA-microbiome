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

# cd /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/
# zcat /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/IBD_level-6_rank_inv_DESEQ_chr19.result.meta_score_imputed_222.MetaScore.assoc.gz | grep 49218060 | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/autosome_IBD_222_plus_top_v.txt
# zcat /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/IBD_level-6_rank_inv_DESEQ_chr19.result.meta_score_imputed_222.MetaScore.assoc.gz | grep 49214274 | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' >> /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/autosome_IBD_222_plus_top_v.txt
# zcat /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/w_duration/IBD_level-6_rank_inv_DESEQ_chr19.result.meta_score.MetaScore.assoc.gz | grep -v "#" | grep -v "POS" | grep -v "NA" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > autosome_Ruminococcus.input
# echo -e "MARKER\tREF\tALT\tN\tPVALUE\tBETA\tSEBETA" > autosome_info_header
# cat autosome_IBD_222_plus_top_v.txt autosome_Ruminococcus.input | sort -n -k1,1 > autosome_Ruminococcus_merged.input
# cat autosome_info_header autosome_Ruminococcus_merged.input > autosome_Ruminococcus_merged_w.input
# rm autosome_Ruminococcus_merged.input
# mv autosome_Ruminococcus_merged_w.input autosome_Ruminococcus_merged.input

# #### ANCOMBC
# cd /kimlab_wd/rhdfyd/IBD_study/IBD_psy/mbQTL/result/checking/g__Ruminococcus
# zcat IBD_Genus_rank_inv_DESEQ_chr19.result.meta_score_imputed_222.MetaScore.assoc.gz | grep 49218060 | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > autosome_IBD_222_plus_top_v.txt
# zcat IBD_Genus_rank_inv_DESEQ_chr19.result.meta_score_imputed_222.MetaScore.assoc.gz | grep 49214274 | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' >> autosome_IBD_222_plus_top_v.txt
# zcat IBD_Genus_rank_inv_DESEQ_chr19.result.meta_score.MetaScore.assoc.gz | grep -v "#" | grep -v "POS" | grep -v "NA" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > autosome_Ruminococcus.input
# echo -e "MARKER\tREF\tALT\tN\tPVALUE\tBETA\tSEBETA" > autosome_info_header
# cat autosome_IBD_222_plus_top_v.txt autosome_Ruminococcus.input | sort -n -k1,1 > autosome_Ruminococcus_merged.input
# cat autosome_info_header autosome_Ruminococcus_merged.input > autosome_Ruminococcus_merged_w.input
# rm autosome_Ruminococcus_merged.input
# mv autosome_Ruminococcus_merged_w.input autosome_Ruminococcus_merged.input


###
# meta-analysis
# /home1/rhdfyd/code/pipeline_QTL/4_1.running_metal.csh /kimlab_wd/rhdfyd/IBD_study/IBD_psy/correlation/mbQTL_input/meta/genus/Ruminococcus_torques Ruminococcus_torques
# paste Ruminococcus_torques_sort.chr_pos Ruminococcus_torques_meta_result.info | cat meta_result_header - > autosome_result_sorted_w_h
# if [ $# -ne 3 ]; then
# 	echo "USAGE: ./0_mbQTL_input mbQTL_result_dir out_dir disease"
# exit 1
# else

# # arguments ##
# WD=$1
# out_dir=$2
# disease=$3

# mkdir ${out_dir}

# # ls $WD | grep k__ > $out_dir/taxa_list
# ls -l $WD | grep "^d" | awk '{print $9}' > $out_dir/taxa_list
# count=(`wc -l $out_dir/taxa_list`)

# for t_count in $(seq 1 $count); do

# taxa_name=`head -n${t_count} $out_dir/taxa_list | tail -n1`
# taxa_dir=`find $WD -name ${taxa_name}`

# mkdir $out_dir/${taxa_name}

# # Rvtest
# zcat chr*.result.meta_score.MetaScore.assoc.gz |  grep -v "#" | grep -v "POS" | grep -v "NA" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$16"\t"$15"\t"(1/$14)}' > autosome_$taxa_name.input


# done

# fi
