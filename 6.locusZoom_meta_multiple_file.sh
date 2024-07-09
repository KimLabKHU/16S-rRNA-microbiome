#!/bin/bash

if [ "$#" -ne 3 ]; then
        echo "USAGE: ./locusZoom_epacts.sh Working_dir, vcf_file_dir (calculate LD), distance to define loci"

exit 1
fi

eval "$(conda shell.bash hook)"
conda activate python2.7

WD=$1
vcf_file_DIR=$2
dist=$3

count=(`wc -l $WD/Temp_taxa_list`)

for t_count in $(seq 1 $count); do

taxa_name=`head -n${t_count} $WD/Temp_taxa_list | tail -n1`

taxa_dir=`find $WD -name ${taxa_name}`

awk '{if ($7 < 5e-08) print}' $taxa_dir/LDSC/first_input_with_rsID | awk '{OFS="\t"}{print $3,$1,$2,$7}' > $taxa_dir/LDSC/mbQTL_signal_5e-8.txt

if [ -s $taxa_dir/LDSC/mbQTL_signal_5e-8.txt ] ; then
flank=$(($dist + 200))

output_dir=${taxa_dir}/locuszoom_result
mkdir $output_dir

Rscript-4.1 summarize.R $taxa_dir/LDSC/mbQTL_signal_5e-8.txt $dist $flank ${output_dir} ${taxa_name} 5e-8

mkdir $taxa_dir/overlapped_w_knownLoci_result
Rscript-4.1 Overlap_checking_for_QTL.R $output_dir/associated_region.summary_${taxa_name}_$dist-5e-8 $taxa_dir/overlapped_w_knownLoci_result

count=(`cat $output_dir/associated_region.summary_${taxa_name}_${dist}-5e-8 | grep -vw "LocusNo" |wc -l`)

for locus_count in $(seq 1 $count); do

Chr_number=`cat $output_dir/associated_region.summary_${taxa_name}_${dist}-5e-8 | grep -vw "LocusNo" | awk '{print $8}' | head -n${locus_count} | tail -n1`
POS_Start=`cat $output_dir/associated_region.summary_${taxa_name}_${dist}-5e-8 | grep -vw "LocusNo" | awk '{print $6}' | head -n${locus_count} | tail -n1`
POS_End=`cat $output_dir/associated_region.summary_${taxa_name}_${dist}-5e-8 | grep -vw "LocusNo" | awk '{print $7}' | head -n${locus_count} | tail -n1`

echo -e "MarkerName\tP-value" > ${output_dir}/header_metal
sed 1d $taxa_dir/LDSC/first_input_with_rsID | awk '{OFS="\t"}{print $3,$7}' > ${output_dir}/metal_locus
cat ${output_dir}/header_metal ${output_dir}/metal_locus > ${output_dir}/metal_locus_final

final_vcf_file=$vcf_file_DIR/chr${Chr_number}.*.vcf.gz

locuszoom \
--metal ${output_dir}/metal_locus_final \
--ld-vcf ${final_vcf_file} \
--build hg19 \
--chr ${Chr_number} \
--start ${POS_Start} --end ${POS_End} \
--rundir ${output_dir} \
--gene-tab gencode

done

else
        echo -e "$taxa_dir/LDSC/mbQTL_signal_5e-8.txt is empty"
        # rm $taxa_dir/LDSC/mbQTL_signal_5e-8.txt
fi

done

Date=`date '+%F  %r' | cut -d ' ' -f 1`
echo "${Date} locusZoom_epacts_all.sh ${WD} ${vcf_file_DIR} ${dist}" > ${WD}/${Date}.6.locusZoom_meta_multiple_file.sh.history

