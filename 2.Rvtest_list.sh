#!/bin/bash 

if [ $# -ne 8 ]; then
	echo "USAGE: ./rvtest.sh vcf_DIR ped_DIR level out_dir chr_start chr_end deseq_result_list rvtest_DIR"
exit 1
else

vcf_dir=$1
ped_dir=$2
level=$3
out_dir=$4
chr_start=$5
chr_end=$6
taxa_list=$7
rvtest_DIR=$8

mkdir $out_dir

cat $taxa_list > $out_dir/temp_list
count=(`wc -l $out_dir/temp_list`)

ped_file=`ls $ped_dir/*_level-${level}_rank_inv_DESEQ.ped | head -n2 | tail -n1`

output_name=`basename $ped_file .ped`

for t_count in $(seq 1 $count); do

taxa_name=`head -n${t_count} $out_dir/temp_list | tail -n1`
mkdir $out_dir/${t_count}_${taxa_name}
final_out=`find $out_dir -name *_${taxa_name}`

for chr_number in $(seq $chr_start $chr_end); do


$rvtest_DIR/rvtest \
--inVcf $vcf_dir/chr${chr_number}_nodup.vcf.gz \
--out $final_out/${output_name}_chr${chr_number}.result.meta_score \
--covar $ped_file --covar-name "Phenotype,Sex,Age,BMI,DAS,Alcohol,Smoking,Sig_OR,PC1,PC2,PC3,PC4,PC5" \
--pheno $ped_file --pheno-name $taxa_name \
--dosage "DS" --meta "score" \
--noweb --numThread 6 --freqUpper 0.995000 --freqLower 0.005000

done
done
mv $out_dir/temp_list $out_dir/taxa_list_${output_name}
Date=`date '+%F  %r' | cut -d ' ' -f 1`
echo "${Date} sh ${vcf_dir} ${ped_dir} ${level} ${out_dir} ${chr_start} ${chr_end} ${taxa_list}" ${rvtest_DIR}" > ${out_dir}/${Date}_2.Rvtest_list.sh

fi
