#!/bin/bash 
 
if [ $# -ne 5 ]; then
	echo "USAGE: ./epacts_input.sh clinic PC taxa_norm OUT_DIR choice (IBD or UC or CD)"
exit 1
else

clinic=$1
PC=$2
taxa_dir=$3
OUT_DIR=$4
choice=$5

mkdir ${OUT_DIR}

for i in $(seq 2 6)
do

if [ $choice == "IBD" ]; then

taxa_norm=${taxa_dir}/level-${i}_vsMF_IBD_nor_count_OR_log2.csv

# Rscript-4.1 /home/rhdfyd/code/epacts/epacts_input.Rscript ${clinic} ${PC} ${taxa_norm} ${OUT_DIR} 1 ${i}

Rscript epacts_input_rank_rinv.Rscript ${clinic} ${PC} ${taxa_norm} ${OUT_DIR} 1 ${i}
sed -i 's/#FID/FID/g' ${OUT_DIR}/IBD_level-${i}_rank_inv_DESEQ.ped

elif [ $choice == "UC" ]; then

taxa_norm=${taxa_dir}/level-${i}_UC_DESEQ_nor_count_filt_HADS.csv

Rscript-4.1 epacts_input_rank_rinv.Rscript ${clinic} ${PC} ${taxa_norm} ${OUT_DIR} 2 ${i}

elif [ $choice == "CD" ]; then

taxa_norm=${taxa_dir}/level-${i}_CD_DESEQ_nor_count_filt_HADS.csv

Rscript-4.1 epacts_input_rank_rinv.Rscript ${clinic} ${PC} ${taxa_norm} ${OUT_DIR} 3 ${i}

else
exit 1
fi
done
Date=`date '+%F  %r' | cut -d ' ' -f 1`
echo "${Date} sh ${clinic} ${PC} ${taxa_dir} ${OUT_DIR} ${choice}" > ${OUT_DIR}/${Date}_1.rvtest_input_pipe.sh
fi