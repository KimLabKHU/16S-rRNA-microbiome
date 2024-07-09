#! /bin/csh -f

if ($#argv < 6) then
        echo "USAGE: ./first_extracting_from_vcf vcf_file_lists(sorted_Chr) outdir loci_list chr_start chr_end samplelist"; exit 1
endif

# arguments ##
set vcflists = $argv[1] # 'sorted' list of imputated vcf files
set RESULTS_DIR = $argv[2] # output directory
set loci_list = $argv[3] # list of collected associated loci
set i = $argv[4] # start line of vcf list
set END = $argv[5] # end line of vcf list
set samplelist=$argv[6]

# outdir ##
mkdir ${RESULTS_DIR}

# making position file for extracting (chr      pos)
# cat ${loci_list} | cut -f 2,3 | sed 1d > ${RESULTS_DIR}/loci_positions_1st
cut -f 2,3 ${loci_list} | sort -g -k1,2 > ${RESULTS_DIR}/loci_positions_1st

# extracting associated loci from vcf files ## 
while ($i <= $END)
echo "Now [$i]-th process"

set input = `head -n${i} ${vcflists} | tail -n1` # selecting i-th input
set output=${RESULTS_DIR}/chr${i}_loci_1st
set DS_input=${RESULTS_DIR}/chr${i}_loci_1st.recode.vcf

# extracting using vcftools (vcf)
vcftools --gzvcf ${input} \
--out ${output} \
--recode \
--keep ${samplelist} \
--positions ${RESULTS_DIR}/loci_positions_1st


# extracting using vcftools (dosage)
vcftools --gzvcf ${DS_input} \
--out ${output} \
--extract-FORMAT-info DS 

if (`grep -v "#" ${RESULTS_DIR}/chr${i}_loci_1st.recode.vcf | wc -l` == 0) then
rm ${RESULTS_DIR}/chr${i}_loci_1st.recode.vcf ${RESULTS_DIR}/chr${i}_loci_1st.DS.FORMAT
endif

@ i++

end
