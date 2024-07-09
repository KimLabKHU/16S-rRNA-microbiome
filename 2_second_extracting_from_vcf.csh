#! /bin/csh -f

if ($#argv < 5) then
        echo "USAGE: ./first_extracting_from_vcf vcf_file_lists(sorted_Chr) outdir loci_list(updated) chr_start chr_end"; exit 1
endif

# arguments ##
set vcflists = $argv[1] # 'sorted' list of imputated vcf files
set RESULTS_DIR = $argv[2] # output directory
set loci_list = $argv[3] # list of collected associated loci
set i = $argv[4] # start line of vcf list
set END = $argv[5] # end line of vcf list

# making position file for extracting (chr      pos)
cat ${loci_list} | cut -f 2,3 | sed 1d > ${RESULTS_DIR}/loci_positions_2nd


# extracting associated loci from vcf files ## 
while ($i <= $END)
echo "Now [$i]-th process"

set input = `head -n${i} ${vcflists} | tail -n1` # selecting i-th input
set output = ${RESULTS_DIR}/chr${i}_loci_2nd

# extracting using vcftools
vcftools --gzvcf ${input} \
--out ${output} \
--recode \
--positions ${RESULTS_DIR}/loci_positions_2nd

# extracting using vcftools (dosage)
vcftools --gzvcf ${output} \
--out ${output} \
--extract-FORMAT-info DS

@ i++
end
