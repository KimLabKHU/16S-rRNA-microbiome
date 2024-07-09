#! /bin/csh -f

if ($#argv < 1) then
        echo "USAGE: ./run_metal WD(WD/input/already_exists) name"; exit 1
endif

# arguments #
set WD = $argv[1] # working directory 
set name = $argv[2] # name

set ex_script_u = metal.ex.script.upper # upper part of script
set ex_script_l = metal.ex.script.lower # lower part of script

find ${WD}/*.input | sed 's/^/PROCESS /g' > ${WD}/temp_input_list # input list for script

set my_script = metal.all.script # final script name

cat ${ex_script_u} ${WD}/temp_input_list ${ex_script_l} > ${WD}/${my_script} # merging into one script

# clearing up #
rm ${WD}/*temp*

# running #
set METAL = metal # directory of metal

${METAL} ${WD}/${my_script} # running with script file
mv METAANALYSIS* ${WD}/ # moving to working directory
mv ${WD}/METAANALYSIS1.TBL ${WD}/${name}.meta

cat ${WD}/${name}.meta | awk '$12 > 226' | grep -v "MarkerName" > ${WD}/${name}_sort.meta ## Filtering SNPs based on sample size
cut -f1 ${WD}/${name}_sort.meta | awk -F':' '{OFS="\t"}{print $1,$2}' > ${WD}/${name}_sort.chr_pos
awk '{OFS="\t"}{print $2,$3,$12,$4,$5,$6}' ${WD}/${name}_sort.meta > ${WD}/${name}_meta_result.info
echo -e "chr\tPOS\tA1\tA2\tN\tBETA\tSE\tP" > ${WD}/meta_result_header
paste ${WD}/${name}_sort.chr_pos ${WD}/${name}_meta_result.info > ${WD}/autosome_result_sorted
paste ${WD}/${name}_sort.chr_pos ${WD}/${name}_meta_result.info | cat ${WD}/meta_result_header - > ${WD}/autosome_result_sorted_w_h

