#! /bin/csh -f

if ($#argv < 3) then
        echo "USAGE: ./merging_data_and_calculating_GRS working_directory(with dosage files) loci_list updated_or_not(Y;updated/N)"; exit 1
endif


# arguments ##
set WD = $argv[1] # working directory
set loci_list = $argv[2] # list of collected associated loci (updated or not)
set update = $argv[3] # updated list or not?


# merging and calculating ##
Rscript GRS_calculation.R ${WD} ${loci_list} ${update}
