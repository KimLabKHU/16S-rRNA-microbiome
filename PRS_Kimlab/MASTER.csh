#! /bin/csh -f

if ($#argv < 3) then
        echo "USAGE: ./MASTER.csh vcflists working_directory loci_list "; exit 1
endif


# arguments ##
set vcflists = $argv[1]
set working_directory = $argv[2] 
set loci_list = $argv[3] 

set i = 1 # first chromosome
set END = 22 # last chromosome


echo "---- first_extracting_from_vcf ----"
./1_first_extracting_from_vcf.csh ${vcflists} ${working_directory} ${loci_list} ${i} ${END}

echo "---- merging_data_and_calculating_GRS ----"
./3_merging_data_and_calculating_GRS.csh ${working_directory} ${loci_list} N
