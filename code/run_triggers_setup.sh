#!/bin/bash

imp_dir=$1
array=$2
sbatch=$3
ebatch=$4
schr=$5 
echr=$6
queue=$7

for batch in $(seq $sbatch $ebatch); do

	echo -e "\n*****************BATCH ${batch}\n"

	## run these two commands only if batch-specific dir isn't already set up
	if [ ! -d ${imp_dir}/batch${batch} ]; then
		echo -e "\t setting up BATCH ${batch} directory\n"
		make batch_directory imp_dir=${imp_dir} batch=$batch
		
		cp ${imp_dir}/batch_sets/batch_set${batch}.txt ${imp_dir}/batch${batch}/keeplists/samps2mask.txt
	fi
	
	make plinkfiles imp_dir=${imp_dir}/batch${batch} project=${array} start_chrom=${schr} end_chrom=${echr} queue=${queue}
	
done
