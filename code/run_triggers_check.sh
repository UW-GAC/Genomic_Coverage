#!/bin/bash

imp_dir=$1
array=$2
sbatch=$3
ebatch=$4
schr=$5
echr=$6

for batch in $(seq $sbatch $ebatch); do

	echo -e "\n*****************BATCH ${batch}\n"

	make check_imputation project=${array} imp_dir=${imp_dir}/batch${batch} start_chrom=${schr} end_chrom=${echr} array_dir=${imp_dir}

done