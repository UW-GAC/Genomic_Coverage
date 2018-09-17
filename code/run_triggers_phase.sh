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

	make shapeit2 imp_dir=${imp_dir}/batch${batch} project=${array} start_chrom=${schr} end_chrom=${echr} queue=${queue}

done