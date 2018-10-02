#!/bin/bash

imp_dir=$1
array=$2
sbatch=$3
ebatch=$4
schr=$5
echr=$6
queue=$7
want_segs=$8
code_dir=$9

for batch in $(seq $sbatch $ebatch); do

	echo -e "\n*****************BATCH ${batch}\n"

	make impute2 imp_dir=${imp_dir}/batch${batch} project=${array} start_chrom=${schr} end_chrom=${echr} queue=${queue} want_segments=${want_segs} code_dir=${code_dir} array_dir=${imp_dir}

done