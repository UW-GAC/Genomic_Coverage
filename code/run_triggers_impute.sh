#!/bin/bash

imp_dir=$1
array=$2
sbatch=$3
ebatch=$4
prephased=$5
schr=$6
echr=$7
queue=$8
want_segs=$9

for batch in $(seq $sbatch $ebatch); do

	echo -e "\n*****************BATCH ${batch}\n"

	make impute2 imp_dir=${imp_dir}/batch${batch} project=${array} prephased=${prephased} start_chrom=${schr} end_chrom=${echr} queue=${queue} want_segments=${want_segs}

done