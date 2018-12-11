#!/bin/bash

imp_dir=$1
schr=$2
echr=$3

for chr in $(seq $schr $echr); do

	# For the chromosome in question, determine how many 5MB-10MB segments
	chunk_map=${imp_dir}/imputation_segments.csv

	numpart=$(awk -F, '$1=='${chr}' {print $0}' ${chunk_map} | sed 's/,/ /g' | sort -nrk 2 | awk 'NR==1 {print $2}')

	echo -e "Checking chr$chr, for $numpart segments\n"

	for set in $(seq 1 $numpart); do	
	
		for panel in AFR AMR EUR EAS SAS; do 
		
			# check that metrics file exists, at non zero size
			
			fn=${imp_dir}/metrics/metrics_${panel}_chr${chr}_set${set}.csv.gz
			
			if [ ! -s ${fn} ]; then
			
				echo -e "\tMissing .metrics for chr${chr}, seg${set}, panel ${panel}"
				echo -e "\t\t${fn}"
			
			fi
			
		done # close panel loop
	done # close segment loop
done # close chr loop