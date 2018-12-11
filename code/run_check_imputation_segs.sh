#!/bin/bash
# After imputation has been completed, check for presence of all expected imputation segments - both gprobs and metrics
# Useful for when imputing large datasets, where some jobs may have gotten killed on the cluster due to exceeding available memory 
# NOTE - Jobs can die even after writing out .gprobs.gz and metrics - check for "Killed" messages in the qsub output files	

# SN 12/11/2018

project=$1
imp_dir=$2
schr=$3
echr=$4
array_dir=$5

for chr in $(seq $schr $echr); do
	
	# For the chromosome in question, determine how many 5MB-10MB segments
	# look at the segmentation definition file
	
	chunk_map=${array_dir}/imputation_segments.csv

#	$(awk '$2=='${chunk_curr}'{print $5}' ${chunk_map})
	
	numpart=$(awk -F, '$1=='${chr}' {print $0}' ${chunk_map} | sed 's/,/ /g' | sort -nrk 2 | awk 'NR==1 {print $2}')

	echo -e "Checking chr$chr, for $numpart segments\n"
	
	for set in $(seq 1 $numpart);do	
		
	# check that gprobs file exists, at non-zero size
	fn=${imp_dir}/imputed/${project}_chr${chr}.set${set}.gprobs.gz
	if [ ! -s ${fn} ]; then
	
		echo -e "\tMissing .gprobs for chr${chr}, seg${set}"
	
	fi
	
	# check that metrics file exists, at non-zero size
	fn=${imp_dir}/metrics/${project}_chr${chr}.set${set}.metrics
	if [ ! -s ${fn} ]; then
	
		echo -e "\tMissing .metrics for chr${chr}, seg${set}"
	
	fi
		
	done # close segment loop
	echo " "
done # close chrom loop

## Report the number of strand flips reported by IMPUTE2. A high number may indicate problems with the plus strand alignment in the study data. Note IMPUTE2 can only recognize and reconcile strand misalignments at strand unambiguous SNPs

echo -e "Reporting on number of strand misalignments that IMPUTE2 detected at strand unambiguous SNPs:"
numFlip=$(egrep "flipped strand due to allele mismatch" ${imp_dir}/imputed/*summary |  awk '{ sum+=$9} END {print sum}')
numTot=$(egrep "flipped strand due to allele mismatch" ${imp_dir}/imputed/*summary |  awk '{ sum+=$12} END {print sum}')
echo -e "\t${numFlip} flips out of ${numTot} total study variants\n"