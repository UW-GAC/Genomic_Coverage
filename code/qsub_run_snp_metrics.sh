#!/bin/bash
# using the generic 'runRscript_array.sh' script to submit the R code - takes care of exporting the R library paths and the R exit code
# SN 1/27/2015

email=$1
queue=$2
imp_dir=$3
project=$4
prephased=$5
schr=$6
echr=$7
sbatch=$8
ebatch=$9
segments=${10}
code_dir=${11}

echo calculating metrics for chr $schr to $echr

# for high memory jobs - specify number of CPUS per job. E.g., if you want only 3 jobs per bigmem node (12 CPUs/node), set this argument to 4 (12 CPUs on a node, 4 CPUS/job -> 3 jobs per node)
# if set to "-1", will not restrict number of jobs per node - i.e., one job per core
cpus_per_job="6"

pe_local_arg=""

if [ "$cpus_per_job" -ne -1 ]; then
	echo -e "\tasking for ${cpus_per_job} CPUs per job\n"
	pe_local_arg="-pe local ${cpus_per_job}"
fi

for chr in $(seq $schr $echr);do
	echo "Starting on chr $chr"
	
	# For the chromosome in question, determine how many 5MB-10MB segments were required to impute
	chunk_map=${imp_dir}/imputation_segments.csv

	numpart=$(awk -F, '$1=='${chr}' {print $0}' ${chunk_map} | sed 's/,/ /g' | sort -nrk 2 | awk 'NR==1 {print $2}')
	
	# if segments=0, then will run all segs as array jobs
	# if segments != 0, will run each segment in arg as separate job
	# sge task ID will = segment #, final argument to R script
	for set in ${segments} ; do

		if [ "$set" -eq 0 ]; then
			echo -e "\tProcessing chr$chr in ${numpart} segments"
			range="1-"${numpart}		
		else
			echo -e "\tstarting segment ${set}"
			range=${set}			
		fi

		# -m e -M $email # turning off email notifications
		qsub -q $queue -N mets_chr${chr} -S /bin/bash -j y -t ${range} \
		${pe_local_arg} ${code_dir}/runRscript_array.sh \
		${code_dir}/run_snp_metrics.R $imp_dir $project $chr $prephased $sbatch $ebatch

		done # close set loop
	echo -e "\n"
done # close chrom loop

