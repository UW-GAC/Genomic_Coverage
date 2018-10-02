#!/bin/bash
# Qsub the run_impute2.sh script

email=$1
queue=$2
project=$3
imp_dir=$4
ref_dir=$5
maf_column=$6
maf_limit=$7
strand_align=$8
schr=${9}
echr=${10}
chunk_size=${11}
buffer=${12}
test=${13}
prephased=${14}
segments=${15}
cpus_per_job=${16}
code_dir=${17}
array_dir=${18}

out_dir="imputed"
task="impute"

echo $PWD

echo "submitting imputation segments for chr $schr to $echr"

if [ $test -eq 1 ]; then
	echo -e "\tJust running a test with 2 samples"
	elif [ $test -eq 0 ] ; then
	echo -e "\tNOTE: about to impute all your samples, is that what you want?"
fi

if [ $prephased -eq 1 ]; then
	echo -e "\tUsing prephased data from SHAPEIT"
	job_name="p-imp"
	elif [ $prephased -eq 0 ] ; then
	echo -e "\tImputing into unphased data"
	job_name="imp"
fi

pe_local_arg=""
if [ "$cpus_per_job" -ne -1 ]; then
	echo -e "\tasking for ${cpus_per_job} CPUs per job\n"
	pe_local_arg="-pe local ${cpus_per_job}"
fi

for chr in $(seq $schr $echr);do
	echo "Starting on chr $chr"

	## determine number of chunks in the chromosome
	chunk_map=$array_dir/imputation_segments.csv
	
	echo -e "\t looking at ${chunk_map} for segment information"

	numpart=$(awk -F, '$1=='${chr}' {print $0}' ${chunk_map} | sed 's/,/ /g' | sort -nrk 2 | awk 'NR==1 {print $2}')

	echo -e "\t\tImputing chr$chr in ${numpart} segments"

	# if segments=0, then will impute all segs as array jobs
	# if segments != 0, will impute each segment in arg as separate job
	for set in ${segments} ; do

		if [ "$set" -eq 0 ]; then
			echo -e "\tImputing chr$chr in ${numpart} segments"
			range="1-"${numpart}		
		else
			echo -e "\tstarting segment ${set}"
			range=${set}			
		fi
			
		# to turn email notifications back on: -m e -M ${email}
		 qsub -q $queue -N ${job_name}${chr}.${project} \
		 -S /bin/bash -j y -t ${range} -cwd $pe_local_arg \
		 $code_dir/run_impute2.sh $project $imp_dir $ref_dir $maf_column $maf_limit $strand_align $buffer $chr $chunk_size $test $prephased $array_dir

	 done # close segment loop
	echo ""  # insert blank line at end of chrom's submission
done # close chrom loop