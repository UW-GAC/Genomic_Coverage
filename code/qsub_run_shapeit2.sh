#!/bin/bash
# Qsub the run_shapeit2.sh script, with the option of using multiple cores

email=$1
queue=$2
project=$3
imp_dir=$4
ref_dir=$5
schr=$6
echr=$7
test=$8
num_cores=$9
graph=${10}
states=${11}
duoHMM=${12}
par=${13}

echo -e "Using SHAPEIT version2 to phase imputation substrate, from chr $schr to $echr\n"

if [ $test -gt 0 ]; then
	echo -e "\tDoing a test run"
	elif [ $test -eq 0 ] ; then
	echo -e "\tNOTE: about to phase all your samples, is that what you want?"
fi

for chr in $(seq $schr $echr);do
	echo "Starting on chr $chr"

	# remove stdout/stderr file
	rm -rf ${imp_dir}/phased/qsub_shapeit_chr${chr}.out
	
	# -pe mpi ${num_cores} # using "-pe local" flag instead, 6/22/15
	
	  qsub -m as -M $email -q $queue -pe local ${num_cores} \
	  -N phase.${project}${chr} -S /bin/bash -j y \
	  -o ${imp_dir}/phased/qsub_shapeit_chr${chr}.out  \
	  $(dirname $(which $0))/run_shapeit2.sh $project $imp_dir $ref_dir $chr $test $num_cores $graph $states $duoHMM $par

done # close chrom loop