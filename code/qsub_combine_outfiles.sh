#!/bin/bash
# Qsub the run_combine_outfiles.sh script
	
email=$1
queue=$2
project=$3
imp_dir=$4
schr=$5
echr=$6

# make directory to hold combined files, unless it already exists
if [ ! -d ${imp_dir}/metrics_combined ]; then
  mkdir  ${imp_dir}/metrics_combined
fi

echo "combining per-segment metrics files for chr $schr to $echr"

# loop through different chroms
for chr in $(seq $schr $echr);do
	echo -e "\tStarting on chr $chr"

	## determine number of chunks in the chromosome
	chunk_map=${imp_dir}/imputation_segments.csv

	numpart=$(awk -F, '$1=='${chr}' {print $0}' ${chunk_map} | sed 's/,/ /g' | sort -nrk 2 | awk 'NR==1 {print $2}')
	
	echo -e "\t\tchr$chr imputed in $numpart segments"

	# remove stdout/stderr file
	rm -rf ${imp_dir}/metrics_combined/qsub_combine_${chr}.out

	qsub -m e -M $email -q $queue -N ${project}.comb.chr${chr} -S /bin/bash -j y \
	-o ${imp_dir}/metrics_combined/qsub_combine_${chr}.out \
	$(dirname $(which $0))/run_combine_outfiles.sh $project $imp_dir $chr $numpart

done # close chrom loop