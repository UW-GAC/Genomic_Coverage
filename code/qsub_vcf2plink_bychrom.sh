#!/bin/bash
# Qsub the creation of per-chrom PLINK files, starting with the binary subject level PLINK

email=$1
queue=$2
project=$3
imp_dir=$4
schr=$5
echr=$6

echo formatting PLINK files for chr $schr to $echr

for chr in $(seq $schr $echr);do
	echo "Starting on chr $chr"

	# remove stdout/stderr file
	rm -rf ${imp_dir}/plinkfiles/qsub_fmt_chr${chr}.out
	
    qsub -M $email -q $queue -N ${project}${chr}.fmt -S /bin/bash -j y \
	-o ${imp_dir}/plinkfiles/qsub_fmt_chr${chr}.out  \
	  $(dirname $(which $0))/run_vcf2plink_bychrom.sh $project $imp_dir $chr

done # close chrom loop
