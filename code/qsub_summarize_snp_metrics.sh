#!/bin/bash
# Qsub plotting metrics, comparing across different arrays
# using the generic 'runRscript_array.sh' script to submit the R code - takes care of exporting the R library paths and the R exit code

# list of arrays concatenated into 1 string, separated by ":"
arrays=$1
schr=$2
echr=$3 
queue=$4
code_dir=$5
output_dir=$6

# parallelize over Phase 3 panels 
for p in AFR AMR EAS EUR SAS; do

	qsub -m e -q $queue -N summarize_${p} -S /bin/bash -j y \
	${code_dir}/runRscript_array.sh \
	${code_dir}/summarize_snp_metrics.R ${arrays} ${p} ${schr} ${echr} ${output_dir}

done # close panel loop
