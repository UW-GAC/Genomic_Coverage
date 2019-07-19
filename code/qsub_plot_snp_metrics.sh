#!/bin/bash
# Qsub plotting metrics, comparing across different arrays
# using the generic 'runRscript_array.sh' script to submit the R code - takes care of exporting the R library paths and the R exit code
# SN 1/27/2015

email="sarahcn@uw.edu"

# list of arrays concatenated into 1 string, separated by ":"
arrays=$1
schr=$2
echr=$3 
queue=$4
code_dir=$5
output_dir=$6
resource_dir=$7

# Runs across all panels (AFR, AMR, EAS, EUR, SAS)
qsub -m e -M $email -q $queue -N plotMets -S /bin/bash -j y \
${code_dir}/runRscript_array.sh \
${code_dir}/plot_snp_metrics.R ${arrays} ${schr} ${echr} ${output_dir} ${resource_dir}


