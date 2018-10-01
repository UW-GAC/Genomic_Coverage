#!/bin/bash

imp_dir=$1
project=$2
schr=$3
echr=$4
queue=$5
sbatch=$6
ebatch=$7
code_dir=$8

echo $project: editing SHAPEIT2 haps files for chr $schr to $echr

# array job
range="$schr-$echr"

# -m e -M $email # turning off email notifications - 552 per array!
qsub  -m e -q $queue -N fixHaps \
-t ${range} -S /bin/bash -j y \
${code_dir}/runRscript_array.sh \
${code_dir}/run_fix_haps.R $imp_dir $project $sbatch $ebatch

