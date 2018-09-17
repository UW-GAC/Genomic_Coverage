#!/bin/bash
# Run IMPUTE2 imputation

# imp dir includes batch# - i.e., MEGArray/batch_iter/batch#

# Define variables
project=$1
imp_dir=$2
ref_dir=$3
maf_column=$4
maf_filt=$5
strand_align=$6
buffer=$7
chr=$8
chunk_size=$9
test=${10}
prephased=${11}

chunk_curr=${SGE_TASK_ID}

echo "For chr$chr, chunk $chunk_curr of size $chunk_size Mb"

chunk_map=../resources/imputation_segments.csv

chunk_start=$(awk -F, '$1=='${chr}' && $2=='${chunk_curr}' {print $5}' ${chunk_map})
chunk_end=$(awk -F, '$1=='${chr}' && $2=='${chunk_curr}' {print $6}' ${chunk_map})

echo -e "\tImputing chr$chr from pos $chunk_start to $chunk_end"
echo -e "\t\tStarted `date` on $HOSTNAME"

## save start time to a variable (timestamp)
START=$(date +%s)

## convert maf filter to integer, if actual maf threshold value
maf_int=$(echo $maf_filt | sed 's/[.]//g')

## Point to current 1000 Genomes Project reference dataset
# 1000G phase 3 file paths
 
ref_hap="${ref_dir}/impute2_fmt/1000GP_Phase3_chr${chr}.hap.gz"
ref_leg="${ref_dir}/annotated_legend/1000GP_Phase3_chr${chr}.legend.gz"

if [ $maf_int -gt 0 ]; then
	echo -e "\t\tfiltering on $maf_column"
	maf_line="-filt_rules_l ${maf_column}<${maf_filt}"
	elif [ $maf_int -eq 0 ]; then	
	echo -e "\t\tno ref maf filters"
	maf_line=""
fi

# TRUE/FALSE flag for aligning strands based on reference MAF
if [ $strand_align -eq 1 ]; then
	echo -e "\t\tInvoking strand alignment"
	strand_line="-align_by_maf_g"
	elif [ $strand_align -eq 0 ]; then
	strand_line=""
fi

# TRUE/FALSE flag for test/debugging run with only 2 samples
if [ $test -eq 1 ]; then
	echo -e "\tStarting test run with 2 samples"
	test_line="-nind 2"
	elif [ $test -eq 0 ]; then
	echo -e "\tNOTE: about to impute all your samples, is that what you want?"
	test_line=""
fi

# Set options for chr X imputation
if [ $chr -eq 23 -a $prephased -eq 0 ]; then
	echo -e "\t\tInvoking additional parameters for chrom X (NOTE: '.sample' file needs 'sex' column)"
	chrX_line="-sample_g ${imp_dir}/genfiles/${project}_chr${chr}.sample -chrX"
	
	elif [ $chr -eq 23 -a $prephased -eq 1 ]; then
	echo "		Invoking additional parameters for chrom X (NOTE: '.sample' file needs 'sex' column)"
	chrX_line="-sample_known_haps_g ${imp_dir}/phased/${project}_chr${chr}.sample.gz -chrX"
	# chrX_line="-chrX"
	
	elif [ $chr != 23 ]; then
	chrX_line=""
fi

# Specify file paths for output files
mets_dir="metrics"
metrics_file="${imp_dir}/${mets_dir}/${project}_chr${chr}.set${chunk_curr}.metrics"
out_dir="imputed"
output_file="${imp_dir}/${out_dir}/${project}_chr${chr}.set${chunk_curr}.gprobs"

# RUN IMPUTE2
#  excluding from the ref panel the test individuals
## syntax to exclude to-be-imputed samples - not needed here as -known_haps_g is already reduced to batch's samples:
# -exclude_samples_g  ${imp_dir}/batch${batch}/batch_sets/samps_exclude_batch${batch}.txt

impute2 -use_prephased_g -m ${ref_dir}/genetic_map_chr${chr}_combined_b37.txt \
-h ${ref_hap} -l ${ref_leg} \
-int ${chunk_start} ${chunk_end} -buffer ${buffer} -allow_large_regions \
-known_haps_g ${imp_dir}/phased/${project}_chr${chr}.haps.gz ${chrX_line} \
-sample_g ${imp_dir}/phased/${project}_chr${chr}.sample.gz \
-os 0 2 -o_gz -exclude_samples_h ${imp_dir}/keeplists/samps2mask.txt \
-sample_h  ${ref_dir}/1000GP_Phase3_mod.sample \
${maf_line} ${strand_line} ${test_line} -o ${output_file} -i ${metrics_file}

echo -e "\n\t\tFinished `date`"

## save end time to a variable, calculate elapsed time (rounded number of hours)
END=$(date +%s)
numhours=$((($END-$START)/3600))

echo "Total time to run was $numhours hours"
 

