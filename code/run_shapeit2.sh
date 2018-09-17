#!/bin/bash
# Use (multi-threaded) SHAPEIT to phase study data, prior to imputation

project=$1
imp_dir=$2
ref_dir=$3
chr=$4
test=$5
num_cores=$6
graph=$7
states=$8
duoHMM=$9
par=${10}

echo -e "Working on $project, chr$chr\n"
echo -e "Started `date` on $HOSTNAME\n"

## save start time to a variable (timestamp)
START=$(date +%s)

## If phasing chrX, add --chrX flag for SHAPEIT2
chrX_line=""
if [ $chr -eq 23 ]; then
	echo -e "Invoking --chrX flag and --noped\n"
	chrX_line="--chrX --noped"	
fi

## If duoHMM=1, add flag to run with known pedigree structure
# allow for any chr besides X (23)
duoHMM_line=""
if [ $duoHMM -eq 1 -a $chr -ne 23 ]; then
	echo -e "Invoking --duoHMM flag\n"
	duohmm_line="--duohmm"
fi

# if duoHMM is set to 1 but we're on chrX, tell user that we're making sure SHAPEIT2 doesn't try to use the family info (as it can't do so correctly)
if [ $duoHMM -eq 1 -a $chr -eq 23 ]; then
	echo -e "You've asked for -duoHMM, but because it doesn't work on chrX we'll avoid problems by invoking '--noped' instead\n"
	# duohmm_line="--noped"
fi

## Run with basic data checks, to output missingness and Mendelian incons by SNP and by sample
if [ $test -eq 1 ]; then
	echo -e "\tRunning data checks; cannot invoke duohmm\n"
	shapeit -check -B ${imp_dir}/plinkfiles/${project}_chr${chr} \
	-M ${ref_dir}/genetic_maps/genetic_map_chr${chr}_combined_b37.txt \
	-T ${num_cores} ${chrX_line} \
	-L ${imp_dir}/phased/shapeit${logtype}_chr${chr}.log
	exit
fi

## Define execute line based on outputting most likely haplotypes, vs. graphs with haplotype uncertainty
if [ $graph -eq 0 ]; then
	echo -e "\tOutputting most likely haplotypes\n"
	# execute_line="-O ${imp_dir}/phased/${project}_chr${chr}"
	execute_line="-O ${imp_dir}/phased/${project}_chr${chr}.haps.gz ${imp_dir}/phased/${project}_chr${chr}.sample.gz"
	test_line=""
	logtypeA=""
elif [ $graph -eq 1 ]; then
	echo -e "\tOutputting graphs with phasing uncertainty\n"
	execute_line="--output-graph ${imp_dir}/phased/${project}_chr${chr}.hgraph.gz" 
	test_line=""
	logtypeA="graph"
fi

## Define test lines based on running pre-phasing on subsets of the data
## Run all samples on first 20 SNPs
if [ $test -eq 2 ]; then
	endpos=$(awk 'NR==21 {print $4}' ${imp_dir}/plinkfiles/${project}_chr${chr}.bim)
	echo -e "\tStarting a test run with all samples, only up to SNP position $endpos\n"
	test_line="--input-to ${endpos} --output-to ${endpos}"
	logtypeB="_test_snps"
fi

## Run first 20 samples on all SNPs
if [ $test -eq 3 ]; then
	echo -e "\tStarting test run with first 10 samples\n"
	awk 'NR<21 {print $2}' ${imp_dir}/plinkfiles/${project}_chr${chr}.fam > ${imp_dir}/phased/samps_include_fortest_chr${chr}.txt
	test_line="--include-ind ${imp_dir}/phased/samps_include_fortest_chr${chr}.txt"
	logtypeB="_test_samples"
fi

# Production run:
shapeit -B ${imp_dir}/plinkfiles/${project}_chr${chr} \
-M ${ref_dir}/genetic_maps/genetic_map_chr${chr}_combined_b37.txt ${execute_line} ${test_line} \
-S ${states} -T ${num_cores} ${chrX_line} ${duohmm_line} \
-L ${imp_dir}/phased/shapeit${logtypeA}${logtypeB}_chr${chr}.log

echo -e "Finished `date`"

# If on chromX and we're processing PAR, separately pre-phase each PAR region
# pre-phase as if autosome: don't use '--chrX' flag; ok to use default effective population size

if [ $par -eq 1 -a $chr -eq 23 ]; then
	echo -e "\nadditionally prephasing PAR1 and PAR2 for pseudo-autosomal imputation\n"
	
	# allow duoHMM on PAR 
	if [ $duoHMM -eq 1 ]; then
		echo -e "\nInvoking --duoHMM flag for PAR1 and PAR2 pre-phasing\n"
		duohmm_line="--duohmm"
	fi
	
	for chr in X_PAR1 X_PAR2
	do
		echo -e "\tpre-phasing ${chr}\n"
		
		## Define execute line based on outputting most likely haplotypes, vs. graphs with haplotype uncertainty
		if [ $graph -eq 0 ]; then
			echo -e "\tOutputting most likely haplotypes\n"
			execute_line="-O ${imp_dir}/phased/${project}_chr${chr}.haps.gz ${imp_dir}/phased/${project}_chr${chr}.sample.gz"
			test_line=""
			logtypeA=""
		elif [ $graph -eq 1 ]; then
			echo -e "\tOutputting graphs with phasing uncertainty\n"
			execute_line="--output-graph ${imp_dir}/phased/${project}_chr${chr}.hgraph.gz" 
			test_line=""
			logtypeA="graph"
		fi
		
		shapeit -B ${imp_dir}/plinkfiles/${project}_chr${chr} \
		-M ${ref_dir}/genetic_maps/genetic_map_chr${chr}_combined_b37.txt ${execute_line} \
		-S ${states} -T ${num_cores} ${duohmm_line} \
		-L ${imp_dir}/phased/shapeit${logtypeA}${logtypeB}_chr${chr}.log
		
		done # close loop through PAR1, PAR2

fi

## save end time to a variable, calculate elapsed time (rounded number of hours)
END=$(date +%s)
numhours=$((($END-$START)/3600))

echo -e "\nTotal time to run was $numhours hours"

## Note on runtime written to log file: The running time given by SHAPEIT is the cumulative time that the program spends using CPU, and not the elapsed running time (i.e. with mult-threading) 