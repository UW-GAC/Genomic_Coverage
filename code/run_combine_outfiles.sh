#!/bin/bash
# combine genomic coverage metrics files

# Define variables
project=$1
imp_dir=$2
chr=$3
numpart=$4

## only one file type - metrics
ft="metrics"

## loop through panels
for p in AFR AMR EAS EUR SAS; do
	
	# Prepare file to hold combined metrics - if it exists, rename it and issue a warning
	# include ancestry group and array name in filename
	comb_file="${imp_dir}/metrics_combined/${project}_metrics_${p}_chr${chr}.csv"
	rm -rf ${comb_file}
	
	if [ -e ${comb_file}.gz ]; then
		echo "WARNING! Combined file exists - renaming"
			mv ${comb_file}.gz ${imp_dir}/metrics_combined/SAVED_${project}_metrics_${p}_chr${chr}.csv.gz
		else
			echo "Starting to combine chr$chr, $ft, $p"
	fi
		
	for chunk in $(seq 1 $numpart); do
		# bypass the header row
		zcat ${imp_dir}/metrics/metrics_${p}_chr${chr}_set${chunk}.csv.gz | awk 'NR>1 {print $0}'  >> $comb_file
	done # close chunk loop
	
	mv $comb_file $comb_file.tmp
	
	## sort by base pair position (should already be sorted, so this is mainly sanity check)
	## need -g due to some bp in scientific notation
	sort -t',' -gk2 $comb_file.tmp > $comb_file.tmp2
	
	## tack on header
	echo "site.id,site.pos,panel.maf,REF,ALT,dos.cor,geno.conc,ma.conc,imp.AA,imp.AB,imp.BB,imp.NA,obs.AA,obs.AB,obs.BB,impute_type" > $comb_file
	
	cat $comb_file.tmp2 >> $comb_file
	
	echo -e "\t\tCounting number of SNPs in combined ${ft} file:"
	wc -l $comb_file
	echo -e "\n"
	
	## gzip output file
	gzip $comb_file	
	
	## remove temporary files
	rm $comb_file.tmp*
	echo -e "\tSee combined files in ${imp_dir}/metrics_combined\n"
	
done # close panel loop