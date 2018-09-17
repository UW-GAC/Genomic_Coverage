#!/bin/bash

# Take SNP array annotation and write out a list of SNP positions to query out of 1000G VCF
# vcftools option 
# --positions <filename>
# Include/exclude a set of sites on the basis of a list of positions in a file. Each line of the input file should contain a (tab-separated) chromosome and position. The file should have a header line. 

# operates on gzipped snp annots
# prints only positions, for use with PLINK2 to format vcf to PLINK

imp_dir=$1
schr=$2
echr=$3
snpannot=$4

echo "Using SNP array annotation file ${snpannot}"

# loop through each chr and make chr-specific keep list

for chr in $(seq $schr $echr);do
	echo "Making list of positions for chr $chr"
	
	# specify output file name, clear if exists
	outfn=${imp_dir}/sitelists/snp_positions_chr${chr}.txt
	rm -rf ${outfn}

	# write header
    echo -e "chr\tposition" >> ${outfn}
		
	# Read in SNP array .csv; skip 7 lines + 1 header row
	if [ $chr -eq 23 ]; then
		zcat ${snpannot} | awk -v OFS="\t" -F,  'NR>8 && $10=="X" {print "X", $11}' | sort -nuk 2 > ${outfn}
	else
		zcat ${snpannot} | awk -v OFS="\t" -F,  'NR>8 && $10=='$chr' {print '$chr', $11}' | sort -nuk 2 > ${outfn}
	
	fi

	echo " "
	
done
