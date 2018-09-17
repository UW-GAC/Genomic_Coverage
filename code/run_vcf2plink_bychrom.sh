#!/bin/bash
# Create a PLINK version of the 1000G samples on array SNPs

project=$1
dir=$2
chr=$3

## Set location of the 1000G Phase 3(reference panel) VCF files
kgdir=../resources/refpanel/vcf

## example tweaks if vcftools is being finicky:
# export PATH=/projects/resources/software/apps/bin:$PATH
# export LD_LIBRARY_PATH=/projects/resources/software/apps/lib
# export PERL5LIB=/projects/resources/software/apps/vcftools/src/perl/

# VCF file name differs for chrX vs autosomes
vcf_fn=${kgdir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz

if [ $chr -eq 23 ]; then
	vcf_fn=${kgdir}/ALL.chrX.phase3_shapeit2_mvncall_integrated.20130502.genotypes.vcf.gz
fi

# Subset the VCF file - reformat as plink 
vcftools --gzvcf ${vcf_fn} \
--keep ${dir}/keeplists/samps2mask.txt --positions ${dir}/sitelists/snp_positions_chr${chr}.txt \
--plink --out ${dir}/plinkfiles_nosex/${project}_chr${chr}

# continue with conversion to binary PLINK (with unique positions)

## Check for duplicate SNPs at a position - write out the variant name at the second instance to exclude from the final PLINK dataset (SHAPEIT2 chokes when >1 SNP at a given map position)
awk 'seen[$4]++ == 1 awk {print $2}'  ${dir}/plinkfiles_nosex/${project}_chr${chr}.map > ${dir}/keeplists/dupsnp_chr${chr}_exclude.txt

## Update the sex in the plink file; write out as binary PLINK files
plink --noweb --file ${dir}/plinkfiles_nosex/${project}_chr${chr} --update-sex ${dir}/batch_sets/samps_updatesex.txt \
 --exclude ${dir}/keeplists/dupsnp_chr${chr}_exclude.txt --make-bed --out  ${dir}/plinkfiles/${project}_chr${chr}
 
 # remove intermediate .ped and .map files
rm ${dir}/plinkfiles_nosex/${project}_chr${chr}.ped
rm ${dir}/plinkfiles_nosex/${project}_chr${chr}.map

# Note - the final PLINK file uses the base pair position in the name field, if no rsID in the source VCF, i.e.:
# [sarahcn@pearson0 plinkfiles]$ egrep -v rs Omni2.5M_chr22.bim | head
# 22      17796921        0       17796921        TTAAC   T
# 22      18640527        0       18640527        0       G
# 22      19965563        0       19965563        0       C
# 22      21804527        0       21804527        A       AG


echo "Done with chr$chr"
