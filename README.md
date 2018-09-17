This is a pipeline to use genome-wide imputation to assess the genomic coverage of select genotyping arrays. The code materials here will likely require user revision and modification to work on local systems, i.e. this is not meant to "out of the box" work from start to finish.

The place to start is in `code/Makefile_Driver.` Following the steps there will guide you through the remaining code and scripts.

# Requirements

1. Local copy of 1000G Phase 3 VCF files
	* Can be downloaded here: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/]([ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/]
1. Local copy of 1000G Phase 3 in IMPUTE2 reference format
	* Can be downloaded here: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
	* The `./resources/annotated_legend` directory here contains custom annotated versions of the legend files, with a column that will allow imputation of variants with at least two copies in any of the super populations. Users will still need to download the reference panel haplotype files from IMPUTE2 (i.e. ".hap.gz")
1. Genetic map files
	* Can be downloaded here: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
1. Computing cluster
1. Array manifest(s) with chromosome and position in build 37 coordinates
1. Software - local, working installations of:
	* IMPUTE2
	* SHAPEIT2
	* PLINK
	* vcftools
	* R

# References

Nelson, S. C., Doheny, K. F., Pugh, E. W., Romm, J. M., Ling, H., Laurie, C. A., … Laurie, C. C. (2013). Imputation-based genomic coverage assessments of current human genotyping arrays. G3: Genes, Genomes, Genetics, 3(10), 1795–1807. https://doi.org/10.1534/g3.113.007161

Nelson, S. C., Romm, J. M., Doheny, K. F., Pugh, E. W., Laurie, C. C. (2017). Imputation-Based Genomic Coverage Assessments of Current Genotyping Arrays: Illumina HumanCore, OmniExpress, Multi-Ethnic global array and sub-arrays, Global Screening Array, Omni2.5M, Omni5M, and Affymetrix UK Biobank. bioRxiv 150219; doi: https://doi.org/10.1101/150219. 

# Items to check/update once a user downloads/forks code

1. Paths
	1. working directory
	1. 1000G Phase 3 VCF
	1. 1000G Phase 3 IMPUTE2 reference files
	
	1. any compute-cluster specific arguments
		* i.e., email for qsub notification
	1. array information (name, manifest location, possible custom script to parse manifest)

	
Sarah Nelson

UW Genetic Analysis Center

sarahcn@uw.edu
