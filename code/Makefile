# Makefile for imputation-based genomic coverage assessment
# Passing dynamic variables in from command line, in Makefile_Driver
# Variables set via command line in Makefile_Driver will override variable assignments below (!)

################################################################ Define variables
## 1) Define system and environment variables

# Code is in same directory as this Makefile
code_dir = $(shell pwd)

# set user email for cluster notifications
email = sarahcn@uw.edu
# email = <YOUR EMAIL HERE>
# set conditions for cluster notification (recommend against requesting notifications for job start and end, as this will generate thousands of emails from IMPUTE2 jobs)
email_notify = as 

# computer cluster name
queue = all.q
start_chrom = 21
end_chrom = 21

# flag to include pseudo-autosomal regions (PAR1 & PAR2) 
# PAR = 1 will run PAR1 and PAR2 with chromosome X (23)
# PAR = 0 (or any value besides 1) will not run PAR1 and PAR2
PAR = 0

## 2) Describe reference panel
## for phase 3
ref_dir = ../resources/refpanel
## column name in legend file to use for filtering imputation targets
maf_column = ma.cnt.gte2.allpanels
## value below which to exclude SNPs in maf_column (usually "1"=passes filter; "0"=fails filter)
maf_limit = 1

## 3) SHAPEIT2 settings (see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#model)
# number of cores for mulithreading
num_cores = 6
graph = 0
## SHAPEIT-specific test flag, where 0=non-test; 1=run "-check" mode to output missingness and Mendelian incons by snp, by sample; 2=test few SNPs on all samples; 3=test few samples on all SNPs
test_shapeit = 0
states = 200
duoHMM = 0

## 4) IMPUTE2 settings
# to process all segments of your chromosome(s), specify "0" for want_segments
# otherwise, specify "1 10 28" etc. to re-run IMPUTE2 on just those segments
# i.e. if certain segments failed due to memory constraints
want_segments = "0"

# for high memory jobs - specify number of CPUS per job. E.g., if you want only 3 jobs per compute node when you have 12 CPUs/node, set this argument to 4 (12 CPUs on a node, 4 CPUS/job -> 3 jobs per node)
# if set to "-1", will not restrict number of jobs per node - i.e., one job per core
# for genomic coverage - restricting to 6 jobs per node (2 cores per job) after some were dying
cpus_per_job = "2"

chunk_size_mb = 5 
buffer_size_kb = 500
prephased = 1
strand_align = 0

# 5) Options to control testing and report writing
## In general, triggers use "test=0" to mean non-test mode, "test>0" to mean test mode.  See individual scripts for detail
test = 0
## when creating tar file - do line and column counts if verbose=1
verbose = 0

################################################################ Triggers
test:
	echo "Hi, I'm working on batch $(batch)"
	
array_directory:
	mkdir $(imp_dir)
	mkdir $(imp_dir)/sitelists/
	mkdir $(imp_dir)/metrics/
	ln -s ../../resources/batch_sets $(imp_dir)/batch_sets
	ln -s ../../resources/imputation_segments.csv $(imp_dir)/imputation_segments.csv
	
# write out a list of chrom/position on selected SNP array (will comprise imputation basis)
print_array_snps:
	$(code_dir)/make_arraySNP_list$(arrayspec).sh $(imp_dir) $(start_chrom) $(end_chrom) $(array_annot) $(chunk_size)
	
batch_directory:
	mkdir $(imp_dir)/batch$(batch)
	# chgrp -R cidr $(imp_dir)/batch$(batch)
	# chmod 2770 $(imp_dir)/batch$(batch)
	mkdir $(imp_dir)/batch$(batch)/keeplists
	mkdir $(imp_dir)/batch$(batch)/plinkfiles
	mkdir $(imp_dir)/batch$(batch)/plinkfiles_nosex
	mkdir $(imp_dir)/batch$(batch)/phased
	mkdir -p $(imp_dir)/batch$(batch)/imputed
	mkdir -p $(imp_dir)/batch$(batch)/metrics/logfiles
	ln -s ../sitelists $(imp_dir)/batch$(batch)/sitelists
	ln -s ../batch_sets $(imp_dir)/batch$(batch)/batch_sets
	
plinkfiles:
	$(code_dir)/qsub_vcf2plink_bychrom.sh $(email) $(queue) $(project) $(imp_dir) $(start_chrom) $(end_chrom)

shapeit2:
	$(code_dir)/qsub_run_shapeit2.sh $(email) $(queue) $(project) $(imp_dir) $(ref_dir) $(start_chrom) $(end_chrom) $(test_shapeit) $(num_cores) $(graph) $(states) $(duoHMM) $(PAR)
	
impute2:
	$(code_dir)/qsub_impute2.sh $(email) $(queue) $(project) $(imp_dir) $(ref_dir) $(maf_column) $(maf_limit) $(strand_align) $(start_chrom) $(end_chrom) $(chunk_size_mb) $(buffer_size_kb) $(test) $(prephased) $(want_segments) $(cpus_per_job) $(code_dir) $(array_dir)

# check for successful completion of imputation jobs
check_imputation:
	$(curr_dir)/run_check_imputation_segs.sh $(project) $(imp_dir) $(start_chrom) $(end_chrom)
	
# Calculate metrics: concordance and correlation between observed and imputed non-array SNPs
calculate_metrics:
	$(curr_dir)/qsub_run_snp_metrics.sh $(email) $(queue) $(imp_dir) $(project) $(prephased) $(start_chrom) $(end_chrom) $(start_batch) $(end_batch) $(want_segments)

# Combine per-segment metrics files into gz-compressed per-chromosome files (for posting to CIDR URL)
combine_outfiles:
	$(curr_dir)/qsub_combine_outfiles.sh $(email) $(queue) $(project) $(imp_dir)  $(start_chrom) $(end_chrom) 



