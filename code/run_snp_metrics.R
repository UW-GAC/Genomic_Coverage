## Calculate metrics: concordance and correlation between observed and imputed non-array SNPs
## Steps:
## (1) loop through batches of 10% total 1000G panel and creates a combined gprobs
## (2) from matrix of combined gprobs, calculate allelic dosage and most likely genotype
## (3) read in observed 1000G data for given chrom; match up samples and sites with imputed data
## separately by continental panel:
## (4) calculate correlation, concordance, and minor allele concordance  
## (5) keep imputation basis (array) SNPs in metrics file as "controls" (i.e., dos, cor should =1)
## (6) write out final csv of metrics, including chr# and panel in the file name

## Running on whole chroms takes prohibitive amount of memory; running on imputation segments instead
## ...Most segs are 5MB; some centromeric and telomeric segs are 10-20 MB, to get enough type 2 sites instead

###### Historical note - explored different ways to read in oberved 1000G data:
# Options:
# 1 - VCF, usign VariantAnnotation. might work to use expanded VCF instead of collapsed VCF - i.e., may "decompose" multi-allelic vars
# 2 - SeqArray GDS - doesn't solve problem of decomposing multi-allelic vars. but there might be SeqVarTools functions that would help. also Stephanie recently wrote an "alleleDosage" fcn that creates sets of matrices with dosages of specified alleles
# 3 - ultimately decided to use IMPUTE2 haps/legend files, with mutli-allelic vars decomposed into bi-allelic, as it most closely matches what we're imputing. use 'scan' to read in haps files more quickly. paste pairs of haps columns together to get dosage of ref allele

########################################################
options(stringsAsFactors = FALSE)

# runRscript_array should be exporting the R library directory, but not working
# so klugeing for time being:
# lib.path <- "/projects/geneva/gcc-fs2/R_packages/library"

library(GWASTools)
library(stringr)
library(readr)

args <- commandArgs(TRUE)
dir.sel <- args[1]
array <- args[2]
chr <- as.integer(args[3])
phased <- as.integer(args[4])
start.batch <- as.integer(args[5])
end.batch <- as.integer(args[6])

# chrX would require further coding, to deal with "-" character in 1000G haps files
if (chr %in% 23) stop("chrX not yet supported -- need to edit code to deal with '-' character in 1000G haps files")

# segment is last argument - the SGE task ID
seg <- as.integer(args[7])

# set to test=1 to run in test mode
test <- 0

# in non-test mode, read all rows
if (test == 0){
  rows.sel <- -1
}

# For testing, set nrows to small number (e.g, 1K or 5K) to get subset of SNPs
if (test == 1) {
  rows.sel <- 1000
  message("running in test mode on first ", rows.sel," rows of data...\n")
}

# list 1000G superpopulations/continental panels
panels <-  c("AFR", "AMR", "EAS", "EUR", "SAS")

# point to directory containing 1000G data
kgdir <-  "../resources/refpanel/"

# point to sample batch assignment file
batch.set <- read.table(file="../resources/batch_sets/batch_assignment.txt", header=TRUE, as.is=TRUE)

# Read in imputation segment info, to determine coordinates of this segment
chunk.fn <- file.path(dir.sel,"imputation_segments.csv")
chunk.map <- read.csv(chunk.fn, header=TRUE)
bp.start <- chunk.map$bp.start[chunk.map$chrom==chr & chunk.map$segment==seg]
bp.end <- chunk.map$bp.end[chunk.map$chrom==chr & chunk.map$segment==seg]

geno.comb <- NULL
dos.comb <- NULL

cat("Starting", date(),"\n")
cat("Processing chr",chr,"on segment",seg,"from",bp.start,"-",bp.end,"\n")

# Loop over sample batches
for (batch in start.batch:end.batch) {
	cat("Reading in batch",batch,"\n")

	# prepare to read in gprobs
	gprobs.fn <- paste(dir.sel,"/batch",batch,"/imputed/",array,"_chr",chr,".set",seg,".gprobs.gz", sep="")
	n.samps <- sum(batch.set$batch==batch)
	cols.type <- c(rep("character",times=2),
	 "numeric",rep("character",times=2),
	 rep("numeric", times=(3*n.samps)))

	# read in gprobs as data frame
	gpr <- read.table(gzfile(gprobs.fn), colClasses=cols.type,
                          comment.char="", nrows=rows.sel)

	# read in metrics
	cols.type <-  c(rep("character",times=2), "numeric",
                        rep("character",times=2), rep("numeric", times=7))
	mets.fn <-  paste(dir.sel,"/batch",batch,"/metrics/",array,"_chr",chr,".set",seg,".metrics", sep="")
	mets <- read.table(mets.fn, header=TRUE, colClasses=cols.type,
                           comment.char="", nrows=rows.sel, as.is=TRUE)

	# check gprobs and metrics dimensions
	if(!allequal(mets[,2],gpr[,2])) {
          stop("Metrics and .gprobs files do not represent same sites")}

	# for troubleshooting, attach batch number to metrics
	assign(paste("mets",batch,sep=""),mets)

  ## IMPUTE2 metrics file now stores alleles
	## save dataframe of SNP annotation (ID, position, type, alleles A and B)
	# snp.save <- mets[,c(2:5,9)]
		 
	# read in sample information - when phased by batch
    if(phased==1) {
		samp.fn <-  paste(dir.sel,"/batch",batch,"/phased/",array,"_chr",chr,".sample.gz", sep="")
		samp <- read.table(gzfile(samp.fn), skip=2)
      	samp.ids <- samp$V2
      } else message("only set up to process pre-phased by batch")
        
	## double check batch composition
	stopifnot(allequal(sort(samp.ids),
					   sort(batch.set$sample[batch.set$batch %in% batch])))

	## using GWASTools internal functions
	## to calculate dosage and most likely genotype (wrt allele A)
	gpr.matx <- as.matrix(gpr[,-(1:5)])
        ## for testing
	# colnames(gpr.matx) <- rep(samp.ids,each=3,len=ncol(gpr.matx))

	dos.matx <- GWASTools:::.probToDosage(gpr.matx)
	geno.matx <- GWASTools:::.probToGenotype(gpr.matx, prob.threshold=0)
	
	# give sample IDs as column names
	colnames(dos.matx) <- colnames(geno.matx) <- samp.ids

	# give site positions+alleleA+alleleB as rownames				
    mets$identify <- str_trim(paste(format(mets$position,scientific=FALSE),
                            mets$a0, mets$a1,sep="-")) 
	rownames(dos.matx) <- rownames(geno.matx) <- mets$identify

    # compare sites against batch 1
	if (batch>1){
		if(!allequal(mets1$position,mets$position) |
		!allequal(mets1$type,mets$type) |
		!allequal(mets1$a0,mets$a0) |
		!allequal(mets1$a1,mets$a1))
		stop("Batches don't have matching SNP information (name, location, alleles)")
	}

	# Combine batches (most likely genotype and allelic dosages)
	geno.comb <- cbind(geno.comb, geno.matx)
	dos.comb <- cbind(dos.comb, dos.matx)

} # close loop over batches

## save vector of imputed sample IDs, in order they appear in imputed results (geno.comb and dos.comb)
stopifnot(allequal(colnames(geno.comb), colnames(dos.comb)))
imp.samps <- colnames(geno.comb)

## read in 1000G sample info (1 row per sample)
kgsamp.fn <- file.path(kgdir,"impute2_fmt/1000GP_Phase3.sample")
kg.samp <- read.table(kgsamp.fn, header=TRUE, as.is=TRUE, na="")
names(kg.samp)[1] <- "sample"
nkg.samp <- nrow(kg.samp)

## change sort order to align with combined dosage and geno matrices
kg.samp.srt <- kg.samp[match(colnames(geno.comb),kg.samp$sample),]
stopifnot(allequal(kg.samp.srt$sample, colnames(geno.comb)))
stopifnot(allequal(kg.samp.srt$sample, colnames(dos.comb)))

## load in IMPUTE2 1000G legend file and add continental MAFs to mets1 df
kg.leg.fn <- paste0(kgdir,"annotated_legend/1000GP_Phase3_chr",chr,".legend.gz")
# trying to speed up with readr function; read_table didn't get delimiters right
kg.leg <- read_delim(file=kg.leg.fn, delim=" ")

## match variants on position, ref, and alt
kg.leg$identify <- str_trim(paste(format(kg.leg$position,scientific=FALSE),
                            kg.leg$a0, kg.leg$a1,sep="-"))
stopifnot(sum(duplicated(kg.leg$identify)) == 0)
kg.cols <- c("id","position","a0","a1","identify","TYPE","ALL.aaf",
             "AFR.maf", "AMR.maf","EAS.maf","EUR.maf", "SAS.maf")
# use most recent batch's mets columns (we checked earlier that these cols are consistent across batches)
mets.cols <- c("identify","position","a0","a1","type")
snp.comb <- merge(mets, kg.leg[,kg.cols],
                  by="identify", all.x=TRUE, all.y=FALSE, sort=FALSE)

## read in 1000G "truth" genotypes from IMPUTE2 haps file
# IMPUTE2 developers have already done the work of decomposing multi-allelic vars in the VCF to bi-allelic vars in the imputation ref data
# added benefit that IMPUTE2 haps/legends files were used to do the imputation and will thus correspond to the imputed "study" (i.e., test subset of 1000G) data

haps.fn <- paste0(kgdir,"impute2_fmt/1000GP_Phase3_chr",chr,".hap.gz")
# system(paste("zcat", haps.fn, " | wc -l"))

# haps files are very large - use "scan" and 
# only read in rows corresponding to given segment 
# variant dimension of legend file matches variant dimension of hap file
# "mets" below is metrics file of last segment
pos.sel <- kg.leg$position >= bp.start & kg.leg$position <= bp.end 
positions <- kg.leg$position[pos.sel]

# in case start and end positions aren't unique
line.start <- min(which(kg.leg$position == min(positions)))
line.end <- max(which(kg.leg$position == max(positions)))
nskip <- line.start - 1

# determine number of values to read - i.e., number of phased haplotype alleles
ngeno <- 2*nkg.samp*sum(pos.sel)

message("scanning in haps file, lines ", line.start, " through ", line.end)
system.time(haps <- scan(haps.fn, skip=nskip, n=ngeno))

# convert to matrix
hap.matx <- matrix(haps, nrow=sum(pos.sel), ncol=2*nkg.samp, byrow=TRUE)

# sanity check using ALL.aaf from legend file
message("checking ALT allele freq: manual calculation from haps vs. legend file col...\n")
frq.chk <- kg.leg[pos.sel,c("identify","ALL.aaf")]

# sum by row - as "1" indicates 1 copy of alt allele
stopifnot(nrow(frq.chk) == nrow(hap.matx))
frq.chk$alt.cnt <- rowSums(hap.matx)
frq.chk$alt.freq <- frq.chk$alt.cnt/(2*nkg.samp)
stopifnot(allequal(round(frq.chk$ALL.aaf,4), round(frq.chk$alt.freq, 4)))

# convert haplotypes to genotypes (1 col/sample): paste every two columns together
# note this step is sloow
kg.matx <- matrix(paste0(hap.matx[,c(TRUE,FALSE)], hap.matx[,c(FALSE, TRUE)]),
           nrow=sum(pos.sel), ncol=nkg.samp)

# label cols with sample names
colnames(kg.matx) <- kg.samp$sample

# label rows with "identify" variable: position + ref + alt
stopifnot(nrow(kg.matx) == sum(pos.sel))
rownames(kg.matx) <- kg.leg$identify[pos.sel]

# reduce to variants in dos.comb, geno.comb
if(sum(!is.element(rownames(dos.comb), rownames(kg.matx)))>0){
  stop("Not all study variants appear in reference data")}

# convert into dosage of REF allele, as REF is the allele counted in our imputed study data
obs <- kg.matx
obs[kg.matx %in% "00"] <- 2
obs[kg.matx %in% c("10", "01")] <- 1
obs[kg.matx %in% "11"] <- 0

# change to numeric
class(obs) <- "numeric"

# match sample and snp dimension between reference and "study" data
vcf.use <- obs[match(rownames(dos.comb), rownames(obs)),
               match(colnames(dos.comb), colnames(obs))]

# save snapshots of the data
message("\nPreview of truth genotypes:\n")
vcf.use[1:4,1:6]

message("\nPreview of imputed dosages:\n")
dos.comb[1:4, 1:6]

message("\nPreview of imputed genotypes:\n")
geno.comb[1:4, 1:6]

# ...verify that imputation and VCF observed data match in sample and site dimensions:
if (!allequal(rownames(dos.comb), rownames(vcf.use))){
  stop("Matrix of observed genotypes is not aligned to imputed data in variant dimension")}

if(!allequal(colnames(dos.comb), colnames(vcf.use))){
  stop("Matrix of observed genotypes is not aligned to imputed data in sample dimension")}

## separate imputed data into continental/super-population panels
## only calculating metrics for variants polymorphic (MAF>0) in given super population
for (p in panels){

	## flag (T/F) columns (samples) to keep
	panel.samps <- kg.samp$sample[kg.samp$GROUP==p] 
	col.flag <- is.element(colnames(dos.comb), panel.samps) 

    message("Creating metrics for ",sum(col.flag), " ", p," samples")
        
	## flag (T/F) rows (sites) to keep (exclude monomorphs) 
	maf.col <- paste(p,"maf",sep=".")
	row.flag <- snp.comb[,maf.col] > 0 # variant (at least one copy of minor allele) in this panel
	snp.mets <- snp.comb[row.flag,]

	## prune matrices of imputed and observed results to
	## ...(1) sites variant in continental panel (rows) and
	## ...(2) samples in continenal panel (columns)
	geno.panel <- geno.comb[row.flag, col.flag]
	dos.panel <- dos.comb[row.flag, col.flag]
	vcf.panel <- vcf.use[row.flag, col.flag]

	# double check alignment of variants and samples
	stopifnot(allequal(rownames(geno.panel), rownames(vcf.panel)))
	stopifnot(allequal(rownames(geno.panel), rownames(dos.panel)))
	stopifnot(allequal(colnames(geno.panel), colnames(vcf.panel)))
	stopifnot(allequal(colnames(geno.panel), colnames(dos.panel)))                 
        
	# create empty data frame to hold metrics
	met <- data.frame(site.id=snp.mets$rs_id,
        site.pos=snp.mets$position.x, panel.maf=snp.mets[,maf.col],
	REF=snp.mets$a0.x, ALT=snp.mets$a1.x,
	dos.cor=NA, geno.conc=NA, ma.conc=NA,
	imp.AA=NA, imp.AB=NA, imp.BB=NA, imp.NA=NA,
	obs.AA=NA, obs.AB=NA, obs.BB=NA,impute_type=snp.mets$type)

	## loop through imputed sites - calculate metrics 
	n.snp <- sum(row.flag)
	for (i in 1:n.snp){

        # report progress
        if(i %% 1e5 ==0) message("\tworking on variant ", prettyNum(i, big.mark=","))
          
        geno.tmp <- geno.panel[i,]
        dos.tmp <- dos.panel[i,]
        vcf.tmp <- vcf.panel[i,]

        met$imp.AA[i] <- sum(geno.tmp==2, na.rm=TRUE)
        met$imp.AB[i] <- sum(geno.tmp==1, na.rm=TRUE)
        met$imp.BB[i] <- sum(geno.tmp==0, na.rm=TRUE)
		
        # if there is a tie between one or more gprobs, then most likely imputed geno = NA
        met$imp.NA[i] <- sum(is.na(geno.tmp))

        met$obs.AA[i] <- sum(vcf.tmp==2)
        met$obs.AB[i] <- sum(vcf.tmp==1)
        met$obs.BB[i] <- sum(vcf.tmp==0)

        # calculate correlation (can always square it later)
        met$dos.cor[i] <- cor(dos.tmp, vcf.tmp, use="pairwise.complete.obs")
        met$geno.conc[i] <- round(mean(geno.tmp==vcf.tmp, na.rm=TRUE),5)

        # A allele freq from observed VCF
        Po <- mean(vcf.tmp)/2

        # Minor allele concordance, between observed and most likely genotype
        vcf.ma.cnt <- vcf.tmp
        geno.ma.cnt <- geno.tmp

        ## if allele A is not minor allele, convert genotypes to count of minor allele
        if (Po>0.5) {
                vcf.ma.cnt[vcf.tmp==2] <- 0
                vcf.ma.cnt[vcf.tmp==0] <- 2

                geno.ma.cnt[geno.tmp==2] <- 0
                geno.ma.cnt[geno.tmp==0] <- 2
        }

        dat <- as.data.frame(cbind(vcf.ma.cnt, geno.ma.cnt))
        names(dat) <- c("obs","imp")

        ## denominator of minor allele concordance -- pairs where either obs or imp has 1 copy of minor allele, non missing observed
        dat$incl <-  (dat$obs>0 | dat$imp>0) & !is.na(dat$imp)
        den.ma <- sum(dat$incl)
        ## numerator of minor allele concordance -- where the two match
        num.ma <-  sum(dat$obs[dat$incl]==dat$imp[dat$incl])
        met$ma.conc[i] <- round(num.ma/den.ma,5)

        } # close loop through variants

	## save as metrics file, where filename includes chrom and continental panel
    message("\twriting csv for ",p, "\n")
	met.fn <- paste(dir.sel,"/metrics/metrics_",p,"_chr",chr,"_set",seg, ".csv",sep="")
	write.csv(file=met.fn, met, quote=FALSE,row.names=FALSE)

    ## gzip file
    system(paste("gzip -f", met.fn))

} # close loop through continental panel

message("\nFinished ", date(),"\n")
