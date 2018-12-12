## Write out csv summary file of metrics by array, by MAF group
## Will be used to plot imputation metrics across different arrays:
## - Using clustered bar charts
## - Separately by test panel set (Phase 3: AFR, AMR, EUR, EAS, SAS)
## - grouping SNPs by MAF bin (<1%, between 1% and 5%, >1%, and >5%)
## - including observed (array SNPs)

options(stringsAsFactors = FALSE)
library(gdsfmt)
library(GWASTools)

args <- commandArgs(TRUE)
array.list <- args[1]
# 1000G panel
p <- args[2]
start.chrom <- as.integer(args[3])
end.chrom <- as.integer(args[4])
out.dir <- args[5]

## to run chroms as specified in Makefile
chrs <- start.chrom:end.chrom
message("\nrunning chroms ", start.chrom," through ", end.chrom, "\n")

# define r2 threshold to use (historically we've used 0.8)
r2.thresh <- 0.8

## define the column name for continental panel MAF
parent.maf <- paste(tolower(p), ".maf",sep="")

## define MAF for 2 copies of minor allele (lower bound for sites we want to include in plots and summaries)
batch.sets <- read.table(file="../resources/batch_sets/batch_assignment.txt",
                         header=TRUE, as.is=TRUE)
panel.dat <- as.data.frame(table(batch.sets$super_pop))
names(panel.dat) <- c("panel","nsamps")
grps <- as.character(panel.dat$panel)

#  calculate MAF when 2 copies of minor allele
panel.dat$maf.2copies <- round(2/(2*panel.dat$nsamps),4)

maf.min <- panel.dat$maf.2copies[panel.dat$panel==p]
message("Two copies of minor allele in ",p,"=",maf.min, " MAF\n")

message("For ",p," 1000G samples\n")

arrays <-  unlist(strsplit(array.list, ":"))
n.arrays <- length(arrays)

# make data frame for array annotations
array.lab <- data.frame(name=arrays)

# define color coding for different 1000G panels
cols <- rainbow(n=length(grps)*length(arrays),v=0.8,s=0.7)

## byrow=TRUE gives (too?) similar color palette for each panel group
cols.mat <- matrix(cols, nrow=length(grps))
rownames(cols.mat) <- grps

array.lab$cols <- cols.mat[rownames(cols.mat)==p]

## Set up MAF categories, for # of categories and MAF boundaries
maf.bins <- data.frame(bin.id=1:4,
                       maf.low=c(maf.min,0.01,0.01,0.05),
                       maf.high=c(0.01,0.05,0.5,0.5))

# construct labels for bins
maf.bins$zero.left.bound <- maf.bins$maf.low==0
maf.bins$right.bound <- maf.bins$maf.high==0.5
maf.bins$lab.0 <- paste(round(maf.bins$maf.low,3)*100, "%<=MAF<",
                        maf.bins$maf.high*100, "%", sep="")
maf.bins$lab.1 <- paste(round(maf.bins$maf.low,3)*100, "%<=MAF", sep="")
maf.bins$lab.2 <-  paste("MAF<=",maf.bins$maf.high*100, "%", sep="")

# Create final label based on lower and upper bounds
maf.bins$label[maf.bins$zero.left.bound] <- maf.bins$lab.2[maf.bins$zero.left.bound] 
maf.bins$label[maf.bins$right.bound] <- maf.bins$lab.1[maf.bins$right.bound]
maf.bins$label[is.na(maf.bins$label)] <- maf.bins$lab.0[is.na(maf.bins$label)]

maf.bins[,c(1:3,9)]

################## Prepare data frame to hold metrics summaries

## Values to report per array (including imputed and observed SNPs)
## Limit to sites with at least TWO copies of minor allele in given panel

# count of total sites (imputed + observed)
array.lab$n.tot <- NA

# count of imputation target sites (as count of observed/imputation basis can be n.tot-n.imputed)
array.lab$n.imputed <- NA

# count of total sites per MAF bin (imputed + observed)
array.lab$n.maf4 <- array.lab$n.maf3 <- array.lab$n.maf2 <- array.lab$n.maf1 <- NA

# count of imputed sites per MAF bin 
array.lab$n.imputed.maf4 <- array.lab$n.imputed.maf3 <- array.lab$n.imputed.maf2 <- array.lab$n.imputed.maf1 <- NA

# count of sites with non/NA dosage r2 (variation in imputed as well as observed data)
array.lab$n.nonmiss.r2.maf4 <- array.lab$n.nonmiss.r2.maf3 <- array.lab$n.nonmiss.r2.maf2 <- array.lab$n.nonmiss.r2.maf1 <- NA

# count of total sites passing r2 threshold, over all MAF bins
# to use as a measure of overall coverage in the line plot with array density on x-axis
array.lab$n.r2pass <- NA

# count of total sites passing r2 threshold  per MAF bin
array.lab$n.r2pass.maf4 <- array.lab$n.r2pass.maf3 <- array.lab$n.r2pass.maf2 <- array.lab$n.r2pass.maf1 <- NA

# mean dosage r2, overall and per MAF bin
array.lab$mnr2 <- NA
array.lab$mnr2.maf4 <- array.lab$mnr2.maf3 <-  array.lab$mnr2.maf2 <- array.lab$mnr2.maf1 <- NA

# mean minor allele concordance, overall and per MAF bin
array.lab$mnMAconc <- NA
array.lab$mnMAconc.maf4 <- array.lab$mnMAconc.maf3 <-  array.lab$mnMAconc.maf2 <- array.lab$mnMAconc.maf1 <- NA

# mean genotype concordance, overall and per MAF bin
array.lab$mnconc <- NA
array.lab$mnconc.maf4 <- array.lab$mnconc.maf3 <-  array.lab$mnconc.maf2 <- array.lab$mnconc.maf1 <- NA

################## Read in metrics, looping through chroms

## Read in imputation segment info, to know how many files to read in per chrom
chunk.fn <- "../resources/imputation_segments.csv"
chunk.all <- read.csv(chunk.fn)
	
## loop through arrays, reading in SNP metrics (metrics files are imputed + observed)
for (a in arrays) {
	message("Processing imputed metrics files for ",a ," array\n")

	mets.allchr <- NULL
	
	## Loop through chroms for array
	for (chr in chrs){

		chunk.map <- chunk.all[chunk.all$chrom %in% chr,]

		message("\tFor chr",chr," in ",max(chunk.map$segment)," segments...\n")

		## loop through segments for chrom
		for (seg in chunk.map$segment){
			met.fn <- paste0("../output/",a,"/metrics/metrics_",p,"_chr",chr,"_set",seg, ".csv.gz")
			met <- read.csv(file=met.fn, as.is=TRUE)
			met$chr <- chr
			met$segment <- seg
			mets.allchr <- rbind(mets.allchr, met)
		}# close segment loop
	} # close chr loop

        # flag "impute_type=0" as imputed
        mets.allchr$imputed <- TRUE
        mets.allchr$imputed[mets.allchr$impute_type > 0] <- FALSE

        # calculate dosage r2 from dosage r
        mets.allchr$dos.r2 <- mets.allchr$dos.cor * mets.allchr$dos.cor

	## retain list of all imputation basis sites, regardless of panel MAF
	## will use to assess how many of array sites are represented in the imputation basis
	all.impbasis.sites <- paste(mets.allchr$chr[!mets.allchr$imputed],
							mets.allchr$site.pos[!mets.allchr$imputed],sep="-")

	## save this list with the array name
	assign(paste(a,".all.impbasis.sites",sep=""),all.impbasis.sites)

	## only keep sites with >= 2 minor alleles in the continental panel (some imputation basis SNPs may be invariant at this point; imputation target SNPs already pruned to those variant in continental panel, when calculating metrics)
	mets <- mets.allchr[round(mets.allchr$panel.maf,4) >= maf.min,]

	## sanity check: print first few MAF values for panel
	mafs.print <- head(sort(unique(round(mets$panel.maf, 4))),4)
	message("Lowest MAF values for sites that will be plotted: ",
                paste(mafs.print, collapse=", "),"\n")

	mets$array <- a

	# create MAF groupings, use T/F flags for MAF group, as groups may not be mutually exclusive
    # create a rounded panel maf
    mets$panel.maf <- round(mets$panel.maf, 4)
        
	mets$maf1 <- mets$panel.maf>=maf.bins$maf.low[1] &
                     mets$panel.maf<maf.bins$maf.high[1] # 2 minor allele copies < MAF <1%
	mets$maf2 <- mets$panel.maf>=maf.bins$maf.low[2] &
		     mets$panel.maf<maf.bins$maf.high[2] # between 1% and 5%
	mets$maf3 <- mets$panel.maf>=maf.bins$maf.low[3] &
		     mets$panel.maf<=maf.bins$maf.high[3] # > 1% MAF
	mets$maf4 <- mets$panel.maf>=maf.bins$maf.low[4] &
		     mets$panel.maf<=maf.bins$maf.high[4] # > 5% MAF

	## save with array name
	# assign(paste("mets",a,sep="."), mets)

	# store info in array.lab data frame
	i <- which(array.lab$name==a)

	array.lab$n.tot[i] <- nrow(mets)
	array.lab$n.imputed[i] <-  sum(mets$imputed)

	## store counts and summary metrics in array.lab data frame
	# count of total sites per MAF bin
	array.lab$n.maf1[i] <- sum(mets$maf1)
	array.lab$n.maf2[i] <- sum(mets$maf2)
	array.lab$n.maf3[i] <- sum(mets$maf3)
	array.lab$n.maf4[i] <- sum(mets$maf4)

	# count of total imputed sites per MAF bin
	array.lab$n.imputed.maf1[i] <- sum(mets$maf1 & mets$imputed)
	array.lab$n.imputed.maf2[i] <- sum(mets$maf2 & mets$imputed)
	array.lab$n.imputed.maf3[i] <- sum(mets$maf3 & mets$imputed)
	array.lab$n.imputed.maf4[i] <- sum(mets$maf4 & mets$imputed)

	# count of total SNPs passing r2 threshold , over all MAF bins
	array.lab$n.r2pass[i] <- sum(mets$dos.r2> r2.thresh, na.rm=TRUE)

	# count of non-missing dos.r2 (some imputed results may lack variation, even though we're restricted to sites with variation in observed genotypes)
	array.lab$n.nonmiss.r2.maf1[i] <- sum(!is.na(mets$dos.r2[mets$maf1]))
	array.lab$n.nonmiss.r2.maf2[i] <- sum(!is.na(mets$dos.r2[mets$maf2]))
	array.lab$n.nonmiss.r2.maf3[i] <- sum(!is.na(mets$dos.r2[mets$maf3]))
	array.lab$n.nonmiss.r2.maf4[i] <- sum(!is.na(mets$dos.r2[mets$maf4]))

	# count of total sites passing r2 threshold, per MAF bin
	array.lab$n.r2pass.maf1[i] <-  sum(mets$dos.r2[mets$maf1] > r2.thresh,na.rm=TRUE)
	array.lab$n.r2pass.maf2[i] <-  sum(mets$dos.r2[mets$maf2] > r2.thresh,na.rm=TRUE)
	array.lab$n.r2pass.maf3[i] <-  sum(mets$dos.r2[mets$maf3] > r2.thresh,na.rm=TRUE)
	array.lab$n.r2pass.maf4[i] <-  sum(mets$dos.r2[mets$maf4] > r2.thresh,na.rm=TRUE)

	# mean dosage r2, overall and per MAF bin
    array.lab$mnr2[i] <- mean(mets$dos.r2, na.rm=TRUE)
	array.lab$mnr2.maf1[i] <- mean(mets$dos.r2[mets$maf1], na.rm=TRUE)
	array.lab$mnr2.maf2[i] <- mean(mets$dos.r2[mets$maf2], na.rm=TRUE)
	array.lab$mnr2.maf3[i] <- mean(mets$dos.r2[mets$maf3], na.rm=TRUE)
	array.lab$mnr2.maf4[i] <- mean(mets$dos.r2[mets$maf4], na.rm=TRUE)

	# mean minor allele concordance, overall and per MAF bin
    array.lab$mnMAconc[i] <- mean(mets$ma.conc,na.rm=TRUE)
	array.lab$mnMAconc.maf1[i] <- mean(mets$ma.conc[mets$maf1],na.rm=TRUE)
	array.lab$mnMAconc.maf2[i] <- mean(mets$ma.conc[mets$maf2],na.rm=TRUE)
	array.lab$mnMAconc.maf3[i] <- mean(mets$ma.conc[mets$maf3],na.rm=TRUE)
	array.lab$mnMAconc.maf4[i] <- mean(mets$ma.conc[mets$maf4],na.rm=TRUE)

	# mean genotype concordance, overall and per MAF bin
    array.lab$mnconc[i] <- mean(mets$geno.conc,na.rm=TRUE)
	array.lab$mnconc.maf1[i] <- mean(mets$geno.conc[mets$maf1],na.rm=TRUE)
	array.lab$mnconc.maf2[i] <- mean(mets$geno.conc[mets$maf2],na.rm=TRUE)
	array.lab$mnconc.maf3[i] <- mean(mets$geno.conc[mets$maf3],na.rm=TRUE)
	array.lab$mnconc.maf4[i] <- mean(mets$geno.conc[mets$maf4],na.rm=TRUE)

	message("Finished reading in metrics files ", date(), "\n\n")

} # close array loop

## save array lab as .csv, where filename indicates both (1) panel and (2) chrs
fn.out <- paste("array_summary_",p,"_chrs",paste(chrs,collapse="."),"_",
                paste(array.lab$name,collapse="-"),".csv",sep="")

## write out metrics summaries
write.csv(array.lab[,setdiff(names(array.lab),"cols")],
          file=file.path(out.dir,fn.out),
		  quote=FALSE, row.names=FALSE)



