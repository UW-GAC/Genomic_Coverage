## Add columns onto IMPUTE2 1000G Phase 3 legend files
## Allows for consideration of minor allele count in > 1 continental panel when defining imputation target variants
# otherwise 'filt_rules_l' argument only excludes in an "OR" fashion - see http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-filt_rules_l

#######################################
rm(list=objects())
options(stringsAsFactors = FALSE, digits=10)

# process command line args
args = commandArgs(trailingOnly=TRUE)
base_dir <- gsub("code","", args[1])
schr <- args[2]
echr <- args[3]

chromosomes <- schr:echr
if(23 %in% chromosomes) chromosomes <- c(chromosomes, "X_PAR1","X_PAR2")

# read in sample information to determine number of samples in each of the super populations
samp_fn <- "../resources/batch_sets/batch_assignment.txt"
samp <- read.table(file=samp_fn, header=TRUE, as.is=TRUE, na.strings="")
dim(samp); head(samp)

length(unique(samp$pop)) 
table(samp$super_pop, useNA='ifany')

## need to link MAF (fields in legend) with minor allele count (desired basis for threshold)
# make a data frame that is one row per group, MAF of 2 copies and then MAF of 4 copies
maf.auto <- as.data.frame(table(samp$super_pop))
names(maf.auto) <- c("group","nsamp")
maf.auto$group <- as.character(maf.auto$group)
# vectorize the calculations
maf.auto$nchr <- 2*maf.auto$nsamp #  for autosomes
maf.auto$maf1copy <- round(1/maf.auto$nchr, 10)
maf.auto$maf2copies <- round(2/maf.auto$nchr, 10)
maf.auto$maf3copies <- round(3/maf.auto$nchr, 10)
maf.auto$maf4copies <- round(4/maf.auto$nchr, 10)
maf.auto

## will use maf.auto for PAR1 and PAR2, where males are diploid

### need alternate maf.dat df for non-PAR chrX
maf.chrx <- maf.auto[,1:2]
male.frq <-  as.data.frame(table(samp$super_pop[samp$sex %in% "male"]))
stopifnot(all.equal(as.character(male.frq$Var1), maf.chrx$group))
maf.chrx$nmale <- male.frq$Freq
# rest are females
maf.chrx$nfemale <- with(maf.chrx, nsamp - nmale)
# females have 2 X chroms, males have 1
maf.chrx$nchr <- with(maf.chrx, 2*nfemale + nmale)
maf.chrx$maf1copy <- round(1/maf.chrx$nchr, 10)
maf.chrx$maf2copies <- round(2/maf.chrx$nchr, 10)
maf.chrx$maf3copies <- round(3/maf.chrx$nchr, 10)
maf.chrx$maf4copies <- round(4/maf.chrx$nchr, 10)
maf.chrx

# collect count by chromosome
chrom.dat <- data.frame(matrix(NA, nrow=length(chromosomes), ncol=5))
names(chrom.dat) <- c("chrom","nvar","n.biallelic.snp",
                     "n.ma.cnt.gte4.allpanels","n.ma.cnt.gte2.allpanels")

chrom.dat$chrom <- chromosomes

message("Note that minor allele frequency in each panel is calculated as frequency of allele that is minor across all samples, so that frequency of consistent allele is given\n\n")

# see for description of legend file columns:
# http://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html

leg_orig <- file.path(base_dir,"resources/refpanel/impute2_fmt")
leg_annot <- file.path(base_dir, "/resources/refpanel/annotated_legend")

for (chrom in chromosomes){
  
  # determine expected chrom count based on autosome vs X chr
  if (is.element(chrom, c(1:22, "X_PAR1","X_PAR2"))) {
    maf.dat <- maf.auto
  } else if (is.element(chrom, 23)) {
    message("adjusting expected MAF values based on males being hemizygous for non-PAR chrX")
    maf.dat <- maf.chrx
  }

  # read in original legend file 
  leg_fn <- paste0(leg_orig, "/1000GP_Phase3_chr",chrom,".legend.gz")
  leg <- read.table(file=leg_fn, header=TRUE, as.is=TRUE)
  nvar <- nrow(leg)
  nsnp <- sum(leg$TYPE=="Biallelic_SNP")
  pct.snp <- round(nsnp/nvar,3)*100
  message(prettyNum(nvar, big.mark=",")," variants on chrom",
          chrom,", ", pct.snp, "% of which are biallelic SNPs\n")

  ## create panel MAF columns from the current AAF (alt allele freq) cols
  panel.cols.idx <- match(c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL"),
                          names(leg))
  names.aaf <- paste(names(leg)[panel.cols.idx],"aaf",sep=".")
  names(leg)[panel.cols.idx] <- names.aaf

  # logical vector to flag where alternate allele != minor allele -- across all samples
  aaf.diff <- leg$ALL.aaf>0.5
  leg$AFR.maf <- ifelse(aaf.diff, 1 - leg$AFR.aaf, leg$AFR.aaf)
  leg$AMR.maf <- ifelse(aaf.diff, 1 - leg$AMR.aaf, leg$AMR.aaf)
  leg$EAS.maf <- ifelse(aaf.diff, 1 - leg$EAS.aaf, leg$EAS.aaf)
  leg$EUR.maf <- ifelse(aaf.diff, 1 - leg$EUR.aaf, leg$EUR.aaf)
  leg$SAS.maf <- ifelse(aaf.diff, 1 - leg$SAS.aaf, leg$SAS.aaf)

  # MAF rounding keeps on causing headaches - convert MAFs to integer allele counts
  leg$afr.cnt <- round(leg$AFR.maf * maf.dat$nchr[maf.dat$group %in% "AFR"], 0)
  leg$amr.cnt <- round(leg$AMR.maf * maf.dat$nchr[maf.dat$group %in% "AMR"], 0)
  leg$eas.cnt <- round(leg$EAS.maf * maf.dat$nchr[maf.dat$group %in% "EAS"], 0)
  leg$eur.cnt <- round(leg$EUR.maf * maf.dat$nchr[maf.dat$group %in% "EUR"], 0)
  leg$sas.cnt <- round(leg$SAS.maf * maf.dat$nchr[maf.dat$group %in% "SAS"], 0)

  # round AAF and MAF columns to 10 decimal places for writing out
  cols.round <- grepl(".maf", names(leg)) | grepl(".aaf", names(leg))
  nms.cols.round <- names(leg)[cols.round]
  leg[,nms.cols.round] <- round(leg[,nms.cols.round],10)

  ### sanity checks
  chk <- unique(leg[,c("afr.cnt","AFR.maf")])
  head(chk[order(chk$afr.cnt),])
  chk <- unique(leg[,c("amr.cnt","AMR.maf")])
  head(chk[order(chk$amr.cnt),])
  chk <- unique(leg[,c("eur.cnt","EUR.maf")])
  head(chk[order(chk$eur.cnt),])
  chk <- unique(leg[,c("eas.cnt","EAS.maf")])
  head(chk[order(chk$eas.cnt),])
  chk <- unique(leg[,c("sas.cnt","SAS.maf")])
  head(chk[order(chk$sas.cnt),])
  # end sanity checks

  # 4 or more copies of the minor allele in any of the 5 panels:
  ma.cnt.gte4.allpanel <- with(leg, afr.cnt >=4 | amr.cnt >=4 |
                               eas.cnt >=4 | eur.cnt >=4 | sas.cnt >=4)

  leg$ma.cnt.gte4.allpanels <- ifelse(ma.cnt.gte4.allpanel, 1, 0)

  # 2 or more copies of the minor allele in any of the 5 panels:
  ma.cnt.gte2.allpanel <-  with(leg, afr.cnt >=2 | amr.cnt >=2 |
                               eas.cnt >=2 | eur.cnt >=2 | sas.cnt >=2)

  leg$ma.cnt.gte2.allpanels <- ifelse(ma.cnt.gte2.allpanel, 1, 0)

  # do some spot checks
  stopifnot(sum(leg$ma.cnt.gte2.allpanels) > sum(leg$ma.cnt.gte4.allpanels))

  # write out annotated legend file
  legout_fn <- paste0(leg_annot, "/1000GP_Phase3_chr",chrom,".legend")

  write.table(leg, file=legout_fn, quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # gzip it
  cmnd <- paste("gzip -f", legout_fn)
  system(cmnd)

  # add totals: total variants; total bi-allelic SNPs; total imp target based on 4 or more copies in any 1 panel; total imp target based on 2 or more copies in any 1 panel
  idx <- which(chrom.dat$chrom==chrom)
  chrom.dat$nvar[idx] <- nvar
  chrom.dat$n.biallelic.snp[idx] <- nsnp
  chrom.dat$ n.ma.cnt.gte4.allpanels[idx] <- sum(ma.cnt.gte4.allpanel)
  chrom.dat$ n.ma.cnt.gte2.allpanels[idx] <- sum(ma.cnt.gte2.allpanel)
  } # end loop through chrom

## write out variant counts
countFn <- paste0(leg_annot,"/ALL_1000G_phase3_impTargets_summarize_chromosomes",schr,"-",echr,".csv")
write.csv(chrom.dat, file=countFn, quote=FALSE, row.names=FALSE)

## Note end time:
cat("All finished, on", format(Sys.time(), "%a %b %d %Y -- %H:%M:%S") ,"\n")
