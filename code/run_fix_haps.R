# For .haps files output by SHAPEIT, edit monomorphic SNPs so that second allele is minor allele, not '0' or major allele repeated

rm(list=objects())
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)

print(args)

dir.sel <- args[1]
array <- args[2]

sbatch <- as.integer(args[3])
ebatch <- as.integer(args[4])
batches <- sbatch:ebatch

# chr as task ID
chr <- as.integer(args[5])

message("Updating .haps files for chromosome",chr," array ", array,"\n")

# loop through batches, to determine both alleles
for (batch in batches) {
  
  # .haps file name
  haps.fn <-  paste(dir.sel,"/batch",batch,"/phased/",array,"_chr",chr,".haps.gz", sep="")
  haps <- read.table(gzfile(haps.fn),as.is=TRUE)
  
  names(haps)[1:5] <- c("snp.id","rs.id","bp","a0","a1")
  haps$mono <- haps$a1==0 | haps$a0 == haps$a1
  
  # count monomorphic
  n.mono <- sum(haps$mono)
  cat("For batch",batch,"on chrom", chr,"there are", n.mono," monomorphic sites\n")
  
  ## Establish first batch as starting point
  if(batch==1){
    haps.ref <- haps[,c("rs.id","bp","a0","a1","mono")]
  } # close if on batch=1
  
  ## for each successive batch, check if we can update the alleles at the monomorphic sites
  ## first verify that haps file matches up in the SNP dimension
  if(batch>1){
    if(!all.equal(haps$rs.id, haps.ref$rs.id)) {stop("Batch ",batch," does not agree in SNP dimension\n")}
    
    # flag sites where currently haps.ref is mono, but current batch can update
    update <- haps.ref$mono & !haps$mono
    
    cat("Batch",batch,"can update alleles at",sum(update),"sites\n")
    
    if(sum(update)>0){
      
      haps.ref$a0[update] <- haps$a0[update]
      haps.ref$a1[update] <- haps$a1[update]
      
      # update the "mono" flag in haps.ref
      haps.ref$mono <- haps.ref$a1==0 | haps.ref$a0 == haps.ref$a1
      
    } # close if in >0 sites to update
  } # close if on batch>1
} # close batch loop

## for any sites still monomorphic in haps.ref - a1 is already set to major allele

## loop through batches again, to write out .haps/gens files with edited a1 at monomorphic sites
## NOTE - must preserve the original 'a1' allele, so that the haplotypes represented by 0/1 will be mapped to the correct allele

for (batch in batches) {
  
  cat("Updating batch",batch,"\n")
  
  # .haps file name
  haps.fn <-  paste(dir.sel,"/batch",batch,"/phased/",array,"_chr",chr,".haps.gz", sep="")
  haps <- read.table(gzfile(haps.fn),as.is=TRUE)
  
  names(haps)[1:5] <- c("snp.id","rs.id","bp","a0","a1")
  
  # isolate the leading columns of current haps file
  haps.dat <- haps[,1:5]
  
  ## identify monomorphic SNPs
  if(sum(haps.dat$a0==0)>0) {stop("batch ",batch," .haps file has 0 in a0 columnn -- unexpected")}
  haps.dat$mono <- haps.dat$a1==0 | haps.dat$a0==haps.dat$a1
  
  # for monomorphs, a0 is the major/fixed allele and a1 is the minor.
  # but verify that this assumption/understanding is correct
  haplotypes.mono <- haps[haps.dat$mono,-c(1:5)]
  sums <- rowSums(haplotypes.mono)
  stopifnot(length(sums) == sum(haps.dat$mono))
  if(max(sums) > 0) {
    stop("haps file violates our assumption that a0 is fixed allele at monomorphs")}
  
  # check that haps.ref matches dimensions of current haps file
  stopifnot(identical(haps.dat$rs.id, haps.ref$rs.id))
  stopifnot(identical(haps.dat$bp, haps.ref$bp))
  
  ## attach the haps.ref alleles
  haps.dat$a0.base <- haps.ref$a0
  haps.dat$a1.base <- haps.ref$a1
  
  ## create slot to hold updated a1 allele
  haps.dat$a1.new <- haps.dat$a1
  
  ## update a1 where monomorphic and a0 matches hap.ref$a0
  fix1 <- haps.dat$mono & haps.dat$a0==haps.dat$a0.base
  haps.dat$a1.new[fix1] <- haps.dat$a1.base[fix1]
  
  ## update a1 where monomorphic and a1 matches hap.ref$a0
  fix2 <- haps.dat$mono & haps.dat$a0==haps.dat$a1.base
  haps.dat$a1.new[fix2] <- haps.dat$a0.base[fix2]
  
  # check that a0 has not been altered
  if(!all.equal(haps.dat$a0, haps$a0)) stop("batch ",batch,"- a0 has been altered")
  
  # update a1 in the .haps file
  haps$a1 <- haps.dat$a1.new
  
  # write out edited haps file
  haps.fn.new <-  paste(dir.sel,"/batch",batch,"/phased/",array,"_chr",chr,".haps", sep="")
  write.table(haps, file=haps.fn.new, row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  
  # rename old file
  haps.fn.old <-  paste(dir.sel,"/batch",batch,"/phased/",array,"_chr",chr,"_premonofix.haps.gz", sep="")
  file.rename(from=haps.fn, to=haps.fn.old)
  
  # gzip new files
  system(paste("gzip",haps.fn.new))
  cat("Finished with batch",batch,"\n")
  
} # close batch loop

