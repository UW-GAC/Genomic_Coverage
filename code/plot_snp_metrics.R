## Plot imputation metrics - modelled after G3 2013 manuscript
## Use csv summary file written our during main "plot_snp_metrics.R"
## run loop through each of the four continental panels
## make barplots by continental panel, plus line plots by density

rm(list=objects())
options(stringsAsFactors = FALSE)

# library(GWASTools)
library(RColorBrewer)
library(plyr)
library(ggplot2)

## set figure output version
## ft <- "eps"
# ft <- "png"
ft <- "pdf"

pdfdim <- 7 # change back to 7 to get default dimensions
pngdim <- 700 # pixels is default unit

args <- commandArgs(TRUE)
array.list <- args[1]
start.chrom <- as.integer(args[2])
end.chrom <- as.integer(args[3])
output_dir <- args[4]
resource_dir <- args[5]
chrs <- start.chrom:end.chrom

## arguments for interactive use:
# array.list <- "ME_global:ME_eur:ME_amr:HumanCore:OmniExpress:Omni2.5M:Omni5M:Affy_UKBio"
# array.list <- "GSA:ME_global"

# manually setting to chrs 1 and 22 (ones initially tested for 1000G Phase 3 coverage)
# chrs <- c(1,22)

arrays <-  unlist(strsplit(array.list, ":"))
n.arrays <-length(arrays)

panels <- grps <- c("AFR","AMR", "EAS", "EUR", "SAS")

# make data frame for array annotations
array.master <- data.frame(name=arrays)

## add  color palette
## colorblind friendly palette from
## http://my.safaribooksonline.com/book/-/9781449363086/12dot-using-colors-in-plots/recipe_colors_palette_discrete_colorblind_html

# luckily there are 8 arrays n our Phase 3 comparison
## color blind palette
## cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
##           "#0072B2", "#D55E00", "#CC79A7")

# defining different palette for 9 array version (adding GSA, May 2017)

cols <- toupper(c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"))

if(length(cols) < n.arrays) stop("need to define more color blind friendly colors given the number of arrays you have")

cb_palette <- cols[1:n.arrays]

## calculate lower bound of first MAF bin for each panel (2 copies of minor allele)
## define MAF for 2 copies of minor allele (lower bound for sites we want to include in plots and summaries)
batch.sets <- read.table(file=file.path(resource_dir,"batch_sets/batch_assignment.txt"), header=TRUE, as.is=TRUE)
panel.dat <- as.data.frame(table(batch.sets$super_pop))
names(panel.dat) <- c("panel","nsamps")
grps <- panel.dat$panel

#  calculate MAF when 2 copies of minor allele, rounded to 4 decimal places
panel.dat$maf.2copies <- round(2/(2*panel.dat$nsamps),4)

###################################### loop through panels
sum <- NULL
for (p in panels){

    # store the maf associated with 2 copies of minor allele in the panel we're working on
    maf.min <- panel.dat$maf.2copies[panel.dat$panel==p]
    message("\nTwo copies of minor allele in ",p,"=",maf.min, " MAF\n")

    ## Set up MAF categories, for # of categories and MAF boundaries
    maf.bins <- data.frame(bin.id=1:4,
                           maf.low=c(maf.min,0.01,0.01,0.05),
                           maf.high=c(0.01,0.05,0.5,0.5))

    # construct labels for bins
    maf.bins$zero.left.bound <- maf.bins$maf.low==0
    maf.bins$right.bound <- maf.bins$maf.high==0.5
    #maf.bins$lab.0 <- paste(round(maf.bins$maf.low,3)*100, "%<MAF<=",maf.bins$maf.high*100, "%", sep="")
    maf.bins$lab.0 <- paste(round(maf.bins$maf.low,3)*100, "%<MAF<=",maf.bins$maf.high*100, "%", sep="")
    maf.bins$lab.1 <- paste(round(maf.bins$maf.low,3)*100, "%<MAF", sep="")
    maf.bins$lab.2 <-  paste("MAF<",maf.bins$maf.high*100, "%", sep="")

    # Create final label based on lower and upper bounds
    maf.bins$label[maf.bins$zero.left.bound] <- maf.bins$lab.2[maf.bins$zero.left.bound] 
    maf.bins$label[maf.bins$right.bound] <- maf.bins$lab.1[maf.bins$right.bound]
    maf.bins$label[is.na(maf.bins$label)] <- maf.bins$lab.0[is.na(maf.bins$label)]

    ## assign bin labels manually to get <= properly
    left1 <- paste(round(maf.bins$maf.low[1],3)*100,"%<",sep="")
    maf1.lbl <-  bquote(.(left1) ~ MAF<="1%")
    maf2.lbl <- bquote(.("1%<") ~ MAF<="5%")
    maf3.lbl <- "MAF>1%"
    maf4.lbl <- "MAF>5%"

    maf.bins[,c(1:3,9)]

    ## read in csv summary of panel metrics
    fn <- paste(output_dir,"/array_summary_",p,"_chrs",paste(chrs,collapse="."),"_", paste(array.master$name,collapse="-"),".csv",sep="")

    ## manually specify where i've renamed 1-22 as "1-22"
    if(start.chrom %in% 1 & end.chrom %in% 22){
      fn <-  paste(output_dir,"/array_summary_",p,"_chrs1-22_", paste(array.master$name,collapse="-"),".csv",sep="")
      }
    
    array.lab <- read.csv(fn)

    ## save with panel name
    array.lab$panel <- p
    sum <- rbind(sum, array.lab)
    # assign(paste(p,"summary",sep="."), array.lab)

    array.lab$cols <- cb_palette

    ################## Make clustered bar charts

    # make df with 1 row per array-maf group - for labeling points
    # loop through MAF groups
    maf.dats <- rep(NULL, times=9)
    for (i in maf.bins$bin.id) {
            maf.dt <- array.lab[,grep(paste("maf",i,sep=""),names(array.lab))]
            maf.dt$name <- array.lab$name
            names(maf.dt) <- sub(paste(".maf",i,sep=""),"",names(maf.dt))
            maf.dt$maf.group <- i
            maf.dats <- rbind(maf.dats, maf.dt)
    }

    ## Bar chart 1: fraction passing dosage r2>0.8 threshold

    maf1.tot <- format(max(maf.dats$n[maf.dats$maf.group==1]),
                       big.mark=",", scientific=FALSE)
    maf2.tot <- format(max(maf.dats$n[maf.dats$maf.group==2]),
                       big.mark=",", scientific=FALSE)
    maf3.tot <- format(max(maf.dats$n[maf.dats$maf.group==3]),
                       big.mark=",", scientific=FALSE)
    maf4.tot <- format(max(maf.dats$n[maf.dats$maf.group==4]),
                       big.mark=",", scientific=FALSE)

    maf.bins$n.snp <- c(maf1.tot, maf2.tot, maf3.tot, maf4.tot)
    maf.bins$label.cnt <- paste(maf.bins$label,"\n",maf.bins$n.snp,"\n", sep="")

    mafs <- maf.bins$label
    rows <- array.lab$name

    ## includes imputation basis SNPs variant in test panel set as having r2=1
    array.lab$pass4 <- array.lab$pass3 <- array.lab$pass2 <- array.lab$pass1 <- NA

    for (i in 1:nrow(array.lab)){
            array.lab$pass1[i] <- array.lab$n.r2pass.maf1[i]/array.lab$n.nonmiss.r2.maf1[i] 
            array.lab$pass2[i] <- array.lab$n.r2pass.maf2[i]/array.lab$n.nonmiss.r2.maf2[i] 
            array.lab$pass3[i] <- array.lab$n.r2pass.maf3[i]/array.lab$n.nonmiss.r2.maf3[i] 
            array.lab$pass4[i] <- array.lab$n.r2pass.maf4[i]/array.lab$n.nonmiss.r2.maf4[i] 
    }
    r2pct <- as.matrix(array.lab[,c("pass1","pass2","pass3","pass4")])
    colnames(r2pct) <- maf.bins$label
    rownames(r2pct) <- array.lab$name

    # echo r2pct to log file (as the metric isn't written to the .csv summary sheets, only its components: numerator and denominator)
    r2pct

    # add counts back in to label
    colnames(r2pct) <- maf.bins$label.cnt

    ## labels for bar chart:
    # below x-axis label of MAF bin: # of SNPs variant in the parent continental panel (the same across all arrays, as it includes both imputation basis + imputation target SNPs)

    ## base file names
    img.fn <- paste(output_dir,"/barplots_",p,"%01d.",ft,sep="")

    if (ft=="png") {
        png(img.fn, width=pngdim, height=pngdim )
    }

    # multiple plots in pdf will end up as separate pages in 1 pdf file
    if (ft=="pdf") {
        pdf(img.fn, width=pdfdim, height=pdfdim)
    }

    ## G3: figures should be 10-20 cm in width and 1-25 cm in height.
    ## height, width args are in inches
    ## 1 cm=0.39 in
    if (ft=="eps") {
      setEPS()
      postscript(img.fn,horizontal=TRUE, paper="special",height=7,width=7.8, onefile=FALSE)
    }

    # setting for all plots
    par(xpd=NA)

    ## set multiple panels with either mfrow or layout
    par(mfrow=c(2,2))

    # layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE),
    #       heights=c(1,1,0.08))

    ## margins: bottom, left, top , right
    # initial margins: mar=c(4,4,0.5,0.5)
    ## par(mar=c(3,4,2.5,1.5),
    ##     oma=c(1.5,2,1,1.5), ## all outer margins default to 0
    ##     xpd=TRUE, 
    ##     mgp=c(2.5,1,0))
    
    axis.text <- 0.8

    # barplot width
    bwid <- 0.8

    # Barchart 1: pct passing r2>0.8
    barplot(r2pct, beside=TRUE, width=bwid,
            #main="Threshold r2",
            #main=paste(p, "Genomic Coverage, threshold r2"),
            xlab="", ylab=bquote(.("Fraction passing") ~ italic(r^2)>=0.8),
            cex.axis=axis.text, cex.lab=axis.text,cex=axis.text-0.3,
       #     names.arg=c(maf1.lbl, maf2.lbl, maf3.lbl,maf4.lbl),
            col=array.lab$cols,# space=c(0,0.5),xlim=c(1.5,33),
            ylim=c(-0.03,1), xpd=FALSE)
    ## set legend params: make legend with 4 columns
    ## legend("topleft", array.lab$name, fill=array.lab$cols,
    ##         bty="n",ncol=4,cex=0.7)

    ## add ref lines at y axis tick marks
    # grid(nx=NA,ny=NULL,col="lightgray") ## grid will extend into margins when xpd=TRUE, can't override
    yaxp <- par("yaxp")
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=2, col = "lightgray", xpd=FALSE)

    ## label subpart
    text(0,1.1,bquote(bold(.("A. Threshold")) ~bold(italic(r))^bold(italic("2"))))

    # Barchart 2: mean dosage r2
    mnr2 <- as.matrix(array.lab[,grep("mnr2.maf",names(array.lab))])
    colnames(mnr2) <-  maf.bins$label.cnt
    rownames(mnr2) <- array.lab$name

    barplot(mnr2, beside=TRUE, width=bwid,
            #main="Mean r2",
            # main=paste(p, "Genomic Coverage, mean r2"),
            xlab="", ylab=expression(Mean~italic(r^2)),
            cex.axis=axis.text, cex.lab=axis.text,cex=axis.text-0.3,
            col=array.lab$cols,#space=c(0,0.5),xlim=c(1.5,33),
            ylim=c(-0.03,1), xpd=FALSE)
    ## legend("topleft", array.lab$name,fill=array.lab$cols,
    ##        bty="n", ncol=4,cex=0.7)
    yaxp <- par("yaxp")
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=2, col = "lightgray", xpd=FALSE)

    ## label subpart
    text(0,1.1,bquote(bold(.("B. Mean")) ~bold(italic(r))^bold(italic("2"))))

    # Barchart 3: mean minor allele concordance  
    mnMAc <- as.matrix(array.lab[,grep("mnMAconc.maf",names(array.lab))])
    colnames(mnMAc) <-  maf.bins$label.cnt
    rownames(mnMAc) <- array.lab$name

    barplot(mnMAc, beside=TRUE, width=bwid,
            #main="Minor Allele Concordance",
            #main=paste(p, "Genomic Coverage, minor allele concordance"),
            xlab="", ylab="Mean MA concordance",
            cex.axis=axis.text, cex.lab=axis.text,cex=axis.text-0.3,
            col=array.lab$cols,#space=c(0,0.5),xlim=c(1.5,33),
            ylim=c(-0.03,1), xpd=FALSE)
    ## legend("topleft", array.lab$name,fill=array.lab$cols,
    ##        bty="n",ncol=4,cex=0.7)
    yaxp <- par("yaxp")
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=2, col = "lightgray", xpd=FALSE)

    ## label subpart
    text(7.5,1.1,"C. Mean minor allele concordance",font=2)

    # all sites have (non-missing) minor allele concordance (at least 1 copy of minor allele in the observed data)

    # Barchart 4: mean genotype concordance  
    mnc <-  as.matrix(array.lab[,grep("mnconc.maf",names(array.lab))])
    colnames(mnc) <-  maf.bins$label.cnt
    rownames(mnc) <- array.lab$name

    barplot(mnc, beside=TRUE, width=bwid,
            #main="Genotype Concordance",
            #main=paste(p, "Genomic Coverage, genotype concordance"),
            xlab="", ylab="Mean concordance",
            cex.axis=axis.text, cex.lab=axis.text,cex=axis.text-0.3,
            col=array.lab$cols,
            ylim=c(min(maf.dats$mnconc)-0.05,1),xpd=FALSE)

    yaxp <- par("yaxp")
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=2, col = "lightgray", xpd=FALSE)

    ## label subpart
    text(4.5,1.01,"D. Mean concordance",font=2)

    dev.off()

    ## print legend separately
    pdf.fn <- file.path(output_dir,"barplots_legend.pdf")
    ## ## ## create legend, centered at bottom (a short 3rd row for plot, created by layout)
    #par(mar=c(.1,.1,.1,.1))
    pdf(pdf.fn)
    par(mar=c(0,0,0,0))
    barplot(0,0, axes=FALSE,col="white",xlab=NA, ylab=NA)
    legend("center", array.lab$name,fill=array.lab$cols,#xjust=1, yjust=0.5,
           bty="o", box.col="white",
           ncol=5,bg="white", cex=0.7)

    # text(0,0,expression("">=))

    dev.off()

    cat("Finished",p,"barplots\n")
} ## close panel loop


######## Make line plots, by MAF group
## calculate fraction passing threshold > 0.8
sum$pct.pass.maf1 <- sum$n.r2pass.maf1/sum$n.nonmiss.r2.maf1
sum$pct.pass.maf2 <- sum$n.r2pass.maf2/sum$n.nonmiss.r2.maf2
sum$pct.pass.maf3 <- sum$n.r2pass.maf3/sum$n.nonmiss.r2.maf3
sum$pct.pass.maf4 <- sum$n.r2pass.maf4/sum$n.nonmiss.r2.maf4

# Define array density, using overall number of assays on the array
array.names <- unlist(strsplit(array.list, ":"))

## read in count of unique positions represented (vs. unique assays)

array.cnts <- read.csv(file.path(resource_dir,"array_1000Gph3_overlap.csv"))

## grep the n.sites columns
#### ideally rework this to calculate array density based on metrics files (not this external source)
cols.keep <- names(array.cnts)[grep("n.sites$",names(array.cnts))]
mat.cnts <- as.matrix(array.cnts[,cols.keep])
## sum up over chrs  
sums <- apply(mat.cnts, 2, sum)

arrays.short <- tolower(gsub(".n.sites", "", colnames(mat.cnts)))
array.dens <- data.frame(array.short=arrays.short,
                         n.uniq.sites=sums)

# manually set full array names for graph labels
attach(array.dens)
array.dens$array[array.short %in% "meg"] <- "ME_global"
array.dens$array[array.short %in% "meamr"] <- "ME_amr"
array.dens$array[array.short %in% "meeur"] <- "ME_eur"
array.dens$array[array.short %in% "hc"] <- "HumanCore"
array.dens$array[array.short %in% "omexp"] <- "OmniExpress"
array.dens$array[array.short %in% "om2"] <- "Omni2.5M"
array.dens$array[array.short %in% "om5"] <- "Omni5M"
array.dens$array[array.short %in% "afbb"] <- "Affy_UKBio"
array.dens$array[array.short %in% "gsa"] <- "GSA"
array.dens$array[array.short %in% "pmda"] <- "Affy_PMDA"
detach(array.dens)

# rank by density (least to most dense)
array.tmp <- array.dens[order(array.dens$n.uniq.sites),]

# reduce to arrays in the plot (if not all 8)
array.srt <- array.tmp[is.element(array.tmp$array, arrays),]
array.srt$density.rank <- 1:nrow(array.srt)
array.srt

## define color scheme - 1 color for each panel
## using subset of color blind friendly color scheme, again from R graphics cookbook:
## http://my.safaribooksonline.com/book/-/9781449363086/12dot-using-colors-in-plots/recipe_colors_palette_discrete_colorblind_html
tmp_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
# need 5 panels for 5 superpopulations in 1000G phase 3
cb_palette5 <- tmp_palette[c(2:4,6,8)]

panels.dat <- data.frame(panel=panels,color=cb_palette5)

## add density rank to metrics
sum$density.rank <- array.srt$density.rank[match(sum$name, array.srt$array)]
## add unique sites to metrics
sum$n.uniq.sites <- array.srt$n.uniq.sites[match(sum$name, array.srt$array)]
unique(sum[order(sum$density.rank), c("name","density.rank","n.uniq.sites")])

## option to scale x axis to array density
x.axis.scale <- "n.uniq.sites"
## to plot by density rank 
# x.axis.scale <- "density.rank"

col.rename <- which(names(sum)==x.axis.scale)
## make generic name for value to plot on y-axis 
names(sum)[col.rename] <- "rank.plot"

# order by ranking
sum.mini <- sum[order(sum$rank.plot),]

## loop through (1) metrics and (2) MAF groups
## unicode for <= is "\u2265"??
maf.groups <- 1:4
mets <- c("pct.pass","mnr2","mnconc","mnMAconc")

for (met in mets) {
  
  message("plotting ", met,"\n")

  met.name <- sub(".","",met,fixed=TRUE)

  ## set image file and graphical parameters
  img.fn <- paste(output_dir,"/lineplots_",met.name,".",ft,sep="")
  if(ft=="pdf") {pdf(img.fn, width=pdfdim, height=pdfdim)}
  if(ft=="eps") {setEPS(); postscript(img.fn)}
  if(ft=="png") {png(img.fn, width=pngdim, height=pngdim)}

  par(mar=c(4.1, 4.1, 8.1, 2.1),
      #oma=c(2,0.5,0.5,0.5),
      oma=c(1,2,0,1),
      xpd=TRUE,
      mgp=c(2.5,1,0)) # space between axis title and axis labels - defaults to 3,1,0
  #par(mgp=c(3,0,0)) # changes distance between axis and axis labels (units)
  layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), heights=c(1,1,0.08))
  
  for (maf in maf.groups) {

    ## get all columns with this metric, in this maf group
    base.cols <- c("name","panel","rank.plot")
    cols1 <- names(sum.mini)[grep(met, names(sum.mini))]
    cols.get <- cols1[grep(paste("maf",maf, sep=""),cols1)]
    sum.use <- sum.mini[,c(base.cols, cols.get)]

    # strip off "maf#" from column names
    names(sum.use) <- sub(paste(".maf",maf, sep=""),"", names(sum.use))
    sum.plot <- sum.use[,c(met,"rank.plot","panel")]
    # give the metric we're plotting a generic name
    names(sum.plot)[1] <- "metric"

    ## get the range of values of this metric from all 4 maf bins
    cols <- names(sum.mini)[grep(met,names(sum.mini))]
    sum.ref <- sum.mini[,cols]
    yrange <- range(sum.ref)

    ## set y-axis and sub title values based on metric
     if (met=="pct.pass") {yaxis.txt <- bquote(.("Fraction passing") ~ italic(r^2)>=0.8)} 
     if (met=="mnr2") {yaxis.txt <- expression(Mean~italic(r^2))} 
     if (met=="mnconc") {yaxis.txt <- "Mean concordance"}
     if (met=="mnMAconc") {yaxis.txt <- "Mean minor allele concordance"}
 
    ## add EUR line first
    p <- "EUR"
  
    plot(sum.plot$rank.plot[sum.plot$panel==p],
         sum.plot$metric[sum.plot$panel==p], main=NA,
         xlab=NA,ylab=yaxis.txt, ylim=yrange,
         xaxt="n",
         col=panels.dat$color[panels.dat$panel==p],pch=20)

    lines(sum.plot$rank.plot[sum.plot$panel==p],
          sum.plot$metric[sum.plot$panel==p],
          col=panels.dat$color[panels.dat$panel==p], lwd=2)

    # legend for figure panel
    # locate two major y tick mark above end of plot
    yaxp <- par("yaxp")
    y.loc <-  yaxp[2]+3*((yaxp[2]-yaxp[1])/yaxp[3])

    if (maf==1) {label <- bquote(bold(.("A. 2 copies")) <= ~ bold(MAF)<=bold("0.01"))}
    if (maf==2) {label <- bquote(bold(.("B. 0.01<")) ~ bold(MAF)<=bold("0.05"))}
    if (maf==3) {label <-  "C. MAF>0.01"}
    if (maf==4) {label <- "D. MAF>0.05"}
                 
    text(-500000,y.loc,label,font=2,cex=1,adj=0)

    ## suppress normal xaxis - create with MB
    axis(1,at=seq(1e6, 4e6, by=1e6),
         tick=TRUE, labels=1:4,cex.axis=0.9,padj=0)
    mtext(side=1,"Array density (M)",line=2,cex=0.6)
    
    ## add ref lines at x values
    abline(v=array.dens$n.uniq.sites,lty=2, col="lightgray", xpd=FALSE)

    #### add array names at top
    # create an x axis with ranks and array names
    x.vals <- unique(sum.plot$rank.plot)
   text(x.vals, rep(1, times=length(x.vals)),array.srt$array,
        cex=0.8,pos=4,adj=0,srt=60, offset=-0.5)
    
    ## add ref lines at y axis tick marks
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=2, col = "lightgray", xpd=FALSE)

    ## loop through remaining panels
    for (p in c("AFR", "AMR", "EAS", "SAS")) {
      points(sum.plot$rank.plot[sum.plot$panel==p],
            sum.plot$metric[sum.plot$panel==p],
            pch=20,col=panels.dat$color[panels.dat$panel==p])

      lines(sum.plot$rank.plot[sum.plot$panel==p],
          sum.plot$metric[sum.plot$panel==p],
          col=panels.dat$color[panels.dat$panel==p], lwd=2)

    } # close panel

  } ## close MAF group loop
  
  ## make legend to map colors to panels
    par(mar=c(.1,.1,.1,.1)) 
    plot(0,0, axes=FALSE,pch="",xlab=NA,ylab=NA)
    legend("bottom",title=NULL,
           panels.dat$panel, fill=panels.dat$color, horiz=TRUE,bty="n")
  
  dev.off()
} ## close metric loop