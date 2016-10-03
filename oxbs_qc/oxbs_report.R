#!/usr/bin/env Rscript --vanilla

docstring<- paste("DESCRIPTION
Produce report for oxBS-Seq kit

ARGUMENTS
1. Input file typically produced by methylation2mpileup
2. Reference file giving position and modification of the cytosines
3. Output directory
4. Prefix for output files.

USAGE
oxbs_report.R <input> <reference> <output_dir> <output_prefix>

SEE ALSO
...", sep= '')

opts<- commandArgs(trailingOnly = TRUE)
if (length(opts) !=  4){
    stop(docstring)
}
input<- opts[1]
reference<- opts[2]
output_dir<- opts[3]
output_prefix<- opts[4]

if (! file.exists(input) ){
    stop(sprintf("Input file %s does not exist", input))
}
if (! file.exists(reference) ){
    stop(sprintf("Reference file %s does not exist", reference))
}

dir.create(output_dir, showWarnings = FALSE)
prefix<- file.path(output_dir, output_prefix)

# Read data
# ---------

if (!file.exists(reference)){
    stop(sprintf('Reference file "%s" not found.', reference))
}
if (!file.exists(input)){
    stop(sprintf('Input file "%s" not found.', reference))
}

## Read reference file and check consistency:
## ------------------------------------------
synthColNames<- c('sequence_name', 'base_position', 'base_iupac', 'short_description') ## Unused col: 'base_modified'
# synthColClasses<- c('character', 'integer', 'character', 'character') ## 'character', 
synth<- read.table(reference, header= TRUE, stringsAsFactors= FALSE, sep= '\t')
## Check you have all the required columns
if (! all(synthColNames %in% names(synth))){
    stop(sprintf("\n\nReference file %s does not have one or more of the following columns:\n%s\n", reference, paste(synthColNames, collapse= ', ')))
}
## Check position is numeric
synth$base_position<- as.integer(synth$base_position)
if ( any( is.na(synth$base_position) ) ){
    stop(sprintf("\n\nNA introduced in column base_position of reference file %s. Are all the values numeric?\n", reference))
}

## Bases: make all uppercase:
synth$base_iupac<- toupper(synth$base_iupac)

synth<- synth[which(synth$base_iupac == 'C' | synth$base_iupac == 'G'),]

## Read bedgraph file:
## -------------------
colClasses<-  c('character', 'integer', 'integer', 'numeric', 'integer', 'integer', 'character', 'character')
colNames<- c('chrom', 'posStart', 'pos', 'pct.met', 'cnt.met', 'tot_reads', 'strand', 'library_id')

bdg<- read.table(input, sep= '\t', stringsAsFactors= FALSE, nrow= 1)
if(ncol(bdg) == 8){
    ## This is for concatenated input where last col is library_id
    bdg<- read.table(input, sep= '\t', stringsAsFactors= FALSE, colClasses= colClasses, col.names= colNames)
} else if(ncol(bdg) == 7){
    ## This is for output straight from mpileup2methylation
    bdg<- read.table(input, sep= '\t', stringsAsFactors= FALSE, colClasses= colClasses[1:7], col.names= colNames[1:7])
    bdg$library_id<- basename(input)
} else {
    stop(sprintf('Invalid input file: Expected 7 or 8 columns, found %s', ncol(bdg)))
}
## bdg$bs<- ifelse(grepl('oxBS', bdg$library_id), 'oxBS', 'BS')

# Make dataframe with all positions in all libs
# ---------------------------------------------
synth_libs<- synth[rep(seq_len(nrow(synth)), length(unique(bdg$library_id))), ]
synth_libs$library_id<- rep(unique(bdg$library_id), each= nrow(synth))
## Merge to have 0 count positions included
bdg<- merge(x= bdg, y= synth_libs, by.x= c('chrom', 'pos', 'library_id'), by.y= c('sequence_name', 'base_position', 'library_id'), sort= FALSE, all.y= TRUE)

## CEGX 2016-10-03: Remove SQfC from the summaries and plots
# Line below keeps fC but turns values into zeroes
# bdg = subset(bdg, bdg$chrom != 'SQfC')
# Line below gets rid of fC mention altogether
bdg = bdg[!(bdg$chrom == 'SQfC'),]

bdg<- bdg[order(bdg$library_id, bdg$chrom, bdg$pos),]
outf<- paste(prefix, '.oxqc.txt', sep= '')
write.table(bdg[,c('chrom', 'pos', 'pct.met', 'cnt.met', 'tot_reads', 'strand', 'base_iupac', 'short_description'),], file= outf, sep= '\t', col.names= TRUE, row.names= FALSE, quote= FALSE)


# Summary of libraries by control sequence
# ----------------------------------------
bdg$mod<- sub('\\+$|-$', '', bdg$short_description, perl= TRUE)
bdg$dummy<- 'all' ## Just used to have later a column for chrom

#bdg.summary<- aggregate(bdg[, c('pct.met', 'cnt.met', 'tot_reads')],
#    by= list(library_id= bdg$library_id, chrom= bdg$chrom, mod= bdg$mod),
#    mean, na.rm= TRUE)

bdg.summary<- aggregate(bdg[, c('cnt.met', 'tot_reads')],
    by= list(library_id= bdg$library_id, chrom= bdg$chrom, mod= bdg$mod),
    sum, na.rm= TRUE)
bdg.summary$pct.met<- round(100 * (bdg.summary$cnt.met / bdg.summary$tot_reads),2)

bdg.summary.ct<- reshape(bdg.summary[, c('library_id', 'pct.met', 'cnt.met', 'tot_reads', 'chrom', 'mod')], 
                      v.names= c('pct.met', 'cnt.met', 'tot_reads'), 
                      idvar= c('chrom', 'mod'), 
                      timevar= 'library_id', direction= 'wide')

bdg.avg<- aggregate(bdg[, c('cnt.met', 'tot_reads')],
    by= list(library_id= bdg$library_id, chrom= bdg$dummy, mod= bdg$mod),
    sum, na.rm= TRUE)
bdg.avg$pct.met<- round( 100 * (bdg.avg$cnt.met / bdg.avg$tot_reads), 2 )

bdg.avg.ct<- reshape(bdg.avg[, c('library_id', 'chrom', 'pct.met', 'cnt.met', 'tot_reads', 'mod')], 
                      v.names= c('pct.met', 'cnt.met', 'tot_reads'), 
                      idvar= c('chrom', 'mod'), 
                      timevar= 'library_id', direction= 'wide')

bdg.summary.ct<- rbind(bdg.summary.ct[order(bdg.summary.ct$chrom, bdg.summary.ct$mod),], bdg.avg.ct)

outfile<- paste(prefix, '.oxqc_summary.txt', sep= '') ## .oxqc_avg.txt
write.table(bdg.summary.ct, file= outfile, sep= '\t', col.names= TRUE, row.names= FALSE, quote= FALSE)

# Plot percent methylated
# -------------------------------------------
libraries<- unique(bdg$library_id)
chroms<- unique(bdg$chrom)

minx<- min(synth$base_position)
maxx<- max(synth$base_position)

pdf(paste(prefix, '.conversion.pdf', sep= ''), height= 24/2.54, width= 17/2.54, pointsize= 9)
for(lib in libraries){
    par(mfrow= c(length(chroms), 1), las= 2, tcl= -0.2, cex.axis= 0.8, mar= c(1,6,3,0), oma= c(2,0,2,1), cex.lab= 1, mgp= c(2,0.5,0), bty= 'n', xaxs= 'i')
    for(chrom in chroms){
        pdata<- bdg[which(bdg$library_id == lib & bdg$chrom == chrom), ]
        pdata<- merge(pdata, cbind(pos= minx:maxx), all.y= T)
        cols<- ifelse(pdata$short_description == '5mC+', 'aquamarine4',
               ifelse(pdata$short_description == '5mC-', 'aquamarine3',
               ifelse(pdata$short_description == '5hmC+', 'darkred',
               ifelse(pdata$short_description == '5hmC-', 'red',
               ifelse(pdata$short_description == '5fC+', 'dodgerblue',
               ifelse(pdata$short_description == '5fC-', 'lightblue',
               ifelse(pdata$short_description == 'C+', 'grey20',
               ifelse(pdata$short_description == 'C-', 'grey30', 'yellow'))))))))
        plot(type= 'n', x= pdata$pos, xlim=c(minx-1, maxx+1) , ylim= c(0, 100), y= pdata$pct.met, main= NA, ylab= chrom, xlab= '', yaxt= 'n', xaxt= 'n')
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= 'grey90', border= 'transparent')
        grid(nx= NA, ny= NULL, col= 'white', lty= 'solid')
        abline(v= seq(1, max(pdata$pos), by= 10), col= 'white', lty= 'solid')
        rect(xleft=pdata$pos - 0.5, xright=pdata$pos + 0.5, ybottom= 0, ytop=pdata$pct.met, col= cols, border= 'grey90')
        axis(side= 2, cex.axis= 1, col= 'black', col.tick= 'black')
        text(labels= pdata$short_description, x= minx:maxx, y= 0, srt= 90, xpd= NA, adj= c(1.2, 0.5), cex= 0.9)
        mtext(side= 3, line= 0.1, text= minx:maxx, at= minx:maxx, cex= 0.4, las= 3)
    }
}
mtext(text= paste(lib, '\nConversion', sep= ''), outer= TRUE, line= -1.5, cex= 1.2, font= 2, las= 1)

# Coverage
# ------------------------------------------------------------------------------
libraries<- unique(bdg$library_id)
chroms<- unique(bdg$chrom)
pdf(paste(prefix, '.coverage.pdf', sep= ''), height= 24/2.54, width= 17/2.54, pointsize= 9)
for(lib in libraries){
    par(mfrow= c(length(chroms), 1), las= 2, tcl= -0.2, cex.axis= 0.8, mar= c(1,6,3,0), oma= c(2,0,2,1), cex.lab= 1, mgp= c(2,0.5,0), bty= 'n', xaxs= 'i')
    for(chrom in chroms){
        pdata<- bdg[which(bdg$library_id == lib & bdg$chrom == chrom), ]
        pdata<- merge(pdata, cbind(pos= minx:maxx), all.y= T)
        cols<- ifelse(pdata$short_description == '5mC+', 'aquamarine4',
               ifelse(pdata$short_description == '5mC-', 'aquamarine3',
               ifelse(pdata$short_description == '5hmC+', 'darkred',
               ifelse(pdata$short_description == '5hmC-', 'red',
               ifelse(pdata$short_description == '5fC+', 'dodgerblue',
               ifelse(pdata$short_description == '5fC-', 'lightblue',
               ifelse(pdata$short_description == 'C+', 'grey20',
               ifelse(pdata$short_description == 'C-', 'grey30', 'yellow'))))))))
        plot(type= 'n', x= pdata$pos, xlim=c(minx-1, maxx+1) , ylim= c(0, max(bdg$tot_reads[which(bdg$library == lib)], na.rm= TRUE)), y= pdata$tot_reads, main= NA, ylab= '', xlab= '', yaxt= 'n', xaxt= 'n')
        mtext(side= 2, line= 4, text= chrom, las= 0, cex= 0.8)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col= 'grey90', border= 'transparent')
        grid(nx= NA, ny= NULL, col= 'white', lty= 'solid')
        abline(v= seq(1, max(pdata$pos), by= 10), col= 'white', lty= 'solid')
        rect(xleft=pdata$pos - 0.5, xright=pdata$pos + 0.5, ybottom= 0, ytop=pdata$tot_reads, col= cols, border= 'grey90')
        axis(side= 2, cex.axis= 1, col= 'black', col.tick= 'black')
        text(labels= pdata$short_description, x= minx:maxx, y= 0, srt= 90, xpd= NA, adj= c(1.2, 0.5), cex= 0.9)
        mtext(side= 3, line= 0.1, text= minx:maxx, at= minx:maxx, cex= 0.4, las= 3)
    }
}
mtext(text= paste(lib, '\nCoverage', sep= ''), outer= TRUE, line= -1.5, cex= 1.2, font= 2, las= 1)
graphics.off()

quit(save= 'no')
