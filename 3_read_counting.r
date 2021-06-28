library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx)

bam.files <- Sys.glob("*.bam") # Function to do wildcard expansion (also known as globbing) on file paths. Returns a character vector of matched file paths

# Input: Given a VCF file from previous step, a list of exons (GRanges object with ranges of all exons within UCSC hg19), and hg19 genome
# Output: This function returns a GRanges object, where each entry is a A/G SNP mutation that occurs within an exon.
filter.vcf <- function( vcf.file, exons, genome ){
    vcf <- readVcf( vcf.file, genome ) # Read in Variant Call Format (VCF) files. vcf is a dataframe
# Filter all SNPs found in dbSNP
    vcf <- vcf[ !info(vcf)$DB ] #$DB is a logical field (TRUE or FALSE). ! changes the FALSE value to TRUE. This grabs all the novel SNPs, i.e. removes those found in dbSNP
# Only keep A/G mutations in forward strand and T/C mutations in reverse strand
    mask1 <- rowRanges( vcf )$REF == 'A' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'G' %in% seq ) ) # sapply() function in R Language takes list, vector or data frame as input and gives output in vector or matrix. It is useful for operations on list objects and returns a list object of same length of original set. Syntax: sapply(X, FUN
    mask2 <- rowRanges( vcf )$REF == 'T' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'C' %in% seq ) )
    gr1 <- rowRanges( vcf )[ mask1 ]
    gr2 <- rowRanges( vcf )[ mask2 ]
    strand(gr1) <- '+' # why hard code strand?
    strand(gr2) <- '-'
    a2g.snp <- c( gr1, gr2 ) # combine the 2 Granges objects. I believe these 2 objects are equivalent? Since an A/G transition will lead to a T/C transition on the reverse strand
# Only keep A/G mutations within exons
    a2g.snp <- a2g.snp[ countOverlaps( a2g.snp, exons ) > 0 ]
    return(a2g.snp)
}

vcf.files <- Sys.glob( "*FinalR.vcf" )
# Load all the exons from UCSC gene models
exbygene <- exonsBy( TxDb.Hsapiens.UCSC.hg19.knownGene, "gene" ) # Returns a large compressed GRanges list of all genes within hg19. Grouped by genes.
l <- mclapply( vcf.files, filter.vcf, exbygene, "hg19", mc.cores = 1 ) #exbygene and hg19 are arguments into filter.vcf function
# mclapply is a parallelized version of lapply, it returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X. It relies on forking and hence is not available on Windows unless mc.cores = 1. 

# all A to G and T to C mutations across all vcfs
names(l) <- basename( vcf.files )

# make into a list of GRanges
grl <- GRangesList(l) 
grl
# Rds files store a single R object. According to R documentation: These functions provide the means to save a single R object to a connection (typically a file) and to restore the object, quite possibly under a different name.
# Clicking onto this file automatically generates this R object, so you don't have to rerun all the functions required to generate this object
saveRDS( grl, file = "a2g_snp_filtered_DCD_only.rds" )

# I want unique SNPs only, since there is a possibility of repeated SNPs across the 9 vcfs. Sorted by coordinate
a2g.snp <- sort( unique( unlist( grl, use.names = F ) ) ) 

# Take bam.files, which is a list of all bam file names, and input into function(bf)
# pileup.res output: "piling up" all the reads at each SNP (across all the bam files) we have obtained previously, into a count.
pileup.res <- mclapply( bam.files, function( bf ){
    bf <- BamFile( bf ) # Use BamFile() to create a reference to a BAM file (and optionally its index). The reference remains open across calls to methods, avoiding costly index re-loading.
    # pileup visits each position in the BAM file, subject to constraints implied by which and flag of scanBamParam. For a given position, all reads overlapping the position that are consistent with constraints dictated by flag of scanBamParam and pileupParam are used for counting. Having "which" indicates to only count the reads for those positions in a2g.snp. Additional flags: we want reads that have mapped mates, is part of a proper pair and are not duplicates. We don't care which strand the read is counted from. Max_depth = maximum number of overlapping alignments considered for each position in the pileup. Min_Base_quality = minimum ‘QUAL’ value for each nucleotide in an alignment. 
    res <- pileup( bf, a2g.snp, 
        scanBamParam=ScanBamParam( flag = scanBamFlag( hasUnmappedMate=F, isProperPair=T, isDuplicate=F ), which = a2g.snp ), 
        pileupParam=PileupParam( distinguish_strands=F, min_base_quality=10, max_depth=1e4 ) )
    return( res ) }, mc.cores = 1 )

class(pileup.res) # a list
pileup.res[[1]] # first bam file

names( pileup.res ) <- basename( bam.files ) # renames according to file name

# Saving pileup object for easier future access
saveRDS( pileup.res, file = "pileup_res_DCD_only.rds" )

# if( F ){
pileup.res <- lapply( pileup.res, function( df ){
    df <- dcast( df, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F )
    df$strand <- as.character( strand( a2g.snp ) )
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    stopifnot( all( snp.id == df$which_label ) )
    rownames( df ) <- df$which_label
    df <- split( df, df$strand ) 
    df$`+` <- data.frame( ref = 'A', alt = 'G', ref.count = df$`+`$A, alt.count = df$`+`$G, row.names = rownames( df$`+` ) )
    df$`-` <- data.frame( ref = 'T', alt = 'C', ref.count = df$`-`$T, alt.count = df$`-`$C, row.names = rownames( df$`-` ) )
    df <- rbind( df$`+`, df$`-` )[ snp.id, ]
    return( df )
} )

# Most simplistic allele annotation - just the gene where the SNP is from
idx <- which( countOverlaps( a2g.snp, exbygene ) == 1 )
a2g.snp.subset <- a2g.snp[idx]
ol  <- findOverlaps( a2g.snp.subset, exbygene )
entrez.id <- names( exbygene )[ subjectHits(ol) ]
gene.symbol <- mget( entrez.id, org.Hs.egSYMBOL, ifnotfound = NA )
stopifnot( all( elementNROWS(gene.symbol) == 1 ) )
gene.symbol <- unlist( gene.symbol )
anno <- DataFrame( entrez.id, gene.symbol )

anno$annotation <- "cds"
utr3bytx <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)
utr5bytx <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)
anno$annotation[countOverlaps( a2g.snp.subset, utr5bytx ) > 0] <- 'utr5'
anno$annotation[countOverlaps( a2g.snp.subset, utr3bytx ) > 0] <- 'utr3'

# Should we go back and check the edits found in dbSNP? Maybe some of them are real after all?
# Add back dbSNP annotation
# dbsnp.gr <- endoapply( grl, function( gr ) gr[ grep( '^rs', names(gr) )] )
# dbsnp.gr <- sort( unique( unlist( dbsnp.gr, use.names = F ) ) )

# anno$dbsnp <- ""
# ol <- findOverlaps( a2g.snp.subset, dbsnp.gr )
# anno$dbsnp[ queryHits( ol ) ] <- names( dbsnp.gr )[ subjectHits( ol ) ]

mcols( a2g.snp.subset ) <- anno


pileup.res <- lapply( pileup.res, function( df, a2g.snp ){
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    df <- df[ snp.id, ]
    df <- cbind( as.data.frame( a2g.snp), df )
    df$end <- NULL
    df$width <- NULL
    return( df )
}, a2g.snp.subset )
anno.df <- pileup.res[[1]][(1:(ncol(pileup.res[[1]])-2))]
pileup.res <- lapply( pileup.res, function( df ) df[-(1:(ncol(pileup.res[[1]])-2))] )
pileup.res <- do.call( "cbind", pileup.res )
pileup.res <- cbind( anno.df, pileup.res )

write.xlsx( pileup.res, file = "snp_counts_dedupped.xlsx" )

# }
